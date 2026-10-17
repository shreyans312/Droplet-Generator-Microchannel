#!/usr/bin/env python3
"""Audit a reconstructed rigid or prescribed-geometry VOF case without calibration.

Flux integration from saved fields is a sampling diagnostic, not a substitute
for a solver-time-step conservation history. Run with /usr/bin/python3 (VTK).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "cases/flexible-rigid"))
from measure import _internal_mesh


def scalar_values(path: Path, patch: str | None = None) -> np.ndarray:
    """Read an OpenFOAM scalar list, respecting the field's binary architecture."""
    data = path.read_bytes()
    if patch is None:
        offset = data.index(b"internalField")
    else:
        offset = data.index(b"boundaryField")
        match = re.search(rb"\b" + patch.encode() + rb"\s*\{", data[offset:])
        if match is None:
            raise ValueError(f"missing patch {patch}: {path}")
        offset += match.end()
    match = re.search(rb"(?:internalField|value)\s+uniform\s+([-+0-9.eE]+)\s*;", data[offset:])
    nonuniform = re.search(rb"(?:internalField|value)\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(", data[offset:])
    if match is not None and (nonuniform is None or match.start() < nonuniform.start()):
        # All uniform boundary flux fields used here are zero; nonzero scalar
        # values need a face count before they can be integrated.
        value = float(match.group(1))
        if patch is not None and value != 0:
            raise ValueError(f"nonzero uniform patch flux requires face count: {path}")
        return np.array([value])
    if nonuniform is None:
        raise ValueError(f"missing scalar list: {path} {patch}")
    count = int(nonuniform.group(1))
    start = offset + nonuniform.end()
    if re.search(rb"format\s+binary", data[:data.index(b"dimensions")]):
        arch = re.search(rb'arch\s+"(LSB|MSB);label=\d+;scalar=(32|64)"', data)
        # Foundation v13 writes native binary without an arch entry. These
        # local cases were produced by the native double-precision solver.
        dtype = np.dtype("=f8") if arch is None else np.dtype(("<" if arch.group(1) == b"LSB" else ">") + ("f8" if arch.group(2) == b"64" else "f4"))
        if data[start + count * dtype.itemsize:start + count * dtype.itemsize + 1] != b")":
            raise ValueError(f"binary scalar width or count mismatch: {path}")
        return np.frombuffer(data, dtype=dtype, count=count, offset=start).astype(float)
    end = data.index(b")", start)
    values = np.fromstring(data[start:end].decode(), sep=" ")
    if len(values) != count:
        raise ValueError(f"scalar count mismatch: {path}")
    return values


def mesh_metrics(case: Path, time: float) -> tuple[np.ndarray, float]:
    marker = case / f"{case.name}.foam"
    marker.touch()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(marker))
    reader.UpdateInformation()
    reader.UpdateTimeStep(time)
    mesh = _internal_mesh(reader.GetOutput())
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(mesh)
    sizes.ComputeVolumeOn()
    sizes.Update()
    volumes = vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Volume")).copy()
    # Cases generated here use orthogonal hexahedral block meshes.
    minimum = math.inf
    for i in range(mesh.GetNumberOfCells()):
        bounds = mesh.GetCell(i).GetBounds()
        minimum = min(minimum, *(bounds[k + 1] - bounds[k] for k in (0, 2, 4)))
    return volumes, minimum


def statistics_for(rows: list[dict[str, str]], metric: str = "axial_diameter_um") -> dict[str, float | int | None]:
    times = [float(r["crossing_time_s"]) for r in rows]
    diameters = [float(r[metric]) for r in rows]
    mean = statistics.mean(diameters) if diameters else None
    sd = statistics.stdev(diameters) if len(diameters) > 1 else None
    return {
        "count": len(rows), "diameter_mean_um": mean, "diameter_sd_um": sd,
        "diameter_cv_percent": 100 * sd / mean if mean and sd is not None else None,
        "frequency_hz": (len(times) - 1) / (times[-1] - times[0]) if len(times) > 1 else None,
    }


def step_conservation(case: Path, flow_rate: float) -> dict | None:
    def history(name: str) -> dict[float, tuple[float, ...]]:
        rows = {}
        paths = [path for directory in (case / "postProcessing").glob(name + "*") for path in directory.glob("*/*.dat")]
        for path in sorted(paths):
            for line in path.read_text().splitlines():
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                values = tuple(map(float, line.split()))
                rows[values[0]] = values[1:]
        return rows
    inventories, fluxes = history("phaseVolume"), history("boundaryFlux")
    times = sorted(set(inventories) & set(fluxes))
    if len(times) < 2:
        return None
    time_array = np.array(times)
    dt = np.diff(time_array)
    volume = np.array([inventories[t][0] for t in times])
    flux = np.array([fluxes[t][0] for t in times])
    residual = np.diff(volume) + dt * flux[1:]
    supplied = flow_rate * (times[-1] - times[0])
    return {
        "samples": len(times), "first_time_s": times[0], "last_time_s": times[-1],
        "inventory_change_nl": float(volume[-1] - volume[0]) / 1e-12,
        "integrated_net_outward_flux_nl": float(dt @ flux[1:]) / 1e-12,
        "residual_percent_of_injected_phase": 100 * float(residual.sum()) / supplied,
        "absolute_accumulated_residual_percent_of_injected_phase": 100 * float(abs(residual).sum()) / supplied,
        "precision_limitation": "Absolute accumulated residual can include per-step rounding noise from inherited low-precision history streams; inspect signed residual and history precision together.",
        "integration": "backward-step phase flux, consistent with Euler; history must cover every solver step",
    }


def audit(case: Path) -> dict:
    metadata = json.loads((case / "case.json").read_text())
    end_match = re.search(r"\bendTime\s+([0-9.eE+-]+)", (case / "system/controlDict").read_text())
    target = float(end_match.group(1))
    times = []
    for entry in case.iterdir():
        try:
            t = float(entry.name)
        except ValueError:
            continue
        if entry.is_dir() and (entry / "alphaPhi.water").exists():
            times.append((t, entry))
    times.sort()
    if not times:
        raise RuntimeError("reconstructed phase-flux fields are required")
    volumes, minimum = mesh_metrics(case, times[-1][0])
    sigma = float(metadata["surface_tension_n_m"])
    capillary = math.sqrt((1630 + 1060) * minimum**3 / (4 * math.pi * sigma)) if sigma > 0 else math.inf
    patches = ("dispersedInlet", "continuousUpper", "continuousLower", "outlet", "walls", "bottom", "top")
    history = []
    for t, directory in times:
        alpha = scalar_values(directory / "alpha.water")
        if len(alpha) != len(volumes):
            raise ValueError("alpha and mesh cell ordering/count mismatch")
        phase_flux = {p: float(scalar_values(directory / "alphaPhi.water", p).sum()) for p in patches}
        total_flux = {p: float(scalar_values(directory / "phi", p).sum()) for p in patches}
        history.append({
            "time_s": t, "inventory_nl": float(alpha @ volumes) / 1e-12,
            "alpha_min": float(alpha.min()), "alpha_max": float(alpha.max()),
            "net_phase_flux_nl_s": sum(phase_flux.values()) / 1e-12,
            "dispersed_inlet_flux_nl_s": phase_flux["dispersedInlet"] / 1e-12,
            "outlet_phase_flux_nl_s": phase_flux["outlet"] / 1e-12,
            "net_total_flux_nl_s": sum(total_flux.values()) / 1e-12,
        })
    ts = np.array([r["time_s"] for r in history])
    flux = np.array([r["net_phase_flux_nl_s"] for r in history])
    integrated = float(np.trapz(flux, ts))
    inventory_change = history[-1]["inventory_nl"] - history[0]["inventory_nl"]
    supplied = float(metadata["flow_rates_m3_s"]["dispersed"]) * (ts[-1] - ts[0]) / 1e-12
    history_path = case / "conservation_snapshot_history.csv"
    with history_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(history[0]))
        writer.writeheader()
        writer.writerows(history)
    dt_values = []
    for log in case.glob("log.foamRun*"):
        dt = None
        opener = gzip.open if log.suffix == ".gz" else open
        with opener(log, "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("deltaT = "):
                    dt = float(line.split()[2])
                elif dt is not None and line.startswith("Time = "):
                    dt_values.append(dt)
    with (case / "droplet_measurements.csv").open(newline="") as handle:
        measurements = list(csv.DictReader(handle))
    primary = sorted([r for r in measurements if float(r["alpha_threshold"]) == .5], key=lambda r: float(r["crossing_time_s"]))
    blocks = []
    analysis_start = float(metadata.get("analysis_start_time_s", .3))
    span = (target - analysis_start) / 3
    for i in range(3):
        lo, hi = analysis_start + i * span, analysis_start + (i+1) * span + 1e-9
        selected = [r for r in primary if lo <= float(r["crossing_time_s"]) < hi]
        blocks.append({"start_s": lo, "end_s": hi, **statistics_for(selected)})
    report = {
        "case": str(case), "completed_output": bool(abs(ts[-1] - target) < 1e-8),
        "target_time_s": target, "last_time_s": float(ts[-1]), "minimum_cell_edge_um": minimum * 1e6,
        "conservative_inviscid_capillary_bound_s": capillary,
        "capillary_bound_is_diagnostic": "Viscosity may relax this bound; time-step convergence and a static-drop benchmark are required before attributing the breakup error to it.",
        "time_steps_max_s": max(dt_values) if dt_values else None,
        "time_steps_median_s": float(np.median(dt_values)) if dt_values else None,
        "time_steps_above_capillary_bound_fraction": float(np.mean(np.array(dt_values) > capillary)) if dt_values else None,
        "snapshot_conservation": {
            "inventory_change_nl": inventory_change, "integrated_net_outward_flux_nl": integrated,
            "residual_nl": inventory_change + integrated,
            "residual_percent_of_injected_phase": 100 * (inventory_change + integrated) / supplied,
            "max_abs_total_flux_nl_s": max(abs(r["net_total_flux_nl_s"]) for r in history),
            "inlet_phase_flux_mean_nl_s": statistics.mean(r["dispersed_inlet_flux_nl_s"] for r in history),
            "sampling_limitation": "Trapezoidal quadrature of saved instantaneous fluxes; use every-step inventory and boundary flux in subsequent cases.",
        },
        "solver_step_conservation": step_conservation(case, float(metadata["flow_rates_m3_s"]["dispersed"])),
        "all_crossings": statistics_for(primary), "last_time_blocks": blocks,
        "post_washout_statistics": {m: statistics_for([r for r in primary if float(r["crossing_time_s"]) >= analysis_start], m) for m in ("axial_diameter_um", "area_equivalent_diameter_um", "maximum_feret_um")},
        "diameter_metrics_after_first_two_crossings": {m: statistics_for(primary[2:], m) for m in ("axial_diameter_um", "area_equivalent_diameter_um", "maximum_feret_um")},
        "threshold_crossings": {str(a): sum(float(r["alpha_threshold"]) == a for r in measurements) for a in (.25, .5, .75)},
        "validation_pass": False,
        "required_checks": ["stationary droplet pressure and parasitic currents", "temporal convergence", "spatial convergence", "steady train statistics", "experimental diameter-method and geometry correspondence"],
    }
    (case / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    plot(case, primary, history, metadata["case"])
    return report


def plot(case: Path, rows: list[dict[str, str]], history: list[dict], geometry_name: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), constrained_layout=True)
    ts = [float(r["crossing_time_s"]) for r in rows]
    for key, label in (("axial_diameter_um", "Axial extent"), ("area_equivalent_diameter_um", "Projected-area equivalent")):
        axes[0].plot(ts, [float(r[key]) for r in rows], "o-", ms=3, label=label)
    paper_diameter, paper_frequency = (88, 14) if geometry_name == "unstretched" else (102, 8)
    axes[0].axhline(paper_diameter, color="black", ls="--", label="Paper diameter (method correspondence pending)")
    axes[0].set_ylabel("Diameter (µm)")
    axes[0].legend(fontsize=8)
    period_rows = [r for r in rows if r["preceding_period_s"]]
    axes[1].plot([float(r["crossing_time_s"]) for r in period_rows], [1 / float(r["preceding_period_s"]) for r in period_rows], "o-", ms=3)
    axes[1].axhline(paper_frequency, color="black", ls="--", label="Paper frequency")
    axes[1].set_ylabel("Crossing frequency (Hz)")
    axes[1].legend(fontsize=8)
    axes[2].plot([r["time_s"] for r in history], [r["inventory_nl"] for r in history])
    axes[2].set_ylabel("Dispersed inventory (nL)")
    for ax in axes:
        ax.set_xlabel("Simulated time (s)")
        ax.grid(True, alpha=.35)
    fig.suptitle("Raw CFD diagnostic — not validated or calibrated")
    fig.savefig(case / "audit.png", dpi=180)
    fig.savefig(case / "audit.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.case.resolve()), indent=2))
