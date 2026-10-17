#!/usr/bin/env python3
"""Compare raw flow-focusing CFD droplet trains to Roshan et al. Fig. 5.

No calibration factors are applied. Optical diameter is provisionally
represented by top-view axial extent, as in the paper's plug-length model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
Q_DISPERSED_M3_S = 30e-9 / 3600
NL = 1e-12
EXPERIMENT = {
    "unstretched": {"diameter_um": 88.0, "diameter_sd_um": 0.33, "frequency_hz": 14.0},
    "stretch5mm": {"diameter_um": 102.0, "diameter_sd_um": 0.92, "frequency_hz": 8.0},
}


def spherical_equivalent_from_rate_and_frequency_um(frequency_hz: float) -> float:
    """Mass-balance diagnostic; not a substitute for projected diameter."""
    volume = Q_DISPERSED_M3_S / frequency_hz
    return (6 * volume / math.pi) ** (1 / 3) * 1e6


def output_time_coverage(case: Path) -> dict[str, object]:
    control = (case / "system/controlDict").read_text()
    match = re.search(r"\bendTime\s+([0-9.eE+-]+)\s*;", control)
    target = float(match.group(1)) if match else None
    times = []
    for path in case.iterdir():
        if path.is_dir():
            try:
                times.append(float(path.name))
            except ValueError:
                pass
    last = max(times) if times else None
    return {
        "last_output_time_s": last,
        "requested_end_time_s": target,
        "reached_requested_end_time": last >= target - 1e-8 if last is not None and target is not None else None,
    }


def evaluate(name: str, mesh: str, tag: str = "robust") -> dict[str, object]:
    case_name = name + (f"_{tag}" if tag else "")
    case = HERE / f"runs-{mesh}" / case_name
    path = case / "droplet_measurements.csv"
    expected = EXPERIMENT[name]
    if not (case / "case.json").exists():
        return {"case": name, "case_directory": case_name, "mesh": mesh, "status": "case_not_generated"}
    metadata = json.loads((case / "case.json").read_text())
    throat_width = float(metadata["geometry_um"]["throat_width_um"])
    cross_spacing = float(metadata.get("cross_stream_cell_um", metadata["nominal_cell_um"]))
    # blockMesh rounds the two half-throat blocks independently.
    throat_cells = 2 * max(1, round((throat_width / 2) / cross_spacing))
    if not path.exists():
        return {"case": name, "case_directory": case_name, "mesh": mesh, "status": "measurement_missing"}
    with path.open(newline="") as handle:
        rows = sorted(
            (row for row in csv.DictReader(handle) if abs(float(row["alpha_threshold"]) - 0.5) < 1e-12),
            key=lambda row: float(row["crossing_time_s"]),
        )
    summary_path = path.with_suffix(".json")
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    counts = summary.get("crossing_counts", {}).get(case_name, {})
    coverage = summary.get("time_coverage", {}).get(case_name) or output_time_coverage(case)
    threshold_robust = all(str(key) in counts for key in (0.25, 0.5, 0.75)) and len({
        counts[str(key)] for key in (0.25, 0.5, 0.75)
    }) == 1
    if len(rows) < 4:
        return {"case": name, "case_directory": case_name, "mesh": mesh, "status": "too_few_droplets", "crossings": len(rows)}
    # Exclude the first two startup crossings; requiring >=2 later values.
    steady = rows[2:]
    axial = [float(row["axial_diameter_um"]) for row in steady]
    periods = [float(b["crossing_time_s"]) - float(a["crossing_time_s"]) for a, b in zip(steady, steady[1:])]
    diameter = statistics.mean(axial)
    frequency = 1 / statistics.mean(periods) if periods else None
    diameter_sd = statistics.stdev(axial) if len(axial) > 1 else 0
    core_rows = [
        row for row in steady
        if row.get("vof_phase_volume_nl") and row.get("volume_pair_quality") == "interpolated"
    ]
    volume_rows = [
        row for row in steady
        if row.get("full_vof_slab_volume_nl") and row.get("slab_separated") == "True"
    ]
    volumes = [float(row["full_vof_slab_volume_nl"]) for row in volume_rows]
    mean_volume = statistics.mean(volumes) if len(volumes) >= 5 else None
    expected_volume = Q_DISPERSED_M3_S / frequency / NL if frequency else None
    dropwise_ratios = [
        float(row["full_vof_slab_volume_nl"]) * NL / (Q_DISPERSED_M3_S * float(row["preceding_period_s"]))
        for row in volume_rows if row is not steady[0] and row.get("preceding_period_s")
        and float(row["preceding_period_s"]) > 0
    ]
    peaks = [float(row["peak_alpha"]) for row in steady if row.get("peak_alpha")]
    min_peak_alpha = min(peaks) if peaks else None
    outlet_core_present = len(peaks) == len(steady) and min_peak_alpha is not None and min_peak_alpha >= 0.9
    balance_ready = (
        len(volume_rows) == len(steady) and len(steady) >= 50
        and threshold_robust and outlet_core_present
        and throat_cells >= 8 and coverage.get("reached_requested_end_time") is True
    )
    frequency_detection_robust = (
        threshold_robust and len(volume_rows) == len(steady)
        and outlet_core_present
        and len(steady) >= 50 and throat_cells >= 8
        and coverage.get("reached_requested_end_time") is True
    )
    if coverage.get("reached_requested_end_time") is False:
        status = "incomplete_run"
    elif throat_cells < 8:
        status = "underresolved_throat"
    elif not threshold_robust:
        status = "threshold_sensitive_detection"
    elif not outlet_core_present:
        status = "diffuse_outlet_interface"
    elif len(volume_rows) != len(steady):
        status = "unpaired_3d_volumes"
    elif len(steady) < 50:
        status = "insufficient_droplets"
    else:
        status = "needs_mesh_and_time_convergence"
    return {
        "case": name,
        "case_directory": case_name,
        "mesh": mesh,
        "status": status,
        "validation_pass": False,
        "throat_cells_across": throat_cells,
        "crossings": len(rows),
        "steady_droplets": len(steady),
        "cfd_axial_diameter_um": diameter,
        "cfd_diameter_sd_um": diameter_sd,
        "cfd_diameter_cv_percent": 100 * diameter_sd / diameter,
        "experiment_diameter_um": expected["diameter_um"],
        "experiment_diameter_sd_um": expected["diameter_sd_um"],
        "diameter_error_percent": 100 * (diameter / expected["diameter_um"] - 1),
        "cfd_frequency_hz": frequency,
        "cfd_frequency_detection_robust": frequency_detection_robust,
        "experiment_frequency_hz": expected["frequency_hz"],
        "frequency_error_percent": 100 * (frequency / expected["frequency_hz"] - 1) if frequency else None,
        "cfd_paired_volume_count": len(volume_rows),
        "cfd_mean_full_vof_slab_volume_nl": mean_volume,
        "cfd_mean_threshold_core_phase_volume_nl": (
            statistics.mean(float(row["vof_phase_volume_nl"]) for row in core_rows) if core_rows else None
        ),
        "cfd_q_over_f_volume_nl": expected_volume,
        "cfd_volume_balance_ratio": mean_volume / expected_volume if balance_ready and mean_volume and expected_volume else None,
        "cfd_mean_dropwise_volume_balance_ratio": statistics.mean(dropwise_ratios) if balance_ready and dropwise_ratios else None,
        "cfd_min_droplet_peak_alpha": min_peak_alpha,
        "cfd_peak_alpha_sample_count": len(peaks),
        "outlet_has_pure_phase_core": outlet_core_present,
        "cfd_crossings_alpha025": counts.get("0.25"),
        "cfd_crossings_alpha05": counts.get("0.5"),
        "cfd_crossings_alpha075": counts.get("0.75"),
        "threshold_robust": threshold_robust,
        "cfd_last_output_time_s": coverage.get("last_output_time_s"),
        "cfd_requested_end_time_s": coverage.get("requested_end_time_s"),
        "cfd_reached_requested_end_time": coverage.get("reached_requested_end_time"),
        "experiment_q_over_f_volume_nl": Q_DISPERSED_M3_S / expected["frequency_hz"] / NL,
        "cfd_q_over_f_spherical_equivalent_um": (
            spherical_equivalent_from_rate_and_frequency_um(frequency) if frequency else None
        ),
        "experiment_q_over_f_spherical_equivalent_um": spherical_equivalent_from_rate_and_frequency_um(
            expected["frequency_hz"]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", choices=("screening", "throat_refined", "focused", "medium", "fine"), default="screening")
    parser.add_argument("--tag", default="robust")
    args = parser.parse_args()
    rows = [evaluate(name, args.mesh, args.tag) for name in EXPERIMENT]
    fields = list(dict.fromkeys(key for row in rows for key in row))
    output = HERE / f"comparison-{args.mesh}{'-' + args.tag if args.tag else ''}.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)
    print(output)


if __name__ == "__main__":
    main()
