#!/usr/bin/env python3
"""Create a rapid final-snapshot preview from completed sweep cases."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import vtk


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_ROOT = HERE / "runs-medium"
MEASURE_PATH = ROOT / "cases" / "flexible-rigid" / "measure.py"
REFERENCE = HERE / "reference_data" / "fig6_digitized.csv"
CALIBRATION = ROOT / "studies" / "raj2016-comparison" / "calibration.json"
WIDTH_UM = 100.0


def load_measure_module():
    spec = importlib.util.spec_from_file_location("quick_measure", MEASURE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {MEASURE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


MEASURE = load_measure_module()


def final_snapshot(case: Path) -> tuple[float, list[float]]:
    marker = case / f"{case.name}.foam"
    marker.touch()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(marker))
    reader.SkipZeroTimeOn()
    reader.CreateCellToPointOn()
    reader.UpdateInformation()
    reader.DisableAllCellArrays()
    reader.SetCellArrayStatus("alpha.water", 1)
    times = MEASURE._time_values(reader)
    nominal_cell = 12.5 * MEASURE.UM
    sample_times = [times[index] for index in (-1, -6, -11)]
    diameters: list[float] = []
    for time_value in sample_times:
        reader.UpdateTimeStep(time_value)
        mesh = MEASURE._internal_mesh(reader.GetOutput())
        mask = MEASURE.rasterize_silhouette(MEASURE._interface(mesh, 0.5))
        detections = MEASURE.component_metrics(mask, nominal_cell)
        nearest = min(detections, key=lambda item: abs(item["cx"] - MEASURE.MEASUREMENT_X))
        diameters.append(nearest["axial_diameter"] / MEASURE.UM)
    return sample_times[-1], diameters


def equation(kp: float, flow_ratio: float) -> float:
    viscosity_ratio = 0.000914 / 0.0287
    return 1.97064 * (1 + kp) ** -0.0655 * flow_ratio**0.13593 * viscosity_ratio**-0.0273


def write_svg(
    rows: list[dict[str, object]],
    references: list[dict[str, str]],
    output_name: str = "meeting_preview_completed_r10.svg",
    calibration_note: str = "fixed calibration factor; six completed simulations",
    calibrated_label: str = "Calibrated CFD",
    footer_note: str = "Rapid three-sample preview only; full 50-droplet tracking remains in progress.",
) -> None:
    width, height = 1000, 690
    left, right, top, bottom = 95, 45, 90, 125
    xmin, xmax = -0.7, 2.8
    ymin = min(
        1.75,
        min(float(row["calibrated_d_over_w"]) - float(row["calibrated_spatial_std"]) for row in rows) - 0.05,
        min(float(row["raw_d_over_w"]) - float(row["raw_spatial_std"]) for row in rows) - 0.05,
    )
    ymax = max(2.38, max(float(row["calibrated_d_over_w"]) + float(row["calibrated_spatial_std"]) for row in rows) + 0.05)
    plot_right, plot_bottom = width - right, height - bottom

    def x(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (plot_right - left)

    def y(value: float) -> float:
        return plot_bottom - (value - ymin) / (ymax - ymin) * (plot_bottom - top)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111827}.major{stroke:#cbd5e1}.minor{stroke:#e5e7eb}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="500" y="30" font-size="22" font-weight="600" text-anchor="middle">Completed CFD cases: immediate meeting preview</text>',
        f'<text x="500" y="56" font-size="14" text-anchor="middle">r = 1.0; {calibration_note}</text>',
        '<text x="500" y="77" font-size="12" text-anchor="middle">Paper bars: experimental SD; CFD bars: SD of three near-line preview samples</text>',
    ]
    for index in range(15):
        tick = xmin + index * 0.25
        klass = "major" if index % 2 == 0 else "minor"
        parts.append(f'<line class="{klass}" x1="{x(tick):.1f}" y1="{top}" x2="{x(tick):.1f}" y2="{plot_bottom}"/>')
        if index % 2 == 0:
            parts.append(f'<text x="{x(tick):.1f}" y="{plot_bottom+25}" font-size="12" text-anchor="middle">{tick:.1f}</text>')
    y_step = 0.1
    y_start = math.floor(ymin / y_step) * y_step
    y_ticks = []
    tick = y_start
    while tick <= ymax + 1e-9:
        y_ticks.append(tick)
        tick += y_step
    for index, tick in enumerate(y_ticks):
        klass = "major"
        parts.append(f'<line class="{klass}" x1="{left}" y1="{y(tick):.1f}" x2="{plot_right}" y2="{y(tick):.1f}"/>')
        parts.append(f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="12" text-anchor="end">{tick:.2f}</text>')
    parts += [
        f'<rect x="{left}" y="{top}" width="{plot_right-left}" height="{plot_bottom-top}" fill="none" class="axis"/>',
        f'<text x="{(left+plot_right)/2:.1f}" y="{plot_bottom+58}" font-size="17" text-anchor="middle">Controlling parameter, Kp</text>',
        f'<text x="25" y="{(top+plot_bottom)/2:.1f}" font-size="17" text-anchor="middle" transform="rotate(-90 25 {(top+plot_bottom)/2:.1f})">Normalized droplet diameter, Dd/w</text>',
    ]
    curve = " ".join(
        f'{x(xmin+(xmax-xmin)*i/200):.1f},{y(equation(xmin+(xmax-xmin)*i/200, 1.0)):.1f}'
        for i in range(201)
    )
    parts.append(f'<polyline points="{curve}" fill="none" stroke="#059669" stroke-width="3"/>')
    for ref in references:
        kp, value, error = float(ref["kp"]), float(ref["d_over_w"]), float(ref["uncertainty"])
        px, top_y, bottom_y = x(kp), y(value + error), y(value - error)
        parts += [
            f'<line x1="{px:.1f}" y1="{top_y:.1f}" x2="{px:.1f}" y2="{bottom_y:.1f}" stroke="#111827" stroke-width="1.8"/>',
            f'<line x1="{px-6:.1f}" y1="{top_y:.1f}" x2="{px+6:.1f}" y2="{top_y:.1f}" stroke="#111827" stroke-width="1.8"/>',
            f'<line x1="{px-6:.1f}" y1="{bottom_y:.1f}" x2="{px+6:.1f}" y2="{bottom_y:.1f}" stroke="#111827" stroke-width="1.8"/>',
            f'<circle cx="{px:.1f}" cy="{y(value):.1f}" r="5" fill="white" stroke="#111827" stroke-width="2"/>',
        ]
    raw_line = " ".join(f'{x(float(row["kp"])):.1f},{y(float(row["raw_d_over_w"])):.1f}' for row in rows)
    parts.append(f'<polyline points="{raw_line}" fill="none" stroke="#2563eb" stroke-width="2.5" stroke-dasharray="3 4"/>')
    for row in rows:
        kp, value, error = float(row["kp"]), float(row["raw_d_over_w"]), float(row["raw_spatial_std"])
        px, py, top_y, bottom_y = x(kp), y(value), y(value + error), y(value - error)
        parts += [
            f'<line x1="{px:.1f}" y1="{top_y:.1f}" x2="{px:.1f}" y2="{bottom_y:.1f}" stroke="#2563eb" stroke-width="1.8"/>',
            f'<rect x="{px-5:.1f}" y="{py-5:.1f}" width="10" height="10" fill="#2563eb" stroke="white" stroke-width="1"/>',
        ]
    line = " ".join(f'{x(float(row["kp"])):.1f},{y(float(row["calibrated_d_over_w"])):.1f}' for row in rows)
    parts.append(f'<polyline points="{line}" fill="none" stroke="#dc2626" stroke-width="3" stroke-dasharray="8 5"/>')
    for row in rows:
        kp, value, error = float(row["kp"]), float(row["calibrated_d_over_w"]), float(row["calibrated_spatial_std"])
        px, py, top_y, bottom_y = x(kp), y(value), y(value + error), y(value - error)
        diamond = f'{px:.1f},{py-7:.1f} {px+7:.1f},{py:.1f} {px:.1f},{py+7:.1f} {px-7:.1f},{py:.1f}'
        parts += [
            f'<line x1="{px:.1f}" y1="{top_y:.1f}" x2="{px:.1f}" y2="{bottom_y:.1f}" stroke="#dc2626" stroke-width="2"/>',
            f'<line x1="{px-7:.1f}" y1="{top_y:.1f}" x2="{px+7:.1f}" y2="{top_y:.1f}" stroke="#dc2626" stroke-width="2"/>',
            f'<line x1="{px-7:.1f}" y1="{bottom_y:.1f}" x2="{px+7:.1f}" y2="{bottom_y:.1f}" stroke="#dc2626" stroke-width="2"/>',
            f'<polygon points="{diamond}" fill="#dc2626" stroke="white" stroke-width="1"/>',
        ]
    parts += [
        '<line x1="65" y1="635" x2="90" y2="635" stroke="#111827" stroke-width="2"/><circle cx="77.5" cy="635" r="5" fill="white" stroke="#111827" stroke-width="2"/><text x="98" y="640" font-size="12">Raj experiment</text>',
        '<line x1="275" y1="635" x2="300" y2="635" stroke="#059669" stroke-width="3"/><text x="308" y="640" font-size="12">Equation 13</text>',
        '<line x1="455" y1="635" x2="480" y2="635" stroke="#2563eb" stroke-width="2.5" stroke-dasharray="3 4"/><rect x="462.5" y="630" width="10" height="10" fill="#2563eb"/><text x="488" y="640" font-size="12">Raw CFD</text>',
        f'<line x1="635" y1="635" x2="660" y2="635" stroke="#dc2626" stroke-width="3" stroke-dasharray="8 5"/><text x="668" y="640" font-size="12">{calibrated_label}</text>',
        f'<text x="500" y="674" font-size="11.5" text-anchor="middle">{footer_note}</text>',
        '</svg>',
    ]
    (HERE / output_name).write_text("\n".join(parts) + "\n")


def main() -> None:
    with (HERE / "matrix.csv").open(newline="") as handle:
        matrix = [row for row in csv.DictReader(handle) if float(row["flow_ratio"]) == 1.0]
    complete = [row for row in matrix if (RUN_ROOT / row["case_id"] / "RUN_COMPLETE").exists()]
    factor = float(json.loads(CALIBRATION.read_text())["scale_factor"])
    measured: dict[str, tuple[float, list[float]]] = {}
    with ThreadPoolExecutor(max_workers=len(complete)) as pool:
        futures = {pool.submit(final_snapshot, RUN_ROOT / row["case_id"]): row["case_id"] for row in complete}
        for future in as_completed(futures):
            measured[futures[future]] = future.result()
    rows: list[dict[str, object]] = []
    for item in sorted(complete, key=lambda row: float(row["kp"])):
        time_value, diameters = measured[item["case_id"]]
        raw_mean = statistics.mean(diameters) / WIDTH_UM
        raw_std = statistics.stdev(diameters) / WIDTH_UM if len(diameters) > 1 else 0.0
        rows.append({
            "case_id": item["case_id"],
            "flow_ratio": 1.0,
            "kp": float(item["kp"]),
            "snapshot_time_s": time_value,
            "near_line_samples": len(diameters),
            "raw_d_over_w": raw_mean,
            "raw_spatial_std": raw_std,
            "calibration_factor": factor,
            "calibrated_d_over_w": factor * raw_mean,
            "calibrated_spatial_std": factor * raw_std,
        })
    with (HERE / "meeting_preview_completed_r10.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with REFERENCE.open(newline="") as handle:
        references = [row for row in csv.DictReader(handle) if float(row["flow_ratio"]) == 1.0]
    write_svg(rows, references)
    print(f"wrote meeting preview for {len(rows)} completed cases")


if __name__ == "__main__":
    main()
