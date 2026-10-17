#!/usr/bin/env python3
"""Write raw Figure 6 CSV, summary, and SVG outputs."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MEASURE = ROOT / "cases" / "flexible-rigid" / "measure.py"
CALIBRATION = ROOT / "studies" / "raj2016-comparison" / "calibration.json"


WIDTH = 100e-6


def measured_axial_diameters(case: Path) -> tuple[list[float], int]:
    output = case / "droplet_measurements.csv"
    subprocess.run(
        [
            "/usr/bin/python3",
            str(MEASURE),
            "--case",
            str(case),
            "--case-name",
            case.name,
            "--output",
            str(output),
            "--threshold",
            "0.5",
        ],
        check=True,
    )
    with output.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if abs(float(row["alpha_threshold"]) - 0.5) < 1e-12
        ]
    return [float(row["axial_diameter_um"]) * 1e-6 for row in rows], sum(
        row["resolution_flag"] != "resolved" for row in rows
    )


def reference(kp: float, ratio: float, viscosity_ratio: float) -> float:
    return 1.97064 * (1 + kp) ** -0.0655 * ratio**0.13593 * viscosity_ratio**-0.0273


def write_svg(rows: list[dict[str, str]], output_tag: str) -> None:
    width, height = 1100, 720
    left, right, top, bottom = 95, 35, 55, 85
    xmin, xmax, ymin, ymax = -0.75, 2.8, 0.9, 2.65
    colors = {0.5: "#111827", 1.0: "#dc2626", 1.5: "#2563eb"}

    def x(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (width - left - right)

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="550" y="30" font-size="22" text-anchor="middle">Figure 6b — raw image-matched axial CFD diameter</text>',
    ]
    for tick in (-0.5, 0, 0.5, 1, 1.5, 2, 2.5):
        parts += [
            f'<line class="grid" x1="{x(tick):.1f}" y1="{top}" x2="{x(tick):.1f}" y2="{height-bottom}"/>',
            f'<text x="{x(tick):.1f}" y="{height-bottom+28}" font-size="15" text-anchor="middle">{tick:g}</text>',
        ]
    for tick in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="15" text-anchor="end">{tick:.1f}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="{(left+width-right)/2}" y="{height-25}" font-size="18" text-anchor="middle">Kp</text>',
        f'<text x="25" y="{(top+height-bottom)/2}" font-size="18" text-anchor="middle" transform="rotate(-90 25 {(top+height-bottom)/2})">Dd / w</text>',
    ]
    for ratio in (0.5, 1.0, 1.5):
        curve = []
        for index in range(121):
            kp = xmin + (xmax - xmin) * index / 120
            curve.append((kp, reference(kp, ratio, 0.03185)))
        parts.append(
            f'<polyline points="{" ".join(f"{x(px):.1f},{y(py):.1f}" for px, py in curve)}" '
            f'fill="none" stroke="{colors[ratio]}" stroke-width="2" opacity="0.7"/>'
        )
        for row in (item for item in rows if float(item["flow_ratio"]) == ratio):
            px = float(row["kp"])
            paper = float(row["paper_d_over_w"])
            uncertainty = float(row["paper_std"])
            parts += [
                f'<line x1="{x(px):.1f}" y1="{y(paper-uncertainty):.1f}" x2="{x(px):.1f}" y2="{y(paper+uncertainty):.1f}" stroke="{colors[ratio]}" stroke-width="2"/>',
                f'<circle cx="{x(px):.1f}" cy="{y(paper):.1f}" r="5" fill="{colors[ratio]}"/>',
            ]
            if row["raw_cfd_d_over_w"]:
                raw = float(row["raw_cfd_d_over_w"])
                parts.append(
                    f'<path d="M {x(px)-7:.1f} {y(raw)-7:.1f} L {x(px)+7:.1f} {y(raw)+7:.1f} '
                    f'M {x(px)-7:.1f} {y(raw)+7:.1f} L {x(px)+7:.1f} {y(raw)-7:.1f}" '
                    f'stroke="{colors[ratio]}" stroke-width="3"/>'
                )
    for index, ratio in enumerate((0.5, 1.0, 1.5)):
        lx = 700 + index * 120
        parts += [
            f'<line x1="{lx}" y1="65" x2="{lx+24}" y2="65" stroke="{colors[ratio]}" stroke-width="3"/>',
            f'<text x="{lx+30}" y="70" font-size="14">r={ratio:g}</text>',
        ]
    parts += [
        '<circle cx="720" cy="95" r="5" fill="#111827"/>',
        '<text x="735" y="100" font-size="14">paper ± SD</text>',
        '<path d="M 895 88 L 909 102 M 895 102 L 909 88" stroke="#111827" stroke-width="3"/>',
        '<text x="918" y="100" font-size="14">raw CFD</text>',
        '</svg>',
    ]
    (HERE / "plots").mkdir(exist_ok=True)
    (HERE / "plots" / f"figure6_raw_{output_tag}.svg").write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", choices=("screening", "medium", "fine"), default="screening")
    parser.add_argument("--completed-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    with (HERE / "matrix.csv").open(newline="") as handle:
        matrix = list(csv.DictReader(handle))
    with (HERE / "reference_data" / "fig6_digitized.csv").open(newline="") as handle:
        references = list(csv.DictReader(handle))
    scale_factor = float(json.loads(CALIBRATION.read_text())["scale_factor"])

    measurements: dict[str, tuple[list[float], int]] = {}
    targets: list[tuple[str, Path]] = []
    for item in matrix:
        case = HERE / f"runs-{args.mesh}" / item["case_id"]
        if not case.exists():
            continue
        if args.completed_only and not (case / "RUN_COMPLETE").exists():
            continue
        targets.append((item["case_id"], case))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(measured_axial_diameters, case): case_id
            for case_id, case in targets
        }
        for future in as_completed(futures):
            case_id = futures[future]
            measurements[case_id] = future.result()
            print(f"{case_id}: measurement complete", flush=True)

    rows = []
    for item in matrix:
        kp = float(item["kp"])
        ratio = float(item["flow_ratio"])
        paper = min(
            (row for row in references if float(row["flow_ratio"]) == ratio),
            key=lambda row: abs(float(row["kp"]) - kp),
        )
        diameters, under_resolved = measurements.get(item["case_id"], ([], 0))
        mean = statistics.mean(diameters) / WIDTH if diameters else None
        std = statistics.stdev(diameters) / WIDTH if len(diameters) > 1 else None
        paper_value = float(paper["d_over_w"])
        equation_value = reference(kp, ratio, float(item["viscosity_ratio"]))
        calibrated_mean = mean * scale_factor if mean is not None else None
        calibrated_std = std * scale_factor if std is not None else None
        flags = []
        if len(diameters) < 50:
            flags.append("insufficient_droplets")
        if under_resolved:
            flags.append("under_resolved")
        quality = "validation_candidate" if not flags else ";".join(flags)
        rows.append({
            "case_id": item["case_id"],
            "kp": item["kp"],
            "flow_ratio": item["flow_ratio"],
            "paper_d_over_w": paper["d_over_w"],
            "paper_std": paper["uncertainty"],
            "equation13_d_over_w": f"{equation_value:.7g}",
            "raw_cfd_d_over_w": "" if mean is None else f"{mean:.7g}",
            "raw_cfd_std": "" if std is None else f"{std:.7g}",
            "calibration_factor": f"{scale_factor:.12g}",
            "calibrated_cfd_d_over_w": "" if calibrated_mean is None else f"{calibrated_mean:.7g}",
            "calibrated_cfd_std": "" if calibrated_std is None else f"{calibrated_std:.7g}",
            "tracked_droplets": str(len(diameters)),
            "under_resolved_droplets": str(under_resolved),
            "raw_error_percent": "" if mean is None else f"{100 * (mean-paper_value) / paper_value:.7g}",
            "calibrated_error_vs_paper_percent": "" if calibrated_mean is None else f"{100 * (calibrated_mean-paper_value) / paper_value:.7g}",
            "calibrated_error_vs_equation13_percent": "" if calibrated_mean is None else f"{100 * (calibrated_mean-equation_value) / equation_value:.7g}",
            "quality_flag": quality,
        })
    output_tag = f"{args.mesh}-partial" if args.completed_only else args.mesh
    with (HERE / f"results-{output_tag}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_svg(rows, output_tag)

    complete = [row for row in rows if row["raw_cfd_d_over_w"]]
    errors = [abs(float(row["raw_error_percent"])) for row in complete]
    calibrated_errors = [abs(float(row["calibrated_error_vs_paper_percent"])) for row in complete]
    summary = [
        "# Raw validation summary",
        "",
        f"- Cases measured: {len(complete)}/{len(rows)}",
        f"- Cases with at least 50 tracked droplets: {sum(int(row['tracked_droplets']) >= 50 for row in rows)}/{len(rows)}",
    ]
    if errors:
        summary += [
            f"- Raw mean absolute error: {statistics.mean(errors):.2f}%",
            f"- Raw maximum absolute error: {max(errors):.2f}%",
            f"- Fixed calibration factor: {scale_factor:.9f}",
            f"- Calibrated mean absolute error vs digitized data: {statistics.mean(calibrated_errors):.2f}%",
            f"- Calibrated maximum absolute error vs digitized data: {max(calibrated_errors):.2f}%",
        ]
    summary += [
        "",
        "The single calibration factor is imported unchanged from the rigid r=0.5, Kp=0 baseline; no pointwise refitting is performed.",
        "",
    ]
    (HERE / f"SUMMARY-{output_tag}.md").write_text("\n".join(summary))
    print(f"Measured {len(complete)}/{len(rows)} cases")


if __name__ == "__main__":
    main()
