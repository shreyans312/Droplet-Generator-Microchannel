#!/usr/bin/env python3
"""Compile retained mesh-convergence diagnostics from completed case snapshots."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
LEVELS = ("screening", "medium", "fine")
CASES = ("rigid", "flexible")
H_UM = {"screening": 25.0, "medium": 12.5, "fine": 8.333333}
NZ = {"screening": 4, "medium": 8, "fine": 12}
FLEXIBLE_MIN_DEPTH_UM = 86.1070356867
METRICS = (
    "axial_diameter_um",
    "transverse_diameter_um",
    "maximum_feret_um",
    "area_equivalent_diameter_um",
)


def values(pattern: str, text: str) -> list[float]:
    return [float(value.rstrip(".")) for value in re.findall(pattern, text, flags=re.MULTILINE)]


def parse_check_mesh(path: Path) -> dict[str, float | int | str]:
    text = path.read_text()
    def one(pattern: str) -> float:
        match = re.search(pattern, text)
        return float(match.group(1).rstrip(".")) if match else math.nan
    return {
        "cell_count": round(one(r"cells:\s+(\d+)")),
        "max_aspect_ratio": one(r"Max aspect ratio = ([+\-\deE.]+)"),
        "min_cell_volume_m3": one(r"Min volume = ([+\-\deE.]+)"),
        "max_cell_volume_m3": one(r"Max volume = ([+\-\deE.]+)"),
        "max_non_orthogonality_deg": one(r"Mesh non-orthogonality Max: ([+\-\deE.]+)"),
        "mean_non_orthogonality_deg": one(r"Mesh non-orthogonality Max: [+\-\deE.]+ average: ([+\-\deE.]+)"),
        "max_skewness": one(r"Max skewness = ([+\-\deE.]+)"),
        "mesh_check": "pass" if "Mesh OK." in text else "fail",
    }


def parse_solver(path: Path) -> dict[str, float | int | str]:
    text = path.read_text()
    co = values(r"Courant Number mean: [+\-\deE.]+ max: ([+\-\deE.]+)", text)
    alpha_co = values(r"Interface Courant Number mean: [+\-\deE.]+ max: ([+\-\deE.]+)", text)
    delta_t = values(r"deltaT = ([+\-\deE.]+)", text)
    alpha_min = values(r"Min\(alpha.water\) = ([+\-\deE.]+)", text)
    alpha_max = values(r"Max\(alpha.water\) = ([+\-\deE.]+)", text)
    cumulative = values(r"cumulative = ([+\-\deE.]+)", text)
    clock = values(r"ClockTime = ([+\-\deE.]+) s", text)
    times = values(r"^Time = ([+\-\deE.]+)s", text)
    warnings = len(re.findall(r"\b(?:WARNING|FOAM Warning)\b", text, flags=re.IGNORECASE))
    return {
        "completed_time_s": max(times) if times else math.nan,
        "time_steps": len(times),
        "min_delta_t_s": min(delta_t) if delta_t else math.nan,
        "max_delta_t_s": max(delta_t) if delta_t else math.nan,
        "max_courant": max(co) if co else math.nan,
        "max_interface_courant": max(alpha_co) if alpha_co else math.nan,
        "min_alpha": min(alpha_min) if alpha_min else math.nan,
        "max_alpha": max(alpha_max) if alpha_max else math.nan,
        "max_abs_cumulative_continuity": max(map(abs, cumulative)) if cumulative else math.nan,
        "clock_time_s": max(clock) if clock else math.nan,
        "warning_count": warnings,
        "solver_status": "complete" if "\nEnd\n" in text or text.rstrip().endswith("End") else "incomplete",
    }


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def mean(rows: list[dict[str, str]], field: str) -> float:
    data = [float(row[field]) for row in rows if row.get(field, "")]
    return statistics.mean(data) if data else math.nan


def apparent_order(coarse: float, medium: float, fine: float) -> float:
    numerator = coarse - medium
    denominator = medium - fine
    if numerator * denominator <= 0 or abs(denominator) < 1e-12:
        return math.nan
    observed = numerator / denominator
    candidates = [index / 100 for index in range(1, 1001)]
    def ratio(order: float) -> float:
        return (25.0**order - 12.5**order) / (12.5**order - 8.333333**order)
    best = min(candidates, key=lambda order: abs(ratio(order) - observed))
    return best if best not in (candidates[0], candidates[-1]) else math.nan


def convergence_label(coarse: float, medium: float, fine: float) -> str:
    if any(math.isnan(value) for value in (coarse, medium, fine)):
        return "insufficient_data"
    first = medium - coarse
    second = fine - medium
    if first * second < 0:
        return "non_monotonic"
    if abs(second) >= abs(first):
        return "monotonic_not_contracting"
    return "monotonic_contracting"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_svg(summary: list[dict[str, object]], output_dir: Path) -> None:
    width, height = 1040, 570
    panels = {"rigid": (80, 500), "flexible": (560, 980)}
    colours = {
        "axial_diameter_um": "#dc2626",
        "transverse_diameter_um": "#2563eb",
        "maximum_feret_um": "#7c3aed",
        "area_equivalent_diameter_um": "#d97706",
    }
    def x(value: float, left: int, right: int) -> float:
        return left + (25.0 - value) / (25.0 - 8.333333) * (right - left)
    def y(value: float) -> float:
        return 485 - (value - 80) / 180 * 360
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="520" y="28" font-size="21" text-anchor="middle">Droplet observables under global mesh refinement</text>',
    ]
    for case_name, (left, right) in panels.items():
        parts.append(f'<text x="{(left+right)/2}" y="57" font-size="18" text-anchor="middle">{case_name.capitalize()}</text>')
        for tick in range(80, 261, 20):
            parts.append(f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{right}" y2="{y(tick):.1f}"/>')
            parts.append(f'<text x="{left-8}" y="{y(tick)+5:.1f}" font-size="12" text-anchor="end">{tick}</text>')
        parts += [
            f'<line class="axis" x1="{left}" y1="485" x2="{right}" y2="485"/>',
            f'<line class="axis" x1="{left}" y1="125" x2="{left}" y2="485"/>',
            f'<text x="{(left+right)/2}" y="525" font-size="14" text-anchor="middle">Nominal cell size (um), coarse to fine</text>',
        ]
        for h in (25.0, 12.5, 8.333333):
            parts.append(f'<text x="{x(h,left,right):.1f}" y="505" font-size="12" text-anchor="middle">{h:.2f}</text>')
        case_rows = {str(row["mesh_level"]): row for row in summary if row["case"] == case_name}
        for metric in METRICS:
            points = []
            for level in LEVELS:
                value = float(case_rows[level][f"mean_{metric}"])
                if not math.isnan(value):
                    points.append((x(H_UM[level], left, right), y(value)))
            if len(points) > 1:
                parts.append('<polyline points="' + ' '.join(f'{px:.1f},{py:.1f}' for px, py in points) + f'" fill="none" stroke="{colours[metric]}" stroke-width="2"/>')
            for px, py in points:
                parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4" fill="{colours[metric]}"/>')
    legend = (("axial_diameter_um", "axial"), ("transverse_diameter_um", "transverse"), ("maximum_feret_um", "max Feret"), ("area_equivalent_diameter_um", "area-equivalent"))
    for index, (metric, label) in enumerate(legend):
        lx = 275 + index * 145
        parts += [f'<line x1="{lx}" y1="553" x2="{lx+25}" y2="553" stroke="{colours[metric]}" stroke-width="3"/>', f'<text x="{lx+31}" y="558" font-size="13">{label}</text>']
    parts += ['<text x="22" y="305" font-size="15" text-anchor="middle" transform="rotate(-90 22 305)">Diameter (um)</text>', '</svg>']
    (output_dir / "mesh_convergence.svg").write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    quality: list[dict[str, object]] = []
    measurements_out: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    means_by_case: dict[str, dict[str, dict[str, float]]] = {case: {} for case in CASES}
    counts_by_case: dict[str, dict[str, int]] = {case: {} for case in CASES}

    for level in LEVELS:
        combined_measurements = args.input / f"{level}_droplet_measurements.csv"
        measurements = (
            read_csv(combined_measurements)
            if combined_measurements.exists()
            else [
                row
                for case_name in CASES
                for row in read_csv(args.input / f"{level}_{case_name}_droplet_measurements.csv")
            ]
        )
        for case_name in CASES:
            check = parse_check_mesh(args.input / f"{level}_{case_name}_checkMesh.log")
            solver = parse_solver(args.input / f"{level}_{case_name}_foamRun.log")
            interface_path = args.input / f"{level}_{case_name}_interface.csv"
            interface = read_csv(interface_path) if interface_path.exists() else []
            depth_min = 100.0 if case_name == "rigid" else FLEXIBLE_MIN_DEPTH_UM
            quality.append({
                "case": case_name,
                "mesh_level": level,
                "nominal_xy_cell_um": H_UM[level],
                "min_vertical_cell_um": depth_min / NZ[level],
                "max_vertical_cell_um": 100.0 / NZ[level],
                "cells_across_junction_width": 100.0 / H_UM[level],
                "cells_across_undeformed_depth": NZ[level],
                **check,
            })
            case_measurements = [row for row in measurements if row["case"] == case_name]
            primary = [row for row in case_measurements if math.isclose(float(row["alpha_threshold"]), 0.5)]
            low = [row for row in case_measurements if math.isclose(float(row["alpha_threshold"]), 0.25)]
            high = [row for row in case_measurements if math.isclose(float(row["alpha_threshold"]), 0.75)]
            primary.sort(key=lambda row: float(row["crossing_time_s"]))
            periods = [float(b["crossing_time_s"]) - float(a["crossing_time_s"]) for a, b in zip(primary, primary[1:])]
            metric_means = {metric: mean(primary, metric) for metric in METRICS}
            threshold_sensitivity: dict[str, float] = {}
            for metric in METRICS:
                centre = metric_means[metric]
                alternatives = [mean(low, metric), mean(high, metric)]
                threshold_sensitivity[metric] = (
                    100 * max(abs(value - centre) for value in alternatives) / centre
                    if centre and all(not math.isnan(value) for value in alternatives)
                    else math.nan
                )
            means_by_case[case_name][level] = metric_means
            counts_by_case[case_name][level] = len(primary)
            for row in case_measurements:
                measurements_out.append({
                    "case": case_name,
                    "mesh_level": level,
                    "alpha_threshold": row["alpha_threshold"],
                    "track_id": row["track_id"],
                    "crossing_time_s": row["crossing_time_s"],
                    **{metric: row[metric] for metric in METRICS},
                    "cells_across_axial": float(row["axial_diameter_um"]) / H_UM[level],
                    "cells_across_transverse": float(row["transverse_diameter_um"]) / H_UM[level],
                    "cells_across_maximum_feret": float(row["maximum_feret_um"]) / H_UM[level],
                    "cells_across_area_equivalent": float(row["area_equivalent_diameter_um"]) / H_UM[level],
                    "legacy_eight_cell_flag": row["resolution_flag"],
                })
            summary.append({
                "case": case_name,
                "mesh_level": level,
                "cell_count": check["cell_count"],
                "nominal_cell_um": H_UM[level],
                "droplets": len(primary),
                "under_resolved_by_eight_cell_rule": sum(row["resolution_flag"] != "resolved" for row in primary),
                "formation_period_estimate_s": statistics.mean(periods) if periods else math.nan,
                **{f"mean_{metric}": value for metric, value in metric_means.items()},
                "mean_cells_across_axial": metric_means["axial_diameter_um"] / H_UM[level],
                "mean_cells_across_transverse": metric_means["transverse_diameter_um"] / H_UM[level],
                "axial_threshold_sensitivity_percent": threshold_sensitivity["axial_diameter_um"],
                "area_equivalent_threshold_sensitivity_percent": threshold_sensitivity["area_equivalent_diameter_um"],
                "mean_effective_interface_band_um": mean(interface, "effective_band_um"),
                "mean_effective_interface_band_cells": mean(interface, "effective_band_cells"),
                **solver,
            })

    convergence: list[dict[str, object]] = []
    for case_name in CASES:
        for metric in METRICS:
            coarse, medium, fine = (means_by_case[case_name][level][metric] for level in LEVELS)
            trend = convergence_label(coarse, medium, fine)
            counts = [counts_by_case[case_name][level] for level in LEVELS]
            convergence.append({
                "case": case_name,
                "observable": metric,
                "screening_mean_um": coarse,
                "medium_mean_um": medium,
                "fine_mean_um": fine,
                "screening_to_medium_change_percent": 100 * (medium - coarse) / coarse,
                "medium_to_fine_change_percent": 100 * (fine - medium) / medium,
                "trend": trend,
                "equal_droplet_counts_across_meshes": len(set(counts)) == 1,
                "apparent_order_if_admissible": (
                    apparent_order(coarse, medium, fine)
                    if trend == "monotonic_contracting" and len(set(counts)) == 1
                    else math.nan
                ),
            })

    write_csv(args.output / "mesh_quality.csv", quality)
    write_csv(args.output / "mesh_measurements.csv", measurements_out)
    write_csv(args.output / "mesh_summary.csv", summary)
    write_csv(args.output / "mesh_convergence.csv", convergence)
    write_svg(summary, args.output)
    print(f"wrote {len(summary)} case/mesh summaries and {len(convergence)} convergence rows")


if __name__ == "__main__":
    main()
