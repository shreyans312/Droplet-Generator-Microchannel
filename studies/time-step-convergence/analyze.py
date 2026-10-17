#!/usr/bin/env python3
"""Compile time-step convergence statistics and publication-ready SVGs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIG6_DATA = ROOT / "studies" / "tjunction-sweep" / "reference_data" / "fig6_digitized.csv"
VARIANTS = ("coarse_40us", "baseline_20us", "fine_10us", "ultrafine_5us")
CONVERGENCE_VARIANTS = ("baseline_20us", "fine_10us", "ultrafine_5us")
CASES = ("rigid", "flexible")
CAP_US = {"coarse_40us": 40.0, "baseline_20us": 20.0, "fine_10us": 10.0, "ultrafine_5us": 5.0}
METRICS = (
    "axial_diameter_um",
    "transverse_diameter_um",
    "maximum_feret_um",
    "area_equivalent_diameter_um",
)
PRIMARY_THRESHOLD = 0.5
SAMPLE_SIZE = 50
CHANNEL_WIDTH_UM = 100.0
FLOW_RATIO = 0.5
VISCOSITY_RATIO = 0.000914 / 0.0287


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def percentile(data: list[float], fraction: float) -> float:
    ordered = sorted(data)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def parse_solver(path: Path, cap_s: float) -> dict[str, float | int | str]:
    patterns = {
        "dt": re.compile(r"^deltaT = ([+\-\deE.]+)"),
        "time": re.compile(r"^Time = ([+\-\deE.]+)s"),
        "co": re.compile(r"^Courant Number mean: [+\-\deE.]+ max: ([+\-\deE.]+)"),
        "aco": re.compile(r"^Interface Courant Number mean: [+\-\deE.]+ max: ([+\-\deE.]+)"),
        "amin": re.compile(r"Min\(alpha.water\) = ([+\-\deE.]+)"),
        "amax": re.compile(r"Max\(alpha.water\) = ([+\-\deE.]+)"),
        "continuity": re.compile(r"cumulative = ([+\-\deE.]+)"),
        "clock": re.compile(r"ClockTime = ([+\-\deE.]+) s"),
    }
    values: dict[str, list[float]] = {name: [] for name in patterns}
    warning_count = 0
    ended = False
    with path.open(errors="replace") as handle:
        for line in handle:
            warning_count += int("warning" in line.lower())
            ended |= line.strip() == "End"
            for name, pattern in patterns.items():
                match = pattern.search(line)
                if match:
                    values[name].append(float(match.group(1).rstrip(".")))
    dt = values["dt"]
    times = values["time"]
    at_cap = sum(value >= 0.99 * cap_s for value in dt)
    return {
        "solver_status": "complete" if ended else "incomplete",
        "completed_time_s": max(times, default=math.nan),
        "time_steps": len(times),
        "min_delta_t_us": min(dt, default=math.nan) * 1e6,
        "mean_delta_t_us": statistics.mean(dt) * 1e6 if dt else math.nan,
        "median_delta_t_us": statistics.median(dt) * 1e6 if dt else math.nan,
        "p05_delta_t_us": percentile(dt, 0.05) * 1e6,
        "p95_delta_t_us": percentile(dt, 0.95) * 1e6,
        "max_delta_t_us": max(dt, default=math.nan) * 1e6,
        "fraction_steps_at_cap": at_cap / len(dt) if dt else math.nan,
        "max_courant": max(values["co"], default=math.nan),
        "max_interface_courant": max(values["aco"], default=math.nan),
        "min_alpha": min(values["amin"], default=math.nan),
        "max_alpha": max(values["amax"], default=math.nan),
        "max_abs_cumulative_continuity": max(map(abs, values["continuity"]), default=math.nan),
        "clock_time_s": max(values["clock"], default=math.nan),
        "warning_count": warning_count,
    }


def lag1_correlation(data: list[float]) -> float:
    if len(data) < 3 or statistics.pstdev(data) == 0:
        return 0.0
    left, right = data[:-1], data[1:]
    mean_left, mean_right = statistics.mean(left), statistics.mean(right)
    numerator = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    denominator = math.sqrt(
        sum((a - mean_left) ** 2 for a in left) * sum((b - mean_right) ** 2 for b in right)
    )
    return numerator / denominator if denominator else 0.0


def t_critical_95(degrees_freedom: float) -> float:
    if degrees_freedom <= 1:
        return 12.706
    table = (
        (2, 4.303), (3, 3.182), (4, 2.776), (5, 2.571),
        (6, 2.447), (8, 2.306), (10, 2.228), (12, 2.179),
        (15, 2.131), (20, 2.086), (25, 2.060), (30, 2.042),
        (40, 2.021), (50, 2.009), (60, 2.000), (120, 1.980), (1e9, 1.960),
    )
    previous_df, previous_value = table[0]
    for current_df, current_value in table[1:]:
        if degrees_freedom <= current_df:
            fraction = (degrees_freedom - previous_df) / (current_df - previous_df)
            return previous_value + fraction * (current_value - previous_value)
        previous_df, previous_value = current_df, current_value
    return 1.960


def sample_stats(data: list[float]) -> dict[str, float | int]:
    count = len(data)
    if count == 0:
        return {"n": 0, "mean": math.nan, "std": math.nan, "lag1": math.nan, "effective_n": 0, "ci_low": math.nan, "ci_high": math.nan}
    mean = statistics.mean(data)
    if count == 1:
        return {"n": 1, "mean": mean, "std": math.nan, "lag1": math.nan, "effective_n": 1, "ci_low": math.nan, "ci_high": math.nan}
    std = statistics.stdev(data)
    rho = max(-0.95, min(0.95, lag1_correlation(data)))
    effective_n = max(2.0, min(float(count), count * (1 - rho) / (1 + rho)))
    half_width = t_critical_95(effective_n - 1) * std / math.sqrt(effective_n)
    return {
        "n": count,
        "mean": mean,
        "std": std,
        "lag1": rho,
        "effective_n": effective_n,
        "ci_low": mean - half_width,
        "ci_high": mean + half_width,
    }


def add_stats(row: dict[str, object], prefix: str, values: list[float]) -> None:
    for key, value in sample_stats(values).items():
        row[f"{prefix}_{key}"] = value


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"no rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def convergence_class(values: list[float]) -> str:
    if len(values) != 3 or any(math.isnan(value) for value in values):
        return "insufficient_data"
    first = values[1] - values[0]
    second = values[2] - values[1]
    if first * second < 0:
        return "non_monotonic"
    if abs(second) >= abs(first):
        return "monotonic_not_contracting"
    return "monotonic_contracting"


def paper_correlation(kp: float) -> float:
    """Raj et al. correlation, returned as Dd/w for this study's r and lambda."""
    return (
        1.97064
        * (1 + kp) ** -0.0655
        * FLOW_RATIO**0.13593
        * VISCOSITY_RATIO**-0.0273
    )


def write_paper_comparison(
    summary: list[dict[str, object]],
    output_dir: Path,
    variant: str = "ultrafine_5us",
    output_stem: str = "paper_comparison",
) -> None:
    references = [
        row for row in read_csv(FIG6_DATA)
        if math.isclose(float(row["flow_ratio"]), FLOW_RATIO)
    ]
    indexed = {(str(row["case"]), str(row["variant"])): row for row in summary}
    case_kp = {"rigid": 0.0, "flexible": 1.695}
    metric_fields = (
        ("axial", "axial"),
        ("transverse", "transverse"),
        ("maximum_feret", "max_feret"),
        ("area_equivalent", "area_equivalent"),
    )
    rows: list[dict[str, object]] = []
    for case_name in CASES:
        kp = case_kp[case_name]
        paper = min(references, key=lambda row: abs(float(row["kp"]) - kp))
        source = indexed[(case_name, variant)]
        for observable, prefix in metric_fields:
            cfd_mean = float(source[f"{prefix}_mean"])
            paper_mean = float(paper["d_over_w"]) * CHANNEL_WIDTH_UM
            rows.append({
                "case": case_name,
                "cfd_variant": variant,
                "cfd_kp": kp,
                "flow_ratio": FLOW_RATIO,
                "observable": observable,
                "cfd_mean_um": cfd_mean,
                "cfd_std_um": source[f"{prefix}_std"],
                "cfd_ci95_low_um": source[f"{prefix}_ci_low"],
                "cfd_ci95_high_um": source[f"{prefix}_ci_high"],
                "distinct_droplets": source["droplets_used"],
                "effective_sample_size": source[f"{prefix}_effective_n"],
                "actual_mean_delta_t_us": source["mean_delta_t_us"],
                "max_courant": source["max_courant"],
                "max_interface_courant": source["max_interface_courant"],
                "paper_nearest_kp": paper["kp"],
                "paper_digitized_mean_um": paper_mean,
                "paper_reported_sd_um": float(paper["uncertainty"]) * CHANNEL_WIDTH_UM,
                "paper_correlation_um_at_cfd_kp": paper_correlation(kp) * CHANNEL_WIDTH_UM,
                "raw_difference_from_digitized_percent": 100 * (cfd_mean - paper_mean) / paper_mean,
                "diagnostic_digitized_to_cfd_ratio_not_applied": paper_mean / cfd_mean,
                "diagnostic_correlation_to_cfd_ratio_not_applied": paper_correlation(kp) * CHANNEL_WIDTH_UM / cfd_mean,
                "paper_data_provenance": paper["source"],
                "comparison_status": "screening_only_spatially_unconverged_observable_not_unique",
            })
    write_csv(output_dir / f"{output_stem}.csv", rows)
    write_paper_svg(rows, references, variant, output_dir / f"{output_stem}.svg")
    calibration_stem = output_stem.replace("paper_comparison", "calibration_diagnostic")
    write_calibration_diagnostic(rows, output_dir / f"{calibration_stem}.csv", output_dir / f"{calibration_stem}.svg")


def write_calibration_diagnostic(
    comparison: list[dict[str, object]],
    csv_output: Path,
    svg_output: Path,
) -> None:
    fields = (
        "case", "cfd_variant", "observable", "cfd_mean_um",
        "paper_digitized_mean_um", "paper_correlation_um_at_cfd_kp",
        "diagnostic_digitized_to_cfd_ratio_not_applied",
        "diagnostic_correlation_to_cfd_ratio_not_applied",
        "comparison_status",
    )
    compact = [{field: row[field] for field in fields} for row in comparison]
    write_csv(csv_output, compact)

    width, height = 1080, 650
    left, right, top, bottom = 90, 35, 85, 135
    ymin, ymax = 1.4, 1.7
    labels = {
        "axial": "axial",
        "transverse": "transverse",
        "maximum_feret": "max Feret",
        "area_equivalent": "area equiv.",
    }

    def x(index: int) -> float:
        return left + (index + 0.5) / len(comparison) * (width - left - right)

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="540" y="30" font-size="22" text-anchor="middle">Post-hoc scale-factor diagnostic — factors are NOT applied</text>',
        '<text x="540" y="55" font-size="14" text-anchor="middle">A valid calibration would require one frozen, transferable model evaluated on independent holdout data</text>',
    ]
    for tick in (1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="13" text-anchor="end">{tick:.2f}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="25" y="{(top+height-bottom)/2}" font-size="16" text-anchor="middle" transform="rotate(-90 25 {(top+height-bottom)/2})">Reference diameter / raw CFD diameter</text>',
        f'<line x1="{x(3)+60:.1f}" y1="{top}" x2="{x(3)+60:.1f}" y2="{height-bottom}" stroke="#9ca3af" stroke-width="2"/>',
    ]
    for index, row in enumerate(comparison):
        digitized = float(row["diagnostic_digitized_to_cfd_ratio_not_applied"])
        correlation = float(row["diagnostic_correlation_to_cfd_ratio_not_applied"])
        px = x(index)
        parts += [
            f'<line x1="{px:.1f}" y1="{y(digitized):.1f}" x2="{px:.1f}" y2="{y(correlation):.1f}" stroke="#9ca3af"/>',
            f'<circle cx="{px:.1f}" cy="{y(digitized):.1f}" r="6" fill="#111827"/>',
            f'<rect x="{px-5:.1f}" y="{y(correlation)-5:.1f}" width="10" height="10" fill="#2563eb"/>',
            f'<text x="{px:.1f}" y="{height-bottom+25}" font-size="12" text-anchor="middle">{labels[str(row["observable"])]}</text>',
            f'<text x="{px:.1f}" y="{height-bottom+44}" font-size="12" text-anchor="middle">{str(row["case"])[0].upper()}</text>',
        ]
    parts += [
        '<circle cx="330" cy="590" r="6" fill="#111827"/>',
        '<text x="344" y="595" font-size="13">digitized Figure 6b / CFD</text>',
        '<rect x="565" y="584" width="11" height="11" fill="#2563eb"/>',
        '<text x="584" y="595" font-size="13">Equation 13 / CFD</text>',
        '<text x="540" y="628" font-size="13" text-anchor="middle">R = rigid; F = flexible. Variation across points disproves a single observable-independent factor.</text>',
        '</svg>',
    ]
    svg_output.write_text("\n".join(parts) + "\n")


def write_paper_svg(
    comparison: list[dict[str, object]],
    references: list[dict[str, str]],
    variant: str,
    output: Path,
) -> None:
    width, height = 1120, 720
    left, right, top, bottom = 100, 35, 70, 125
    xmin, xmax, ymin, ymax = -0.5, 2.7, 0.9, 2.2
    colours = {
        "axial": "#dc2626",
        "transverse": "#2563eb",
        "maximum_feret": "#7c3aed",
        "area_equivalent": "#d97706",
    }
    offsets = {"axial": -0.045, "transverse": -0.015, "maximum_feret": 0.015, "area_equivalent": 0.045}

    def x(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (width - left - right)

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="560" y="30" font-size="21" text-anchor="middle">Raw CFD observables vs digitized Raj et al. Figure 6b (r = 0.5)</text>',
        f'<text x="560" y="54" font-size="13" text-anchor="middle">{variant}: no fitted calibration; spatial convergence and observable equivalence are not established</text>',
    ]
    for tick in (-0.5, 0, 0.5, 1, 1.5, 2, 2.5):
        parts += [
            f'<line class="grid" x1="{x(tick):.1f}" y1="{top}" x2="{x(tick):.1f}" y2="{height-bottom}"/>',
            f'<text x="{x(tick):.1f}" y="{height-bottom+27}" font-size="13" text-anchor="middle">{tick:g}</text>',
        ]
    for tick in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="13" text-anchor="end">{tick:.1f}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="{(left+width-right)/2}" y="{height-82}" font-size="16" text-anchor="middle">Kp</text>',
        f'<text x="25" y="{(top+height-bottom)/2}" font-size="17" text-anchor="middle" transform="rotate(-90 25 {(top+height-bottom)/2})">Diameter / 100 um channel width</text>',
    ]
    curve = []
    for index in range(161):
        kp = xmin + (xmax - xmin) * index / 160
        curve.append((kp, paper_correlation(kp)))
    parts.append(
        '<polyline points="' + ' '.join(f'{x(px):.1f},{y(py):.1f}' for px, py in curve)
        + '" fill="none" stroke="#4b5563" stroke-width="2" stroke-dasharray="7 5"/>'
    )
    for row in references:
        kp = float(row["kp"])
        mean = float(row["d_over_w"])
        sd = float(row["uncertainty"])
        parts += [
            f'<line x1="{x(kp):.1f}" y1="{y(mean-sd):.1f}" x2="{x(kp):.1f}" y2="{y(mean+sd):.1f}" stroke="#111827" stroke-width="2"/>',
            f'<circle cx="{x(kp):.1f}" cy="{y(mean):.1f}" r="5" fill="#111827"/>',
        ]
    for row in comparison:
        observable = str(row["observable"])
        kp = float(row["cfd_kp"]) + offsets[observable]
        mean = float(row["cfd_mean_um"]) / CHANNEL_WIDTH_UM
        low = float(row["cfd_ci95_low_um"]) / CHANNEL_WIDTH_UM
        high = float(row["cfd_ci95_high_um"]) / CHANNEL_WIDTH_UM
        colour = colours[observable]
        parts += [
            f'<line x1="{x(kp):.1f}" y1="{y(low):.1f}" x2="{x(kp):.1f}" y2="{y(high):.1f}" stroke="{colour}" stroke-width="2"/>',
            f'<path d="M {x(kp)-6:.1f} {y(mean)-6:.1f} L {x(kp)+6:.1f} {y(mean)+6:.1f} M {x(kp)-6:.1f} {y(mean)+6:.1f} L {x(kp)+6:.1f} {y(mean)-6:.1f}" stroke="{colour}" stroke-width="3"/>',
        ]
    legend = (("axial", "CFD axial"), ("transverse", "CFD transverse"), ("maximum_feret", "CFD max Feret"), ("area_equivalent", "CFD area-equivalent"))
    for index, (observable, label) in enumerate(legend):
        lx = 320 + index * 185
        colour = colours[observable]
        parts += [
            f'<path d="M {lx} 653 L {lx+12} 665 M {lx} 665 L {lx+12} 653" stroke="{colour}" stroke-width="3"/>',
            f'<text x="{lx+19}" y="664" font-size="12">{label}</text>',
        ]
    parts += [
        '<circle cx="105" cy="659" r="5" fill="#111827"/>',
        '<text x="118" y="664" font-size="12">paper mean +/- SD</text>',
        '<line x1="105" y1="688" x2="133" y2="688" stroke="#4b5563" stroke-width="2" stroke-dasharray="7 5"/>',
        '<text x="141" y="693" font-size="12">paper correlation</text>',
        '</svg>',
    ]
    output.write_text("\n".join(parts) + "\n")


def write_svg(summary: list[dict[str, object]], output: Path) -> None:
    width, height = 1140, 760
    panel_left = {"rigid": 80, "flexible": 600}
    colours = {
        "axial": "#dc2626", "transverse": "#2563eb",
        "max_feret": "#7c3aed", "area_equivalent": "#d97706",
    }
    rows = {(str(row["case"]), str(row["variant"])): row for row in summary}

    def x(dt: float, left: int) -> float:
        return left + (40 - dt) / 35 * 430

    def y_diameter(value: float) -> float:
        return 330 - (value - 90) / 70 * 230

    def y_period(value: float) -> float:
        return 660 - (value - 0.035) / 0.035 * 220

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="570" y="28" font-size="21" text-anchor="middle">Time-step sensitivity on the fixed medium mesh</text>',
    ]
    for case_name in CASES:
        left, right = panel_left[case_name], panel_left[case_name] + 430
        parts.append(f'<text x="{(left+right)/2}" y="60" font-size="18" text-anchor="middle">{case_name.capitalize()}</text>')
        for tick in (90, 100, 110, 120, 130, 140, 150, 160):
            parts.append(f'<line class="grid" x1="{left}" y1="{y_diameter(tick):.1f}" x2="{right}" y2="{y_diameter(tick):.1f}"/>')
            parts.append(f'<text x="{left-8}" y="{y_diameter(tick)+4:.1f}" font-size="11" text-anchor="end">{tick}</text>')
        parts += [f'<line class="axis" x1="{left}" y1="330" x2="{right}" y2="330"/>', f'<line class="axis" x1="{left}" y1="100" x2="{left}" y2="330"/>']
        for metric, label in (("axial", "axial"), ("transverse", "transverse"), ("max_feret", "max Feret"), ("area_equivalent", "area-equivalent")):
            points = []
            for variant in VARIANTS:
                row = rows[(case_name, variant)]
                px = x(float(row["mean_delta_t_us"]), left)
                mean = float(row[f"{metric}_mean"])
                low = float(row[f"{metric}_ci_low"])
                high = float(row[f"{metric}_ci_high"])
                py = y_diameter(mean)
                points.append((px, py))
                parts += [
                    f'<line x1="{px:.1f}" y1="{y_diameter(low):.1f}" x2="{px:.1f}" y2="{y_diameter(high):.1f}" stroke="{colours[metric]}"/>',
                    f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4" fill="{colours[metric]}"/>',
                ]
            parts.append('<polyline points="' + ' '.join(f'{px:.1f},{py:.1f}' for px, py in points) + f'" fill="none" stroke="{colours[metric]}" stroke-width="2"/>')
        for tick in (0.035, 0.045, 0.055, 0.065, 0.070):
            parts.append(f'<line class="grid" x1="{left}" y1="{y_period(tick):.1f}" x2="{right}" y2="{y_period(tick):.1f}"/>')
            parts.append(f'<text x="{left-8}" y="{y_period(tick)+4:.1f}" font-size="11" text-anchor="end">{tick:.3f}</text>')
        parts += [f'<line class="axis" x1="{left}" y1="660" x2="{right}" y2="660"/>', f'<line class="axis" x1="{left}" y1="440" x2="{left}" y2="660"/>']
        period_points = []
        for variant in VARIANTS:
            row = rows[(case_name, variant)]
            px = x(float(row["mean_delta_t_us"]), left)
            mean = float(row["period_mean"])
            low = float(row["period_ci_low"])
            high = float(row["period_ci_high"])
            py = y_period(mean)
            period_points.append((px, py))
            parts += [f'<line x1="{px:.1f}" y1="{y_period(low):.1f}" x2="{px:.1f}" y2="{y_period(high):.1f}" stroke="#059669"/>', f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4" fill="#059669"/>']
        parts.append('<polyline points="' + ' '.join(f'{px:.1f},{py:.1f}' for px, py in period_points) + '" fill="none" stroke="#059669" stroke-width="2"/>')
        for cap in (40, 20, 10, 5):
            parts.append(f'<text x="{x(cap,left):.1f}" y="682" font-size="11" text-anchor="middle">{cap}</text>')
        parts.append(f'<text x="{(left+right)/2}" y="705" font-size="13" text-anchor="middle">Actual mean accepted time step (us)</text>')
    parts += [
        '<text x="20" y="215" font-size="14" text-anchor="middle" transform="rotate(-90 20 215)">Diameter (um)</text>',
        '<text x="20" y="550" font-size="14" text-anchor="middle" transform="rotate(-90 20 550)">Period (s)</text>',
    ]
    for index, (metric, label) in enumerate((("axial", "axial"), ("transverse", "transverse"), ("max_feret", "max Feret"), ("area_equivalent", "area-equivalent"))):
        left = 245 + index * 150
        parts += [f'<line x1="{left}" y1="735" x2="{left+24}" y2="735" stroke="{colours[metric]}" stroke-width="3"/>', f'<text x="{left+30}" y="740" font-size="12">{label}</text>']
    parts += ['<line x1="845" y1="735" x2="869" y2="735" stroke="#059669" stroke-width="3"/>', '<text x="875" y="740" font-size="12">period</text>', '</svg>']
    output.write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=HERE / "retained")
    parser.add_argument("--logs", type=Path, default=HERE / "runs")
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--paper-only", action="store_true")
    parser.add_argument("--paper-variant", choices=VARIANTS, default="ultrafine_5us")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.paper_only:
        existing = read_csv(args.output / "time_step_summary.csv")
        stem = f"paper_comparison_{args.paper_variant}_provisional"
        write_paper_comparison(existing, args.output, args.paper_variant, stem)
        return
    manifest = json.loads((args.input / "manifest.json").read_text())
    raw_out: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []

    for variant in VARIANTS:
        cap_s = float(manifest["variants"][variant])
        for case_name in CASES:
            measurements = read_csv(args.input / f"{variant}_{case_name}_measurements.csv")
            for item in measurements:
                raw_out.append({"variant": variant, "max_delta_t_us": CAP_US[variant], **item})
            primary = [row for row in measurements if math.isclose(float(row["alpha_threshold"]), PRIMARY_THRESHOLD)]
            primary.sort(key=lambda row: float(row["crossing_time_s"]))
            selected = primary[:SAMPLE_SIZE]
            row: dict[str, object] = {
                "case": case_name,
                "variant": variant,
                "configured_max_delta_t_us": CAP_US[variant],
                "total_primary_crossings": len(primary),
                "droplets_used": len(selected),
                "target_droplets_met": "yes" if len(selected) == SAMPLE_SIZE else "no",
                "under_resolved_count": sum(item["resolution_flag"] != "resolved" for item in selected),
            }
            row.update(parse_solver(args.logs / f"{variant}_{case_name}_foamRun.log", cap_s))
            for metric, prefix in (
                ("axial_diameter_um", "axial"),
                ("transverse_diameter_um", "transverse"),
                ("maximum_feret_um", "max_feret"),
                ("area_equivalent_diameter_um", "area_equivalent"),
            ):
                add_stats(row, prefix, [float(item[metric]) for item in selected])
            crossing_times = [float(item["crossing_time_s"]) for item in selected]
            periods = [later - earlier for earlier, later in zip(crossing_times, crossing_times[1:])]
            add_stats(row, "period", periods)
            period_low = float(row["period_ci_low"])
            period_high = float(row["period_ci_high"])
            row["frequency_hz"] = 1 / float(row["period_mean"]) if periods else math.nan
            row["frequency_ci_low_hz"] = 1 / period_high if period_high > 0 else math.nan
            row["frequency_ci_high_hz"] = 1 / period_low if period_low > 0 else math.nan
            if len(periods) >= 4:
                middle = len(periods) // 2
                first = statistics.mean(periods[:middle])
                second = statistics.mean(periods[middle:])
                row["period_first_half_s"] = first
                row["period_second_half_s"] = second
                row["period_drift_percent"] = 100 * (second - first) / first
            else:
                row["period_first_half_s"] = math.nan
                row["period_second_half_s"] = math.nan
                row["period_drift_percent"] = math.nan
            summary.append(row)

    write_csv(args.output / "time_step_measurements.csv", raw_out)
    write_csv(args.output / "time_step_summary.csv", summary)

    convergence: list[dict[str, object]] = []
    indexed = {(str(row["case"]), str(row["variant"])): row for row in summary}
    quantities = (
        ("axial", "um"), ("transverse", "um"), ("max_feret", "um"),
        ("area_equivalent", "um"), ("period", "s"), ("frequency", "hz"),
    )
    for case_name in CASES:
        for quantity, unit in quantities:
            field = "frequency_hz" if quantity == "frequency" else f"{quantity}_mean"
            values = [float(indexed[(case_name, variant)][field]) for variant in CONVERGENCE_VARIANTS]
            baseline_to_fine = 100 * (values[1] - values[0]) / values[0]
            fine_to_ultrafine = 100 * (values[2] - values[1]) / values[1]
            convergence.append({
                "case": case_name,
                "quantity": quantity,
                "unit": unit,
                "baseline_20us_value": values[0],
                "fine_10us_value": values[1],
                "ultrafine_5us_value": values[2],
                "baseline_to_fine_change_percent": baseline_to_fine,
                "fine_to_ultrafine_change_percent": fine_to_ultrafine,
                "trend": convergence_class(values),
            })
    write_csv(args.output / "time_step_convergence.csv", convergence)
    write_svg(summary, args.output / "time_step_convergence.svg")
    write_paper_comparison(summary, args.output)


if __name__ == "__main__":
    main()
