#!/usr/bin/env python3
"""Generate absolute and rigid-normalized Raj et al. Figure 6 comparisons."""

from __future__ import annotations

import csv
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUMMARY = ROOT / "studies" / "time-step-convergence" / "time_step_summary.csv"
REFERENCE = ROOT / "studies" / "tjunction-sweep" / "reference_data" / "fig6_digitized.csv"
SWEEP_RESULTS = ROOT / "studies" / "tjunction-sweep" / "results-medium.csv"
VARIANT = "fine_10us"
WIDTH_UM = 100.0
FLOW_RATIO = 0.5
VISCOSITY_RATIO = 0.000914 / 0.0287
CASE_KP = {"rigid": 0.0, "flexible": 1.695}
METRICS = (
    ("axial", "axial"),
    ("transverse", "transverse"),
    ("maximum_feret", "max_feret"),
    ("area_equivalent", "area_equivalent"),
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def correlation(kp: float, flow_ratio: float = FLOW_RATIO) -> float:
    return (
        1.97064
        * (1 + kp) ** -0.0655
        * flow_ratio**0.13593
        * VISCOSITY_RATIO**-0.0273
        * WIDTH_UM
    )


def build_rows() -> list[dict[str, object]]:
    summary = {
        row["case"]: row
        for row in read_csv(SUMMARY)
        if row["variant"] == VARIANT
    }
    references = [
        row for row in read_csv(REFERENCE)
        if math.isclose(float(row["flow_ratio"]), FLOW_RATIO)
    ]
    rows: list[dict[str, object]] = []
    for case_name in ("rigid", "flexible"):
        kp = CASE_KP[case_name]
        paper = min(references, key=lambda row: abs(float(row["kp"]) - kp))
        paper_mean = float(paper["d_over_w"]) * WIDTH_UM
        rows.extend((
            {
                "case": case_name,
                "series": "paper_digitized",
                "observable": "paper_droplet_size",
                "kp": kp,
                "diameter_um": paper_mean,
                "uncertainty_type": "experimental_standard_deviation",
                "uncertainty_um": float(paper["uncertainty"]) * WIDTH_UM,
                "source": "Raj et al. Figure 6b vector extraction",
            },
            {
                "case": case_name,
                "series": "paper_equation_13",
                "observable": "paper_droplet_size_correlation",
                "kp": kp,
                "diameter_um": correlation(kp),
                "uncertainty_type": "not_published",
                "uncertainty_um": "",
                "source": "Raj et al. Equation 13; R2=0.97",
            },
        ))
        current = summary[case_name]
        for observable, prefix in METRICS:
            rows.append({
                "case": case_name,
                "series": f"cfd_{observable}",
                "observable": observable,
                "kp": kp,
                "diameter_um": current[f"{prefix}_mean"],
                "uncertainty_type": "lag1_adjusted_95_percent_sampling_interval",
                "uncertainty_um": (
                    float(current[f"{prefix}_ci_high"])
                    - float(current[f"{prefix}_ci_low"])
                ) / 2,
                "source": "current uncalibrated interface-based CFD estimator",
            })
    return rows


def absolute_svg(rows: list[dict[str, object]], output: Path) -> None:
    width, height = 980, 650
    left, right, top, bottom = 90, 35, 70, 105
    ymin, ymax = 100.0, 210.0
    case_x = {"rigid": 320.0, "flexible": 700.0}
    series = (
        "paper_digitized", "paper_equation_13", "cfd_axial",
        "cfd_transverse", "cfd_maximum_feret", "cfd_area_equivalent",
    )
    offsets = dict(zip(series, (-75, -45, -15, 15, 45, 75)))
    colours = {
        "paper_digitized": "#111827", "paper_equation_13": "#059669",
        "cfd_axial": "#dc2626", "cfd_transverse": "#2563eb",
        "cfd_maximum_feret": "#7c3aed", "cfd_area_equivalent": "#d97706",
    }
    labels = {
        "paper_digitized": "paper mean +/- SD", "paper_equation_13": "paper Eq. 13",
        "cfd_axial": "CFD axial", "cfd_transverse": "CFD transverse",
        "cfd_maximum_feret": "CFD max Feret", "cfd_area_equivalent": "CFD area-equivalent",
    }

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    indexed = {(str(row["case"]), str(row["series"])): row for row in rows}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="490" y="28" font-size="21" text-anchor="middle">Raj et al. Figure 6b: absolute droplet-size comparison</text>',
        '<text x="490" y="52" font-size="13" text-anchor="middle">Raw CFD; no fitted diameter multiplier</text>',
    ]
    for tick in range(100, 211, 10):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-10}" y="{y(tick)+5:.1f}" font-size="12" text-anchor="end">{tick}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="24" y="{(top+height-bottom)/2}" font-size="16" text-anchor="middle" transform="rotate(-90 24 {(top+height-bottom)/2})">Diameter (um)</text>',
    ]
    for case_name, cx in case_x.items():
        for name in series:
            row = indexed[(case_name, name)]
            value = float(row["diameter_um"])
            px = cx + offsets[name]
            uncertainty = float(row["uncertainty_um"]) if row["uncertainty_um"] != "" else 0.0
            if uncertainty:
                parts.append(f'<line x1="{px}" y1="{y(value-uncertainty):.1f}" x2="{px}" y2="{y(value+uncertainty):.1f}" stroke="{colours[name]}" stroke-width="2"/>')
            parts.append(f'<circle cx="{px}" cy="{y(value):.1f}" r="6" fill="{colours[name]}"/>')
        parts.append(f'<text x="{cx}" y="{height-bottom+28}" font-size="16" text-anchor="middle">{case_name.capitalize()}</text>')
    for index, name in enumerate(series):
        lx = 105 + (index % 3) * 285
        ly = 590 + (index // 3) * 26
        parts += [
            f'<circle cx="{lx}" cy="{ly}" r="5" fill="{colours[name]}"/>',
            f'<text x="{lx+13}" y="{ly+5}" font-size="12">{labels[name]}</text>',
        ]
    parts += ['</svg>']
    output.write_text("\n".join(parts) + "\n")


def normalized_outputs(rows: list[dict[str, object]]) -> None:
    indexed = {(str(row["case"]), str(row["series"])): row for row in rows}
    series = (
        "paper_digitized", "paper_equation_13", "cfd_axial",
        "cfd_transverse", "cfd_maximum_feret", "cfd_area_equivalent",
    )
    normalized: list[dict[str, object]] = []
    for name in series:
        rigid = float(indexed[("rigid", name)]["diameter_um"])
        flexible = float(indexed[("flexible", name)]["diameter_um"])
        normalized.append({
            "series": name,
            "rigid_diameter_um": rigid,
            "flexible_diameter_um": flexible,
            "flexible_over_rigid": flexible / rigid,
            "size_reduction_percent": 100 * (1 - flexible / rigid),
            "normalization_status": "within_series_rigid_baseline; no fitted multiplier",
        })
    write_csv(HERE / "figure6_rigid_normalized.csv", normalized)

    width, height = 980, 620
    left, right, top, bottom = 85, 30, 75, 130
    ymin, ymax = 0.88, 1.00

    def x(index: int) -> float:
        return left + (index + 0.5) / len(normalized) * (width - left - right)

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    labels = ("Paper data", "Paper Eq. 13", "CFD axial", "CFD transverse", "CFD max Feret", "CFD area equiv.")
    colours = ("#111827", "#059669", "#dc2626", "#2563eb", "#7c3aed", "#d97706")
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="490" y="28" font-size="21" text-anchor="middle">Flexible response normalized by each series rigid baseline</text>',
        '<text x="490" y="52" font-size="13" text-anchor="middle">Tests relative size reduction without calibrating CFD magnitude</text>',
    ]
    for tick in (0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.00):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-10}" y="{y(tick)+5:.1f}" font-size="13" text-anchor="end">{tick:.2f}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="24" y="{(top+height-bottom)/2}" font-size="16" text-anchor="middle" transform="rotate(-90 24 {(top+height-bottom)/2})">Flexible / rigid diameter</text>',
    ]
    for index, (row, label, colour) in enumerate(zip(normalized, labels, colours)):
        value = float(row["flexible_over_rigid"])
        px = x(index)
        parts += [
            f'<rect x="{px-31:.1f}" y="{y(value):.1f}" width="62" height="{y(ymin)-y(value):.1f}" fill="{colour}" opacity="0.88"/>',
            f'<text x="{px:.1f}" y="{y(value)-9:.1f}" font-size="13" text-anchor="middle">{value:.3f}</text>',
            f'<text x="{px:.1f}" y="{height-bottom+25}" font-size="12" text-anchor="middle">{label}</text>',
            f'<text x="{px:.1f}" y="{height-bottom+45}" font-size="11" text-anchor="middle">{float(row["size_reduction_percent"]):.1f}% reduction</text>',
        ]
    parts += [
        '<text x="490" y="585" font-size="12" text-anchor="middle">Normalization removes absolute scale and therefore cannot establish absolute-diameter validation.</text>',
        '</svg>',
    ]
    (HERE / "figure6_rigid_normalized.svg").write_text("\n".join(parts) + "\n")


def calibrated_outputs(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Fit one rigid axial scale factor and evaluate flexible as holdout."""
    indexed = {(str(row["case"]), str(row["series"])): row for row in rows}
    summary = {
        row["case"]: row
        for row in read_csv(SUMMARY)
        if row["variant"] == VARIANT
    }
    rigid_raw = float(indexed[("rigid", "cfd_axial")]["diameter_um"])
    rigid_target = float(indexed[("rigid", "paper_equation_13")]["diameter_um"])
    factor = rigid_target / rigid_raw
    calibrated: list[dict[str, object]] = []
    for case_name in ("rigid", "flexible"):
        kp = CASE_KP[case_name]
        raw = float(indexed[(case_name, "cfd_axial")]["diameter_um"])
        equation = float(indexed[(case_name, "paper_equation_13")]["diameter_um"])
        digitized = float(indexed[(case_name, "paper_digitized")]["diameter_um"])
        value = factor * raw
        raw_standard_deviation = float(summary[case_name]["axial_std"])
        raw_ci_half_width = (
            float(summary[case_name]["axial_ci_high"])
            - float(summary[case_name]["axial_ci_low"])
        ) / 2
        calibrated.append({
            "case": case_name,
            "role": "calibration_fit" if case_name == "rigid" else "independent_kp_holdout",
            "flow_ratio": FLOW_RATIO,
            "kp": kp,
            "raw_cfd_axial_um": raw,
            "rigid_baseline_scale_factor": factor,
            "calibrated_cfd_axial_um": value,
            "calibrated_cfd_standard_deviation_um": factor * raw_standard_deviation,
            "calibrated_cfd_95_ci_half_width_um": factor * raw_ci_half_width,
            "paper_equation_13_um": equation,
            "paper_digitized_um": digitized,
            "error_vs_equation_13_percent": 100 * (value - equation) / equation,
            "error_vs_digitized_percent": 100 * (value - digitized) / digitized,
        })
    write_csv(HERE / "figure6_calibrated.csv", calibrated)

    rigid = calibrated[0]
    flexible = calibrated[1]
    cfd_reduction = 100 * (
        1 - float(flexible["calibrated_cfd_axial_um"]) / float(rigid["calibrated_cfd_axial_um"])
    )
    equation_reduction = 100 * (
        1 - float(flexible["paper_equation_13_um"]) / float(rigid["paper_equation_13_um"])
    )
    digitized_reduction = 100 * (
        1 - float(flexible["paper_digitized_um"]) / float(rigid["paper_digitized_um"])
    )
    metadata = {
        "scale_factor": factor,
        "fit": "Raj Equation 13 rigid point, Kp=0, divided by fine_10us rigid CFD axial mean",
        "holdout": "flexible point at Kp=1.695",
        "holdout_error_vs_equation_13_percent": flexible["error_vs_equation_13_percent"],
        "calibrated_cfd_size_reduction_percent": cfd_reduction,
        "equation_13_size_reduction_percent": equation_reduction,
        "digitized_size_reduction_percent": digitized_reduction,
        "trend_difference_vs_equation_13_percentage_points": cfd_reduction - equation_reduction,
    }
    import json

    (HERE / "calibration.json").write_text(json.dumps(metadata, indent=2) + "\n")

    reference_points = sorted(
        (
            row
            for row in read_csv(REFERENCE)
            if math.isclose(float(row["flow_ratio"]), FLOW_RATIO)
        ),
        key=lambda row: float(row["kp"]),
    )
    width, height = 900, 640
    left, right, top, bottom = 95, 40, 85, 95
    xmin, xmax, ymin, ymax = -0.4, 2.6, 175.0, 207.0

    def x(kp: float) -> float:
        return left + (kp - xmin) / (xmax - xmin) * (width - left - right)

    def y(value: float) -> float:
        return height - bottom - (value - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="450" y="30" font-size="21" text-anchor="middle">Raj et al. Figure 6b: calibrated CFD trend</text>',
        f'<text x="450" y="54" font-size="13" text-anchor="middle">One frozen rigid-baseline factor = {factor:.6f}; flexible Kp point is holdout</text>',
    ]
    for tick in range(175, 208, 5):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="{width-right}" y2="{y(tick):.1f}"/>',
            f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="13" text-anchor="end">{tick}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="{(left+width-right)/2}" y="590" font-size="16" text-anchor="middle">Kp</text>',
        f'<text x="25" y="{(top+height-bottom)/2}" font-size="16" text-anchor="middle" transform="rotate(-90 25 {(top+height-bottom)/2})">Droplet diameter (um)</text>',
    ]
    for tick in (-0.4, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
        parts += [
            f'<line class="grid" x1="{x(tick):.1f}" y1="{top}" x2="{x(tick):.1f}" y2="{height-bottom}"/>',
            f'<text x="{x(tick):.1f}" y="{height-bottom+26}" font-size="12" text-anchor="middle">{tick:g}</text>',
        ]

    digitized_points = " ".join(
        f'{x(float(row["kp"])):.1f},{y(float(row["d_over_w"]) * WIDTH_UM):.1f}'
        for row in reference_points
    )
    parts.append(f'<polyline points="{digitized_points}" fill="none" stroke="#111827" stroke-width="2" opacity="0.75"/>')
    for row in reference_points:
        kp = float(row["kp"])
        value = float(row["d_over_w"]) * WIDTH_UM
        uncertainty = float(row["uncertainty"]) * WIDTH_UM
        cap = 5.0
        parts += [
            f'<line x1="{x(kp):.1f}" y1="{y(value-uncertainty):.1f}" x2="{x(kp):.1f}" y2="{y(value+uncertainty):.1f}" stroke="#111827" stroke-width="1.5"/>',
            f'<line x1="{x(kp)-cap:.1f}" y1="{y(value-uncertainty):.1f}" x2="{x(kp)+cap:.1f}" y2="{y(value-uncertainty):.1f}" stroke="#111827" stroke-width="1.5"/>',
            f'<line x1="{x(kp)-cap:.1f}" y1="{y(value+uncertainty):.1f}" x2="{x(kp)+cap:.1f}" y2="{y(value+uncertainty):.1f}" stroke="#111827" stroke-width="1.5"/>',
            f'<circle cx="{x(kp):.1f}" cy="{y(value):.1f}" r="4.5" fill="#111827"/>',
        ]

    equation_points = " ".join(
        f'{x(xmin + (xmax-xmin)*index/160):.1f},{y(correlation(xmin + (xmax-xmin)*index/160)):.1f}'
        for index in range(161)
    )
    parts.append(f'<polyline points="{equation_points}" fill="none" stroke="#059669" stroke-width="3"/>')

    cfd_points = " ".join(
        f'{x(float(row["kp"])):.1f},{y(float(row["calibrated_cfd_axial_um"])):.1f}'
        for row in calibrated
    )
    parts.append(f'<polyline points="{cfd_points}" fill="none" stroke="#dc2626" stroke-width="3" stroke-dasharray="8 5"/>')
    for row in calibrated:
        px = x(float(row["kp"]))
        value = float(row["calibrated_cfd_axial_um"])
        uncertainty = float(row["calibrated_cfd_standard_deviation_um"])
        py = y(value)
        parts += [
            f'<line x1="{px:.1f}" y1="{y(value-uncertainty):.1f}" x2="{px:.1f}" y2="{y(value+uncertainty):.1f}" stroke="#dc2626" stroke-width="2.2"/>',
            f'<line x1="{px-7:.1f}" y1="{y(value-uncertainty):.1f}" x2="{px+7:.1f}" y2="{y(value-uncertainty):.1f}" stroke="#dc2626" stroke-width="2.2"/>',
            f'<line x1="{px-7:.1f}" y1="{y(value+uncertainty):.1f}" x2="{px+7:.1f}" y2="{y(value+uncertainty):.1f}" stroke="#dc2626" stroke-width="2.2"/>',
            f'<circle cx="{px:.1f}" cy="{py:.1f}" r="7" fill="#dc2626"/>',
            f'<text x="{px:.1f}" y="{py-13:.1f}" font-size="12" text-anchor="middle">CFD {value:.2f}</text>',
        ]

    legend = (
        ("Paper digitized +/- SD", "#111827"),
        ("Paper Equation 13", "#059669"),
        ("Calibrated CFD axial +/- droplet SD", "#dc2626"),
    )
    for index, (label, colour) in enumerate(legend):
        lx = 105 + index * 250
        parts += [
            f'<line x1="{lx}" y1="605" x2="{lx+26}" y2="605" stroke="{colour}" stroke-width="3"/>',
            f'<text x="{lx+34}" y="610" font-size="12">{label}</text>',
        ]
    parts += [
        f'<text x="450" y="75" font-size="12" text-anchor="middle">Flexible holdout error vs Equation 13: {float(flexible["error_vs_equation_13_percent"]):+.2f}%</text>',
        '</svg>',
    ]
    (HERE / "figure6_calibrated.svg").write_text("\n".join(parts) + "\n")
    return calibrated


def completed_sweep_rows() -> list[dict[str, object]]:
    """Load the 18-point sweep when analysis has produced complete results."""
    if not SWEEP_RESULTS.exists():
        return []
    rows = read_csv(SWEEP_RESULTS)
    if len(rows) != 18 or any(not row.get("calibrated_cfd_d_over_w") for row in rows):
        return []
    return [
        {
            "case": row["case_id"],
            "role": "fixed_factor_sweep",
            "flow_ratio": float(row["flow_ratio"]),
            "kp": float(row["kp"]),
            "calibrated_cfd_axial_um": float(row["calibrated_cfd_d_over_w"]) * WIDTH_UM,
            "calibrated_cfd_standard_deviation_um": (
                float(row["calibrated_cfd_std"]) * WIDTH_UM
                if row["calibrated_cfd_std"] else 0.0
            ),
            "calibrated_cfd_95_ci_half_width_um": "",
            "error_vs_digitized_percent": float(row["calibrated_error_vs_paper_percent"]),
            "error_vs_equation_13_percent": float(row["calibrated_error_vs_equation13_percent"]),
        }
        for row in rows
    ]


def publication_outputs(calibrated: list[dict[str, object]]) -> None:
    """Create a paper-style all-ratio plot and its complete numerical table."""
    references = sorted(
        read_csv(REFERENCE),
        key=lambda row: (float(row["flow_ratio"]), float(row["kp"])),
    )
    publication_rows: list[dict[str, object]] = []
    for row in references:
        flow_ratio = float(row["flow_ratio"])
        kp = float(row["kp"])
        paper_ratio = float(row["d_over_w"])
        paper_uncertainty = float(row["uncertainty"])
        equation_ratio = correlation(kp, flow_ratio) / WIDTH_UM
        common = {
            "flow_ratio": flow_ratio,
            "kp": kp,
            "mean_95_ci_half_width_d_over_w": "",
            "mean_95_ci_half_width_um": "",
        }
        publication_rows.extend((
            {
                **common,
                "series": "raj_experiment",
                "d_over_w": paper_ratio,
                "diameter_um": paper_ratio * WIDTH_UM,
                "uncertainty_type": "experimental_standard_deviation",
                "uncertainty_d_over_w": paper_uncertainty,
                "uncertainty_um": paper_uncertainty * WIDTH_UM,
                "percent_error_vs_digitized": "",
            },
            {
                **common,
                "series": "raj_equation_13",
                "d_over_w": equation_ratio,
                "diameter_um": equation_ratio * WIDTH_UM,
                "uncertainty_type": "not_published",
                "uncertainty_d_over_w": "",
                "uncertainty_um": "",
                "percent_error_vs_digitized": 100 * (equation_ratio - paper_ratio) / paper_ratio,
            },
        ))

    for row in calibrated:
        flow_ratio = float(row.get("flow_ratio", FLOW_RATIO))
        kp = float(row["kp"])
        value_um = float(row["calibrated_cfd_axial_um"])
        standard_deviation_um = float(row["calibrated_cfd_standard_deviation_um"])
        raw_ci = row.get("calibrated_cfd_95_ci_half_width_um", "")
        ci_half_width_um = float(raw_ci) if raw_ci != "" else None
        ratio_references = [
            item
            for item in references
            if math.isclose(float(item["flow_ratio"]), flow_ratio)
        ]
        paper = min(ratio_references, key=lambda item: abs(float(item["kp"]) - kp))
        paper_ratio = float(paper["d_over_w"])
        publication_rows.append({
            "flow_ratio": flow_ratio,
            "kp": kp,
            "series": "calibrated_cfd_axial",
            "d_over_w": value_um / WIDTH_UM,
            "diameter_um": value_um,
            "uncertainty_type": "droplet_standard_deviation_n_50",
            "uncertainty_d_over_w": standard_deviation_um / WIDTH_UM,
            "uncertainty_um": standard_deviation_um,
            "mean_95_ci_half_width_d_over_w": "" if ci_half_width_um is None else ci_half_width_um / WIDTH_UM,
            "mean_95_ci_half_width_um": "" if ci_half_width_um is None else ci_half_width_um,
            "percent_error_vs_digitized": 100 * (value_um / WIDTH_UM - paper_ratio) / paper_ratio,
        })
    write_csv(HERE / "figure6_all_series.csv", publication_rows)

    width, height = 1200, 790
    left, right, top, bottom = 105, 45, 105, 140
    xmin, xmax, ymin, ymax = -0.7, 2.8, 1.68, 2.58
    plot_right = width - right
    plot_bottom = height - bottom

    def x(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (plot_right - left)

    def y(value: float) -> float:
        return plot_bottom - (value - ymin) / (ymax - ymin) * (plot_bottom - top)

    colours = {0.5: "#0072B2", 1.0: "#D55E00", 1.5: "#009E73"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111827}.axis{stroke:#111827;stroke-width:2}.major{stroke:#cbd5e1;stroke-width:1}.minor{stroke:#e5e7eb;stroke-width:.8}.tick{font-size:14px}.legend{font-size:13px}</style>',
        f'<clipPath id="plot"><rect x="{left}" y="{top}" width="{plot_right-left}" height="{plot_bottom-top}"/></clipPath>',
        '<text x="600" y="32" font-size="23" font-weight="600" text-anchor="middle">Droplet size versus membrane-control parameter</text>',
        '<text x="600" y="58" font-size="14" text-anchor="middle">Raj et al. Figure 6b data, Equation 13, and fixed-factor calibrated CFD</text>',
        '<text x="600" y="80" font-size="12.5" text-anchor="middle">Raj markers: +/-1 experimental SD; CFD diamonds: +/-1 droplet SD (n = 50)</text>',
    ]

    for index in range(15):
        tick = xmin + index * 0.25
        klass = "major" if index % 2 == 0 else "minor"
        parts.append(f'<line class="{klass}" x1="{x(tick):.1f}" y1="{top}" x2="{x(tick):.1f}" y2="{plot_bottom}"/>')
        if index % 2 == 0:
            parts.append(f'<text class="tick" x="{x(tick):.1f}" y="{plot_bottom+26}" text-anchor="middle">{tick:.1f}</text>')
    for index in range(19):
        tick = ymin + index * 0.05
        klass = "major" if index % 2 == 0 else "minor"
        parts.append(f'<line class="{klass}" x1="{left}" y1="{y(tick):.1f}" x2="{plot_right}" y2="{y(tick):.1f}"/>')
        if index % 2 == 0:
            parts.append(f'<text class="tick" x="{left-13}" y="{y(tick)+5:.1f}" text-anchor="end">{tick:.2f}</text>')
    parts += [
        f'<rect x="{left}" y="{top}" width="{plot_right-left}" height="{plot_bottom-top}" fill="none" class="axis"/>',
        f'<text x="{(left+plot_right)/2:.1f}" y="{plot_bottom+49}" font-size="18" text-anchor="middle">Controlling parameter, K<tspan baseline-shift="sub" font-size="13">p</tspan></text>',
        f'<text x="31" y="{(top+plot_bottom)/2:.1f}" font-size="18" text-anchor="middle" transform="rotate(-90 31 {(top+plot_bottom)/2:.1f})">Normalized droplet diameter, D<tspan baseline-shift="sub" font-size="13">d</tspan>/w</text>',
        '<g clip-path="url(#plot)">',
    ]

    for flow_ratio in (0.5, 1.0, 1.5):
        colour = colours[flow_ratio]
        equation_points = " ".join(
            f'{x(xmin + (xmax-xmin)*index/240):.1f},{y(correlation(xmin + (xmax-xmin)*index/240, flow_ratio)/WIDTH_UM):.1f}'
            for index in range(241)
        )
        parts.append(f'<polyline points="{equation_points}" fill="none" stroke="{colour}" stroke-width="3"/>')
        for row in references:
            if not math.isclose(float(row["flow_ratio"]), flow_ratio):
                continue
            kp = float(row["kp"])
            value = float(row["d_over_w"])
            uncertainty = float(row["uncertainty"])
            px = x(kp)
            top_y = y(value + uncertainty)
            bottom_y = y(value - uncertainty)
            parts += [
                f'<line x1="{px:.1f}" y1="{top_y:.1f}" x2="{px:.1f}" y2="{bottom_y:.1f}" stroke="{colour}" stroke-width="1.8"/>',
                f'<line x1="{px-6:.1f}" y1="{top_y:.1f}" x2="{px+6:.1f}" y2="{top_y:.1f}" stroke="{colour}" stroke-width="1.8"/>',
                f'<line x1="{px-6:.1f}" y1="{bottom_y:.1f}" x2="{px+6:.1f}" y2="{bottom_y:.1f}" stroke="{colour}" stroke-width="1.8"/>',
                f'<circle cx="{px:.1f}" cy="{y(value):.1f}" r="5.5" fill="white" stroke="{colour}" stroke-width="2.5"/>',
            ]

    completed_ratios: list[float] = []
    for flow_ratio in (0.5, 1.0, 1.5):
        cfd_rows = sorted(
            (
                row
                for row in calibrated
                if math.isclose(float(row.get("flow_ratio", FLOW_RATIO)), flow_ratio)
            ),
            key=lambda row: float(row["kp"]),
        )
        if not cfd_rows:
            continue
        completed_ratios.append(flow_ratio)
        colour = colours[flow_ratio]
        cfd_line = " ".join(
            f'{x(float(row["kp"])):.1f},{y(float(row["calibrated_cfd_axial_um"])/WIDTH_UM):.1f}'
            for row in cfd_rows
        )
        parts.append(f'<polyline points="{cfd_line}" fill="none" stroke="{colour}" stroke-width="3" stroke-dasharray="9 6"/>')
        for row in cfd_rows:
            kp = float(row["kp"])
            value = float(row["calibrated_cfd_axial_um"]) / WIDTH_UM
            uncertainty = float(row["calibrated_cfd_standard_deviation_um"]) / WIDTH_UM
            px, py = x(kp), y(value)
            top_y, bottom_y = y(value + uncertainty), y(value - uncertainty)
            points = f'{px:.1f},{py-8:.1f} {px+8:.1f},{py:.1f} {px:.1f},{py+8:.1f} {px-8:.1f},{py:.1f}'
            parts += [
                f'<line x1="{px:.1f}" y1="{top_y:.1f}" x2="{px:.1f}" y2="{bottom_y:.1f}" stroke="{colour}" stroke-width="2.5"/>',
                f'<line x1="{px-8:.1f}" y1="{top_y:.1f}" x2="{px+8:.1f}" y2="{top_y:.1f}" stroke="{colour}" stroke-width="2.5"/>',
                f'<line x1="{px-8:.1f}" y1="{bottom_y:.1f}" x2="{px+8:.1f}" y2="{bottom_y:.1f}" stroke="{colour}" stroke-width="2.5"/>',
                f'<polygon points="{points}" fill="{colour}" stroke="white" stroke-width="1.2"/>',
            ]
    parts.append('</g>')

    legend_y = 731
    for index, flow_ratio in enumerate((0.5, 1.0, 1.5)):
        lx = 110 + index * 245
        colour = colours[flow_ratio]
        parts += [
            f'<line x1="{lx}" y1="{legend_y}" x2="{lx+35}" y2="{legend_y}" stroke="{colour}" stroke-width="3"/>',
            f'<circle cx="{lx+17.5}" cy="{legend_y}" r="5" fill="white" stroke="{colour}" stroke-width="2"/>',
            f'<text class="legend" x="{lx+45}" y="{legend_y+5}">Raj data + Eq. 13, r = {flow_ratio:.1f}</text>',
        ]
    lx = 855
    parts += [
        f'<line x1="{lx}" y1="{legend_y}" x2="{lx+35}" y2="{legend_y}" stroke="#374151" stroke-width="3" stroke-dasharray="8 5"/>',
        f'<polygon points="{lx+17.5},{legend_y-7} {lx+24.5},{legend_y} {lx+17.5},{legend_y+7} {lx+10.5},{legend_y}" fill="#374151"/>',
        f'<text class="legend" x="{lx+45}" y="{legend_y+5}">Calibrated CFD +/- SD</text>',
        f'<text x="600" y="766" font-size="12.5" text-anchor="middle">CFD points present: {len(calibrated)}/18; completed flow ratios: {", ".join(f"r = {ratio:.1f}" for ratio in completed_ratios)}.</text>',
        '</svg>',
    ]
    (HERE / "figure6_all_series.svg").write_text("\n".join(parts) + "\n")


def main() -> None:
    rows = build_rows()
    write_csv(HERE / "figure6_absolute.csv", rows)
    absolute_svg(rows, HERE / "figure6_absolute.svg")
    normalized_outputs(rows)
    calibrated = calibrated_outputs(rows)
    sweep = completed_sweep_rows()
    publication_outputs(sweep or calibrated)


if __name__ == "__main__":
    main()
