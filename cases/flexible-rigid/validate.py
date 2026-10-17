#!/usr/bin/env python3
"""Compare raw, image-matched droplet measurements with the paper."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
WIDTH = 100e-6
FLOW_RATIO = 0.5
VISCOSITY_RATIO = 0.000914 / 0.0287
PRIMARY_THRESHOLD = 0.5


def reference_diameter(kp: float) -> float:
    return (
        1.97064
        * (1 + kp) ** -0.0655
        * FLOW_RATIO**0.13593
        * VISCOSITY_RATIO**-0.0273
        * WIDTH
    )


def _mean(rows: list[dict[str, str]], field: str) -> float | None:
    values = [float(row[field]) for row in rows if row.get(field, "")]
    return statistics.mean(values) if values else None


def _stdev(rows: list[dict[str, str]], field: str) -> float | None:
    values = [float(row[field]) for row in rows if row.get(field, "")]
    return statistics.stdev(values) if len(values) > 1 else None


def analyze_measurements() -> list[dict[str, object]]:
    path = ROOT / "droplet_measurements.csv"
    if not path.exists():
        raise SystemExit("droplet_measurements.csv is missing; run /usr/bin/python3 measure.py")
    with path.open(newline="") as handle:
        measurements = list(csv.DictReader(handle))

    results: list[dict[str, object]] = []
    for case_name, kp in (("rigid", 0.0), ("flexible", 1.695)):
        case_rows = [row for row in measurements if row["case"] == case_name]
        primary = [row for row in case_rows if math.isclose(float(row["alpha_threshold"]), PRIMARY_THRESHOLD)]
        low = [row for row in case_rows if math.isclose(float(row["alpha_threshold"]), 0.25)]
        high = [row for row in case_rows if math.isclose(float(row["alpha_threshold"]), 0.75)]
        axial = _mean(primary, "axial_diameter_um")
        transverse = _mean(primary, "transverse_diameter_um")
        area_equivalent = _mean(primary, "area_equivalent_diameter_um")
        maximum_feret = _mean(primary, "maximum_feret_um")
        reference_um = reference_diameter(kp) * 1e6
        error = None if axial is None else 100 * (axial - reference_um) / reference_um
        low_axial = _mean(low, "axial_diameter_um")
        high_axial = _mean(high, "axial_diameter_um")
        sensitivity_values = [value for value in (low_axial, high_axial) if value is not None]
        sensitivity = None
        if axial and sensitivity_values:
            sensitivity = 100 * max(abs(value - axial) for value in sensitivity_values) / axial
        aspect_ratio = None if not axial or not transverse else axial / transverse
        axial_std = _stdev(primary, "axial_diameter_um")
        flags: list[str] = []
        if len(primary) < 50:
            flags.append("insufficient_droplets")
        if any(row["resolution_flag"] != "resolved" for row in primary):
            flags.append("under_resolved")
        if aspect_ratio is not None and not 0.9 <= aspect_ratio <= 1.1:
            flags.append("shape_mismatch")
        if sensitivity is not None and sensitivity > 5:
            flags.append("threshold_sensitive")
        if not primary:
            flags.append("no_crossings")
        results.append(
            {
                "case": case_name,
                "Kp": kp,
                "reference_diameter_um": reference_um,
                "raw_axial_diameter_um": "" if axial is None else axial,
                "raw_transverse_diameter_um": "" if transverse is None else transverse,
                "raw_area_equivalent_diameter_um": "" if area_equivalent is None else area_equivalent,
                "raw_maximum_feret_um": "" if maximum_feret is None else maximum_feret,
                "axial_error_percent": "" if error is None else error,
                "axial_std_um": "" if axial_std is None else axial_std,
                "mean_axial_to_transverse_ratio": "" if aspect_ratio is None else aspect_ratio,
                "alpha_0.25_axial_mean_um": "" if low_axial is None else low_axial,
                "alpha_0.75_axial_mean_um": "" if high_axial is None else high_axial,
                "threshold_sensitivity_percent": "" if sensitivity is None else sensitivity,
                "tracked_droplets": len(primary),
                "under_resolved_droplets": sum(row["resolution_flag"] != "resolved" for row in primary),
                "quality_flag": "validation_candidate" if not flags else ";".join(flags),
            }
        )
    return results


def write_results(rows: list[dict[str, object]]) -> None:
    with (ROOT / "validation_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_svg(rows: list[dict[str, object]]) -> None:
    width, height = 760, 500
    left, top, bottom = 90, 55, 75
    ymax = 240.0

    def y(value: float) -> float:
        return height - bottom - value / ymax * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.grid{stroke:#d1d5db}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="380" y="30" font-size="21" text-anchor="middle">Image-matched region-B diameter</text>',
    ]
    for tick in range(0, 241, 20):
        parts += [
            f'<line class="grid" x1="{left}" y1="{y(tick):.1f}" x2="730" y2="{y(tick):.1f}"/>',
            f'<text x="{left-12}" y="{y(tick)+5:.1f}" font-size="14" text-anchor="end">{tick}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="730" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
        f'<text x="24" y="{(top+height-bottom)/2}" font-size="17" text-anchor="middle" transform="rotate(-90 24 {(top+height-bottom)/2})">Diameter (um)</text>',
    ]
    for x, row in zip((270, 550), rows):
        paper = float(row["reference_diameter_um"])
        parts += [
            f'<circle cx="{x}" cy="{y(paper):.1f}" r="7" fill="#2563eb"/>',
            f'<text x="{x}" y="{height-bottom+28}" font-size="16" text-anchor="middle">{row["case"]}</text>',
        ]
        if str(row["raw_axial_diameter_um"]):
            axial = float(row["raw_axial_diameter_um"])
            parts.append(f'<path d="M {x-8} {y(axial)-8:.1f} L {x+8} {y(axial)+8:.1f} M {x-8} {y(axial)+8:.1f} L {x+8} {y(axial)-8:.1f}" stroke="#dc2626" stroke-width="4"/>')
        if str(row["raw_area_equivalent_diameter_um"]):
            area = float(row["raw_area_equivalent_diameter_um"])
            parts.append(f'<rect x="{x-6}" y="{y(area)-6:.1f}" width="12" height="12" fill="#d97706"/>')
    parts += [
        '<circle cx="155" cy="468" r="6" fill="#2563eb"/>',
        '<text x="167" y="473" font-size="14">paper</text>',
        '<path d="M 270 462 L 282 474 M 270 474 L 282 462" stroke="#dc2626" stroke-width="3"/>',
        '<text x="289" y="473" font-size="14">axial</text>',
        '<rect x="390" y="462" width="12" height="12" fill="#d97706"/>',
        '<text x="409" y="473" font-size="14">area-equivalent</text>',
        '</svg>',
    ]
    (ROOT / "validation_plot.svg").write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot-existing", action="store_true")
    args = parser.parse_args()
    if args.plot_existing:
        with (ROOT / "validation_results.csv").open(newline="") as handle:
            rows: list[dict[str, object]] = list(csv.DictReader(handle))
    else:
        rows = analyze_measurements()
        write_results(rows)
    write_svg(rows)
    for row in rows:
        print(
            f'{row["case"]}: paper={float(row["reference_diameter_um"]):.2f} um, '
            f'axial={row["raw_axial_diameter_um"]}, area-equivalent={row["raw_area_equivalent_diameter_um"]}, '
            f'n={row["tracked_droplets"]}, {row["quality_flag"]}'
        )


if __name__ == "__main__":
    main()
