#!/usr/bin/env /usr/bin/python3
"""Fail-closed verification report for the existing Raj Figure 6 CFD sweep.

The gate intentionally cannot certify experimental validation while the
Figure 6 device mask and optical diameter definition are undocumented.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

from volume_audit import audit


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def baseline_volumes(case: Path, time_value: float) -> tuple[float | None, float | None]:
    rows = [
        row for row in audit(case, time_value)
        if row["alpha_cutoff"] == 0.5 and not row["attached_to_inlet"]
    ]
    upstream = [
        float(row["water_volume_pL"])
        for row in rows
        if 400 <= float(row["x_centroid_um"]) <= 600
        and float(row["water_volume_pL"]) > 500
    ]
    downstream = [
        float(row["water_volume_pL"])
        for row in rows
        if 1200 <= float(row["x_centroid_um"]) <= 1750
        and row["fully_in_region_b"]
        and float(row["water_volume_pL"]) > 500
    ]
    return (
        statistics.median(upstream) if upstream else None,
        statistics.median(downstream) if downstream else None,
    )


def render(mesh: str, time_value: float) -> str:
    with (HERE / f"results-{mesh}.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    with (ROOT / "studies/mesh-convergence/mesh_convergence.csv").open(newline="") as handle:
        convergence = list(csv.DictReader(handle))

    sampled = sum(int(row["tracked_droplets"]) >= 50 for row in rows)
    raw_errors = [abs(float(row["raw_error_percent"])) for row in rows]
    mesh_pass = all(
        row["trend"] not in {"non_monotonic", "insufficient_data"}
        and math.isfinite(float(row["fine_mean_um"]))
        for row in convergence
    )
    lines = [
        f"# Figure 6 publication-readiness gate ({mesh} mesh)",
        "",
        "**Status: FAIL — do not present these CFD points as validated predictions.**",
        "",
        f"- Droplet-count criterion: {sampled}/{len(rows)} conditions reach 50 tracked droplets.",
        f"- Raw diameter comparison: {statistics.mean(raw_errors):.2f}% mean absolute error; "
        f"{max(raw_errors):.2f}% maximum error.",
        f"- Existing three-level mesh study: {'PASS' if mesh_pass else 'FAIL'}; "
        "the fine rigid case has no crossing and flexible droplet size changes non-monotonically.",
        "- Experimental geometry/observable match: FAIL; the paper does not give a complete "
        "Figure 6 T-junction mask or a precise optical-diameter algorithm.",
        "- Flexible-wall physics: the current mesh has an imposed static profile, "
        "not pressure-coupled two-way fluid–structure interaction.",
        "",
        f"## Direct 3-D VOF volume check at t = {time_value:g} s",
        "",
        "Volumes are sums of alpha.water × CFD cell volume for *alpha ≥ 0.5* "
        "connected cores. The core excludes some interfacial water, so values "
        "are conservative lower estimates. Upstream is x = 400–600 µm; "
        "downstream is x = 1200–1750 µm. Rows are near Kp = 0.",
        "",
        "| Flow ratio | Upstream core (pL) | Downstream core (pL) | Down/up | Interpretation |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    baselines = [
        (0.5, "fig6_r050_km00010"),
        (1.0, "fig6_r100_km00010"),
        (1.5, "fig6_r150_km00010"),
    ]
    for ratio, name in baselines:
        upstream, downstream = baseline_volumes(HERE / f"runs-{mesh}" / name, time_value)
        if upstream is None or downstream is None:
            lines.append(f"| {ratio:g} | unavailable | unavailable | — | inconclusive snapshot |")
            continue
        volume_ratio = downstream / upstream
        interpretation = (
            "approximately conserved"
            if abs(volume_ratio - 1) <= 0.2
            else "topology/merger investigation required"
        )
        lines.append(
            f"| {ratio:g} | {upstream:.0f} | {downstream:.0f} | "
            f"{volume_ratio:.2f} | {interpretation} |"
        )
    lines += [
        "",
        "The r = 1.0 and 1.5 increases are visible even at alpha = 0.5; they "
        "cannot be repaired by changing a plotted diameter or calibration. "
        "A connectivity change can arise from physical coalescence or from "
        "under-resolved oil films/VOF bridging. The present runs do not "
        "distinguish these mechanisms.",
        "",
        "## Required before a publication claim",
        "",
        "1. Obtain or measure the experimental T-junction geometry and optical "
        "diameter rule, including where in region B it was read.",
        "2. Resolve the r-dependent connectivity change with a local interface/film "
        "refinement and time-step study; compare upstream generation volume, "
        "downstream volume, and outlet dispersed-phase flux.",
        "3. Re-run a geometry-matched rigid baseline and then pressure-coupled "
        "flexible case; demonstrate mesh/time-step independence and ≥50 "
        "steady droplets per condition. Report raw errors before any fit.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", choices=("medium",), default="medium")
    parser.add_argument("--time", type=float, default=1.5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = render(args.mesh, args.time)
    if args.output:
        args.output.write_text(result)
        print(args.output)
    else:
        print(result)


if __name__ == "__main__":
    main()
