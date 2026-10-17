#!/usr/bin/env python3
"""Fit and plot one r=1.0 multiplier for the urgent meeting preview."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "meeting_preview_completed_r10.csv"
REFERENCE = HERE / "reference_data" / "fig6_digitized.csv"
OUTPUT = HERE / "meeting_preview_refit_r10.csv"


def load_preview_module():
    path = HERE / "quick_completed_preview.py"
    spec = importlib.util.spec_from_file_location("meeting_preview", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    with SOURCE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    with REFERENCE.open(newline="") as handle:
        references = [row for row in csv.DictReader(handle) if float(row["flow_ratio"]) == 1.0]
    pairs: list[tuple[float, float]] = []
    paper_by_case: dict[str, float] = {}
    for row in rows:
        kp = float(row["kp"])
        paper = min(references, key=lambda item: abs(float(item["kp"]) - kp))
        raw = float(row["raw_d_over_w"])
        target = float(paper["d_over_w"])
        pairs.append((raw, target))
        paper_by_case[row["case_id"]] = target
    factor = sum(raw * target for raw, target in pairs) / sum(raw**2 for raw, _ in pairs)
    fitted: list[dict[str, object]] = []
    for row in rows:
        raw = float(row["raw_d_over_w"])
        raw_std = float(row["raw_spatial_std"])
        target = paper_by_case[row["case_id"]]
        value = factor * raw
        fitted.append({
            **row,
            "calibrated_d_over_w": value,
            "calibrated_spatial_std": factor * raw_std,
            "meeting_fit_factor": factor,
            "error_vs_paper_percent": 100 * (value - target) / target,
            "fit_role": "in_sample_r1_global_least_squares",
        })
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fitted[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(fitted)
    preview = load_preview_module()
    preview.write_svg(
        fitted,
        references,
        output_name="meeting_preview_refit_r10.svg",
        calibration_note=f"in-sample r = 1.0 fit factor = {factor:.6f}",
        calibrated_label="Refitted CFD",
        footer_note="In-sample presentation fit to these six points; not an independent validation result.",
    )
    mean_error = sum(abs(float(row["error_vs_paper_percent"])) for row in fitted) / len(fitted)
    print(f"factor={factor:.12f}; in-sample MAE={mean_error:.3f}%")


if __name__ == "__main__":
    main()
