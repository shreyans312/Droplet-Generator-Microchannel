#!/usr/bin/env python3
"""Write corrected alpha=0.5 silhouette components at every saved time."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path


def load_measure(path: Path):
    spec = importlib.util.spec_from_file_location("droplet_measure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--measure-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    measure = load_measure(args.measure_script.resolve())
    times, frames, _, _ = measure.case_frames(args.case.resolve(), 0.5)
    rows = []
    for time_value, detections in zip(times, frames):
        for component_id, item in enumerate(detections, start=1):
            rows.append(
                {
                    "time_s": time_value,
                    "component_id": component_id,
                    "centroid_x_um": item["cx"] / measure.UM,
                    "centroid_y_um": item["cy"] / measure.UM,
                    "axial_diameter_um": item["axial_diameter"] / measure.UM,
                    "transverse_diameter_um": item["transverse_diameter"] / measure.UM,
                    "maximum_feret_um": item["maximum_feret"] / measure.UM,
                    "area_equivalent_diameter_um": item["area_equivalent_diameter"] / measure.UM,
                }
            )
    fields = list(rows[0]) if rows else [
        "time_s", "component_id", "centroid_x_um", "centroid_y_um",
        "axial_diameter_um", "transverse_diameter_um", "maximum_feret_um",
        "area_equivalent_diameter_um",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} component samples")


if __name__ == "__main__":
    main()
