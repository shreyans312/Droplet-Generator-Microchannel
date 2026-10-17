#!/usr/bin/env python3
"""Prepare fixed-medium-mesh cases for the time-step convergence study."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "cases" / "flexible-rigid"
RUNS = HERE / "runs"
END_TIME = 3.0
WRITE_INTERVAL = 0.005
VARIANTS = {
    "coarse_40us": 40e-6,
    "baseline_20us": 20e-6,
    "fine_10us": 10e-6,
    "ultrafine_5us": 5e-6,
}


def replace_entry(text: str, key: str, value: str) -> str:
    pattern = rf"(?m)^(\s*{re.escape(key)}\s+)[^;]+;"
    replaced, count = re.subn(pattern, rf"\g<1>{value};", text, count=1)
    if count != 1:
        raise RuntimeError(f"expected one {key} entry, found {count}")
    return replaced


def load_mesh_generator():
    path = SOURCE / "mesh_generator.py"
    spec = importlib.util.spec_from_file_location("mesh_generator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def copy_case(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True)
    shutil.copytree(source / "0", destination / "0")
    shutil.copytree(
        source / "constant",
        destination / "constant",
        ignore=shutil.ignore_patterns("polyMesh"),
    )
    shutil.copytree(source / "system", destination / "system")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing generated runs directory",
    )
    parser.add_argument(
        "--variant",
        choices=tuple(VARIANTS),
        help="prepare one missing variant without replacing completed runs",
    )
    args = parser.parse_args()
    if RUNS.exists():
        if args.force:
            shutil.rmtree(RUNS)
        elif not args.variant:
            raise SystemExit(f"{RUNS} exists; use --force to replace generated cases")

    mesh_generator = load_mesh_generator()
    source_control = (SOURCE / "common" / "controlDict").read_text()
    expected = {
        "adjustTimeStep": "yes",
        "maxCo": "0.25",
        "maxAlphaCo": "0.15",
        "maxDeltaT": "2e-5",
    }
    for key, value in expected.items():
        if not re.search(rf"(?m)^\s*{key}\s+{re.escape(value)}\s*;", source_control):
            raise RuntimeError(f"active baseline no longer has {key}={value}")

    manifest = {
        "purpose": "time-step convergence; not experimental validation",
        "mesh_level": "medium",
        "end_time_s": END_TIME,
        "write_interval_s": WRITE_INTERVAL,
        "initial_delta_t_s": 1e-6,
        "adaptive_controls": {"maxCo": 0.25, "maxAlphaCo": 0.15},
        "variants": VARIANTS,
        "notes": [
            "Only maxDeltaT changes between temporal variants.",
            "Binary compressed output changes storage only, not solver arithmetic.",
            "Actual accepted time steps must be read from solver logs.",
        ],
    }
    RUNS.mkdir(parents=True, exist_ok=True)
    (RUNS / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    selected = {args.variant: VARIANTS[args.variant]} if args.variant else VARIANTS
    for variant, cap in selected.items():
        if (RUNS / variant).exists():
            raise SystemExit(f"generated variant already exists: {RUNS / variant}")
        common = RUNS / variant / "common"
        shutil.copytree(SOURCE / "common", common)
        control = source_control
        for key, value in (
            ("endTime", f"{END_TIME:g}"),
            ("writeInterval", f"{WRITE_INTERVAL:g}"),
            ("maxDeltaT", f"{cap:.8g}"),
            ("writeFormat", "binary"),
            ("writeCompression", "on"),
        ):
            control = replace_entry(control, key, value)
        (common / "controlDict").write_text(control)

        for case_name, kp in (("rigid", 0.0), ("flexible", 1.695)):
            case = RUNS / variant / case_name
            copy_case(SOURCE / case_name, case)
            block_mesh, geometry = mesh_generator.generate_block_mesh(kp, "medium")
            (case / "system" / "blockMeshDict").write_text(block_mesh)
            (case / "geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
            print(f"prepared {variant}/{case_name}: cap={cap:.8g} s, {geometry['estimated_cells']} cells")


if __name__ == "__main__":
    main()
