#!/usr/bin/env python3
"""Generate raw Figure 6 validation cases."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "cases" / "flexible-rigid"
sys.path.insert(0, str(SOURCE))

from mesh_generator import generate_block_mesh  # noqa: E402


WIDTH = 100e-6
BASE_HEIGHT = 100e-6
CONTINUOUS_MU = 0.0287
SIGMA = 0.003501
WRITE_INTERVAL = 0.01
MAX_DELTA_T = 1e-5
END_TIME_BY_FLOW_RATIO = {0.5: 3.3, 1.0: 2.5, 1.5: 2.0}


def include_path(text: str) -> str:
    return text.replace("../../common/", "../../../../../cases/flexible-rigid/common/")


def replace_entry(text: str, key: str, value: str) -> str:
    pattern = rf"(?m)^(\s*{re.escape(key)}\s+)[^;]+;"
    replaced, count = re.subn(pattern, rf"\g<1>{value};", text, count=1)
    if count != 1:
        raise RuntimeError(f"expected one {key} entry, found {count}")
    return replaced


def load_rows(group: str) -> list[dict[str, str]]:
    with (HERE / "matrix.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if group == "all":
        return rows
    if group in {"fig6_r05", "fig6_r10", "fig6_r15"}:
        ratio = {"fig6_r05": "0.5", "fig6_r10": "1.0", "fig6_r15": "1.5"}[group]
        return [row for row in rows if row["flow_ratio"] == ratio]
    requested = set(group.split(","))
    return [row for row in rows if row["case_id"] in requested]


def build(row: dict[str, str], overwrite: bool, run_root: Path, mesh_level: str, contact_angle: float) -> None:
    case = run_root / row["case_id"]
    if case.exists() and overwrite:
        shutil.rmtree(case)
    if case.exists():
        return
    for directory in (case / "0", case / "constant", case / "system"):
        directory.mkdir(parents=True, exist_ok=True)

    kp = float(row["kp"])
    ratio = float(row["flow_ratio"])
    ca = float(row["capillary_number"])
    mu_d = float(row["dispersed_dynamic_viscosity_pa_s"])
    rho_d = float(row["dispersed_density_kg_m3"])
    qc = ca * SIGMA / CONTINUOUS_MU * WIDTH * BASE_HEIGHT
    qd = ratio * qc

    alpha = (SOURCE / "common" / "alpha.water.initial").read_text()
    alpha = re.sub(r"theta0\s+[^;]+;", f"theta0 {contact_angle:.8g};", alpha)
    (case / "0" / "alpha.water").write_text(alpha)
    (case / "0" / "p_rgh").write_text(
        "FoamFile { format ascii; class volScalarField; object p_rgh; }\n"
        '#include "../../../../../cases/flexible-rigid/common/p_rgh"\n'
    )
    velocity = (SOURCE / "rigid" / "0" / "U").read_text()
    velocity = re.sub(r"Qc\s+[^;]+;", f"Qc {qc:.12g};", velocity)
    velocity = re.sub(r"Qd\s+[^;]+;", f"Qd {qd:.12g};", velocity)
    (case / "0" / "U").write_text(include_path(velocity))

    for name in ("g", "momentumTransport", "phaseProperties", "physicalProperties.oil"):
        (case / "constant" / name).write_text(include_path((SOURCE / "rigid" / "constant" / name).read_text()))
    water = (SOURCE / "rigid" / "constant" / "physicalProperties.water").read_text()
    water = re.sub(r"nu\s+[^;]+;", f"nu {mu_d / rho_d:.12g};", water)
    (case / "constant" / "physicalProperties.water").write_text(include_path(water))

    block_mesh, geometry = generate_block_mesh(kp, mesh_level)
    (case / "system" / "blockMeshDict").write_text(block_mesh)
    for name in ("fvSchemes", "fvSolution", "setFieldsDict"):
        (case / "system" / name).write_text(include_path((SOURCE / "rigid" / "system" / name).read_text()))
    end_time = END_TIME_BY_FLOW_RATIO[ratio]
    control = "FoamFile { format ascii; class dictionary; object controlDict; }\n" + (SOURCE / "common" / "controlDict").read_text()
    for key, value in (
        ("endTime", f"{end_time:g}"),
        ("writeInterval", f"{WRITE_INTERVAL:g}"),
        ("maxDeltaT", f"{MAX_DELTA_T:.8g}"),
        ("writeFormat", "binary"),
        ("writeCompression", "on"),
    ):
        control = replace_entry(control, key, value)
    (case / "system" / "controlDict").write_text(control)

    metadata = dict(row)
    metadata.update(
        continuous_flow_rate_m3_s=qc,
        dispersed_flow_rate_m3_s=qd,
        contact_angle_degrees=contact_angle,
        end_time_s=end_time,
        write_interval_s=WRITE_INTERVAL,
        maximum_time_step_s=MAX_DELTA_T,
        output_format="binary_compressed",
        geometry=geometry,
    )
    (case / "case.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", default="fig6_r05")
    parser.add_argument("--mesh", choices=("screening", "medium", "fine"), default="screening")
    parser.add_argument("--contact-angle", type=float, default=147.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    selected = load_rows(args.group)
    run_root = HERE / f"runs-{args.mesh}"
    run_root.mkdir(exist_ok=True)
    for row in selected:
        build(row, args.overwrite, run_root, args.mesh, args.contact_angle)
    print(f"Generated {len(selected)} cases in {run_root}")


if __name__ == "__main__":
    main()
