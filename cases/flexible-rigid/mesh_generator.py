#!/usr/bin/env python3
"""Generate a T-junction, widened measurement region, and static wall profile."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
UM = 1e-6
WIDTH = 100 * UM
BASE_HEIGHT = 100 * UM
REGION_B_WIDTH = 200 * UM
MEMBRANE_THICKNESS = 22.5 * UM
POISSON_RATIO = 0.5
BOTTOM = -50 * UM


@dataclass(frozen=True)
class FootprintBlock:
    x0: float
    x1: float
    y0: float
    y1: float


def average_deflection(kp: float) -> float:
    if kp == 0:
        return 0.0
    magnitude = (
        3 * WIDTH * MEMBRANE_THICKNESS**2 * abs(kp)
        / (64 * (1 + POISSON_RATIO))
    ) ** (1 / 3)
    return math.copysign((2 / 3) * magnitude, kp)


def axial_envelope(x: float) -> float:
    if -150 * UM <= x <= 0:
        return 0.5 * (1 + math.cos(math.pi * x / (150 * UM)))
    if 0 < x <= 450 * UM:
        return 0.5 * (1 + math.cos(math.pi * x / (450 * UM)))
    return 0.0


def top_height(x: float, y: float, kp: float) -> float:
    if abs(y) > WIDTH / 2:
        return BOTTOM + BASE_HEIGHT
    transverse = 1 - (2 * y / WIDTH) ** 2
    maximum = 1.5 * average_deflection(kp)
    return BOTTOM + BASE_HEIGHT - maximum * axial_envelope(x) * transverse


def intervals(values: list[float]) -> list[tuple[float, float]]:
    return list(zip(values, values[1:]))


def footprints() -> list[FootprintBlock]:
    main_x = [-600, -300, -150, -50, -25, 0, 25, 50, 150, 300, 450, 600]
    main_y = [-50, -25, 0, 25, 50]
    branch_x = [-50, -25, 0, 25, 50]
    branch_y = [50, 100, 150, 300, 650]
    region_x = [600, 800, 1200, 1900]
    region_y = [-100, -50, -25, 0, 25, 50, 100]
    blocks = [
        FootprintBlock(x0 * UM, x1 * UM, y0 * UM, y1 * UM)
        for x0, x1 in intervals(main_x)
        for y0, y1 in intervals(main_y)
    ]
    blocks += [
        FootprintBlock(x0 * UM, x1 * UM, y0 * UM, y1 * UM)
        for x0, x1 in intervals(branch_x)
        for y0, y1 in intervals(branch_y)
    ]
    blocks += [
        FootprintBlock(x0 * UM, x1 * UM, y0 * UM, y1 * UM)
        for x0, x1 in intervals(region_x)
        for y0, y1 in intervals(region_y)
    ]
    return blocks


def generate_block_mesh(kp: float, mesh_level: str) -> tuple[str, dict[str, float | int | str]]:
    target = {"screening": 25 * UM, "medium": 12.5 * UM, "fine": 8.333333 * UM}[mesh_level]
    nz = {"screening": 4, "medium": 8, "fine": 12}[mesh_level]
    vertices: list[tuple[float, float, float]] = []
    vertex_ids: dict[tuple[float, float, float], int] = {}
    block_rows: list[tuple[tuple[int, ...], int, int, int]] = []
    face_occurrences: dict[frozenset[int], list[tuple[tuple[int, ...], str]]] = defaultdict(list)

    def vertex(x: float, y: float, z: float) -> int:
        key = (round(x, 13), round(y, 13), round(z, 13))
        if key not in vertex_ids:
            vertex_ids[key] = len(vertices)
            vertices.append(key)
        return vertex_ids[key]

    cell_count = 0
    for block in footprints():
        z00 = top_height(block.x0, block.y0, kp)
        z10 = top_height(block.x1, block.y0, kp)
        z11 = top_height(block.x1, block.y1, kp)
        z01 = top_height(block.x0, block.y1, kp)
        ids = (
            vertex(block.x0, block.y0, BOTTOM),
            vertex(block.x1, block.y0, BOTTOM),
            vertex(block.x1, block.y1, BOTTOM),
            vertex(block.x0, block.y1, BOTTOM),
            vertex(block.x0, block.y0, z00),
            vertex(block.x1, block.y0, z10),
            vertex(block.x1, block.y1, z11),
            vertex(block.x0, block.y1, z01),
        )
        nx = max(1, round((block.x1 - block.x0) / target))
        ny = max(1, round((block.y1 - block.y0) / target))
        block_rows.append((ids, nx, ny, nz))
        cell_count += nx * ny * nz
        faces = {
            "xmin": (ids[0], ids[4], ids[7], ids[3]),
            "xmax": (ids[1], ids[2], ids[6], ids[5]),
            "ymin": (ids[0], ids[1], ids[5], ids[4]),
            "ymax": (ids[3], ids[7], ids[6], ids[2]),
            "bottom": (ids[0], ids[3], ids[2], ids[1]),
            "top": (ids[4], ids[5], ids[6], ids[7]),
        }
        for face_type, face in faces.items():
            face_occurrences[frozenset(face)].append((face, face_type))

    patches: dict[str, list[tuple[int, ...]]] = defaultdict(list)
    for records in face_occurrences.values():
        if len(records) != 1:
            continue
        face, face_type = records[0]
        coordinates = [vertices[index] for index in face]
        if face_type == "bottom":
            patch = "bottom"
        elif face_type == "top":
            patch = "membrane"
        elif all(abs(point[0] + 600 * UM) < 1e-12 for point in coordinates):
            patch = "oilInlet"
        elif all(abs(point[1] - 650 * UM) < 1e-12 for point in coordinates):
            patch = "waterInlet"
        elif all(abs(point[0] - 1900 * UM) < 1e-12 for point in coordinates):
            patch = "outlet"
        else:
            patch = "walls"
        patches[patch].append(face)

    lines = [
        "FoamFile { format ascii; class dictionary; object blockMeshDict; }",
        "convertToMeters 1;",
        "",
        "vertices",
        "(",
    ]
    lines += [f"    ({x:.12g} {y:.12g} {z:.12g})" for x, y, z in vertices]
    lines += [");", "", "blocks", "("]
    lines += [
        f"    hex ({' '.join(map(str, ids))}) ({nx} {ny} {depth}) simpleGrading (1 1 1)"
        for ids, nx, ny, depth in block_rows
    ]
    lines += [");", "", "edges ();", "", "boundary", "("]
    patch_types = {"oilInlet": "patch", "waterInlet": "patch", "outlet": "patch"}
    for name in ("oilInlet", "waterInlet", "outlet", "walls", "bottom", "membrane"):
        lines += [f"    {name}", "    {", f"        type {patch_types.get(name, 'wall')};", "        faces", "        ("]
        lines += [f"            ({' '.join(map(str, face))})" for face in patches[name]]
        lines += ["        );", "    }"]
    lines += [");", "", "mergePatchPairs ();", ""]

    avg = average_deflection(kp)
    metadata: dict[str, float | int | str] = {
        "kp": kp,
        "mesh_level": mesh_level,
        "estimated_cells": cell_count,
        "junction_width_m": WIDTH,
        "undeformed_height_m": BASE_HEIGHT,
        "region_b_width_m": REGION_B_WIDTH,
        "region_b_start_m": 600 * UM,
        "measurement_line_m": 1200 * UM,
        "average_junction_deflection_m": avg,
        "maximum_junction_deflection_m": 1.5 * avg,
        "wall_model": "static analytical parabolic profile; not two-way FSI",
        "region_b_source": "width inferred from Figure 6 image scale; axial lengths assumed",
    }
    return "\n".join(lines), metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=("rigid", "flexible", "all"), default="all")
    parser.add_argument("--mesh", choices=("screening", "medium", "fine"), default="screening")
    args = parser.parse_args()
    selected = ("rigid", "flexible") if args.case == "all" else (args.case,)
    for case_name in selected:
        kp = 0.0 if case_name == "rigid" else 1.695
        block_mesh, metadata = generate_block_mesh(kp, args.mesh)
        case = HERE / case_name
        (case / "system" / "blockMeshDict").write_text(block_mesh)
        (case / "geometry.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"{case_name}: {metadata['estimated_cells']} cells")


if __name__ == "__main__":
    main()
