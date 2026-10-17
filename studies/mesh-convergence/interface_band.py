#!/usr/bin/env python3
"""Measure the mixed-cell VOF band in region B at droplet crossing times."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import vtk


X_MIN = 700e-6
X_MAX = 1900e-6
LOW_ALPHA = 0.01
HIGH_ALPHA = 0.99


def internal_mesh(output: vtk.vtkDataObject) -> vtk.vtkDataSet:
    if isinstance(output, vtk.vtkDataSet):
        return output
    if not isinstance(output, vtk.vtkMultiBlockDataSet):
        raise RuntimeError(f"unsupported reader output: {output.GetClassName()}")
    for index in range(output.GetNumberOfBlocks()):
        block = output.GetBlock(index)
        if block is None:
            continue
        metadata = output.GetMetaData(index)
        name = metadata.Get(vtk.vtkCompositeDataSet.NAME()) if metadata else None
        if name == "internalMesh" and isinstance(block, vtk.vtkDataSet):
            return block
        try:
            return internal_mesh(block)
        except RuntimeError:
            pass
    raise RuntimeError("internalMesh was not found")


def interface_area(mesh: vtk.vtkDataSet) -> float:
    converted = vtk.vtkCellDataToPointData()
    converted.SetInputData(mesh)
    converted.PassCellDataOn()
    converted.Update()
    contour = vtk.vtkContourFilter()
    contour.SetInputData(converted.GetOutput())
    contour.SetInputArrayToProcess(
        0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "alpha.water"
    )
    contour.SetValue(0, 0.5)
    contour.Update()
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(contour.GetOutputPort())
    triangles.Update()
    poly = triangles.GetOutput()
    ids = vtk.vtkIdList()
    polygons = poly.GetPolys()
    polygons.InitTraversal()
    area = 0.0
    while polygons.GetNextCell(ids):
        if ids.GetNumberOfIds() != 3:
            continue
        points = np.asarray([poly.GetPoint(ids.GetId(i)) for i in range(3)])
        centroid_x = float(points[:, 0].mean())
        if X_MIN <= centroid_x <= X_MAX:
            area += 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))
    return area


def band_metrics(mesh: vtk.vtkDataSet, nominal_cell: float) -> dict[str, float | int]:
    alpha_array = mesh.GetCellData().GetArray("alpha.water")
    if alpha_array is None:
        raise RuntimeError("cell alpha.water was not loaded")
    alpha = np.asarray([alpha_array.GetTuple1(i) for i in range(mesh.GetNumberOfCells())])

    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(mesh)
    sizes.SetComputeLength(False)
    sizes.SetComputeArea(False)
    sizes.SetComputeVolume(True)
    sizes.Update()
    sized = sizes.GetOutput()
    volume_array = sized.GetCellData().GetArray("Volume")
    volumes = np.asarray([volume_array.GetTuple1(i) for i in range(sized.GetNumberOfCells())])

    centre_filter = vtk.vtkCellCenters()
    centre_filter.SetInputData(mesh)
    centre_filter.Update()
    centres = centre_filter.GetOutput().GetPoints()
    x = np.asarray([centres.GetPoint(i)[0] for i in range(centres.GetNumberOfPoints())])
    region = (x >= X_MIN) & (x <= X_MAX)
    mixed = region & (alpha > LOW_ALPHA) & (alpha < HIGH_ALPHA)
    core_mixed = region & (alpha > 0.1) & (alpha < 0.9)
    area = interface_area(mesh)
    mixed_volume = float(volumes[mixed].sum())
    thickness = mixed_volume / area if area > 0 else math.nan
    return {
        "region_b_cells": int(region.sum()),
        "mixed_cells_0.01_0.99": int(mixed.sum()),
        "mixed_cells_0.1_0.9": int(core_mixed.sum()),
        "mixed_cell_volume_m3": mixed_volume,
        "alpha_0.5_area_m2": area,
        "effective_band_um": thickness * 1e6,
        "effective_band_cells": thickness / nominal_cell,
        "region_alpha_min": float(alpha[region].min()),
        "region_alpha_max": float(alpha[region].max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--mesh-level", choices=("screening", "medium", "fine"), required=True)
    parser.add_argument("--measurements", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.measurements.open(newline="") as handle:
        measurements = [
            row
            for row in csv.DictReader(handle)
            if row["case"] == args.case_name and math.isclose(float(row["alpha_threshold"]), 0.5)
        ]
    nominal = {"screening": 25e-6, "medium": 12.5e-6, "fine": 8.333333e-6}[args.mesh_level]
    marker = args.case / f"{args.case.name}.foam"
    marker_existed = marker.exists()
    marker.touch()
    try:
        reader = vtk.vtkOpenFOAMReader()
        reader.SetFileName(str(marker))
        reader.SkipZeroTimeOn()
        reader.CreateCellToPointOff()
        reader.UpdateInformation()
        reader.DisableAllCellArrays()
        reader.SetCellArrayStatus("alpha.water", 1)
        available = [reader.GetTimeValues().GetValue(i) for i in range(reader.GetTimeValues().GetNumberOfValues())]
        requests: list[dict[str, str]] = measurements or [
            {"track_id": "", "crossing_time_s": str(max(available))}
        ]
        rows: list[dict[str, object]] = []
        cache: dict[float, dict[str, float | int]] = {}
        for measurement in requests:
            crossing = float(measurement["crossing_time_s"])
            sampled = min(available, key=lambda value: abs(value - crossing))
            if sampled not in cache:
                reader.UpdateTimeStep(sampled)
                cache[sampled] = band_metrics(internal_mesh(reader.GetOutput()), nominal)
            rows.append(
                {
                    "case": args.case_name,
                    "mesh_level": args.mesh_level,
                    "track_id": measurement["track_id"],
                    "sample_reason": "crossing" if measurements else "latest_time_no_crossing",
                    "crossing_time_s": crossing,
                    "sampled_time_s": sampled,
                    **cache[sampled],
                }
            )
    finally:
        if not marker_existed:
            marker.unlink(missing_ok=True)

    fields = [
        "case", "mesh_level", "track_id", "sample_reason", "crossing_time_s", "sampled_time_s",
        "region_b_cells", "mixed_cells_0.01_0.99", "mixed_cells_0.1_0.9",
        "mixed_cell_volume_m3", "alpha_0.5_area_m2", "effective_band_um",
        "effective_band_cells", "region_alpha_min", "region_alpha_max",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{args.case_name} {args.mesh_level}: {len(rows)} interface samples")


if __name__ == "__main__":
    main()
