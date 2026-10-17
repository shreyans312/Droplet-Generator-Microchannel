#!/usr/bin/env /usr/bin/python3
"""Independent 3-D VOF droplet-volume and connectivity diagnostic.

This does not infer droplet volume from a top-view diameter. It integrates
alpha.water times the actual CFD cell volume for each connected component.
Two alpha cutoffs are reported so that an interfacial bridge or diffuse tail
cannot silently masquerade as a well-defined isolated droplet.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


UM = 1e-6
PL = 1e-15
REGION_START = 700 * UM
REGION_END = 1900 * UM
MIN_REPORT_VOLUME = 20 * PL


def internal_mesh(case: Path, requested_time: float) -> tuple[vtk.vtkDataSet, float]:
    marker = case / f"{case.name}.foam"
    if not marker.exists():
        raise FileNotFoundError(f"missing OpenFOAM reader marker: {marker}")
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(marker))
    reader.SkipZeroTimeOn()
    reader.UpdateInformation()
    reader.DisableAllCellArrays()
    reader.SetCellArrayStatus("alpha.water", 1)
    available = reader.GetTimeValues()
    times = [available.GetValue(index) for index in range(available.GetNumberOfValues())]
    if not times:
        raise RuntimeError(f"no saved fields in {case}")
    actual_time = min(times, key=lambda value: abs(value - requested_time))
    reader.UpdateTimeStep(actual_time)
    output = reader.GetOutput()
    if output.GetNumberOfBlocks() != 1:
        raise RuntimeError(f"expected one internal-mesh block, found {output.GetNumberOfBlocks()}")
    mesh = output.GetBlock(0)
    if not isinstance(mesh, vtk.vtkDataSet) or mesh.GetCellData().GetArray("alpha.water") is None:
        raise RuntimeError("internal mesh or cell-centred alpha.water is unavailable")
    return mesh, actual_time


def components(mesh: vtk.vtkDataSet, cutoff: float) -> list[dict[str, object]]:
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(mesh)
    threshold.SetInputArrayToProcess(
        0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "alpha.water"
    )
    threshold.SetLowerThreshold(cutoff)
    threshold.SetUpperThreshold(1.0 + 1e-9)
    threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)

    connected = vtk.vtkConnectivityFilter()
    connected.SetInputConnection(threshold.GetOutputPort())
    connected.SetExtractionModeToAllRegions()
    connected.ColorRegionsOn()

    sized = vtk.vtkCellSizeFilter()
    sized.SetInputConnection(connected.GetOutputPort())
    sized.ComputeVolumeOn()
    sized.Update()
    output = sized.GetOutput()
    if output.GetNumberOfCells() == 0:
        return []

    center_filter = vtk.vtkCellCenters()
    center_filter.SetInputData(output)
    center_filter.Update()
    centers = vtk_to_numpy(center_filter.GetOutput().GetPoints().GetData())
    cell_data = output.GetCellData()
    labels = vtk_to_numpy(cell_data.GetArray("RegionId"))
    alpha = vtk_to_numpy(cell_data.GetArray("alpha.water"))
    volumes = vtk_to_numpy(cell_data.GetArray("Volume"))
    result = []
    for label in np.unique(labels):
        selected = labels == label
        xyz = centers[selected]
        water_volume = float(np.dot(alpha[selected], volumes[selected]))
        if water_volume < MIN_REPORT_VOLUME:
            continue
        weights = alpha[selected] * volumes[selected]
        x_centroid = float(np.dot(xyz[:, 0], weights) / water_volume)
        x_min, x_max = float(xyz[:, 0].min()), float(xyz[:, 0].max())
        y_max = float(xyz[:, 1].max())
        # The widened outlet itself reaches y=+100 um. Only the side-arm
        # excursion beyond it identifies a component still linked to inlet.
        attached_to_inlet = y_max > 100 * UM + 1e-9
        fully_in_region_b = (
            x_min > REGION_START + 1e-9
            and x_max < REGION_END - 1e-9
            and not attached_to_inlet
        )
        result.append(
            {
                "alpha_cutoff": cutoff,
                "component_id": int(label),
                "water_volume_pL": water_volume / PL,
                "x_centroid_um": x_centroid / UM,
                "x_min_um": x_min / UM,
                "x_max_um": x_max / UM,
                "y_max_um": y_max / UM,
                "cell_count": int(selected.sum()),
                "attached_to_inlet": attached_to_inlet,
                "fully_in_region_b": fully_in_region_b,
            }
        )
    return sorted(result, key=lambda item: float(item["x_centroid_um"]))


def audit(case: Path, time_value: float) -> list[dict[str, object]]:
    mesh, actual_time = internal_mesh(case, time_value)
    metadata = json.loads((case / "case.json").read_text())
    rows = []
    for cutoff in (0.01, 0.5):
        for item in components(mesh, cutoff):
            rows.append(
                {
                    "case": case.name,
                    "time_s": actual_time,
                    "flow_ratio": metadata["flow_ratio"],
                    "kp": metadata["kp"],
                    **item,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--time", type=float, default=1.5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = audit(args.case.resolve(), args.time)
    fields = [
        "case", "time_s", "flow_ratio", "kp", "alpha_cutoff", "component_id",
        "water_volume_pL", "x_centroid_um", "x_min_um", "x_max_um", "y_max_um",
        "cell_count", "attached_to_inlet", "fully_in_region_b",
    ]
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    else:
        writer = csv.DictWriter(__import__("sys").stdout, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
