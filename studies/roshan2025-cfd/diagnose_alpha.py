#!/usr/bin/env python3
"""Inspect outlet phase-field resolution without assuming a droplet shape."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


MEASURE = Path(__file__).resolve().parents[2] / "cases/flexible-rigid"
sys.path.insert(0, str(MEASURE))
from measure import _internal_mesh  # noqa: E402


def inspect(case: Path, target_time: float) -> None:
    marker = case / f"{case.name}.foam"
    marker.touch()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(marker))
    reader.SkipZeroTimeOff()
    reader.CreateCellToPointOn()
    reader.UpdateInformation()
    reader.DisableAllCellArrays()
    reader.SetCellArrayStatus("alpha.water", 1)
    reader.SetCellArrayStatus("U", 1)
    values = reader.GetTimeValues()
    times = [values.GetValue(i) for i in range(values.GetNumberOfValues())]
    if not times:
        raise RuntimeError("no saved times")
    time = min(times, key=lambda value: abs(value - target_time))
    reader.UpdateTimeStep(time)
    mesh = _internal_mesh(reader.GetOutput())
    alpha_array = mesh.GetCellData().GetArray("alpha.water")
    if alpha_array is None:
        raise RuntimeError("cell alpha.water field missing")
    alpha = vtk_to_numpy(alpha_array)
    velocity_array = mesh.GetCellData().GetArray("U")
    velocity = vtk_to_numpy(velocity_array) if velocity_array is not None else None
    if velocity_array is not None:
        speeds = np.linalg.norm(velocity, axis=1)
        print(f"speed range: {speeds.min():.6g} to {speeds.max():.6g} m/s")
    centers = vtk.vtkCellCenters()
    centers.SetInputData(mesh)
    centers.Update()
    xyz = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
    outlet = (xyz[:, 0] > 450e-6) & (xyz[:, 0] < 990e-6)
    sample = alpha[outlet]
    print(f"{case.name} t={time:.6g}s outlet_cells={len(sample)}")
    print(f"alpha range: {sample.min():.6g} to {sample.max():.6g}")
    for label, lo, hi in (("dispersed inlet", -400, -80), ("junction", -80, 0), ("throat", 0, 412)):
        region = (xyz[:, 0] >= lo * 1e-6) & (xyz[:, 0] < hi * 1e-6)
        print(f"{label}: max={alpha[region].max():.6g}, alpha>0.9={np.count_nonzero(alpha[region] > 0.9)}")
    for threshold in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        print(f"  cells alpha >= {threshold:g}: {np.count_nonzero(sample >= threshold)}")
    for x_um in (500, 550, 600, 650, 750, 850, 950):
        station = np.abs(xyz[:, 0] - x_um * 1e-6) <= 4e-6
        subset = alpha[station]
        if len(subset):
            mean_ux = f", mean_Ux={velocity[station, 0].mean():.5g}m/s" if velocity is not None else ""
            core = station & (alpha > 0.5)
            if np.any(core):
                z_span = (xyz[core, 2].max() - xyz[core, 2].min()) * 1e6
                depth = f", core_z_center_span={z_span:.4g}um"
            else:
                depth = ""
            print(f"  x={x_um}um: max={subset.max():.5g}, alpha>0.5={np.count_nonzero(core)}{depth}{mean_ux}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--time", type=float, default=0.4)
    args = parser.parse_args()
    inspect(args.case.resolve(), args.time)
