#!/usr/bin/env python3
"""Measure projected VOF silhouettes and paired 3-D phase volumes.

An alpha isosurface is only a provisional optical proxy. Threshold dependence,
phase separation, and mesh/time convergence must be checked before validation.

Run this script with the system Python because Ubuntu's VTK bindings are
installed there:

    /usr/bin/python3 measure.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import deque
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parent
UM = 1e-6
X_MIN = 700 * UM
X_MAX = 1900 * UM
Y_MIN = -120 * UM
Y_MAX = 120 * UM
PIXEL_SIZE = 1 * UM
MEASUREMENT_X = 1200 * UM
THRESHOLDS = (0.25, 0.5, 0.75)
PRIMARY_THRESHOLD = 0.5
MAX_TRACK_STEP = 250 * UM
MAX_MISSING_FRAMES = 2
MIN_DIAMETER_CELLS = 2.0
RESOLVED_DIAMETER_CELLS = 8.0

NX = round((X_MAX - X_MIN) / PIXEL_SIZE)
NY = round((Y_MAX - Y_MIN) / PIXEL_SIZE)


def _projected_pixel(point: tuple[float, float, float]) -> tuple[float, float]:
    return ((point[0] - X_MIN) / PIXEL_SIZE, (point[1] - Y_MIN) / PIXEL_SIZE)


def _draw_edge(mask: np.ndarray, p0: tuple[float, float], p1: tuple[float, float]) -> None:
    length = max(abs(p1[0] - p0[0]), abs(p1[1] - p0[1]))
    count = max(2, math.ceil(2 * length) + 1)
    us = np.linspace(p0[0], p1[0], count)
    vs = np.linspace(p0[1], p1[1], count)
    columns = np.floor(us).astype(int)
    rows = np.floor(vs).astype(int)
    valid = (columns >= 0) & (columns < NX) & (rows >= 0) & (rows < NY)
    mask[rows[valid], columns[valid]] = True


def _fill_triangle(mask: np.ndarray, points: list[tuple[float, float]]) -> None:
    min_column = max(0, math.floor(min(point[0] for point in points)))
    max_column = min(NX - 1, math.floor(max(point[0] for point in points)))
    min_row = max(0, math.floor(min(point[1] for point in points)))
    max_row = min(NY - 1, math.floor(max(point[1] for point in points)))
    if min_column > max_column or min_row > max_row:
        return

    x0, y0 = points[0]
    x1, y1 = points[1]
    x2, y2 = points[2]
    denominator = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    if abs(denominator) < 1e-10:
        return
    columns = np.arange(min_column, max_column + 1) + 0.5
    rows = np.arange(min_row, max_row + 1) + 0.5
    xx, yy = np.meshgrid(columns, rows)
    a = ((y1 - y2) * (xx - x2) + (x2 - x1) * (yy - y2)) / denominator
    b = ((y2 - y0) * (xx - x2) + (x0 - x2) * (yy - y2)) / denominator
    c = 1 - a - b
    inside = (a >= -1e-9) & (b >= -1e-9) & (c >= -1e-9)
    mask[min_row : max_row + 1, min_column : max_column + 1] |= inside


def _fill_closed_outlines(mask: np.ndarray) -> np.ndarray:
    """Fill projected closed interfaces using an exterior flood fill."""
    thick = mask.copy()
    for row_shift in (-1, 0, 1):
        for column_shift in (-1, 0, 1):
            if row_shift == column_shift == 0:
                continue
            source_rows = slice(max(0, -row_shift), min(NY, NY - row_shift))
            source_columns = slice(max(0, -column_shift), min(NX, NX - column_shift))
            target_rows = slice(max(0, row_shift), min(NY, NY + row_shift))
            target_columns = slice(max(0, column_shift), min(NX, NX + column_shift))
            thick[target_rows, target_columns] |= mask[source_rows, source_columns]

    outside = np.zeros_like(thick)
    queue: deque[tuple[int, int]] = deque()
    for column in range(NX):
        for row in (0, NY - 1):
            if not thick[row, column] and not outside[row, column]:
                outside[row, column] = True
                queue.append((row, column))
    for row in range(NY):
        for column in (0, NX - 1):
            if not thick[row, column] and not outside[row, column]:
                outside[row, column] = True
                queue.append((row, column))
    while queue:
        row, column = queue.popleft()
        for next_row, next_column in (
            (row - 1, column),
            (row + 1, column),
            (row, column - 1),
            (row, column + 1),
        ):
            if (
                0 <= next_row < NY
                and 0 <= next_column < NX
                and not thick[next_row, next_column]
                and not outside[next_row, next_column]
            ):
                outside[next_row, next_column] = True
                queue.append((next_row, next_column))
    filled = ~outside
    # The one-pixel dilation above only closes raster gaps. Remove that added
    # layer again so the measured contour remains centred on alpha=threshold.
    eroded = filled.copy()
    for row_shift in (-1, 0, 1):
        for column_shift in (-1, 0, 1):
            if row_shift == column_shift == 0:
                continue
            shifted = np.zeros_like(filled)
            source_rows = slice(max(0, -row_shift), min(NY, NY - row_shift))
            source_columns = slice(max(0, -column_shift), min(NX, NX - column_shift))
            target_rows = slice(max(0, row_shift), min(NY, NY + row_shift))
            target_columns = slice(max(0, column_shift), min(NX, NX + column_shift))
            shifted[target_rows, target_columns] = filled[source_rows, source_columns]
            eroded &= shifted
    return eroded


def rasterize_silhouette(surface: vtk.vtkPolyData) -> np.ndarray:
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputData(surface)
    triangles.Update()
    poly = triangles.GetOutput()
    mask = np.zeros((NY, NX), dtype=bool)
    identifiers = vtk.vtkIdList()
    polygons = poly.GetPolys()
    polygons.InitTraversal()
    while polygons.GetNextCell(identifiers):
        if identifiers.GetNumberOfIds() != 3:
            continue
        points = [
            _projected_pixel(poly.GetPoint(identifiers.GetId(index)))
            for index in range(3)
        ]
        _fill_triangle(mask, points)
        _draw_edge(mask, points[0], points[1])
        _draw_edge(mask, points[1], points[2])
        _draw_edge(mask, points[2], points[0])
    return _fill_closed_outlines(mask) if mask.any() else mask


def _convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    points = sorted(set(points))
    if len(points) <= 1:
        return points

    def cross(origin: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _maximum_feret(boundary: list[tuple[int, int]]) -> float:
    hull = _convex_hull(boundary)
    maximum_squared = 0.0
    for index, point in enumerate(hull):
        for other in hull[index + 1 :]:
            maximum_squared = max(
                maximum_squared,
                (point[0] - other[0]) ** 2 + (point[1] - other[1]) ** 2,
            )
    return math.sqrt(maximum_squared) * PIXEL_SIZE


def component_metrics(mask: np.ndarray, nominal_cell: float) -> list[dict[str, float]]:
    visited = np.zeros_like(mask)
    found: list[dict[str, float]] = []
    for seed_row, seed_column in np.argwhere(mask):
        if visited[seed_row, seed_column]:
            continue
        visited[seed_row, seed_column] = True
        queue = deque([(int(seed_row), int(seed_column))])
        pixels: list[tuple[int, int]] = []
        touches_edge = False
        while queue:
            row, column = queue.popleft()
            pixels.append((row, column))
            touches_edge |= row in (0, NY - 1) or column in (0, NX - 1)
            for next_row, next_column in (
                (row - 1, column),
                (row + 1, column),
                (row, column - 1),
                (row, column + 1),
            ):
                if (
                    0 <= next_row < NY
                    and 0 <= next_column < NX
                    and mask[next_row, next_column]
                    and not visited[next_row, next_column]
                ):
                    visited[next_row, next_column] = True
                    queue.append((next_row, next_column))
        if touches_edge:
            continue
        rows = np.fromiter((pixel[0] for pixel in pixels), dtype=int)
        columns = np.fromiter((pixel[1] for pixel in pixels), dtype=int)
        area = len(pixels) * PIXEL_SIZE**2
        area_diameter = math.sqrt(4 * area / math.pi)
        if area_diameter < MIN_DIAMETER_CELLS * nominal_cell:
            continue
        boundary: list[tuple[int, int]] = []
        perimeter_edges = 0
        pixel_set = set(pixels)
        for row, column in pixels:
            exposed = 0
            for neighbor in ((row - 1, column), (row + 1, column), (row, column - 1), (row, column + 1)):
                if neighbor not in pixel_set:
                    exposed += 1
            if exposed:
                boundary.append((column, row))
                perimeter_edges += exposed
        perimeter = perimeter_edges * PIXEL_SIZE
        axial = (columns.max() - columns.min() + 1) * PIXEL_SIZE
        transverse = (rows.max() - rows.min() + 1) * PIXEL_SIZE
        found.append(
            {
                "cx": X_MIN + (float(columns.mean()) + 0.5) * PIXEL_SIZE,
                "cy": Y_MIN + (float(rows.mean()) + 0.5) * PIXEL_SIZE,
                "axial_diameter": axial,
                "transverse_diameter": transverse,
                "area_equivalent_diameter": area_diameter,
                "maximum_feret": _maximum_feret(boundary),
                "circularity": 4 * math.pi * area / perimeter**2 if perimeter else 0.0,
                "pixel_count": float(len(pixels)),
                "cells_across": axial / nominal_cell,
            }
        )
    return sorted(found, key=lambda item: item["cx"])


def _internal_mesh(output: vtk.vtkDataObject) -> vtk.vtkDataSet:
    if isinstance(output, vtk.vtkDataSet):
        return output
    if not isinstance(output, vtk.vtkMultiBlockDataSet):
        raise RuntimeError(f"unsupported OpenFOAM reader output: {output.GetClassName()}")
    for index in range(output.GetNumberOfBlocks()):
        block = output.GetBlock(index)
        if block is None:
            continue
        metadata = output.GetMetaData(index)
        name = metadata.Get(vtk.vtkCompositeDataSet.NAME()) if metadata else None
        if name == "internalMesh" and isinstance(block, vtk.vtkDataSet):
            return block
        try:
            return _internal_mesh(block)
        except RuntimeError:
            pass
    raise RuntimeError("internalMesh was not found in the OpenFOAM case")


def _interface(mesh: vtk.vtkDataSet, threshold: float) -> vtk.vtkPolyData:
    if mesh.GetPointData().GetArray("alpha.water") is None:
        converter = vtk.vtkCellDataToPointData()
        converter.SetInputData(mesh)
        converter.PassCellDataOn()
        converter.Update()
        mesh = converter.GetOutput()
    contour = vtk.vtkContourFilter()
    contour.SetInputData(mesh)
    contour.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "alpha.water")
    contour.SetValue(0, threshold)
    contour.Update()
    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(contour.GetOutput())
    geometry.Update()
    return geometry.GetOutput()


def _volume_components(mesh: vtk.vtkDataSet, threshold: float) -> list[dict[str, float]]:
    """Connected 3-D water regions, retaining cellwise VOF phase volume.

    Connectivity is defined by cells with alpha >= threshold. The reported
    phase volume is sum(alpha * cell volume) in that component; it is therefore
    threshold-dependent and is not a substitute for a convergence study.
    """
    selected = vtk.vtkThreshold()
    selected.SetInputData(mesh)
    selected.SetInputArrayToProcess(
        0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "alpha.water"
    )
    selected.SetLowerThreshold(threshold)
    selected.SetUpperThreshold(1.0 + 1e-6)
    selected.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
    selected.Update()
    if not selected.GetOutput().GetNumberOfCells():
        return []

    connected = vtk.vtkConnectivityFilter()
    connected.SetInputConnection(selected.GetOutputPort())
    connected.SetExtractionModeToAllRegions()
    connected.ColorRegionsOn()
    connected.Update()
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputConnection(connected.GetOutputPort())
    sizes.ComputeVolumeOn()
    sizes.Update()
    volume_mesh = sizes.GetOutput()
    cells = volume_mesh.GetCellData()
    region_ids = vtk_to_numpy(cells.GetArray("RegionId")).astype(int)
    alpha = np.clip(vtk_to_numpy(cells.GetArray("alpha.water")), 0.0, 1.0)
    cell_volumes = vtk_to_numpy(cells.GetArray("Volume"))
    centers = vtk.vtkCellCenters()
    centers.SetInputData(volume_mesh)
    centers.Update()
    xyz = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
    found = []
    for region_id in np.unique(region_ids):
        member = region_ids == region_id
        points = xyz[member]
        weights = alpha[member] * cell_volumes[member]
        phase_volume = float(weights.sum())
        if phase_volume <= 0:
            continue
        centroid = np.average(points, axis=0, weights=weights)
        found.append({
            "cx": float(centroid[0]),
            "cy": float(centroid[1]),
            "volume": phase_volume,
            "cell_volume": float(cell_volumes[member].sum()),
            "peak_alpha": float(alpha[member].max()),
            "x_min": float(points[:, 0].min()),
            "x_max": float(points[:, 0].max()),
            "y_min": float(points[:, 1].min()),
            "y_max": float(points[:, 1].max()),
        })
    return found


def _pair_volumes(
    detections: list[dict[str, float]],
    regions: list[dict[str, float]],
    nominal_cell: float,
) -> None:
    """Attach one complete 3-D component to each projected silhouette."""
    candidates = []
    for detection_index, detection in enumerate(detections):
        detection["volume"] = math.nan
        detection["cell_volume"] = math.nan
        detection["peak_alpha"] = math.nan
        for region_index, region in enumerate(regions):
            # A component clipped by the optical window or still connected to
            # an inlet is not a complete downstream droplet.
            if region["x_min"] <= X_MIN + nominal_cell or region["x_max"] >= X_MAX - nominal_cell:
                continue
            dx = abs(region["cx"] - detection["cx"])
            dy = abs(region["cy"] - detection["cy"])
            if dx <= detection["axial_diameter"] / 2 + nominal_cell and dy <= detection["transverse_diameter"] / 2 + nominal_cell:
                candidates.append((math.hypot(dx, dy), detection_index, region_index))
    paired_detections: set[int] = set()
    paired_regions: set[int] = set()
    for _, detection_index, region_index in sorted(candidates):
        if detection_index in paired_detections or region_index in paired_regions:
            continue
        detections[detection_index]["volume"] = regions[region_index]["volume"]
        detections[detection_index]["cell_volume"] = regions[region_index]["cell_volume"]
        detections[detection_index]["peak_alpha"] = regions[region_index]["peak_alpha"]
        paired_detections.add(detection_index)
        paired_regions.add(region_index)


def _slab_phase_volumes(
    mesh: vtk.vtkDataSet,
    detections: list[dict[str, float]],
    nominal_cell: float,
) -> None:
    """Integrate *all* alpha in axial bins separated by droplet midgaps.

    This avoids discarding the diffuse interface. A bin is usable as a
    per-droplet volume only if both midgap cross-sections are nearly dry.
    """
    if not detections:
        return
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(mesh)
    sizes.ComputeVolumeOn()
    sizes.Update()
    sized = sizes.GetOutput()
    alpha = np.clip(vtk_to_numpy(sized.GetCellData().GetArray("alpha.water")), 0.0, 1.0)
    volumes = vtk_to_numpy(sized.GetCellData().GetArray("Volume"))
    centers = vtk.vtkCellCenters()
    centers.SetInputData(sized)
    centers.Update()
    xyz = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
    x, y = xyz[:, 0], xyz[:, 1]
    in_window = (y >= Y_MIN) & (y <= Y_MAX)
    ordered = sorted(detections, key=lambda item: item["cx"])
    for index, detection in enumerate(ordered):
        left = (ordered[index - 1]["cx"] + detection["cx"]) / 2 if index else X_MIN
        right = (detection["cx"] + ordered[index + 1]["cx"]) / 2 if index + 1 < len(ordered) else X_MAX
        slab = in_window & (x >= left) & (x < right)
        detection["slab_volume"] = float(np.sum(alpha[slab] * volumes[slab]))
        edge_peaks = [
            float(np.max(alpha[in_window & (abs(x - edge) <= nominal_cell)]))
            if np.any(in_window & (abs(x - edge) <= nominal_cell)) else math.nan
            for edge in (left, right)
        ]
        detection["slab_edge_peak_alpha"] = max(edge_peaks)
        inside = (
            detection["cx"] - detection["axial_diameter"] / 2 > X_MIN + nominal_cell
            and detection["cx"] + detection["axial_diameter"] / 2 < X_MAX - nominal_cell
        )
        detection["slab_separated"] = float(inside and max(edge_peaks) < 0.1)


def _time_values(reader: vtk.vtkOpenFOAMReader) -> list[float]:
    values = reader.GetTimeValues()
    return [values.GetValue(index) for index in range(values.GetNumberOfValues())]


def _geometry(case: Path) -> dict[str, object]:
    if (case / "geometry.json").exists():
        return json.loads((case / "geometry.json").read_text())
    metadata = json.loads((case / "case.json").read_text())
    return metadata.get("geometry", metadata)


def configure_window(
    x_min_um: float,
    x_max_um: float,
    y_min_um: float,
    y_max_um: float,
    measurement_x_um: float,
) -> None:
    """Set an optical field of view without changing the underlying mesh."""
    global X_MIN, X_MAX, Y_MIN, Y_MAX, MEASUREMENT_X, NX, NY
    if not (x_min_um < measurement_x_um < x_max_um and y_min_um < y_max_um):
        raise ValueError("measurement line must lie inside a nonempty field of view")
    X_MIN, X_MAX = x_min_um * UM, x_max_um * UM
    Y_MIN, Y_MAX = y_min_um * UM, y_max_um * UM
    MEASUREMENT_X = measurement_x_um * UM
    NX = round((X_MAX - X_MIN) / PIXEL_SIZE)
    NY = round((Y_MAX - Y_MIN) / PIXEL_SIZE)


def write_overlay(case_name: str, time_value: float, mask: np.ndarray) -> None:
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="300" viewBox="0 0 {NX} {NY+60}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="15" y="24" font-family="Arial,sans-serif" font-size="16">{case_name}: alpha.water = 0.5, t = {time_value:.5g} s</text>',
        f'<g transform="translate(0 40) scale(1 -1) translate(0 {-NY})">',
    ]
    for row in range(NY):
        columns = np.flatnonzero(mask[row])
        if not len(columns):
            continue
        start = previous = int(columns[0])
        for column in map(int, columns[1:]):
            if column != previous + 1:
                parts.append(f'<rect x="{start}" y="{row}" width="{previous-start+1}" height="1" fill="#2563eb"/>')
                start = column
            previous = column
        parts.append(f'<rect x="{start}" y="{row}" width="{previous-start+1}" height="1" fill="#2563eb"/>')
    measurement_column = (MEASUREMENT_X - X_MIN) / PIXEL_SIZE
    parts += [
        '</g>',
        f'<line x1="{measurement_column}" y1="40" x2="{measurement_column}" y2="{NY+40}" stroke="#dc2626" stroke-width="2"/>',
        f'<text x="{measurement_column+7}" y="{NY+57}" font-family="Arial,sans-serif" font-size="13">x = {MEASUREMENT_X/UM:g} um</text>',
        '</svg>',
    ]
    (ROOT / f"{case_name}_measurement_overlay.svg").write_text("\n".join(parts) + "\n")


def case_frames(
    case: Path, threshold: float
) -> tuple[list[float], list[list[dict[str, float]]], float, tuple[float, np.ndarray] | None]:
    geometry = _geometry(case)
    nominal_cell = (
        float(geometry["nominal_cell_um"]) * UM
        if "nominal_cell_um" in geometry
        else {
            "screening": 25 * UM,
            "medium": 12.5 * UM,
            "fine": 8.333333 * UM,
        }[str(geometry["mesh_level"])]
    )
    marker = case / f"{case.name}.foam"
    marker.touch()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(marker))
    reader.SkipZeroTimeOn()
    reader.CreateCellToPointOn()
    reader.UpdateInformation()
    reader.DisableAllCellArrays()
    reader.SetCellArrayStatus("alpha.water", 1)
    times = _time_values(reader)
    frames: list[list[dict[str, float]]] = []
    diagnostic: tuple[float, np.ndarray] | None = None
    diagnostic_distance = math.inf
    for time_value in times:
        reader.UpdateTimeStep(time_value)
        mesh = _internal_mesh(reader.GetOutput())
        mask = rasterize_silhouette(_interface(mesh, threshold))
        detections = component_metrics(mask, nominal_cell)
        _pair_volumes(detections, _volume_components(mesh, threshold), nominal_cell)
        _slab_phase_volumes(mesh, detections, nominal_cell)
        frames.append(detections)
        if math.isclose(threshold, PRIMARY_THRESHOLD):
            for detection in detections:
                distance = abs(detection["cx"] - MEASUREMENT_X)
                if distance < diagnostic_distance:
                    diagnostic_distance = distance
                    diagnostic = (time_value, mask.copy())
    return times, frames, nominal_cell, diagnostic


def track_crossings(times: list[float], frames: list[list[dict[str, float]]]) -> list[dict[str, float]]:
    tracks: list[dict[str, object]] = []
    measurements: list[dict[str, float]] = []
    metric_names = (
        "cx",
        "cy",
        "axial_diameter",
        "transverse_diameter",
        "area_equivalent_diameter",
        "maximum_feret",
        "circularity",
        "pixel_count",
        "cells_across",
        "volume",
        "cell_volume",
        "peak_alpha",
        "slab_volume",
        "slab_edge_peak_alpha",
        "slab_separated",
    )
    for frame_index, (time_value, detections) in enumerate(zip(times, frames)):
        active = [
            index
            for index, track in enumerate(tracks)
            if frame_index - int(track["last_frame"]) <= MAX_MISSING_FRAMES
        ]
        candidates: list[tuple[float, int, int]] = []
        for detection_index, detection in enumerate(detections):
            for track_index in active:
                previous = tracks[track_index]["last"]
                assert isinstance(previous, dict)
                backward = float(previous["cx"]) - detection["cx"]
                distance = math.hypot(detection["cx"] - float(previous["cx"]), detection["cy"] - float(previous["cy"]))
                if backward <= 50 * UM and distance <= MAX_TRACK_STEP:
                    candidates.append((distance, track_index, detection_index))
        matched_tracks: set[int] = set()
        matched_detections: set[int] = set()
        for _, track_index, detection_index in sorted(candidates):
            if track_index in matched_tracks or detection_index in matched_detections:
                continue
            track = tracks[track_index]
            previous = track["last"]
            assert isinstance(previous, dict)
            current = detections[detection_index]
            if not bool(track["measured"]) and previous["cx"] < MEASUREMENT_X <= current["cx"]:
                fraction = (MEASUREMENT_X - previous["cx"]) / (current["cx"] - previous["cx"])
                result = {
                    name: float(previous[name]) + fraction * (float(current[name]) - float(previous[name]))
                    for name in metric_names if name not in (
                        "volume", "cell_volume", "peak_alpha", "slab_volume",
                        "slab_edge_peak_alpha", "slab_separated",
                    )
                }
                volume_quality = "interpolated"
                for name in ("volume", "cell_volume"):
                    before, after = float(previous[name]), float(current[name])
                    if math.isfinite(before) and math.isfinite(after):
                        result[name] = before + fraction * (after - before)
                    elif math.isfinite(before) or math.isfinite(after):
                        result[name] = before if math.isfinite(before) else after
                        volume_quality = "one_sided"
                    else:
                        result[name] = math.nan
                        volume_quality = "unmatched"
                result["volume_pair_quality"] = volume_quality
                peaks = [float(value["peak_alpha"]) for value in (previous, current)]
                result["peak_alpha"] = min(value for value in peaks if math.isfinite(value)) if any(map(math.isfinite, peaks)) else math.nan
                result["slab_volume"] = (
                    float(previous["slab_volume"]) + fraction *
                    (float(current["slab_volume"]) - float(previous["slab_volume"]))
                )
                result["slab_edge_peak_alpha"] = max(
                    float(previous["slab_edge_peak_alpha"]),
                    float(current["slab_edge_peak_alpha"]),
                )
                result["slab_separated"] = min(
                    float(previous["slab_separated"]),
                    float(current["slab_separated"]),
                )
                result["time"] = float(track["last_time"]) + fraction * (time_value - float(track["last_time"]))
                result["track_id"] = float(track_index + 1)
                measurements.append(result)
                track["measured"] = True
            track["last"] = current
            track["last_time"] = time_value
            track["last_frame"] = frame_index
            matched_tracks.add(track_index)
            matched_detections.add(detection_index)
        for detection_index, detection in enumerate(detections):
            if detection_index not in matched_detections:
                tracks.append(
                    {
                        "last": detection,
                        "last_time": time_value,
                        "last_frame": frame_index,
                        "measured": detection["cx"] >= MEASUREMENT_X,
                    }
                )
    return measurements


def measure_cases(
    cases: list[tuple[str, Path]],
    output: Path,
    write_overlays: bool,
    thresholds: tuple[float, ...] = THRESHOLDS,
) -> None:
    rows: list[dict[str, object]] = []
    crossing_counts: dict[str, dict[str, int]] = {}
    time_coverage: dict[str, dict[str, object]] = {}
    for case_name, case in cases:
        crossing_counts[case_name] = {}
        for threshold in thresholds:
            print(f"{case_name}: reconstructing alpha.water={threshold:g} silhouettes", flush=True)
            times, frames, nominal_cell, diagnostic = case_frames(case, threshold)
            if case_name not in time_coverage:
                control = (case / "system/controlDict").read_text()
                match = re.search(r"\bendTime\s+([0-9.eE+-]+)\s*;", control)
                target_time = float(match.group(1)) if match else None
                last_time = max(times) if times else None
                time_coverage[case_name] = {
                    "first_output_time_s": min(times) if times else None,
                    "last_output_time_s": last_time,
                    "requested_end_time_s": target_time,
                    "reached_requested_end_time": (
                        last_time >= target_time - 1e-8 if last_time is not None and target_time is not None else None
                    ),
                }
            measurements = track_crossings(times, frames)
            measurements.sort(key=lambda item: item["time"])
            for previous, item in zip(measurements, measurements[1:]):
                item["preceding_period_s"] = item["time"] - previous["time"]
            crossing_counts[case_name][f"{threshold:g}"] = len(measurements)
            if write_overlays and diagnostic is not None:
                write_overlay(case_name, diagnostic[0], diagnostic[1])
            for item in measurements:
                rows.append(
                    {
                        "case": case_name,
                        "alpha_threshold": threshold,
                        "track_id": round(item["track_id"]),
                        "crossing_time_s": item["time"],
                        "preceding_period_s": item.get("preceding_period_s", ""),
                        "instantaneous_frequency_hz": 1 / item["preceding_period_s"] if item.get("preceding_period_s", 0) else "",
                        "centroid_x_um": item["cx"] / UM,
                        "centroid_y_um": item["cy"] / UM,
                        "axial_diameter_um": item["axial_diameter"] / UM,
                        "transverse_diameter_um": item["transverse_diameter"] / UM,
                        "maximum_feret_um": item["maximum_feret"] / UM,
                        "area_equivalent_diameter_um": item["area_equivalent_diameter"] / UM,
                        "circularity": item["circularity"],
                        "silhouette_pixels": round(item["pixel_count"]),
                        "nominal_cells_across_axial_diameter": item["cells_across"],
                        "resolution_flag": "resolved" if item["cells_across"] >= RESOLVED_DIAMETER_CELLS else "under_resolved",
                        "nominal_cell_um": nominal_cell / UM,
                        "vof_phase_volume_nl": item["volume"] / 1e-12 if math.isfinite(item["volume"]) else "",
                        "threshold_cell_volume_nl": item["cell_volume"] / 1e-12 if math.isfinite(item["cell_volume"]) else "",
                        "volume_pair_quality": item["volume_pair_quality"],
                        "peak_alpha": item["peak_alpha"] if math.isfinite(item["peak_alpha"]) else "",
                        "full_vof_slab_volume_nl": item["slab_volume"] / 1e-12,
                        "slab_edge_peak_alpha": item["slab_edge_peak_alpha"],
                        "slab_separated": bool(item["slab_separated"]),
                    }
                )
            print(f"{case_name}: threshold {threshold:g}, {len(measurements)} line crossings", flush=True)
    fields = list(rows[0]) if rows else [
        "case",
        "alpha_threshold",
        "track_id",
        "crossing_time_s",
        "preceding_period_s",
        "instantaneous_frequency_hz",
        "centroid_x_um",
        "centroid_y_um",
        "axial_diameter_um",
        "transverse_diameter_um",
        "maximum_feret_um",
        "area_equivalent_diameter_um",
        "circularity",
        "silhouette_pixels",
        "nominal_cells_across_axial_diameter",
        "resolution_flag",
        "nominal_cell_um",
        "vof_phase_volume_nl",
        "threshold_cell_volume_nl",
        "volume_pair_quality",
        "peak_alpha",
        "full_vof_slab_volume_nl",
        "slab_edge_peak_alpha",
        "slab_separated",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    output.with_suffix(".json").write_text(json.dumps({
        "alpha_thresholds": list(thresholds),
        "crossing_counts": crossing_counts,
        "time_coverage": time_coverage,
        "measurement_x_um": MEASUREMENT_X / UM,
        "window_um": [X_MIN / UM, X_MAX / UM, Y_MIN / UM, Y_MAX / UM],
        "volume_definition": {
            "full_vof_slab_volume_nl": "sum(alpha.water * cell_volume) over all cells between adjacent projected-droplet centroid midgaps",
            "vof_phase_volume_nl": "sum(alpha.water * cell_volume) in connected cells with alpha.water >= alpha_threshold; diagnostic core only",
            "slab_separated": "both bin-edge cross sections peak below alpha.water=0.1 and droplet is wholly inside the view",
        },
    }, indent=2) + "\n")


def _synthetic_surface(radius: float, scale_x: float, scale_y: float) -> vtk.vtkPolyData:
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(160)
    sphere.SetPhiResolution(80)
    transform = vtk.vtkTransform()
    transform.Translate(1200 * UM, 0, 0)
    transform.Scale(scale_x, scale_y, 0.5)
    transformed = vtk.vtkTransformPolyDataFilter()
    transformed.SetInputConnection(sphere.GetOutputPort())
    transformed.SetTransform(transform)
    transformed.Update()
    return transformed.GetOutput()


def self_test() -> None:
    tests = (
        ("circle", _synthetic_surface(100 * UM, 1.0, 1.0), 200.0, 200.0, 200.0),
        ("ellipse", _synthetic_surface(100 * UM, 1.25, 0.75), 250.0, 150.0, math.sqrt(250 * 150)),
    )
    for name, surface, expected_x, expected_y, expected_area in tests:
        values = component_metrics(rasterize_silhouette(surface), 10 * UM)
        if len(values) != 1:
            raise RuntimeError(f"{name}: expected one silhouette, found {len(values)}")
        value = values[0]
        actual = (
            value["axial_diameter"] / UM,
            value["transverse_diameter"] / UM,
            value["area_equivalent_diameter"] / UM,
        )
        expected = (expected_x, expected_y, expected_area)
        errors = [100 * abs(observed - target) / target for observed, target in zip(actual, expected)]
        print(f"{name}: axial={actual[0]:.2f} um, transverse={actual[1]:.2f} um, area-equivalent={actual[2]:.2f} um")
        if max(errors) > 2.0:
            raise RuntimeError(f"{name}: raster measurement error exceeds 2%: {errors}")

    # Two separated VOF cells must retain their individual phase volumes;
    # simply turning either cell into a geometric full-volume droplet is wrong.
    image = vtk.vtkImageData()
    image.SetDimensions(5, 2, 2)
    image.SetSpacing(UM, UM, UM)
    alpha = vtk.vtkDoubleArray()
    alpha.SetName("alpha.water")
    for value in (0.8, 0.0, 0.6, 0.0):
        alpha.InsertNextValue(value)
    image.GetCellData().AddArray(alpha)
    components = sorted(_volume_components(image, 0.5), key=lambda item: item["cx"])
    if len(components) != 2 or any(
        not math.isclose(component["volume"], expected * UM**3, rel_tol=1e-6)
        for component, expected in zip(components, (0.8, 0.6))
    ):
        raise RuntimeError(f"VOF component volume test failed: {components}")
    configure_window(0, 4, 0, 1, 1)
    detections = [
        {"cx": 0.5 * UM, "axial_diameter": UM},
        {"cx": 2.5 * UM, "axial_diameter": UM},
    ]
    _slab_phase_volumes(image, detections, UM / 4)
    if any(
        not math.isclose(detection["slab_volume"], expected * UM**3, rel_tol=1e-6)
        for detection, expected in zip(detections, (0.8, 0.6))
    ):
        raise RuntimeError(f"full-VOF slab volume test failed: {detections}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--case", type=Path)
    parser.add_argument("--case-name")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--x-min-um", type=float, default=X_MIN / UM)
    parser.add_argument("--x-max-um", type=float, default=X_MAX / UM)
    parser.add_argument("--y-min-um", type=float, default=Y_MIN / UM)
    parser.add_argument("--y-max-um", type=float, default=Y_MAX / UM)
    parser.add_argument("--measurement-x-um", type=float, default=MEASUREMENT_X / UM)
    parser.add_argument(
        "--threshold",
        type=float,
        action="append",
        help="measure only the requested alpha threshold; may be repeated",
    )
    args = parser.parse_args()
    configure_window(
        args.x_min_um, args.x_max_um, args.y_min_um, args.y_max_um,
        args.measurement_x_um,
    )
    if args.self_test:
        self_test()
    elif args.case:
        case = args.case.resolve()
        measure_cases(
            [(args.case_name or case.name, case)],
            (args.output or case / "droplet_measurements.csv").resolve(),
            False,
            tuple(args.threshold) if args.threshold else THRESHOLDS,
        )
    else:
        measure_cases(
            [(name, ROOT / name) for name in ("rigid", "flexible")],
            ROOT / "droplet_measurements.csv",
            True,
        )


if __name__ == "__main__":
    main()
