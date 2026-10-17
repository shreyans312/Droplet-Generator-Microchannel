#!/usr/bin/env python3
"""Stationary 3-D Laplace-drop benchmark using the generator's two fluids.

No experimental droplet size is used to tune this verification case.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RADIUS = 40e-6
SIGMA = .00345
LAPLACE = 2 * SIGMA / RADIUS


def generate(smoothed: bool = False, sharpened: bool = False, geometric: bool = False) -> None:
    from generate import GEOMETRIES, fields, foam_header
    cases = (("coarse_mplic_smoothed", 20, 2.8e-6), ("fine_mplic_smoothed", 40, .99e-6)) if geometric else (("coarse_sharpened", 20, 2.8e-6), ("fine_sharpened", 40, .99e-6)) if sharpened else (("coarse_smoothed", 20, 2.8e-6), ("fine_smoothed", 40, .99e-6)) if smoothed else (("coarse_uncapped", 20, 1e-4), ("coarse_capped", 20, 2.8e-6), ("fine_capped", 40, .99e-6))
    for name, cells, maximum in cases:
        case = HERE / "static-drop" / name
        if case.exists():
            raise FileExistsError(case)
        content = fields(GEOMETRIES["unstretched"], .003, SIGMA, -200, False)
        # Keep only material properties and numerics from the production case.
        content = {k: v for k, v in content.items() if k.startswith("constant/")}
        content["system/blockMeshDict"] = foam_header("blockMeshDict") + f"""convertToMeters 1e-6;
vertices ((-80 -80 -80) (80 -80 -80) (80 80 -80) (-80 80 -80)
          (-80 -80 80) (80 -80 80) (80 80 80) (-80 80 80));
blocks (hex (0 1 2 3 4 5 6 7) ({cells} {cells} {cells}) simpleGrading (1 1 1));
edges ();
boundary (walls {{ type wall; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }});
mergePatchPairs ();
"""
        for field, cls, dim, initial, wall in (
            ("alpha.water", "volScalarField", "[]", "0", "zeroGradient"),
            ("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "0", "zeroGradient"),
            ("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", "fixedValue; value uniform (0 0 0)"),
        ):
            content[f"0/{field}"] = foam_header(field, cls) + f"dimensions {dim};\ninternalField uniform {initial};\nboundaryField {{ walls {{ type {wall}; }} }}\n"
        content["system/controlDict"] = foam_header("controlDict") + f"""solver incompressibleVoF;
startFrom startTime; startTime 0; stopAt endTime; endTime .003;
deltaT 1e-6; adjustTimeStep yes; maxCo .2; maxAlphaCo .1; maxDeltaT {maximum:.12g};
writeControl adjustableRunTime; writeInterval .00025; purgeWrite 0;
writeFormat binary; writeCompression off; timePrecision 10; runTimeModifiable true;
functions {{}}
"""
        baseline = HERE / "runs-focused/unstretched_robust/system"
        for f in ("fvSchemes", "fvSolution"):
            content[f"system/{f}"] = (baseline / f).read_text()
        if smoothed:
            library = HERE / "lib/libcurvatureSmooth.so"
            content["system/controlDict"] += f'libs ("{library}");\n'
            content["system/fvSchemes"] = content["system/fvSchemes"].replace("gradSchemes { default Gauss linear; }", "gradSchemes { default Gauss linear; nHat curvatureSmooth 2 0.5; }")
        if sharpened:
            content["system/controlDict"] = content["system/controlDict"].replace("solver incompressibleVoF;", "solver incompressibleVoFCurvature;")
            content["system/controlDict"] += "surfaceTensionSharpening 0.5;\n"
        if geometric:
            original = "Gauss interfaceCompression vanLeer 1"
            if original not in content["system/fvSchemes"]:
                raise RuntimeError("baseline advection scheme not found")
            content["system/fvSchemes"] = content["system/fvSchemes"].replace(original, "Gauss MPLIC")
        content["system/fvSolution"] = content["system/fvSolution"].replace("momentumPredictor yes;", "momentumPredictor yes;\n    pRefCell 0; pRefValue 0;")
        content["system/decomposeParDict"] = foam_header("decomposeParDict") + "numberOfSubdomains 8; method scotch;\n"
        content["system/setFieldsDict"] = foam_header("setFieldsDict") + f"""defaultValues {{ alpha.water 0; p_rgh 0; }}
zones {{ drop {{ type sphere; centre (0 0 0); radius {RADIUS}; values {{ alpha.water 1; p_rgh {LAPLACE}; }} }} }}
"""
        for relative, data in content.items():
            path = case / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(data)
        (case / "benchmark.json").write_text(json.dumps({"radius_m": RADIUS, "sigma_n_m": SIGMA, "laplace_pressure_pa": LAPLACE, "edge_m": 160e-6/cells, "max_delta_t_s": maximum, "curvature": "smoothed_2_passes_weight_0.5" if smoothed else "stock_csf", "force_colour_sharpening": 0.5 if sharpened else 0.0, "advection": "MPLIC" if geometric else "compression"}, indent=2) + "\n")
        print(case)


def analyse() -> None:
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    sys.path.insert(0, str(ROOT / "cases/flexible-rigid"))
    from measure import _internal_mesh
    results = []
    for case in sorted((HERE / "static-drop").iterdir()):
        if not case.is_dir():
            continue
        metadata = json.loads((case / "benchmark.json").read_text())
        if not (case / "0.003/U").exists():
            results.append({"case": case.name, "completed": False})
            continue
        marker = case / f"{case.name}.foam"
        marker.touch()
        reader = vtk.vtkOpenFOAMReader()
        reader.SetFileName(str(marker))
        reader.SkipZeroTimeOff()
        reader.UpdateInformation()
        tvalues = reader.GetTimeValues()
        times = [tvalues.GetValue(i) for i in range(tvalues.GetNumberOfValues())]
        if max(times) < .003 - 1e-10:
            results.append({"case": case.name, "completed": False, "last_reconstructed_time_s": max(times)})
            continue
        for t in (0, max(times)):
            reader.UpdateTimeStep(t)
            mesh = _internal_mesh(reader.GetOutput())
            size = vtk.vtkCellSizeFilter()
            size.SetInputData(mesh)
            size.ComputeVolumeOn()
            size.Update()
            volume = vtk_to_numpy(size.GetOutput().GetCellData().GetArray("Volume"))
            alpha = vtk_to_numpy(mesh.GetCellData().GetArray("alpha.water"))
            p = vtk_to_numpy(mesh.GetCellData().GetArray("p_rgh"))
            u = vtk_to_numpy(mesh.GetCellData().GetArray("U"))
            core, outside = alpha > .99, alpha < .01
            delta_p = float(np.average(p[core], weights=volume[core]) - np.average(p[outside], weights=volume[outside]))
            speed = float(np.linalg.norm(u, axis=1).max())
            phase_volume = float(alpha @ volume)
            equivalent_radius = (3 * phase_volume / (4 * math.pi))**(1/3)
            results.append({"case": case.name, "completed": True, "time_s": t, "phase_volume_m3": phase_volume, "volume_equivalent_radius_m": equivalent_radius, "pressure_jump_pa": delta_p, "laplace_nominal_pa": LAPLACE, "pressure_error_percent": 100*(delta_p/LAPLACE-1), "pressure_error_relative_volume_radius_percent": 100*(delta_p/(2*SIGMA/equivalent_radius)-1), "maximum_speed_m_s": speed, "parasitic_capillary_number": .042*speed/SIGMA, **metadata})
    (HERE / "static-drop/results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("generate", "generate-smoothed", "generate-sharpened", "generate-mplic", "analyse"))
    a = p.parse_args()
    if a.mode.startswith("generate"):
        generate(a.mode != "generate", a.mode == "generate-sharpened", a.mode == "generate-mplic")
    else:
        analyse()
