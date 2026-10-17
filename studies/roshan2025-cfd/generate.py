#!/usr/bin/env python3
"""Generate a 3-D, three-inlet VOF flow-focuser for Roshan et al. (2025).

Paper DOI 10.1016/j.snb.2024.136617. Baseline dimensions are from section
3.1; 5-mm dimensions are approximate readings of experimental Fig. 3a.
The imposed deformed geometry represents the stabilised stretched device,
not a transient fluid-solid coupling calculation.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
UM = 1e-6
H = 137 * UM
X_INLET = -400 * UM
X_BRANCH_0 = -80 * UM
X_BRANCH_1 = -40 * UM
X_THROAT = 0.0
X_EXPAND = 412 * UM
X_OUTLET = 1012 * UM
SIDE_LENGTH = 200 * UM
Q_CONTINUOUS = 90e-9 / 3600  # 90 microlitres/hour, both side arms together
Q_DISPERSED = 30e-9 / 3600
RHO_DISPERSED = 1630.0
RHO_CONTINUOUS = 1060.0


@dataclass(frozen=True)
class Geometry:
    dispersed_width_um: float
    continuous_width_um: float
    throat_width_um: float
    outlet_width_um: float
    height_um: float = 137.0


GEOMETRIES = {
    "unstretched": Geometry(63.0, 40.0, 29.0, 90.0),
    "stretch5mm": Geometry(123.0, 40.0, 70.0, 153.0),
}


def footprint(geometry: Geometry) -> list[tuple[float, float, float, float]]:
    d = geometry.dispersed_width_um * UM / 2
    c = geometry.throat_width_um * UM / 2
    o = geometry.outlet_width_um * UM / 2
    if not (0 < c < d < o):
        raise ValueError("expected throat < dispersed inlet < outlet width")
    if not math.isclose((X_BRANCH_1 - X_BRANCH_0) / UM, geometry.continuous_width_um):
        raise ValueError("side-arm opening does not match continuous-channel width")
    blocks = []
    for x0, x1 in zip((X_INLET, X_BRANCH_0, X_BRANCH_1), (X_BRANCH_0, X_BRANCH_1, X_THROAT)):
        for y0, y1 in zip((-d, -c, 0, c), (-c, 0, c, d)):
            blocks.append((x0, x1, y0, y1))
    for y0, y1 in ((-d - SIDE_LENGTH, -d), (d, d + SIDE_LENGTH)):
        blocks.append((X_BRANCH_0, X_BRANCH_1, y0, y1))
    for y0, y1 in ((-c, 0), (0, c)):
        blocks.append((X_THROAT, X_EXPAND, y0, y1))
    for y0, y1 in zip((-o, -c, 0, c), (-c, 0, c, o)):
        blocks.append((X_EXPAND, X_OUTLET, y0, y1))
    return blocks


def mesh_dict(geometry: Geometry, level: str) -> tuple[str, int]:
    x_spacing, y_spacing, nz = {
        "screening": (8 * UM, 8 * UM, 6),
        "throat_refined": (8 * UM, 4 * UM, 6),
        "focused": (8 * UM, 4 * UM, 12),
        "medium": (4 * UM, 4 * UM, 12),
        "fine": (2 * UM, 2 * UM, 24),
    }[level]
    vertices = []
    vertex_ids = {}
    hexes = []
    faces = defaultdict(list)
    cell_count = 0

    def vertex(x: float, y: float, z: float) -> int:
        key = (round(x, 13), round(y, 13), round(z, 13))
        if key not in vertex_ids:
            vertex_ids[key] = len(vertices)
            vertices.append(key)
        return vertex_ids[key]

    for x0, x1, y0, y1 in footprint(geometry):
        z0, z1 = -geometry.height_um * UM / 2, geometry.height_um * UM / 2
        ids = (
            vertex(x0, y0, z0), vertex(x1, y0, z0),
            vertex(x1, y1, z0), vertex(x0, y1, z0),
            vertex(x0, y0, z1), vertex(x1, y0, z1),
            vertex(x1, y1, z1), vertex(x0, y1, z1),
        )
        nx, ny = max(1, round((x1 - x0) / x_spacing)), max(1, round((y1 - y0) / y_spacing))
        hexes.append((ids, nx, ny, nz))
        cell_count += nx * ny * nz
        for kind, face in {
            "xmin": (ids[0], ids[4], ids[7], ids[3]),
            "xmax": (ids[1], ids[2], ids[6], ids[5]),
            "ymin": (ids[0], ids[1], ids[5], ids[4]),
            "ymax": (ids[3], ids[7], ids[6], ids[2]),
            "bottom": (ids[0], ids[3], ids[2], ids[1]),
            "top": (ids[4], ids[5], ids[6], ids[7]),
        }.items():
            faces[frozenset(face)].append((face, kind))

    patches = defaultdict(list)
    top_y = geometry.dispersed_width_um * UM / 2 + SIDE_LENGTH
    for records in faces.values():
        if len(records) == 2:
            continue
        if len(records) != 1:
            raise RuntimeError("nonmanifold footprint face")
        face, kind = records[0]
        coords = [vertices[index] for index in face]
        if kind in ("bottom", "top"):
            patch = kind
        elif all(abs(p[0] - X_INLET) < 1e-12 for p in coords):
            patch = "dispersedInlet"
        elif all(abs(p[1] - top_y) < 1e-12 for p in coords):
            patch = "continuousUpper"
        elif all(abs(p[1] + top_y) < 1e-12 for p in coords):
            patch = "continuousLower"
        elif all(abs(p[0] - X_OUTLET) < 1e-12 for p in coords):
            patch = "outlet"
        else:
            patch = "walls"
        patches[patch].append(face)

    lines = [
        "FoamFile { format ascii; class dictionary; object blockMeshDict; }",
        "convertToMeters 1;", "vertices", "(",
    ]
    lines += [f"    ({x:.12g} {y:.12g} {z:.12g})" for x, y, z in vertices]
    lines += [");", "blocks", "("]
    lines += [
        f"    hex ({' '.join(map(str, ids))}) ({nx} {ny} {depth}) simpleGrading (1 1 1)"
        for ids, nx, ny, depth in hexes
    ]
    lines += [");", "edges ();", "boundary", "("]
    for name in ("dispersedInlet", "continuousUpper", "continuousLower", "outlet", "walls", "bottom", "top"):
        patch_type = "patch" if name in ("dispersedInlet", "continuousUpper", "continuousLower", "outlet") else "wall"
        lines += [f"    {name}", "    {", f"        type {patch_type};", "        faces", "        ("]
        lines += [f"            ({' '.join(map(str, face))})" for face in patches[name]]
        lines += ["        );", "    }"]
    lines += [");", "mergePatchPairs ();", ""]
    return "\n".join(lines), cell_count


def foam_header(name: str, field_class: str = "dictionary") -> str:
    return f"FoamFile {{ format ascii; class {field_class}; object {name}; }}\n"


def fields(
    geometry: Geometry, end_time: float, sigma: float,
    initial_water_end_x_um: float, homogeneous: bool,
    max_delta_t: float = 1e-4, write_interval: float = .005,
) -> dict[str, str]:
    inlets = ("dispersedInlet", "continuousUpper", "continuousLower")
    walls = ("walls", "bottom", "top")
    alpha = foam_header("alpha.water", "volScalarField") + "dimensions [];\ninternalField uniform 0;\nboundaryField\n{\n"
    for name in inlets:
        value = 1 if name == "dispersedInlet" and not homogeneous else 0
        alpha += f"    {name} {{ type fixedValue; value uniform {value}; }}\n"
    alpha += "    outlet { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
    for name in walls:
        alpha += f"    {name} {{ type contactAngle; theta0 180; limit gradient; value uniform 0; }}\n"
    alpha += "}\n"

    p = foam_header("p_rgh", "volScalarField") + "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n"
    for name in (*inlets, *walls):
        p += f"    {name} {{ type fixedFluxPressure; value uniform 0; }}\n"
    p += "    outlet { type fixedValue; value uniform 0; }\n}\n"

    U = foam_header("U", "volVectorField") + "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField\n{\n"
    for name, rate in (("dispersedInlet", Q_DISPERSED), ("continuousUpper", Q_CONTINUOUS / 2), ("continuousLower", Q_CONTINUOUS / 2)):
        U += f"    {name} {{ type flowRateInletVelocity; volumetricFlowRate constant {rate:.12g}; value uniform (0 0 0); }}\n"
    U += "    outlet { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
    for name in walls:
        U += f"    {name} {{ type noSlip; }}\n"
    U += "}\n"

    d = geometry.dispersed_width_um * UM / 2
    half_height = geometry.height_um * UM / 2
    initial = foam_header("setFieldsDict") + (
        "defaultValues { alpha.water 0; }\n"
        "zones\n{\n    dispersedBranch\n    {\n        type box;\n"
        f"        box ({X_INLET:.12g} {-d:.12g} {-half_height:.12g}) ({initial_water_end_x_um*UM:.12g} {d:.12g} {half_height:.12g});\n"
        f"        values {{ alpha.water {0 if homogeneous else 1}; }}\n    }}\n}}\n"
    )
    control = foam_header("controlDict") + (
        "solver incompressibleVoF;\nstartFrom startTime;\nstartTime 0;\n"
        f"stopAt endTime;\nendTime {end_time:.8g};\ndeltaT 1e-6;\n"
        f"writeControl adjustableRunTime;\nwriteInterval {write_interval:.12g};\npurgeWrite 0;\n"
        "writeFormat binary;\nwriteCompression off;\nwritePrecision 12;\ntimeFormat general;\ntimePrecision 12;\n"
        "runTimeModifiable true;\nadjustTimeStep yes;\nmaxCo 0.2;\nmaxAlphaCo 0.1;\n"
        f"maxDeltaT {max_delta_t:.12g};\n"
        "functions { phaseVolume { type volFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "fields (alpha.water); operation volIntegrate; cellZone all; writeFields false; log false; "
        "executeControl timeStep; executeInterval 1; writeControl timeStep; writeInterval 1; } "
        "boundaryFlux { type surfaceFieldValue; libs (\"libfieldFunctionObjects.so\"); "
        "patches (dispersedInlet continuousUpper continuousLower outlet walls bottom top); "
        "fields (alphaPhi.water phi); operation sum; writeFields false; log false; "
        "executeControl timeStep; executeInterval 1; writeControl timeStep; writeInterval 1; } }\n"
    )
    return {
        "0/alpha.water": alpha, "0/p_rgh": p, "0/U": U,
        "constant/phaseProperties": foam_header("phaseProperties") + f"phases (water oil);\nsigma {sigma:.12g};\n",
        "constant/physicalProperties.water": foam_header("physicalProperties.water") + "viscosityModel constant;\nnu 8.03680981595e-7;\nrho 1630;\n",
        "constant/physicalProperties.oil": foam_header("physicalProperties.oil") + "viscosityModel constant;\nnu 3.96226415094e-5;\nrho 1060;\n",
        "constant/g": foam_header("g", "uniformDimensionedVectorField") + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n",
        "constant/momentumTransport": foam_header("momentumTransport") + "simulationType laminar;\n",
        "system/setFieldsDict": initial,
        "system/controlDict": control,
    }


def generate(
    name: str, mesh: str, end_time: float, sigma: float = 0.00345,
    initial_water_end_x_um: float = -200, tag: str = "", homogeneous: bool = False,
    robust_numerics: bool = True, alpha_scheme: str = "compression",
    start_from_latest: bool = False, velocity_solver: str = "pbicgstab",
    max_delta_t: float | None = None, capillary_safety: float = .5,
    write_interval: float = .005, surface_force: str = "stock",
) -> Path:
    geometry = GEOMETRIES[name]
    case_name = name if sigma else f"{name}_zero_sigma_diagnostic"
    case_name += "_homogeneous_diagnostic" if homogeneous else ""
    case_name += "_robust" if robust_numerics else ""
    case_name += f"_{tag}" if tag else ""
    case = HERE / f"runs-{mesh}" / case_name
    if case.exists():
        raise FileExistsError(f"case already exists; refusing to overwrite {case}")
    mesh_text, cell_count = mesh_dict(geometry, mesh)
    source = HERE.parents[1] / "cases/flexible-rigid/common"
    if not (-400 < initial_water_end_x_um < -80):
        raise ValueError("initial water front must be upstream of both continuous-phase arms")
    if not (0 < capillary_safety <= 1) or write_interval <= 0:
        raise ValueError("capillary safety must be in (0,1] and write interval positive")
    x_spacing, y_spacing, nz = {
        "screening": (8 * UM, 8 * UM, 6), "throat_refined": (8 * UM, 4 * UM, 6),
        "focused": (8 * UM, 4 * UM, 12), "medium": (4 * UM, 4 * UM, 12),
        "fine": (2 * UM, 2 * UM, 24),
    }[mesh]
    minimum_edge = geometry.height_um * UM / nz
    for x0, x1, y0, y1 in footprint(geometry):
        minimum_edge = min(minimum_edge,
            (x1-x0) / max(1, round((x1-x0)/x_spacing)),
            (y1-y0) / max(1, round((y1-y0)/y_spacing)))
    # Conservative inviscid CSF capillary-wave criterion. Viscosity may relax
    # stability but does not remove the need for temporal convergence.
    capillary_bound = math.sqrt((RHO_CONTINUOUS + RHO_DISPERSED) * minimum_edge**3 / (4 * math.pi * sigma)) if sigma else None
    if max_delta_t is None:
        max_delta_t = capillary_safety * capillary_bound if capillary_bound else 1e-4
    if max_delta_t <= 0:
        raise ValueError("maxDeltaT must be positive")
    content = fields(geometry, end_time, sigma, initial_water_end_x_um, homogeneous, max_delta_t, write_interval)
    if start_from_latest:
        content["system/controlDict"] = content["system/controlDict"].replace(
            "startFrom startTime;", "startFrom latestTime;"
        )
    content["system/blockMeshDict"] = mesh_text
    content["system/decomposeParDict"] = foam_header("decomposeParDict") + "numberOfSubdomains 8;\nmethod scotch;\n"
    for key in ("fvSchemes", "fvSolution"):
        content[f"system/{key}"] = foam_header(key) + (source / key).read_text()
    if velocity_solver not in ("smooth", "pbicgstab"):
        raise ValueError(f"unsupported velocity solver: {velocity_solver}")
    if velocity_solver == "pbicgstab":
        original = "U\n    {\n        solver          smoothSolver;\n        smoother        symGaussSeidel;"
        replacement = "U\n    {\n        solver          PBiCGStab;\n        preconditioner  DILU;"
        if original not in content["system/fvSolution"]:
            raise RuntimeError("expected baseline velocity solver is missing")
        content["system/fvSolution"] = content["system/fvSolution"].replace(original, replacement)
    if alpha_scheme not in ("compression", "mplic"):
        raise ValueError(f"unsupported alpha advection scheme: {alpha_scheme}")
    if alpha_scheme == "mplic":
        original = "Gauss interfaceCompression vanLeer 1"
        replacement = "Gauss MPLIC"
        if original not in content["system/fvSchemes"]:
            raise RuntimeError("expected baseline alpha advection scheme is missing")
        content["system/fvSchemes"] = content["system/fvSchemes"].replace(original, replacement)
    if robust_numerics:
        content["system/fvSchemes"] = content["system/fvSchemes"].replace(
            "Gauss linearUpwind grad(U)", "Gauss upwind"
        )
        content["system/fvSolution"] = content["system/fvSolution"].replace(
            "nOuterCorrectors  1", "nOuterCorrectors  3"
        ).replace("nCorrectors       2", "nCorrectors       3")
    content["system/fvSolution"] = content["system/fvSolution"].replace(
        "\n}\nPIMPLE",
        "\n    Phi { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }\n}\nPIMPLE",
    )
    content["system/fvSolution"] += "\npotentialFlow { nNonOrthogonalCorrectors 0; }\n"
    if surface_force not in ("stock", "smoothed", "sharpened"):
        raise ValueError(f"unsupported surface force: {surface_force}")
    if surface_force != "stock":
        content["system/controlDict"] += f'\nlibs ("{HERE / "lib/libcurvatureSmooth.so"}");\n'
        original_grad = "gradSchemes { default Gauss linear; }"
        if original_grad not in content["system/fvSchemes"]:
            raise RuntimeError("expected gradient scheme is missing")
        content["system/fvSchemes"] = content["system/fvSchemes"].replace(original_grad, "gradSchemes { default Gauss linear; nHat curvatureSmooth 2 0.5; }")
    if surface_force == "sharpened":
        content["system/controlDict"] = content["system/controlDict"].replace("solver incompressibleVoF;", "solver incompressibleVoFCurvature;")
        content["system/controlDict"] += "surfaceTensionSharpening 0.5;\n"
    for relative, value in content.items():
        path = case / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)
    metadata = {
        "paper_doi": "10.1016/j.snb.2024.136617",
        "case": name, "mesh_level": mesh,
        "nominal_cell_um": {"screening": 8, "throat_refined": 8, "focused": 8, "medium": 4, "fine": 2}[mesh],
        "cross_stream_cell_um": {"screening": 8, "throat_refined": 4, "focused": 4, "medium": 4, "fine": 2}[mesh],
        "cell_count": cell_count,
        "geometry_um": asdict(geometry), "throat_length_um": 412.0,
        "flow_rates_m3_s": {"continuous_total": Q_CONTINUOUS, "dispersed": Q_DISPERSED},
        "surface_tension_n_m": sigma,
        "minimum_cell_edge_um": minimum_edge / UM,
        "max_delta_t_s": max_delta_t,
        "capillary_bound_s": capillary_bound,
        "capillary_safety": capillary_safety,
        "write_interval_s": write_interval,
        "surface_force": surface_force,
        "surface_force_coefficients": {"curvature_passes": 2, "smoothing_weight": .5, "force_colour_sharpening": .5 if surface_force == "sharpened" else 0} if surface_force != "stock" else {},
        "conservation_history": "phase inventory and total boundary phase/mixture flux at every solver step",
        "homogeneous_diagnostic": homogeneous,
        "robust_numerics": robust_numerics,
        "alpha_advection_scheme": alpha_scheme,
        "velocity_linear_solver": velocity_solver,
        "start_from_latest": start_from_latest,
        "initial_water_end_x_um": initial_water_end_x_um,
        "continuous_phase": "TMPTMA (OpenFOAM phase alias: oil)",
        "dispersed_phase": "HFE-7500 (OpenFOAM phase alias: water)",
        "model": "3-D VOF on measured quasi-static unstretched/stretched geometry; not transient FSI",
        "provenance": "baseline widths section 3.1; 5-mm widths approximate experimental Fig. 3a readings",
        "assumptions": "near-junction connectivity matches Fig. 1; straight inlet profiles idealise the tapered/flared entrances in Fig. 3; inlet/outlet lengths are approximate; throat length held fixed after stretch",
    }
    (case / "case.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return case


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", choices=("screening", "throat_refined", "focused", "medium", "fine"), default="focused")
    parser.add_argument("--case", choices=("unstretched", "stretch5mm", "both"), default="both")
    parser.add_argument("--end-time", type=float, default=1.2)
    parser.add_argument("--zero-sigma-diagnostic", action="store_true")
    parser.add_argument("--initial-water-end-x-um", type=float, default=-200)
    parser.add_argument("--tag", default="")
    parser.add_argument("--alpha-scheme", choices=("compression", "mplic"), default="compression")
    parser.add_argument("--velocity-solver", choices=("smooth", "pbicgstab"), default="pbicgstab")
    parser.add_argument("--start-from-latest", action="store_true")
    parser.add_argument("--max-delta-t", type=float, default=None, help="explicit override for time-step sensitivity; default is capillary-limited")
    parser.add_argument("--capillary-safety", type=float, default=.5)
    parser.add_argument("--write-interval", type=float, default=.005)
    parser.add_argument("--surface-force", choices=("stock", "smoothed", "sharpened"), default="stock", help="experimental verified-in-static-drop force discretisation; not diameter calibration")
    parser.add_argument("--homogeneous-diagnostic", action="store_true")
    numerical = parser.add_mutually_exclusive_group()
    numerical.add_argument("--robust-numerics", dest="robust_numerics", action="store_true", default=True)
    numerical.add_argument("--legacy-numerics-diagnostic", dest="robust_numerics", action="store_false")
    args = parser.parse_args()
    names = GEOMETRIES if args.case == "both" else (args.case,)
    for name in names:
        path = generate(
            name, args.mesh, args.end_time,
            0.0 if args.zero_sigma_diagnostic else 0.00345,
            args.initial_water_end_x_um, args.tag, args.homogeneous_diagnostic,
            args.robust_numerics, args.alpha_scheme, args.start_from_latest,
            args.velocity_solver,
            args.max_delta_t, args.capillary_safety, args.write_interval, args.surface_force,
        )
        print(f"generated {path}")


if __name__ == "__main__":
    main()
