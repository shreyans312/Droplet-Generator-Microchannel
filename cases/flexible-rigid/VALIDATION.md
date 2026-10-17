# Droplet measurement and numerical verification

The cases are not validated against Raj et al. The present results are a
numerical-resolution study of the uncalibrated CFD model.

## Experimental observable

Raj et al. state that region B is widened to measure the undeformed droplet
diameter from microscope images. They do not publish the segmentation rule,
measurement location, or whether diameter is an axial span, fitted-circle
diameter, Feret diameter, or projected-area-equivalent diameter. Figure 6(a)
shows flow-aligned droplets during breakup, while the device description and
Figure 8 describe measurement after relaxation in region B.

The exact Figure 6 observable therefore cannot be established uniquely. The
workflow retains axial, transverse, maximum-Feret, and area-equivalent
diameters. The area-equivalent value

$$
D_A=\sqrt{\frac{4A_{silhouette}}{\pi}}
$$

is the rotation-independent operational proxy for a circular top-view
diameter. Axial and transverse spans remain necessary shape diagnostics. No
observable is selected because it agrees better with the paper.

## Interface reconstruction

For every saved time and `alpha.water` level 0.25, 0.5, and 0.75,
`measure.py`:

1. reads the OpenFOAM phase fraction with `vtkOpenFOAMReader`;
2. interpolates the cell fraction to mesh points;
3. reconstructs the selected iso-interface;
4. projects it onto a 1 um/pixel orthographic top view;
5. tracks connected silhouettes; and
6. records a silhouette when its centroid crosses `x = 1.2 mm`.

The estimator recovers a synthetic 200 um circle as 200 by 200 um and a
synthetic 250 by 150 um ellipse as 250 by 150 um, with area-equivalent errors
below 1%. The experimental equivalence of `x = 1.2 mm` has not been
established.

## Mesh definitions and refinement coverage

All meshes use the same multiblock footprint and static wall geometry. The
refinement is global rather than interface-following.

| Mesh | Nominal x-y spacing | Depth layers | Cells per case | Cells across 100 um width/depth |
|---|---:|---:|---:|---:|
| Screening | 25.00 um | 4 | 2,816 | 4 / 4 |
| Medium | 12.50 um | 8 | 22,528 | 8 / 8 |
| Fine | 8.33 um | 12 | 76,032 | 12 / 12 |

Screening to medium is a factor-two refinement. Medium to fine is a factor-1.5
refinement and is not node-nested. The meshes are systematic global
refinements, but the generator does not preferentially refine:

- the T-junction;
- the dispersed-phase neck;
- the moving VOF interface;
- the widened-region entrance; or
- the measurement region around `x = 1.2 mm`.

The extra block divisions near the junction only describe topology; their cell
size remains the level's global spacing. The measurement line coincides with a
block boundary but has no additional resolution. A paper-correlated diameter
of 185-197 um would span approximately 7.4-7.9, 14.8-15.8, and 22.2-23.6
nominal cells on these grids. Those values are reference dimensions, not CFD
results or acceptance targets.

The rigid grids are orthogonal. For flexible screening/medium/fine grids,
maximum non-orthogonality is 17.30/21.31/21.75 degrees, maximum skewness is
0.357/0.387/0.397, and maximum aspect ratio is 1.381/1.420/1.433. All meshes
pass `checkMesh`. The deformation therefore adds moderate, acceptable mesh
distortion, but its maximum quality metrics do not improve with refinement.
The minimum vertical cell sizes in the deformed junction are 21.53, 10.76,
and 7.18 um.

## Interface resolution

The solver uses algebraic VOF with MULES, two alpha correctors, two subcycles,
and `Gauss interfaceCompression vanLeer 1`. A VOF interface has no physical
numerical thickness. As a diagnostic, the study reports

$$
\delta_{mix}=\frac{V_{0.01<\alpha<0.99}}{A_{\alpha=0.5}}
$$

in region B at crossing times. This is the volume of mixed cells divided by
the reconstructed interface area, not a material-interface thickness.

| Case | Screening | Medium | Fine |
|---|---:|---:|---:|
| Rigid | 55.42 um (2.22 cells) | 21.25 um (1.70 cells) | 17.86 um (2.14 cells) at final slug |
| Flexible | 53.67 um (2.15 cells) | 22.47 um (1.80 cells) | 15.16 um (1.82 cells) |

The physical band narrows under refinement, while it remains approximately
two cells thick. Screening droplets span only about 5-6 cells, so the mixed
band consumes a large fraction of their projected size. Axial threshold
sensitivity decreases from 23.60%/18.87% for rigid/flexible screening to
11.92%/12.11% on medium. Fine flexible sensitivity is 4.23%; fine rigid has no
droplet for this calculation.

## Measured results

The table reports alpha=0.5 values from the corrected silhouette estimator.
The period is estimated from successive measurement-line crossings and is not
a statistically converged formation frequency.

| Case | Mesh | Drops | Axial | Transverse | Max Feret | Area equivalent | Period estimate |
|---|---|---:|---:|---:|---:|---:|---:|
| Rigid | Screening | 2 | 126.50 um | 123.00 um | 125.59 um | 123.21 um | 0.05340 s |
| Rigid | Medium | 3 | 124.33 um | 124.00 um | 126.28 um | 125.05 um | 0.05010 s |
| Rigid | Fine | 0 | - | - | - | - | - |
| Flexible | Screening | 1 | 158.55 um | 162.00 um | 163.09 um | 160.97 um | - |
| Flexible | Medium | 3 | 118.32 um | 119.33 um | 120.60 um | 119.49 um | 0.04491 s |
| Flexible | Fine | 1 | 233.90 um | 200.00 um | 233.52 um | 218.51 um | - |

Fine rigid does not produce a detached droplet train by 0.36 s. Its only
region-B component grows from 252 by 191 um at 0.20 s to a 591 by 200 um slug
at 0.36 s; its centroid reaches only 1.191 mm. Fine flexible produces one
channel-spanning component. These are changes in breakup topology, not small
diameter corrections.

## Convergence assessment

For rigid screening to medium, axial, transverse, maximum-Feret, and
area-equivalent changes are -1.71%, +0.81%, +0.55%, and +1.50%. This apparent
closeness cannot establish convergence because the fine grid does not produce
a measurable droplet.

For flexible screening to medium, the changes are -25.38%, -26.34%, -26.05%,
and -25.77%. Medium to fine changes are +97.69%, +67.61%, +93.64%, and
+82.88%. Every flexible observable is non-monotonic, and the droplet counts
differ across meshes.

No observable demonstrates three-grid asymptotic convergence. An apparent
order and Richardson extrapolation are therefore not reported. The current
CFD predictions are not numerically verified well enough for a physical
comparison with Figure 6.

## Solver stability

| Case | Mesh | Time steps | Max Co | Max interface Co | Minimum alpha |
|---|---|---:|---:|---:|---:|
| Rigid | Screening | 18,011 | 0.017 | 0.017 | -0.00020 |
| Rigid | Medium | 18,111 | 0.331 | 0.217 | -0.00050 |
| Rigid | Fine | 55,541 | 0.514 | 0.279 | -0.00431 |
| Flexible | Screening | 18,011 | 0.025 | 0.025 | -0.00021 |
| Flexible | Medium | 18,085 | 0.332 | 0.252 | -0.00111 |
| Flexible | Fine | 82,168 | 0.491 | 0.338 | -0.00706 |

The configured limits are `maxCo = 0.25` and `maxAlphaCo = 0.15`. Medium and
fine runs experience transient overshoots before adaptive time-step reduction.
All runs complete without solver warnings and cumulative continuity errors
remain below 7e-9, but the fine runs require much smaller time steps and show
larger phase-fraction undershoots. Time-step convergence with the corrected
estimator is still required. Fine cases were run on three MPI ranks per case
for computational practicality; the mesh and solver controls were unchanged.

## Eight-cell criterion

The existing eight-cell flag is a useful warning, not a scientifically
defensible convergence criterion. It does not account for the two-cell mixed
band, neck resolution, curvature accuracy, wall films, contact-angle
discretization, or topology change. Medium droplets span about 9.5-10 cells
and pass the flag but do not agree with fine behavior. Fine flexible spans
24-28 cells and still produces a different breakup regime. Acceptance must be
based on grid convergence of the measured observables and formation topology.

## Remaining limitations and next numerical step

- Only zero to three droplets are available per case and mesh; startup and
  sampling effects are not separated.
- The spatial study is not independent of time-step error because fine-grid
  capillary events trigger a different adaptive history.
- The region-B geometry and measurement location remain inferred assumptions.
- The flexible wall is a fixed deformed geometry, not two-way FSI.
- The contact angle and other physical assumptions remain outside this
  numerical-resolution study.

The smallest defensible next change is a new mesh family with controlled local
resolution along the dispersed inlet, T-junction, neck, downstream interface
path, and measurement region, followed by a separate time-step study. Cell
targets must be chosen from demonstrated convergence, not agreement with
Figure 6. Physical parameters must remain unchanged during that work.

## Retained numerical-verification outputs

- [`mesh_quality.csv`](../../studies/mesh-convergence/mesh_quality.csv)
- [`mesh_measurements.csv`](../../studies/mesh-convergence/mesh_measurements.csv)
- [`mesh_summary.csv`](../../studies/mesh-convergence/mesh_summary.csv)
- [`mesh_convergence.csv`](../../studies/mesh-convergence/mesh_convergence.csv)
- [`mesh_convergence.svg`](../../studies/mesh-convergence/mesh_convergence.svg)
- [`fine_rigid_component_trace.csv`](../../studies/mesh-convergence/fine_rigid_component_trace.csv)

## Primary references

- [Raj et al.](https://doi.org/10.1007/s10404-016-1768-4)
- [OpenFOAM 13 `incompressibleVoF`](https://cpp.openfoam.org/v13/classFoam_1_1solvers_1_1incompressibleVoF.html)
- [VOF simulations of microfluidic T-junctions](https://doi.org/10.1063/1.4978801)
- [Evaluation of the OpenFOAM interface-capturing method](https://doi.org/10.1088/1749-4699/5/1/014016)
