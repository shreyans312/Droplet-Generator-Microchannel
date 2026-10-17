# Flexible and rigid microchannel validation

This repository contains OpenFOAM 13 cases for a raw comparison with Figure 6
of Raj et al. (2016), *Droplet generation in a microchannel with a
controllable deformable wall*.

## Layout

```text
cases/flexible-rigid/       Rigid and spatially deformed cases
studies/tjunction-sweep/    Exact Figure 6 matrix and raw analysis
studies/time-step-convergence/ 50-droplet numerical study
studies/raj2016-comparison/ Calibrated and raw Figure 6 comparisons
references/                 Retained paper
legacy/original-case/       Original repository case
```

Generated meshes, time directories and parameter-run folders are ignored by
Git and can be recreated.

For the project's rigid-versus-flexible trend objective, use
[`studies/raj2016-comparison`](studies/raj2016-comparison/README.md). It compares
50-droplet CFD statistics with Raj et al. Figure 6b and Equation 13. The
calibrated view uses one frozen rigid-baseline multiplier and evaluates the
flexible `Kp = 1.695` point as a holdout.

## Run the two cases

```bash
cd cases/flexible-rigid
./Allrun --parallel --mesh screening
```

Mesh choices are `screening`, `medium`, and `fine`, with approximately 4, 8,
and 12 cells through the undeformed depth. Check meshes without running the
transients with:

```bash
./Allrun --mesh-only --mesh medium
```

The model uses fixed inlet flow rates across `Kp`, a parabolic static wall
profile from the paper equations, and a widened region B. Droplet measurement
now reconstructs the VOF interface with VTK and processes a calibrated top-view
silhouette. Each droplet is measured once as it crosses `x = 1.2 mm`.

The primary reported quantity is flow-direction silhouette diameter. The CSV
also retains transverse diameter, maximum Feret diameter, area-equivalent
diameter, shape ratio, resolution, and sensitivity to `alpha.water` thresholds
of 0.25, 0.5, and 0.75.

Raw CFD values are always retained. The separate calibrated graph applies one
factor, `1.600009`, fitted at rigid `Kp = 0`; it is not refitted for the
flexible point.

## Validation limits

The paper does not publish the complete T-junction axial dimensions, region-B
dimensions, contact angle, or raw Figure 6 data. Region-B width is inferred as
approximately 200 um from the embedded scale bar; its axial placement is an
explicit assumption. The wall is a static analytical profile, not two-way
FSI. These limitations prevent a claim of exact reproduction.

The paper measured more than 50 droplets per condition on at least three
devices. Results with fewer than 50 independently tracked droplets are marked
`insufficient_droplets`.

The retained fine-time-step dataset contains 50 droplets per case. Raw axial
CFD diameter decreases from 123.15 um at `Kp = 0` to 116.63 um at
`Kp = 1.695`, a 5.30% reduction. Raj Equation 13 predicts a 6.29% reduction.
After the single rigid calibration, the flexible holdout is 186.60 um versus
184.65 um from Equation 13, an error of +1.06%.

Regenerate the comparison with:

```bash
cd studies/raj2016-comparison
/usr/bin/python3 generate.py
```

See the [calibrated graph](studies/raj2016-comparison/figure6_calibrated.svg),
[calibrated values](studies/raj2016-comparison/figure6_calibrated.csv), and
[calibration metadata](studies/raj2016-comparison/calibration.json). Raw and
rigid-normalized graphs remain available for audit.

## ParaView

From the repository root:

```bash
paraFoam -case cases/flexible-rigid/rigid
paraFoam -case cases/flexible-rigid/flexible
```

Display `alpha.water` and use an isosurface value of `0.5`.

## Clean

```bash
cd cases/flexible-rigid
./Allclean
```

This removes generated meshes, time directories, logs, and temporary
measurement fields. The small validation CSV is retained as the run record.
