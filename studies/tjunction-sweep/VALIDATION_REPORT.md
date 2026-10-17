# Figure 6 validation protocol

## Reference provenance

The target is Figure 6b of Raj et al. (2016), DOI
`10.1007/s10404-016-1768-4`. The paper reports more than 50 droplets per
condition, experiments on at least three devices, and standard-deviation error
bars. Author raw data and supplementary geometry files were not found on the
publisher record.

The 18 plotted markers and error bars were extracted from the PDF vector
artwork and retained in `reference_data/fig6_digitized.csv`.

## Geometry and flow corrections

- Inlet flow rates remain fixed when `Kp` changes.
- The membrane follows the paper's parabolic transverse displacement and a
  localized static axial envelope.
- The measurement outlet widens from 100 um to 200 um. The latter is inferred
  from Figure 6's embedded 200 um scale bar.
- VOF iso-interfaces are converted into calibrated top-view silhouettes.
- Droplets are measured once as their centroids cross `x = 1.2 mm`.
- Axial, transverse, maximum-Feret, and area-equivalent diameters are retained.
- Alpha levels 0.25, 0.5, and 0.75 quantify interface-threshold sensitivity.
- The empirical factor previously used to scale CFD diameters has been removed.

Every generated case stores the numerical geometry and assumptions in
`case.json`. The active cases store them in `geometry.json`.

## What remains uncertain

The paper does not publish the exact Figure 6 T-junction depth, inlet lengths,
region-B axial location, full widened-region dimensions, wall contact angle,
or spatial controlling-channel pressure. The current 100 um undeformed depth,
axial envelope, and region-B position therefore remain assumptions.

The deformation is a static analytical geometry. Stock OpenFOAM 13 in this
repository has no verified fluid-solid coupling library, and the paper lacks
the boundary data required for a defensible two-way FSI reconstruction.

## Publication threshold

A case is only a validation candidate after:

1. screening, medium, and fine results demonstrate mesh convergence;
2. a second maximum Courant number demonstrates time-step convergence;
3. at least 50 independently tracked droplets are measured;
4. raw errors and simulation dispersion overlap the experimental uncertainty
   to an agreed tolerance.

No calibrated result is accepted as validation.

## Current execution status

The two active 2,816-cell screening cases completed to 0.36 s. Fixed-line
tracking found two rigid droplets with a mean axial diameter of 126.50 um and
one flexible droplet at 158.55 um. Errors against the paper correlation are
-35.80% and -14.13%.

All measurements are under-resolved. Moving the interface threshold from 0.25
to 0.75 changes the rigid mean from 155.72 to 96.64 um and the flexible result
from 188.47 to 130.82 um. The previously reported time-step comparison used the
discarded cell-column estimator and is no longer treated as evidence.

All three mesh levels pass `checkMesh`:

| Level | Cells | Maximum flexible skewness |
|---|---:|---:|
| Screening | 2,816 | 0.357 |
| Medium | 22,528 | 0.387 |
| Fine | 76,032 | 0.397 |

Medium and fine transients remain unexecuted, neither screening case has 50
droplets, and the screening interface is highly threshold-sensitive. The
present result fails the publication threshold.
