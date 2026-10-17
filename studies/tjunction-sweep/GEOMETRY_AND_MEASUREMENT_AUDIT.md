# Raj et al. Figure 6: geometry and droplet-measurement audit

Audit date: 2026-09-29. This is a diagnosis of the existing medium-mesh CFD
results, not a recalibration or a claim of validated agreement.

## What the original article establishes

Source: A. Raj, R. Halder, P. Sajeesh, and A. K. Sen (2016),
*Droplet generation in a microchannel with a controllable deformable wall*,
Microfluidics and Nanofluidics 20:102, DOI
[10.1007/s10404-016-1768-4](https://doi.org/10.1007/s10404-016-1768-4),
especially Figs. 1, 3, 6, 8 and Sections 2, 4.1, 4.3, 5.3, 5.4.
The local full text is `../../references/deformable-microchannel-study.pdf`.

- Figure 6 uses a **T-junction**. The water/glycerol dispersed phase enters
  through the side arm; mineral oil plus 5 wt% Span 85 is the carrier phase.
  A separate control-fluid channel acts through a thin PDMS membrane above
  region A, where droplets form. The droplet channel widens in region B,
  beyond the control channel, to measure the droplet's *undeformed diameter*.
- The membrane thickness is reported as 22.5 ± 2 µm. Section 4.1 calls the
  fabricated control channel 1000 µm wide and the "focussing channel" 100 µm
  wide. The stated **100 × 100 µm junction** in Section 5.4 belongs to the
  *flow-focusing* demonstration of Figure 8; it is not an explicit dimensioned
  drawing of the Figure 6 T-junction.
- The paper shows a top-view optical micrograph with a 200 µm scale bar in
  Figure 6a. It describes high-speed-camera imaging and says the Figure 6
  bars are experimental standard deviations from >50 droplets per condition
  on at least three devices. Its explicit ImageJ procedure concerns
  *fluorescence membrane-deflection measurement*, not droplet segmentation.
- The paper does **not** specify the Figure 6 T-junction's complete dimensioned
  mask (both inlet widths, exact junction depth, region-B transition/length),
  the precise Figure 6 droplet-image segmentation/diameter statistic, the
  location within region B where diameter was read, or the contact angle.
  Figure 6a can support only approximate image-derived dimensions. No
  author-supplied supplementary CAD or sizing protocol was found via the
  publisher record. These omissions prevent an exact geometric reproduction.

The paper's statement that droplets relax to a "spherical shape" in region B
(Section 5.4) does not define a three-dimensional volume-equivalent diameter.
The observable shown is a top-view image; the actual Figure 6 diameter
algorithm remains unknown.

## What this CFD actually does

`../../cases/flexible-rigid/mesh_generator.py` builds a 100 µm-wide,
100 µm-deep junction with a 100 µm-wide side arm. Region B is 200 µm wide
(estimated from Figure 6a's scale bar), starts 600 µm after the junction,
and ends at 1900 µm; the measurement line is at 1200 µm. The axial extent of
the pressure-induced membrane shape is an assumed cosine envelope rather than
a measured pressure/deflection field. `generate_cases.py` assigns 147°
contact angle and infers absolute inlet flow rates from Ca, assumed junction
area, viscosity, and interfacial tension. These parameters are not a complete
author-reported Figure 6 device specification.

`../../cases/flexible-rigid/measure.py` contours `alpha.water = 0.5`, projects
the interface onto the x–y plane, rasterizes at 1 µm/pixel, labels connected
silhouettes, and tracks their centroids between saved frames. It reports the
*axial bounding-box span* when a centroid crosses x = 1200 µm. This is not
volume divided by length, and the current Figure 6 comparison does not use a
volume-equivalent diameter. The code also records projected area-equivalent
diameter, transverse span, and circularity, but the validation CSV uses only
axial span. Saved frames are 0.01 s apart; no per-droplet 3-D volume or
explicit split/merge flag is recorded. Tracks are greedily associated by
centroid proximity and can therefore be ambiguous near a merger or breakup.

## Checks against the saved 18-case medium-mesh sweep

Data: `results-medium.csv`, individual `runs-medium/*/droplet_measurements.csv`,
`runs-medium/*/case.json`, and
`runs-medium/*/postProcessing/phaseVolume/0/volFieldValue.dat`.
The summary below uses all recorded crossings for the scale-factor ratios;
dropping the first two transient crossings has negligible effect on the
non-anomalous cases. Numbers are rounded and are diagnostic, not fit targets.

| Flow ratio r | Mean raw axial error vs Figure 6 | Mean raw projected-area error | Median paper / CFD axial | Median paper / CFD area | Count reaching paper's >50-droplet criterion |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 37.0% | 36.5% | 1.60 | 1.59 | 6/6 |
| 1.0 | 13.2% | 13.9% | 1.16 | 1.17 | 0/6 |
| 1.5 | 13.8% | 12.4% | 1.10 | 1.14 | 0/6 |

Thus a change from axial span to projected area-equivalent diameter cannot
explain the r = 0.5 offset. A single global multiplier calibrated at r = 0.5
is contradicted by the other flow ratios. It must not be presented as
independent validation. The six r = 1.0 cases have 36–47 tracks; the six
r = 1.5 cases have 25–48, below the paper's stated sampling level.

The r = 1.5, Kp ≈ 0 case is particularly anomalous: Figure 6 gives
`Dd/w = 2.320`; the CFD gives raw axial `Dd/w = 3.25` including startup,
where neighboring r = 1.5 cases are about 1.9–2.3. Its first recorded
crossings are 226 and 606 µm, followed by a long run near 317 µm with an
approximately 75 ms crossing period. The near-zero-Kp r = 1.0 and r = 0.5
periods are approximately 62 and 50 ms. This might reflect physical
coalescence/changed breakup or a tracking/segmentation failure; the saved
crossings alone do not distinguish them.

A subsequent direct 3-D connected-component audit (see `volume_audit.py` and
`PUBLICATION_READINESS.md`) rules out a *purely 2-D tracking* explanation for
the near-zero-Kp ratio-dependent growth. At t = 1.5 s, alpha ≥ 0.5 connected
core volume changes from about 917 to 896 pL between upstream and downstream
at r = 0.5, 1125 to 2305 pL at r = 1.0, and 996 to 4205 pL at r = 1.5.
The approximately 2× and 4× downstream changes also recur in adjacent saved
times. This is a change in the simulated phase connectivity. The existing
fields do not establish whether it is physical coalescence or unphysical
VOF bridging across an under-resolved oil film; neither can be corrected by
a plotted diameter multiplier.

The integrated `alpha.water` volume starts at 6000 pL in each case. Its
initial rate of increase before dispersed phase reaches the outlet matches
the imposed dispersed-phase inlet rate to about 0.02%, and the late total
volume remains bounded rather than growing indefinitely. This is a useful
inlet/whole-domain sanity check, **not** a per-droplet mass-conservation
measurement: the integral includes the inlet branch and all droplets, and
the present logs do not record a separated-droplet volume or outlet flux.
As a conditional mass-balance estimate, `Qd × median steady crossing period`
is about 975 pL (r = 0.5), 2437 pL (r = 1.0), and 4393 pL (the abnormal
r = 1.5 case) at Kp ≈ 0. These are *not directly measured droplet volumes*;
the estimate assumes one generated droplet per crossing and no dispersed
phase loss or coalescence between inlet and station.

## What to resolve before changing calibration

1. Obtain the authors' Figure 6 T-junction photomask/dimensioned drawing:
   carrier and dispersed inlet widths and depths, junction depth, actual
   region-B width/transition/length, and the membrane/control-channel overlap.
   Ask separately for the raw Figure 6 droplet-image sizing rule (axial,
   transverse, area-equivalent, fitted circle, or other), measurement station,
   and any image-processing settings.
2. Resolve the downstream connectivity change with local interface/film
   refinement and a time-step study. The new `volume_audit.py` measures
   connected-core **3-D alpha integrals** at a snapshot, but a full
   pinch-off/split/merge history and outlet dispersed-phase flux are still
   required. Do not infer volume from the whole-domain `phaseVolume` trace.
3. After geometry and sizing are matched, run controlled sensitivity and
   mesh/time-step convergence for junction depth, region-B width, contact
   angle, wall shape, and measurement station. Keep raw CFD and any
   calibration clearly separate. A physical baseline needs the same
   predeclared metric across all r and Kp, including the failed cases.
