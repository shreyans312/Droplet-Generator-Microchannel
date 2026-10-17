# Figure 6 publication-readiness gate (medium mesh)

**Status: FAIL — do not present these CFD points as validated predictions.**

- Droplet-count criterion: 6/18 conditions reach 50 tracked droplets.
- Raw diameter comparison: 21.34% mean absolute error; 40.13% maximum error.
- Existing three-level mesh study: FAIL; the fine rigid case has no crossing and flexible droplet size changes non-monotonically.
- Experimental geometry/observable match: FAIL; the paper does not give a complete Figure 6 T-junction mask or a precise optical-diameter algorithm.
- Flexible-wall physics: the current mesh has an imposed static profile, not pressure-coupled two-way fluid–structure interaction.

## Direct 3-D VOF volume check at t = 1.5 s

Volumes are sums of alpha.water × CFD cell volume for *alpha ≥ 0.5* connected cores. The core excludes some interfacial water, so values are conservative lower estimates. Upstream is x = 400–600 µm; downstream is x = 1200–1750 µm. Rows are near Kp = 0.

| Flow ratio | Upstream core (pL) | Downstream core (pL) | Down/up | Interpretation |
| --- | ---: | ---: | ---: | --- |
| 0.5 | 917 | 896 | 0.98 | approximately conserved |
| 1 | 1125 | 2305 | 2.05 | topology/merger investigation required |
| 1.5 | 996 | 4205 | 4.22 | topology/merger investigation required |

The r = 1.0 and 1.5 increases are visible even at alpha = 0.5; they cannot be repaired by changing a plotted diameter or calibration. A connectivity change can arise from physical coalescence or from under-resolved oil films/VOF bridging. The present runs do not distinguish these mechanisms.

## Required before a publication claim

1. Obtain or measure the experimental T-junction geometry and optical diameter rule, including where in region B it was read.
2. Resolve the r-dependent connectivity change with a local interface/film refinement and time-step study; compare upstream generation volume, downstream volume, and outlet dispersed-phase flux.
3. Re-run a geometry-matched rigid baseline and then pressure-coupled flexible case; demonstrate mesh/time-step independence and ≥50 steady droplets per condition. Report raw errors before any fit.
