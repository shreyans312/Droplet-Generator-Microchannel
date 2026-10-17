# Raj et al. Figure 6 comparison

These plots use the current 50-droplet, interface-based CFD measurements and
all 18 digitized Figure 6b markers at `r = 0.5`, `1.0`, and `1.5`. The
publication plot reproduces all three Raj et al. datasets and Equation 13
curves on one gridded axis. The two available calibrated CFD points are shown
only at `r = 0.5`; the plot does not interpolate CFD values for the uncomputed
flow ratios.

The calibrated plot fits one multiplier to the rigid `Kp = 0` axial CFD value
against Raj et al. Equation 13, freezes it, and evaluates the flexible
`Kp = 1.695` point as a holdout. The factor is `1.600009`. The calibrated CFD
diameter falls by 5.30%, compared with 6.29% from Equation 13, and the flexible
holdout error is +1.06%. Calibration does not change the CFD trend; it corrects
the absolute imaging/measurement scale.

Error bars are defined explicitly. Raj et al. markers use the published
experimental standard deviation extracted from the vector PDF. CFD markers
use the standard deviation across 50 measured droplets so that the visual
comparison uses the same dispersion statistic. The lag-1-adjusted 95% sampling
interval of the CFD mean is also retained in the CSV, but is not substituted
for the standard-deviation bars.

The paper does not uniquely define whether its imaged diameter corresponds to
the CFD axial, transverse, maximum-Feret, or area-equivalent observable, so all
four remain visible. The CFD mesh is not spatially converged; these are
comparison plots, not a validation claim.

```bash
/usr/bin/python3 generate.py
```

The main publication outputs are `figure6_all_series.svg` and
`figure6_all_series.csv`. The focused `r = 0.5` outputs are
`figure6_calibrated.svg` and `figure6_calibrated.csv`. `calibration.json` and
the raw and normalized plots remain available for audit.
