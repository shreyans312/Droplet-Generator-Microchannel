# Figure 6 calibrated CFD sweep

The matrix contains the 18 markers extracted from Figure 6b at `r = 0.5, 1.0,
1.5`, `Ca = 0.032`, and `lambda = 0.03185`. Values and standard-deviation bars
are in `reference_data/fig6_digitized.csv`.

## Workflow

```bash
source /opt/openfoam13/etc/bashrc
python3 generate_cases.py --group all --mesh medium --overwrite
python3 run_sweep.py --group all --mesh medium --jobs 12 --prepare-only
python3 run_sweep.py --group all --mesh medium --jobs 12
python3 analyze.py --mesh medium
python3 ../raj2016-comparison/generate.py
```

For an unattended run, `finalize_when_complete.py` waits for all 18 completion
markers and then performs the last two analysis/plotting commands.

Outputs are `results-medium.csv`, `SUMMARY-medium.md`,
`plots/figure6_raw_medium.svg`, and the publication outputs in
`../raj2016-comparison/`. Run directories are ignored by Git.

The analysis uses fixed flow rates across `Kp`, the analytical spatial wall
profile, widened region B, and the shared VTK silhouette measurement. The
production matrix uses the same medium mesh and 10 us maximum time step as the
two cases used to establish the frozen calibration. Saved fields are binary
and compressed at 0.01 s intervals. The run driver preserves every
`alpha.water` snapshot required for silhouette tracking, one complete restart
snapshot, and compressed solver logs.

Analysis writes both raw and calibrated dimensionless axial diameters. The
single factor in `../raj2016-comparison/calibration.json` is applied unchanged
to every point; no pointwise or flow-ratio-specific fitting is permitted.

`screening`, `medium`, and `fine` use about 4, 8, and 12 cells through the
undeformed depth. A defensible result requires mesh and time-step convergence
and at least 50 tracked droplets per condition.

## Publication-readiness audit

The current sweep **fails** publication-readiness checks. Run
`/usr/bin/python3 publication_gate.py --output PUBLICATION_READINESS.md` to
regenerate the fail-closed report. `volume_audit.py` independently measures
connected 3-D water-core volumes from saved VOF fields at two alpha cutoffs.
At near-zero Kp, downstream core volume is approximately twice upstream at
`r = 1.0` and four times upstream at `r = 1.5`, while it is approximately
unchanged at `r = 0.5`. These topology changes, missing experimental geometry,
and failed fine-mesh convergence must be resolved before claiming validation;
the fixed diameter multiplier does not resolve them. See
`GEOMETRY_AND_MEASUREMENT_AUDIT.md` for paper and model provenance.
