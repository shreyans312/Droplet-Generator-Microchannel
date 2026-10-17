# Time-step convergence study

This study tests temporal sensitivity on the fixed 22,528-cell medium mesh. It
does not establish spatial convergence or experimental validation.

Only the adaptive `maxDeltaT` cap changes: 40, 20, 10, and 5 us. All variants
retain `deltaT = 1 us`, `adjustTimeStep = yes`, `maxCo = 0.25`,
`maxAlphaCo = 0.15`, the Euler time scheme, and identical MULES, geometry, and
physical settings. All use the same 5 ms saved-field interval. The accepted
time-step distributions are parsed from the solver logs; the configured caps
are not treated as the effective time steps.

The attempted 40 us level is retained as a Courant-limited sensitivity case,
not as part of the formal refinement sequence: its actual accepted steps are
far below the cap. The systematic convergence assessment therefore uses
20, 10, and 5 us.

The reported period is the interval between successive tracked droplets at
the `x = 1.2 mm` measurement line. In a stationary one-in/one-out droplet
train this is a proxy for the formation period, but it is not a direct
measurement of the pinch-off instant. Linear interpolation between the 5 ms
snapshots refines the crossing-time estimate; it does not create additional
resolved breakup dynamics.

The solver runs are generated under the ignored `runs/` directory. Compact
logs and measurements are staged under `retained/` for analysis and removed
after the final CSV/SVG products and documentation have been produced.
During long runs, `prune_outputs.py` removes only auxiliary fields from older
snapshots. It retains every `alpha.water` field used by `measure.py` and keeps
the newest full snapshot restartable until the solve is complete.

```bash
/usr/bin/python3 prepare.py
source /opt/openfoam13/etc/bashrc
blockMesh -case runs/baseline_20us/rigid
setFields -case runs/baseline_20us/rigid
foamRun -case runs/baseline_20us/rigid
```

Measurements use the current interface/silhouette estimator:

```bash
/usr/bin/python3 ../../cases/flexible-rigid/measure.py \
  --case runs/baseline_20us/rigid \
  --case-name rigid \
  --output retained/baseline_20us_rigid_measurements.csv
/usr/bin/python3 analyze.py
```

The reported confidence intervals account approximately for lag-one serial
correlation. They are numerical repeatability descriptors, not experimental
uncertainty bounds.

No formal temporal order or Richardson extrapolation is inferred from these
adaptive-step calculations: each level is a distribution of accepted time
steps, and the coarsest attempted cap is Courant-limited.
