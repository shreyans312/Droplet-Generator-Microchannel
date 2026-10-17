# Rigid and spatially deformed cases

Both cases use a 100 um junction, a 200 um widened measurement region, and the
same inlet conditions:

| Quantity | Value |
|---|---:|
| Continuous flow rate | 2.34213 uL/min |
| Dispersed flow rate | 1.17107 uL/min |
| Continuous viscosity | 0.0287 Pa s |
| Dispersed viscosity | 0.000914 Pa s |
| Interfacial tension | 0.003501 N/m |
| Contact angle | 147 degrees |

The rigid case uses `Kp = 0`. The flexible case uses `Kp = 1.695`, giving an
analytical average junction deflection of 9.26 um and maximum deflection of
13.89 um. The deformation varies axially and parabolically across the channel.
It is static one-way geometry, not coupled FSI.

Run:

```bash
./Allrun --parallel --mesh screening
```

Check only:

```bash
./Allrun --mesh-only --mesh medium
```

`mesh_generator.py` writes each `blockMeshDict` and `geometry.json`.
`measure.py` reconstructs VOF interfaces, generates calibrated top-view
silhouettes, and tracks their fixed-line crossings. It requires Ubuntu's VTK
bindings and is run with `/usr/bin/python3`. `Allrun` handles this automatically
and first checks the estimator against synthetic circle and ellipse surfaces.

Individual measurements are stored in `droplet_measurements.csv`. The summary
is in `validation_results.csv` and [`validation_plot.svg`](validation_plot.svg).
The [`rigid`](rigid_measurement_overlay.svg) and
[`flexible`](flexible_measurement_overlay.svg) overlays make segmentation
visually auditable. No empirical scaling is applied.
