# Mesh-convergence study

This study checks numerical resolution of the rigid and statically deformed
cases without changing their physical parameters. Screening, medium, and fine
meshes use the same block topology and nominal spacings of 25, 12.5, and 8.33
um. The refinement is global; it does not follow the VOF interface.

`interface_band.py` measures the volume of mixed cells divided by the
alpha=0.5 surface area in region B. This is a numerical-interface diagnostic,
not a physical interface thickness. `analyze.py` combines this value with
`checkMesh`, solver-log, and corrected silhouette measurements.

The retained CSV and SVG files are numerical-verification records. OpenFOAM
meshes, time directories, and logs are intentionally excluded.

## Result

The mesh sequence is not converged. Rigid screening and medium diameters are
close, but the fine case forms a growing channel-spanning slug and has no
measurement-line crossing. Flexible diameter changes by about -25% from
screening to medium and then by +68% to +98%, depending on the definition,
from medium to fine. The flexible sequence is non-monotonic and also changes
droplet count and topology.

The mixed-cell interface band remains approximately two cells thick on all
three meshes. Medium and fine cases also show transient Courant-limit
overshoots and increasingly negative alpha undershoots. No apparent order or
Richardson extrapolation is admissible.

Outputs:

- `mesh_quality.csv`: cell sizes, counts, and `checkMesh` metrics;
- `mesh_measurements.csv`: every corrected silhouette measurement and alpha
  threshold;
- `mesh_summary.csv`: measurement, interface, period, and solver-log summary;
- `mesh_convergence.csv`: changes and trend classification;
- `mesh_convergence.svg`: the four diameter definitions across mesh levels;
- `fine_rigid_component_trace.csv`: evidence for the unresolved breakup
  topology on the rigid fine mesh.
