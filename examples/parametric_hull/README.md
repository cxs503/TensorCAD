# TensorCAD first geometry example: parametric hull

This headless example creates a simple symmetric sampled hull surface from principal dimensions. It exports an OBJ triangle mesh, longitudinal breadth offsets (CSV), and metadata (JSON).

## Run from the repository root

```bash
python -m pip install -e ".[test]"
pytest
tensorcad-hull --length 12 --beam 2.4 --draft 0.9 --nx 61 --nz 17 --output results/hull-demo
```

Without installing the CLI entry point, run:

```bash
python -m tensorcad.cli --output results/hull-demo
```

## Outputs

- `hull.obj`: sampled triangular surface mesh.
- `offsets.csv`: longitudinal station and overall breadth.
- `metadata.json`: input parameters, units, coordinate frame, mesh counts, basic validation results, and known engineering limitations.

## Scope and limitations

This is an initial software integration benchmark, not a production ship hull. The analytic breadth distribution is intentionally simple. Basic validation checks dimensions, finite coordinates, counts, and index ranges; it does not prove watertightness, surface fairness, hydrostatic accuracy, stability, resistance, or class-rule compliance. STEP/B-rep export and CAD-kernel integration are follow-up work after the minimal API and CI are stable.
