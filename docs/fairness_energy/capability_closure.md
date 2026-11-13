# Curvature field and smoothing energy: entry, integration and evidence

The physics and the measured numbers live in [`fairness_energy.md`](../fairness_energy.md)
and [`fairness_energy/report.md`](report.md). This note records the **delivery surface**
added in the m17465 rework round, because the previous delivery scored L1: the module
existed and was tested, but nothing in the product reached it.

## What was missing, and what now exists

| level | requirement | where it now lives |
|---|---|---|
| L2 | non-placeholder implementation | `tensorcad/fairness_energy.py` (585 lines: two independent curvature estimators, three integrands, exactly-rounded accumulation, crease split), `tensorcad/wigley.py::analytic_flank_energy` (the non-mesh reference) |
| L2 | entry point | `tensorcad/fairness_cli.py` → `python -m tensorcad.fairness_cli {sphere,wigley,hull}`; package console script `tensorcad-hull fairness …` via `tensorcad/cli.py`; public API `from tensorcad import smoothing_energy, fairness_energy_report, principal_curvature_field, …` |
| L2 | run artifacts | `runs/m17465/` produced by `runs/m17465/repro_cmd.sh`: `benchmark/{result,gates_verified}.json` + figure, three CLI records, the package-entry record, the three-case cross-check, `pytest.log`, `manifest.json` (sha256 of every artifact) |
| L3 | integration callers | `tensorcad/geometry_audit.py::audit_hull_geometry` now computes the three energies (stage `smoothing_energy`) and carries them in `GeometryAuditReport.to_dict()`; `examples/geometry_audit/run_audit.py` and `tools/generate_gates_m5479.py` therefore report the number on the real audit path. `tensorcad/cli.py` dispatches the `fairness` subcommand, same pattern as `fleet` |
| L4 | benchmark ref vs got | `benchmarks/verified/curvature_fairness/` — `result.json` holds every `reference / measured / error / tolerance` row; gates G1–G6 pass, diagnostic G7 is reported with its failing row |
| L4 | repro command | `bash runs/m17465/repro_cmd.sh` (≈ 2 min on the control plane: benchmark 50 s, cross-check 10 s, tests 60 s) |

## The numbers the benchmark adjudicates

| gate | reference | measured | error |
|---|---|---|---|
| G1 sphere `k`, R = 0.1 … 10 m | `1/R` | 10 … 0.1 1/m | ≤ 2.2e-7 1/m absolute (clause: < 1e-6) |
| G2 sphere `E_bending` / `E_willmore` / `E_gauss`, field path, level 24 | 8π / 4π / 4π | 25.1563 / 12.5781 / 12.5781 | +0.0937 % |
| G3 same on a conditioned mesh, level 32 | 8π / 4π / 4π | 26.5736 / 12.7631 / 12.2395 | +5.73 % / +1.57 % / −2.60 % |
| G4 Wigley smooth `E_bending`, `E_willmore` at 121×37 | 21.5379 / 5.38447 | 20.0849 / 5.0363 | −6.75 % / −6.47 % |
| G5 Wigley Richardson-extrapolated `E_bending` | 21.5379 | 21.5503 | +0.058 % |
| G6 invariance (scale / order / 1.6 km translation) | 0 | 0 / 2.3e-13 / 1.1e-11 | round-off |
| G7 (diagnostic) three-case cross-check | exact field | sphere 1.58 % H, 0.048 % K; ellipsoid 12.83 % H, **22.43 % K** | ellipsoid K breaks 15 % |

## Honest limitations

* **`meshmagick` is WITHHELD.** `python -c 'import meshmagick'` → `ModuleNotFoundError`;
  `pip download meshmagick --no-deps` → no version on the internal index
  (`lpypi.loongnix.cn`); `git+https` clones are refused by the sandbox network policy.
  The acceptance clause's three-case comparison therefore runs against
  `trimesh`'s normal-cycle measures and the closed forms, with the substitution
  recorded in `cross_check.provenance` in `result.json` and in the CLI record. No
  document in this repository claims a `meshmagick` agreement.
* **Ellipsoid Gaussian curvature on a coarse mesh: 22.4 % median** deviation from the
  exact field (G7). The cotangent `H` and the angle-defect `K` are only consistent in
  the refinement limit; on a 24-level marching-tet ellipsoid `H² − K` goes negative at
  some vertices, which the field reports as `nan` and counts
  (`negative_discriminant_vertices`) instead of clipping. The field path (G2) and the
  conditioned sphere (G3) show the estimators are correct where the mesh is fair.
* **Wigley mesh-path convergence is first order and from below** (−20.3 % at 41×13,
  −10.2 % at 81×25, −6.7 % at 121×37). That is the crease filter dropping an `O(h)`
  band, not an error in the integral: the Richardson-extrapolated value is +0.058 %
  from the analytic flank.
* **The crease bucket diverges by design.** Keel, stem and deck edge carry a curvature
  delta in the continuum limit, so a pointwise estimate grows like `1/h`; the split is
  reported (`crease_area_fraction`) rather than folded into a converged-looking total.
* `tests/test_geometry_audit_m5479.py::TestC5GatesArtifact` fails in a fresh checkout
  because `artifacts/m5479_geometry_audit/gates_verified.json` is not committed — a
  pre-existing gap unrelated to this work (`tools/generate_gates_m5479.py` regenerates
  it).

## Reproduce

```bash
bash runs/m17465/repro_cmd.sh
```

`runs/m17465/manifest.json` records the commit, the environment and the sha256 of every
artifact, so the committed bytes are checkable rather than decorative.