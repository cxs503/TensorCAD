# Cap adaptation: dropping profile columns where a hull closes

> Node: `TensorCAD`, self-start seat m13763 (verifier). Code:
> `tensorcad/parametric_hull.py::tessellate_surface(cap_adaptation=…)`,
> `_adaptive_profile_rows`. Evidence: `examples/parametric_hull_cap/run_cap_adaptation.py`
> → `examples/parametric_hull_cap/outputs/{results.json,cap_adaptation.png}`; tests
> `tests/test_cap_adaptation.py`. Every number below is written by that script into
> `results.json` (no hand-copied values); the script exits non-zero if a gate fails.

## 1. The gap this closes

`docs/parametric_hull/README.md` §6.2 recorded an open limitation: the DARPA SUBOFF
bare hull mesh "does not meet criterion (e)" — min scaled Jacobian **0.13**, max aspect
ratio **7.9** — and attributed it to the tensor-product grid turning the polar caps
into stretched fans, concluding that "fixing it needs an unstructured/polar-cap
remesher — out of scope here".

The diagnosis was right, the conclusion was pessimistic. The two caps *do* need fewer
circumferential points, but they do not need a new mesh data structure: the columns can
be dropped **inside** the existing structured builder, and the zipper triangles fall
out of the quad loop for free (see §3).

## 2. Why refining `dx` does not work (measured, not argued)

Structured grid, SUBOFF bare hull, `n_v = 16`:

| `n_u` | `u_spacing` | faces | min scaled J | max aspect | mean aspect |
|---:|---|---:|---:|---:|---:|
| 160 | uniform_x | 9480 | 0.1267 | 7.88 | 2.20 |
| 320 | uniform_x | 19 080 | 0.1562 | 6.38 | 3.56 |
| 640 | uniform_x | 38 280 | 0.1273 | 7.85 | 6.70 |
| 160 | cosine | 9480 | 0.0927 | 10.74 | 1.90 |
| 640 | cosine | 38 280 | 0.0675 | 14.73 | 5.15 |

The Jacobian **does not converge** and cosine clustering towards the tips makes it
*worse*. Anatomy of the worst triangle at the stern (unstructured probe, reproduced in
`work/probe_triangle_anatomy.py` during this node): its shortest edge is
**circumferential** (4.4 mm) while the two long edges are **meridional** (34.5 mm) — the
ring at `r = 0.018 m` still carries 30 points around, i.e. arc `2πr/30 = 3.8 mm` against
`dx = 27 mm`. Refining `dx` shortens the meridional edge and leaves the circumferential
one alone, so the aspect ratio gets *worse*, and that is exactly what the table shows.
The binding constraint is the column count, not the station density.

## 3. Method

`_adaptive_profile_rows` computes, for each station `i`:

```
L_i     polyline length of the section profile (all n_v samples)
h_i     local axial spacing = |x_{i+1} - x_{i-1}| / 2   (one-sided at the ends)
s_i     = clip(round(L_i / (cap_aspect_target * h_i)), 2, n_v - 1)   profiles segments
v_j     = round(j * s_i / (n_v - 1)) / s_i                           j = 0 … n_v-1
```

so the circumferential cell `L_i / s_i` matches the axial one up to
`cap_aspect_target` (default `1.0`). Two properties make this a *small* change:

* the samples stay on the **same** `n_v` index grid, so the existing quad loop,
  `mirror_y` seam handling and waterplane lid are untouched; index `0` and `n_v-1` are
  never moved (`v = 0` and `v = 1` are the keel and deck seam, both on the centre
  plane, which is also why the waterplane area is bit-identical);
* coincident vertices appear in adjacent rows, and the pipeline's existing `weld` +
  `remove_degenerate_faces` (already needed for the degenerate stem column) collapses
  them, which *is* the zipper between a ring of `s_i` and a ring of `s_{i+1}` points.

The option is **off by default** (`cap_adaptation: bool = False`), so every existing
call site and every previously reported number is unchanged.

## 4. Results — SUBOFF bare hull, `n_u = 160`, `n_v = 16`

| variant | faces | min scaled J | max aspect | mean aspect | watertight | manifold | degenerate | volume [m³] |
|---|---:|---:|---:|---:|---|---|---:|---:|
| structured (baseline) | 9480 | 0.1267 | 7.88 | 2.20 | yes | yes | 0 | 0.693885 |
| cap-adapted | 8904 | **0.3251** | **3.03** | 2.03 | yes | yes | 0 | 0.693770 |

Acceptance gates (all evaluated by the example script, `results.json.gates`):

| gate | measured | verdict |
|---|---|---|
| min scaled Jacobian > 0.30 (criterion (e)) | 0.3251 | PASS |
| max aspect ratio < 5.0 (criterion (e)) | 3.03 | PASS |
| watertight and edge-manifold | yes / yes | PASS |
| no degenerate faces | 0 | PASS |
| volume drift vs baseline < 0.1 % | 0.017 % | PASS |
| Wigley: Jacobian not degraded | 0.365 → 0.595 | PASS |
| Wigley: closed-volume drift < 0.5 % | −0.37 % | PASS |

The cap adaptation also **improves the Wigley hull** (min J 0.365 → 0.595, aspect
2.74 → 1.68) because its stem is a cap too, at a closed-volume cost of −0.37 %; the
waterplane area is bit-identical (`ha.waterplane_area == hb.waterplane_area`,
`tests/test_cap_adaptation.py::test_wigley_quality_improves_and_waterplane_is_unchanged`).

## 5. What it costs (the honest ledger)

* **Chord sagitta where columns are dropped.** 7 rings are coarsened at
  `cap_aspect_target = 1.0`; the worst single-ring chord deviation from the exact
  lofted profile is **2.62 mm = 1.03 % of the hull radius** at `x = 3.90 m` (in the
  stern taper). The section polygons are coarsened *only* where the ring is already
  finer than the axial spacing, so the surface error stays below the axial
  discretisation error it is trading against.
* **Volume.** SUBOFF drift 0.017 %; Wigley −0.37 % (its stem runs to a point, so more
  of the hull is coarsened). On a hull where a fine stem matters, keep
  `cap_adaptation=False` for the hydrostatics run and use the adapted mesh for the
  solver.
* **The `-1.93 %` SUBOFF volume error against the analytic bare-hull volume is
  pre-existing** (present in the structured baseline, a cap-resolution property at
  `n_u = 160`); this node changes it by 0.017 %, it does not fix it. The SUBOFF *radius*
  criterion is evaluated on the lofted surface, never on the mesh polygon, so it is
  unaffected.
* **Not a general remesher.** The rule is a per-station profile coarsening; it cannot
  fix a *closed* transom corner, a knuckle, or an anisotropic box — those are not
  caps. It assumes the profile direction runs bottom-to-top on each section, which is
  the `CrossSection` contract.

## 6. Reproduction

```bash
python examples/parametric_hull_cap/run_cap_adaptation.py     # writes outputs/, exits 1 on a gate
python -m pytest tests/test_cap_adaptation.py -q              # 7 tests
```

## 7. Error verdict

Every key quantity above is a machine gate with the threshold stated in the table; the
headline acceptance (`min J > 0.30`, `aspect < 5`) is exceeded by the adapted mesh with
a margin of 8 % and 39 % respectively, and the two operations that could be "gamed" by
moving geometry are bounded separately: volume drift 0.017 % (gate 0.1 %) and chord
sagitta 1.03 % of the radius. Against the analytic bare-hull volume the residual is
1.93 % and is **not** reduced by this change — it is a baseline property, reported as
such rather than absorbed into a "pass".