# Hull geometry quality audit pipeline

Mission 5283 (`TensorCAD`). Code: `tensorcad/geometry_audit.py`, `tensorcad/exchange.py`,
`examples/geometry_audit/run_audit.py`, tests: `tests/test_geometry_audit.py`.
Evidence: `artifacts/m5283_geometry_audit/` (machine gate `gates_verified.json`,
per-case JSON, three figures, `pytest.txt`).

## 1. Case: what is being audited and why

A hull that renders correctly is not a hull a solver can use. Three failure
families are invisible in a viewer and fatal downstream:

* **topology** — a seam gap between two patches, an overlapping shell, a flipped
  normal, a detached solid. The model has a surface but not a boundary of a
  volume, so any displacement integrated from it is fiction;
* **fairness** — a curvature discontinuity (a dent, an unintended knuckle, a bad
  patch join) that no amount of downstream meshing repairs and that shows up later
  as a pressure spike or a structural hard spot;
* **hydrostatics** — the lines plan and the 3-D model disagree. Shipyards check
  the two against each other precisely because each can be individually plausible.

This node builds that check as a machine gate, orders the stages so a failure stops
the stages that depend on it, and drives every criterion with both a positive and a
negative case. Two standard hulls carry the quantitative evidence.

| Case | Form | Principal dimensions | Mesh | Source of the mesh |
|---|---|---|---|---|
| Wigley III | `y = (B/2)(1-\|2xi-1\|^2)(1-(1-zeta)^2)` | L 10 m, B 1.0 m, T 0.625 m, Cb = 4/9 = 0.4444 | 161 x 41 stations x levels, 12 961 vertices, 25 918 triangles | `tensorcad.wigley.build_mesh` |
| Series 60 Cb=0.70 | `y = (B/2)(1-\|2xi-1\|^n) zeta^(1/m)`, n = 2.509, m = 46.7 | L 10 m, B 1.5 m, T 0.6 m | 161 x 41, 12 961 vertices, 25 918 triangles | `tensorcad.series60.build_series60_mesh` |

The Series 60 parent-form exponents come from the published `Cp` values of the
series (Todd 1951); the closed forms `Cb = n*m/((n+1)(m+1))`, `Cp = n/(n+1)`,
`Cm = m/(m+1)` are what pin them, and they are used here as the literature
reference. Wigley III's coefficients (`Cb = Cp = Cm = Cwp = 2/3`-family, `Cb = 4/9`)
are the classical benchmark values, and its wetted area is obtained by
Gauss-Legendre quadrature of the analytic surface (`wigley.analytic_wetted_surface_area`),
because that integral has no elementary primitive.

## 2. Method

### 2.1 Stage order (and why it is enforced)

```
import (STEP facets / patchset / mesh)
  -> watertightness ....... hard stop: displacement of an open shell is not a number about the ship
  -> fairness ............. recorded, does not stop: a sharp bilge does not invalidate the volume
  -> declared L/B/T/Cb .... hard stop: if the model is not the ship it claims, its volume is not evidence
  -> hydrostatic consistency (offsets table vs mesh, < 1 %)
```

`audit_hull_geometry` returns `stages_run`, so a rejected report says exactly how
far it got (a gap defect yields `('watertightness',)` and no hydrostatics at all).

### 2.2 Topology

* Free (one-face) edges are grouped into connected components. A component that
  closes on itself is a ring, and its *opening* is the largest chord between
  vertices at least a third of the ring apart; an open chain's opening is the
  distance between its ends.
* If there is more than one component, the defect size is the **shortest distance
  between components**: that is the distance the surface must move to close. The
  gate fires above `0.001 L` (a 1 cm opening on a 10 m model).
* Independent findings: non-manifold edges (shared by > 2 faces), duplicate
  facets, zero-area facets, winding conflicts (an edge traversed in the same
  direction by both incident faces, with the flipped area reported), signed
  volume <= 0 (inward orientation), and the number of connected shells against the
  expected count.

### 2.3 Fairness

* Discrete Gaussian curvature `K = (2*pi - sum theta_i) / A_mixed` per interior
  vertex, with Meyer's mixed (Voronoi/barycentric) area. Validated against the
  analytic sphere: `K = 1/R^2` to within 5 % at 24 x 48 resolution.
* For every edge whose two faces meet at less than 45 degrees, the relative jump
  `|K_u - K_v| / max(0.5(|K_u|+|K_v|), 0.1 * mean|K|)` is compared with
  **0.5 (50 %)**. The ratio saturates at 2.0 (a sign reversal), which is why the
  reported numbers are bounded by 2.
* Designed creases (dihedral > 45 degrees: deck edge, transom knuckle, stem) are
  excluded together with two edge rings around them, and the count is reported
  (`crease_edges`, `crease_vertices`). Without that exclusion the deck edge of a
  closed solid is an infinite curvature break and every hull fails.
* Flagged edges are joined into connected regions; a region counts as a defect if
  it covers at least `1e-4` of the surface.

### 2.4 Hydrostatics: two independent integrations

* **Offsets table (lines plan).** Half-breadths on a 201 x 51 station x waterline
  grid, integrated by composite Simpson: `V = 2 int A_s(x) dx`, `A_s = 2 int y dz`;
  waterplane `2 int y(x,T) dx`; wetted area `2 int int sqrt(1+y_x^2+y_z^2) dx dz`
  from the table's corner slopes. The table is filled either from *analytic*
  half-breadths (independent of the mesh) or digitised off the mesh by cutting
  each station with `bonjean.section_loops` (shares vertices, not quadrature --
  the ordinary lines-plan-versus-model check).
* **Mesh.** Divergence theorem on the triangles (`hydrostatics.compute_hydrostatics`),
  wetted area by triangle area sum.
* Gate: `|table - mesh| / mesh < 1 %` on volume, wetted area, `Cb` (plus `Awp`,
  `Cp` reported alongside).

### 2.5 Exchange seam

| Route | Status | Evidence |
|---|---|---|
| faceted STEP (`FACETED_BREP` over `POLY_LOOP`), write + read | implemented, round trip exact | REAL |
| neutral `.patchset.json` (patch table + triangles) | implemented | REAL |
| analytic B-rep STEP (`ADVANCED_FACE`, `B_SPLINE_SURFACE`, ...) | refused with the blocking entity names | NOT_SUPPORTED |
| IGES | refused (`ExchangeUnavailable`, needs OpenCASCADE) | NOT_SUPPORTED |
| third-party CAD validation of the written STEP | not available in this environment | WITHHELD |

The refusals are deliberate: tessellating a B-rep by guesswork would audit the
tessellator instead of the model, and a watertightness verdict built on that is
worse than no verdict. Coordinates are written in millimetres (STEP convention)
with the scale recorded, so the round trip is a x1000 then /1000 with a residual
of one ulp (measured 0 deviation on this mesh).

### 2.6 Injectors (the negative half of every criterion)

`inject_seam_offset` (duplicates the seam vertices, snaps the cut onto a vertex
plane, slides the aft patch away by `delta`), `inject_hole`, `inject_flip_faces`,
`inject_duplicate_faces`, `inject_cb_offset` (a waterline-vanishing reweighting
that moves `Cb` by a factor while leaving L, B and T untouched -- scaling the beam
would trip the beam check first and prove nothing about the `Cb` gate) and
`inject_bump` (Gaussian displacement along vertex normals).

### 2.7 Toolchain and reproduction

Python 3.11, standard library for the kernels, `numpy`/`matplotlib` only for the
figures. `pytest` for the criteria.

```bash
python examples/geometry_audit/run_audit.py          # writes artifacts/m5283_geometry_audit/
python -m pytest tests/test_geometry_audit.py -q     # 33 passed
python -m pytest tests/ -q                           # 1271 passed, 14 skipped, 2 failed
```

The two failures in the full suite are **pre-existing and unrelated**:
`tests/test_box_gz_area.py` reads `artifacts/m4873/box_gz_areas.json`, which has
never been committed (`git log -- artifacts/m4873` is empty), so it fails on a clean
checkout of the parent commit as well. Nothing in this node touches that artifact or
`box_gz`; the failure is reported here rather than silenced. Wall clock: 928 s for the
whole suite, 17 s for the audit module alone.

## 3. Results

### 3.1 Machine gates

| Gate | Criterion | Measured | Verdict |
|---|---|---|---|
| C1 | seam gap > 0.001 L detected, recall >= 0.95 | 24/24 injected seams (0.10 %-2.0 % of L), **recall 1.000**, measured/injected ratio median 1.000 | PASS |
| C2 | region with \|dK\|/K > 50 % marked, recall >= 0.95 | 9 of 12 bumps crossed the criterion (spatial jump 0.56-2.00) and all 9 were marked, **recall 1.000**; 3 sub-criterion bumps (0.19, 0.43, 0.48) correctly not marked; clean-hull false-positive area **0.0 %** | PASS |
| C3 | volume / wetted area / Cb: table vs mesh < 1 % | worst **0.55 %** (Series 60 volume), Wigley 0.020 % / 0.003 % / 0.020 % | PASS |
| C4 | +5 % Cb offset refused with the deviation named | 4/4 negative cases refused; `Cb +4.98 %` against a 2 % band, reported as `block_coefficient_mismatch` with `realised Cb 0.466302 against declared 0.444444: +4.92 %` | PASS |
| C5 | Wigley III and Series 60 within 2 % of the closed forms | Wigley V -0.020 %, S -0.002 %, Cb -0.020 %; Series 60 V/Cb **-1.14 %** | PASS |
| C6 | every criterion pytest-ised, `gates_verified.json` written | `33 passed`, gate file self-check (all 6 ids, `all_passed`) | PASS |

Full JSON: `artifacts/m5283_geometry_audit/gates_verified.json`,
`c1_gap_recall.json`, `c2_curvature_recall.json`, `c4_negative_cases.json`.

### 3.2 Hydrostatics of the two standard hulls

| Quantity | Wigley III mesh | Wigley closed form | rel. | Series 60 mesh | Series 60 closed form | rel. |
|---|---|---|---|---|---|---|
| Displacement volume [m^3] | 2.777235 | 2.777778 | -0.020 % | 6.228259 | 6.300000 | -1.14 % |
| Wetted area [m^2] | 14.878795 | 14.879063 | -0.002 % | 21.896353 | -- (no closed form) | -- |
| Block coefficient | 0.444358 | 0.444444 | -0.020 % | 0.692029 | 0.700000 | -1.14 % |
| Waterplane area [m^2] | 6.666406 | 6.666667 | -0.004 % | 10.724510 | 10.724790 | -0.003 % |
| Prismatic coefficient | 0.666641 | 0.666667 | -0.004 % | 0.714967 | 0.714966 | +0.000 % |

The Series 60 residual is the resolved-ness of its steep `zeta^(1/m)` section: the
volume deficit falls 3.68 % -> 2.14 % -> 1.14 % at 41x13 -> 81x21 -> 161x41, i.e. it
is a discretisation error of the *loft*, not of the integration.

| Quantity | Wigley table vs mesh | Series 60 table vs mesh | Gate |
|---|---|---|---|
| Displacement volume | +0.0195 % | +0.552 % | 1 % |
| Wetted area | -0.0027 % | -0.167 % | 1 % |
| Block coefficient | +0.0195 % | +0.552 % | 1 % |
| Waterplane area | +0.0039 % | +0.0046 % | 1 % |
| Prismatic coefficient | +0.0039 % | +0.0046 % | 1 % |

Digitising the table off the mesh instead of from the analytic half-breadths gives
Wigley: volume +0.0156 %, wetted area -0.0104 % -- a check of the quadrature rather
than of the geometry, and the record says which of the two was used.

### 3.3 Negative cases

| Case | Verdict | Reported defects | Stage that stopped it |
|---|---|---|---|
| seam gap 1 % L | rejected | `seam_gap` (0.100 m = 1.00 % L), `component_count` | watertightness |
| 40 flipped facets | rejected | `normal_inconsistency` (120 conflicting edges) | watertightness |
| 24 duplicated facets | rejected | `duplicate_facets`, `nonmanifold_edge` | watertightness |
| Cb +5 % (L/B/T unchanged) | rejected | `block_coefficient_mismatch` +4.98 % | declared geometry |
| 0.05 % L gap (below tolerance) | rejected | no gap defect raised; `component_count` (the shells are separate) | watertightness |

### 3.4 Exchange round trip

| Measurement | Value |
|---|---|
| faceted STEP written for the defected hull | 5 136 985 bytes, 115 430 entities, 25 918 faces |
| vertex deviation after write/read | 0.0 m (mm scale, 17 significant digits) |
| gap found after re-import of the STEP file | yes, measured 1.000 % L (injected 1.0 % L) |
| patchset JSON re-import | gap found, patch ids preserved |
| IGES / B-rep STEP | refused with an actionable message (`ExchangeUnavailable` / `ExchangeError` naming the blocking entities) |

### 3.5 Figures

`artifacts/m5283_geometry_audit/fig_c1_gap_recall.png` -- measured gap against
injected gap for the 24 seams, with the 0.001 L gate; the right panel is the
distribution of measured/injected ratios (all at 1.0, i.e. the detector measures
the defect rather than a proxy).
`artifacts/m5283_geometry_audit/fig_c2_curvature.png` -- Gaussian curvature field
on the Wigley hull (left) and with an injected bump (right); the flagged jump
region is circled where the field breaks.
`artifacts/m5283_geometry_audit/fig_c3_c5_hydrostatics.png` -- the two
consistency routes against the 1 % band (left) and the standard hulls against
their closed forms against the 2 % band (right).

## 4. What the pipeline found (a real defect, fixed here)

Running the audit on the standard hulls immediately produced a 95 % volume error
for **Series 60**, and the natural reading -- "the grid is too coarse" -- is what
`docs/series60.md` used to say. It was wrong. The vertical taper was
`zeta**m` with `m = 24...160` instead of `zeta**(1/m)`, so every section collapsed
to a sliver over the bottom 90 % of the draft:

| grid | Cb=0.60 volume error (inverted exponent) | Cb=0.60 volume error (fixed) |
|---|---|---|
| 41 x 13 | -95.32 % | -3.68 % |
| 81 x 21 | -95.70 % | -2.14 % |
| 161 x 41 | -95.78 % | -1.04 % |
| 241 x 101 | -95.81 % | -- |

Refinement cannot remove it, which is exactly how a formula error is told apart
from a discretisation error. The fix is one exponent (`tensorcad/series60.py`),
the closed form for `VCB` was corrected with it (`(m+1)/(2m+1)`, not
`(m+1)/(2(m+2))`: 0.481 T against 0.510 T at m=24, a 6.5 % error that no
"0 < VCB < T" sanity check catches), and both are pinned by new regression tests,
including one that reproduces the broken exponent so the failure mode cannot come
back silently. `docs/series60.md` carries the corrected formulas and the measured
numbers.

## 5. Notes, pitfalls and applicability

1. **The cut must be snapped to a mesh plane.** Splitting a loft between two rings
   tears it into a staircase, and then the smallest separation between the pieces
   is `|h - delta|` (one station spacing minus the offset) -- measured 7.8 mm for
   a 54.7 mm injection. That is a real defect signature but not a ground truth, so
   `inject_seam_offset(snap=True)` cuts on a vertex plane and the measured gap is
   exactly the injected one. The two cases where an unsnapped cut slipped below the
   0.1 % tolerance are the reason this paragraph exists.
2. **Designed creases must be excluded, not smoothed.** The deck edge of a closed
   solid is a 90-degree break; without the dihedral filter and a two-ring dilation
   every hull is "unfair". `crease_edges` (560 on Wigley, 726 on Series 60) and
   `crease_vertices` keep the exclusion visible.
3. **The relative jump is ill-conditioned where `K -> 0`.** A developable bottom
   next to a curved bilge produces an unbounded *relative* jump however smooth the
   surface is. The denominator floor (`0.1 * mean|K|`) damps this; it does not
   remove it, which is why the Series 60 **Cb = 0.60** parent form at 161 x 41 is
   rejected by the fairness stage (3 regions, 9.7 % of the area, along the keel
   line): that form has a *zero-radius bilge* and a curvature singularity there, so
   the verdict is physically correct -- a production hull needs a finite bilge
   radius. Cb = 0.70 at the same resolution passes (max jump 0.47 against the 0.5
   threshold).
4. **The measure of a small hole matters.** A hole punched through a symmetric
   hull yields two free boundaries, so it is classified `seam_gap` (the pieces are
   apart) rather than `open_boundary`. The kind is informational; the gate acts on
   size.
5. **The criterion is sensitive.** A 6 mm swell over a 0.4 m patch already changes
   `K` by more than 50 % on this hull; only amplitudes below ~3 mm stay under the
   threshold (test `test_small_jump_is_not_marked`). That is the criterion doing
   what it says, not a false positive -- but a user wanting "big dents only" needs
   a larger threshold or an absolute curvature-change filter.
6. **What the exchange seam does not do.** Only faceted STEP is read and written.
   Analytic B-reps and IGES are refused (OpenCASCADE is not installed here), and
   the written STEP file has been validated by an internal round trip only -- no
   third-party CAD kernel ever read it. Interoperability with a real CAD system is
   therefore *unverified*, and this document says so rather than implying it.
7. **Geometry restrictions of the offsets route.** The table assumes a symmetric
   hull whose half-breadth is single valued in `z` below the waterline. A bulbous
   bow, a flare with overhang, or an asymmetric hull needs the mesh route alone (or
   a real loft); the table route would silently under-describe the surface.
8. **Cost.** The full audit (topology + curvature + declared + consistency) is
   7.8 s on a 25 918-triangle hull and 2.3 s on 6 558 triangles, in pure Python;
   the curvature pass dominates (a dihedral angle for every edge, in a loop).

## 6. Error judgement against the baselines

Explicit answers to "is it < 3 %?" for every key quantity with an external
reference (the closed forms and the analytic surface integral; the consistency
gate is an internal two-route check and is listed separately):

| Quantity | Baseline | Error | < 3 %? |
|---|---|---|---|
| Wigley III displacement volume | closed form 2.777778 m^3 | -0.020 % | yes (< 2 % required) |
| Wigley III wetted area | analytic surface integral 14.879063 m^2 | -0.002 % | yes |
| Wigley III block coefficient | 4/9 = 0.444444 | -0.020 % | yes |
| Series 60 Cb=0.70 displacement volume | parent-form closed form 6.300000 m^3 | -1.14 % | yes (< 2 % required) |
| Series 60 Cb=0.70 block coefficient | 0.700000 | -1.14 % | yes |
| Series 60 wetted area | **no closed form exists** | -- | n/a (carried by C3 instead) |
| Wigley volume, table vs mesh | mesh integral 2.777235 m^3 | +0.020 % | yes (< 1 % required) |
| Series 60 volume, table vs mesh | mesh integral 6.228259 m^3 | +0.552 % | yes (< 1 % required) |
| Injected seam gap, 24 cases | injected `delta` | median ratio 1.000, all 24 detected | yes |
| Injected bump regions, 9 cases above threshold | injected centre | all 9 marked, 0 false positives on the clean hull | yes |
| Series 60 (pre-fix) volume | closed form | **-95.8 %** | **no** -- the defect this node found and fixed |

Honest limits, stated rather than hidden: the fairness criterion is
resolution-dependent and its verdict on a form with a curvature singularity
(Series 60 Cb = 0.60, zero-radius bilge) is "rejected"; the STEP route is faceted
only and has no third-party validation; IGES and analytic B-reps are not supported
in this environment; wetted area for the Series 60 parent form has no closed form,
so that one quantity rests on the internal consistency gate.