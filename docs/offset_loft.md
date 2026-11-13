# Offset table -> NURBS loft (`tensorcad.offset_loft`, `tensorcad.offset_forms`)

`GAP-CAD-T1-OFFSET-LOFT`: a ship's geometry arrives as an **offset table** (half-breadths on
a station × waterline grid) and has to become a **surface**. The acceptance is written as a
comparison, not a promise: *the lofted surface's displacement and waterplane area must match
the table's own Simpson integrals to better than 0.5 %, and the dependence on the number of
stations must be measured.*

## Pipeline

```
ShipOffsetTable  ──sections()──▶  CrossSection[]  ──loft_nurbs_surface──▶  NURBSSurface
   (ordinates,                       (y, z per station,                     (one tensor-product
    provenance,                       x = const)                             patch, cubic × cubic)
    sha256)
        │                                                                        │
        ├── table_integrals()  ── Simpson over waterlines, then stations         ├── OffsetLoft.volume()/waterplane_area()
        │   (the baseline: no surface at all)                                    │   (Gauss-Legendre per knot span on the
        │                                                                        │    surface's own basis functions)
        └────────────────────────────────────────────────────────────────────────┴── OffsetLoft.mesh_integrals()
                                                                                     (tessellate + compute_hydrostatics)
```

* **Baseline** (`table_integrals`): sectional areas `2 ∫ y dz` by composite Simpson over the
  waterlines, then `V = ∫ A dx`, `LCB`, `VCB` by Simpson over the stations, waterplane
  `2 ∫ y(x, T) dx` by Simpson over the stations. Uniform grids use the classical 1/3 rule
  with a single 3/8 panel when the interval count is odd; **clustered (non-uniform)
  waterline grids** -- what a real table looks like -- use the unequal-spacing Simpson rule
  (exact for quadratics, `_quadratic_panel_weights`).
* **NURBS route** (`OffsetLoft`): with unit weights the patch is
  `S(u,v) = Σᵢ Σⱼ Bᵢ(u) Bⱼ(v) Pᵢⱼ`, so the section at `u` is the `v`-spline with control
  points `cⱼ(u) = Σᵢ Bᵢ(u) Pᵢⱼ`, the half-section area is the Green's-theorem integral
  `A(u) = ∫ y dz = Σⱼₖ cⱼʸ cₖᶻ M1[j][k]` with `M1[j][k] = ∫ Bⱼ Bₖ' dv` precomputed once,
  and `V = 2 ∫ A(u) dx/du du`. Every factor is polynomial on a knot span, so 12-point
  Gauss per span is exact to round-off. `volume_direct()` / `waterplane_area_direct()` keep
  the slow `NURBSSurface.derivatives` route as an independent check
  (agreement `1e-16` relative on a 13 × 9 table).
* **Mesh route** (`mesh_integrals`): the patch is tessellated, mirrored and closed with a
  waterplane lid, then integrated with `parametric_hull.compute_hydrostatics` -- this is the
  "loft 网格体积" half of the acceptance.

## Provenance and integrity

Every table carries `OffsetProvenance(source, method, url, licence, official, notes)` and a
`sha256` of the canonical ordinates. `write_offset_bundle` writes `<name>.csv` +
`<name>.provenance.json` (with the CSV digest as well); `read_offset_bundle` verifies both,
so an edited ordinate is rejected instead of silently changing the displacement. `official`
is `False` for every table landed so far: they are documented reconstructions of published
hull forms, not the issuing body's own tables.

## Landed tables (`docs/offset_loft/tables/`)

| table | source | grid |
|---|---|---|
| `wigley_iii` | Wigley parabolic hull, `Cb = 4/9`, `Cwp = 2/3`, `VCB = 5T/8` exactly | 21 × 17, uniform |
| `series60_cb060` | Series 60 parent form (Todd 1951), `Cp = 0.625`, `n = 5/3`, `m = 24` | 21 × 17, `p = 3` |
| `kvlcc2_class` | KVLCC2-class tanker calibrated to `L = 320 m, B = 58 m, T = 20.8 m, Cb = 0.8098, LCB = -3.48 % L` | 41 × 17, `p = 3` |

None of the ordinates is typed by hand: they are evaluated at run time from the closed
forms in `tensorcad.geomquality.baselines`, and each table's digest is recorded in
`docs/offset_loft/results.json`.

## Results (`docs/offset_loft/results.json`, `report.md`, `plot.png`)

Worst deviation of the loft from the table Simpson integrals at the landed resolutions:
**0.24 %** (NURBS quadrature, tanker waterplane) and **0.13 %** (loft mesh), against the
0.5 % requirement. Wigley is exact to round-off (`0.0000 %` vs the closed form), which pins
the method; the other two agree with their own closed forms to 0.005 % / 0.20 %.

Two findings worth keeping:

1. **Station count matters only where the plan form has sharp ends.** The full-scale
   `Cb = 0.81` plan form needs ~21 stations before the loft is within 0.5 % of an
   81-station table; 11 stations give 1.6 %. The Wigley is exact at 11.
2. **A uniform waterline grid cannot represent a power-law section.** For `y ~ zeta^(1/m)`
   the first cell above the keel carries much of the section area: with 17 *uniform*
   waterlines the Series 60 table is 1.77 % below its own closed form, and no amount of
   station refinement helps. Clustering the waterlines (`z_j = T (j/(n-1))^p`, `p = 3`,
   first waterline 0.15 mm at model scale) brings the loft to 0.005 % of the closed form.

## Limitations

* Bare hull at zero trim, no appendages, deck or superstructure; the design draft must be a
  tabulated waterline.
* The station curve must be single-valued in `z`; a non-monotone section (bulbous bow,
  overhanging flare) is rejected with a named reason instead of being mis-integrated.
* The fast integration path assumes unit weights and planar sections; a rational patch or a
  section whose `x` varies along `v` is rejected (`OffsetTableError`).
* The KVLCC2-class table is a **class hull**: `Cb`/`LCB` agreement is a consistency check of
  the representation, not an independent KVLCC2 validation. The official KRISO offset table
  is not in this repository and was not reachable from the build sandbox
  (`github.com`/`raw.githubusercontent.com` blocked; only `api.github.com` and PyPI), so no
  real ship's published offset table was landed. The loader exists to take one: drop the CSV
  plus its provenance sidecar in, and `read_offset_bundle` verifies the digest.

## Reproduce

```bash
python examples/offset_loft/run_example.py            # writes docs/offset_loft/*
python -m pytest tests/test_offset_loft.py -q         # 18 tests
```