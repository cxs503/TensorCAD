# m18488 — baseline fleet refresh: authoritative references for 4 hulls, and the two real blockers

Case: `GAP-CAD-T1-BASELINE-FLEET` (standard-hull authority geometry set, 6 hulls).
Subject worktree: branch `agent/m18488`. Machine-readable copies of everything below:
`benchmarks/verified/baseline_fleet/result.json` (benchmark) and
`benchmarks/verified/baseline_fleet/fleet_specs.json` (per-hull provenance: page
SHA-256, matched snippet, scale basis, publisher-table corrections).

## 1. What the case is

The acceptance asks for six standard hulls (KCS, KVLCC2, ONRT, JBC, DTC, Series 60)
as *watertight, manifold* hull meshes whose displacement agrees with the publisher's
own figure to < 1 %, each with a provenance manifest (source URL, licence state,
SHA-256, verification conclusion). Before this run the tree could rebuild **one**
hull (KCS) from the NMRI T2015 archive; the other five were `WITHHELD` /
`NOT_SUPPORTED`, three of them for reasons that turn out to be pipeline defects
rather than missing data.

## 2. Method (what was changed)

1. **Reproducibility repair of the acquisition/build pair.** `fetch_sources.py`
   cached publisher pages as `raw/<key>.html`, while `build_fleet.py` looked for
   `raw/pages/<key>.html`; the page keys also differed (`t2015_kcs` vs `kcs`). A fresh
   clone therefore produced *no buildable hull at all* ("run fetch_sources.py first")
   even with every byte present. The fetcher now writes `raw/pages/<key>.html`, the
   builder accepts both layouts (`locate_page`), and page → ship mapping is explicit
   (`ARCHIVES[ship][2]`). `page_text()` now unescapes HTML entities, which the old
   `&nabla;` pattern depended on by accident; the JBC patterns no longer assume that no
   MathJax markup sits between a label and its unit, and the stale `\(L_{PP}\)` regex
   (which never matched: it required literal backslashes) is replaced by a
   label-anchored pattern.
2. **A fourth-order reference source.** The SIMMAN 2020 *Ship Data* page
   (`https://simman2020.kr/contents/Ship_data_summary.php`, cached with SHA-256 in the
   registry) carries Tables 1–3: full-scale principal particulars for KVLCC2, KCS and
   ONRT — including the **model scale** and the **displacement**. This is the page that
   the previous run never located, which is why KVLCC2 had "no published table" and
   ONRT had "no scale factor obtained". Each ship's row block is sliced out before the
   patterns run, so a field can never be read from a neighbouring table.
3. **Explicit model → full-scale policy** (`GEOMETRY_SCALE_POLICY` +
   `ShipSpec.geometry_scale`/`scale_basis`): `full` (KCS, KVLCC2), `published`
   (ONRT: 48.9355 read from the publisher's own `Scale 1 48.9355 48.9355` row),
   `inferred` (JBC: no published model scale, so the point cloud is normalised to the
   published L_PP **and the manifest says so** — for that hull the volume gate is a
   block-coefficient gate, recorded as `volume_gate_basis`).
4. **Publisher-table correction, machine-checked.** The summary page prints KVLCC2's
   `Bwl (m) = 580`. The publisher's own definition `Cb = ∇/(L_PP·B_WL·T)` with
   ∇ = 312 622 m³, Cb = 0.8098, L_PP = 320 m, T = 20.8 m gives 58.0001 m, i.e. a
   decimal typo. The assembler corrects a discrepancy only when it is clean to within
   1 % and records both values plus the page SHA-256 in
   `fleet_specs.json → evidence.KVLCC2.corrections`. No number is typed in by hand.
5. **Frame guard.** Workshop files do not share an origin (one patch set places the
   *waterline* at z = 0). `build_ship` now detects a hull lying entirely below a plane
   within one draft of its deepest point and shifts the keel to z = 0; the shift is
   reported as `source.frame_shift_z`. A translation cannot change the volume.

## 3. Results

Reference values are machine-extracted from cached publisher bytes (page SHA-256 +
snippet stored per field). `artifacts/baseline_fleet/` holds the rebuilt meshes,
per-hull manifests and figures; `--samples-per-metre 6`, dx = 0.5 m, dz = 0.1 m,
240 ring points.

| hull | published ∇ (m³) | mesh ∇ (m³) | rel. error | watertight | manifold | verdict |
|---|---|---|---|---|---|---|
| KCS | 52 030 | 52 187.98 | **+0.304 %** | yes | yes | **VERIFIED** |
| JBC | 178 369.9 | 145 769.05 | −18.28 % | yes | yes | WITHHELD |
| ONRT | 8 500.78 (from published Cb·Lwl·Bwl·T) | 92.95 | −98.91 % | yes | yes | WITHHELD |
| KVLCC2 | 312 622 | — | — | — | — | WITHHELD (build error) |
| DTC / Series 60 | — | — | — | — | — | NOT_SUPPORTED (no official geometry located) |

Publisher cross-checks that are independent of the mesh: KCS Cb 0.6505 vs 0.6516
(+0.17 %), LCB −1.48 % L_PP (exact); KVLCC2 ∇/(L_PP·B_WL·T) reproduces the printed
Cb 0.8098 with the corrected 58.0 m beam (0.8099).

### 3.1 Error verdict, per key quantity

| quantity | criterion | result | < 1 % ? |
|---|---|---|---|
| KCS displacement | vs published 52 030 m³ | +0.304 % | **yes** |
| KCS block coefficient | vs published 0.6505 | +0.17 % | **yes** |
| JBC displacement | vs published 178 369.9 m³ | −18.28 % | **no** |
| ONRT displacement | vs published 8 500.78 m³ | −98.91 % | **no** (badly wrong: geometry not in one frame) |
| KVLCC2 displacement | vs published 312 622 m³ | not computed | **no** (build error) |

Honest summary: the anchor is **not** complete. One hull is verified; three hulls now
have an authoritative, machine-traced reference and a diagnosed, documented failure
mode instead of "no reference available".

## 4. Root causes of the three remaining failures (diagnosed, not fixed)

1. **ONRT — the patch set is not in one global frame.** The reader documents
   "no transformation-matrix … processing" (`tensorcad/iges.py`), i.e. IGES entity 124
   (transformation matrix) is ignored. The ONRT file is an assembly of appended hull +
   deckhouse + rudders whose members carry local frames: the raw point cloud spans
   1.006 m in x, not the 3.147 m of the model hull, so even the published 48.9355
   scale factor cannot rescue it (extent 42.5 m vs published 154 m, ratio 0.276).
   Fix = apply entity 124 when present (IGES reader), then rebuild.
2. **JBC — the field loses ≈18 % of the volume.** Units and the inferred scale are
   right (x-extent normalised to 280 m), the mesh is watertight, but the reconstructed
   block coefficient is ≈ 0.70 against the published 0.858. The sampled breadth never
   reaches the published B_WL = 45 m, i.e. part of the hull's breadth is being dropped
   (the raw control net spans only ≈ 42.7 m after scaling). Prime suspects, in order:
   the single patch quarantined by the control-net layout guard (1 of 45), and the
   `max|y|` binning of the half-breadth field on a hull whose breadth grows steeply at
   the bilge. Fix = resolve the quarantined patch's layout and re-measure breadth
   against the published 45 m before claiming a volume verdict.
3. **KVLCC2 — `field produced fewer than two usable waterline slices`.** The file is
   full scale in millimetres and its z-axis origin does not match either convention
   (keel at z = 0 or waterline at z = 0): raw z spans −28.07 m … +0.055 m, which is
   neither. The frame guard therefore shifts or does not shift without producing
   waterline slices. Fix = read the z orientation from the data (channel/deck extent vs
   draft) rather than from a single sign test, and verify against the published
   T = 20.8 m.

## 5. Scope limits

* DTC and Series 60 are still `NOT_SUPPORTED`: no official IGES/offset table was
  located this run, and the acceptance forbids substituting a placeholder. The analytic
  `tensorcad/series60.py` generator is *not* an authority source for the case.
* No bytes of the workshop archives are redistributed (neither publisher states a
  licence); the repository stores URL + SHA-256 + HTTP metadata and re-materialises the
  bytes with `fetch_sources.py`.
* `artifacts/` is not committed; `benchmarks/verified/baseline_fleet/` carries the
  machine-readable result and the provenance record so the numbers stay auditable.
* Regression proof: `python -m pytest tests/test_iges_baseline.py tests/test_fleet_gate.py -q`
  → 11 passed, 1 xfailed (the xfail is the pre-existing "DTC has no source" case).

## 6. Reproduce

```bash
python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet   # network, records SHA-256
python tools/baseline_fleet/build_fleet.py  --out artifacts/baseline_fleet --samples-per-metre 6
python -m pytest tests/test_iges_baseline.py tests/test_fleet_gate.py -q
```

## 7. Next steps (ordered by value per unit effort)

1. Implement IGES entity 124 in `tensorcad/iges.py` (submatrix + translation per
   surface) and add a unit test with a synthetic two-frame file; rebuild ONRT. This
   unblocks the acceptance item that is currently *furthest* from passing while being
   the most nearly mechanical fix.
2. Resolve the JBC quarantined patch and audit the breadth reach against published
   B_WL = 45 m before the volume gate.
3. KVLCC2 z-frame detection as described in §4.3, gated on matching the published
   T = 20.8 m.
4. DTC geometry: only worth attempting with a citable official source.