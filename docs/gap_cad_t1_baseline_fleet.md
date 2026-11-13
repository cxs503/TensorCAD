# GAP-CAD-T1-BASELINE-FLEET — baseline-fleet intake gate (mission m17335)

**Repository** TensorCAD · **Branch** `agent/m17335` · **Anchor** GAP-CAD-T1-BASELINE-FLEET
**New code** `tensorcad/fleet/` (`manifest.py`, `build.py`, `gates.py`, `controls.py`, `run.py`, `cli.py`)
· **Tests** `tests/test_fleet_baseline.py` (21 passed)
· **Entry points** `python -m tensorcad.fleet.cli --out-dir runs/m17335_fleet` and
`tensorcad.cli.main(["fleet", ...])` · **Artefacts** `runs/m17335_fleet/`

## 1. Case

The ship-platform geometry library ships 8 STL of which 6 are 304-triangle placeholders, KVLCC1 and
KVLCC2 are byte-identical, and `/generate` only ever emits boxes: hull-form development is blocked on
*input*. The owner asked for the six standard hulls (KCS/KVLCC2/ONRT/JBC/DTC/Series 60) from official
sources, watertight and manifold, volume within 1 % of published displacement, each with a manifest
(source URL / license / sha256 / verdict), plus a replacement plan for the placeholder STLs.

**What this environment could not do, stated up front**: there is no outbound network in this
container (HTTP requests are blocked by policy), so no official IGES or offset file could be
downloaded, and no database role was available for the TK/negative-ledger queries. The *acquisition*
half of the anchor is therefore **WITHHELD**, not delivered. What is delivered and machine-checked is
the half that stops the defect from recurring: a provenance-carrying fleet registry, three geometry
ingest routes (form-family loft, offset table, STL), a gate set that rejects placeholders and
duplicate files by name and number, and a run harness that records every verdict.

## 2. Method

1. **Registry** (`manifest.py`, `data/fleet_manifest.json`): one record per hull with L/B/T/Cb/Cm/LCB,
   a displacement-volume reference that declares its own basis (`literature_value_recalled` vs
   `closure`) and a provenance block (source, URL, license, retrieval, verified). The manifest's own
   sha256 travels into every artefact.
2. **Routes** (`build.py`): `series60_form` (this repository's Series-60 section law, only inside its
   validated Cb window 0.55–0.85), `offset_table` (loft of stations × waterline half-breadths, closing
   by construction: shared keel/deck vertices, fan only where the station has breadth) and
   `imported_stl` (binary/ASCII STL, welded then oriented). A fourth law, `parametric_sections`, solves a
   Beta section-area distribution to the record's exact Cb, Cm and LCB — used for the hulls the Series-60
   law refuses (KVLCC2, JBC, ONRT).
3. **Gates** (`gates.py`): G0 provenance complete, G1 not a placeholder (304-triangle signature, zero
   volume, axis-aligned box), G2 no byte-identical geometry in the fleet, G3 watertight (0 boundary
   edges), G4 manifold + consistently oriented + no zero-area/duplicate faces, G5 single closed shell
   (Euler 2, positive volume), G6 declared particulars internally coherent, G7 volume vs reference < 1 %,
   G8 mesh spans the declared L/B/T. Every gate returns the number it decided on and its tolerance.
4. **Independent volume integration**: each admitted mesh is handed to
   `tensorcad.hydrostatics.compute_hydrostatics` (a different code path from the builder); the two volume
   integrations agree to < 1e-6 relative.
5. **Grid convergence**: each ship is rebuilt on three grids (n/4, n/2, n) and the volume error is fitted
   for observed order `p` plus a Richardson estimate, so a discretisation deficit cannot be hidden by
   tuning a tolerance.
6. **Negative controls** (`controls.py`): the 304-triangle placeholder, a byte-identical hull pair, a Cb
   clipped into the Series-60 window, and a record whose displacement contradicts its own particulars.
   All four must be rejected, and the run records the reason.

Reproduce (single command):

```
python -m tensorcad.fleet.cli --out-dir runs/m17335_fleet --nx 481 --nz 97
python -m pytest tests/test_fleet_baseline.py -q
```

## 3. Results (numbers injected programmatically from `runs/m17335_fleet/evidence.json`)

Production grid: nx = 481, nz = 97; volume tolerance 0.01.

| ship | route | reference basis | V_ref [m³] | V_mesh [m³] | rel. error | tris | boundary edges | Euler | Richardson err | order p |
|---|---|---|---|---|---|---|---|---|---|---|
| KCS | series60_form | literature_value_recalled | 52030.0 | 51850.7 | 0.3445% | 185278 | 0 | 2 | 0.0958% | 1.07 |
| KVLCC2 | parametric_sections | literature_value_recalled | 312622.0 | 310946.1 | 0.5361% | 185278 | 0 | 2 | 0.0049% | 1.06 |
| ONRT | parametric_sections | literature_value_recalled | 8518.1 | 8493.5 | 0.2887% | 185278 | 0 | 2 | 0.2036% | 1.38 |
| JBC | parametric_sections | closure | 178378.2 | 177417.2 | 0.5387% | 185278 | 0 | 2 | 0.0038% | 1.06 |
| DTC | series60_form | closure | 155834.4 | 155407.2 | 0.2741% | 185278 | 0 | 2 | 0.1713% | 1.06 |
| Series60 | series60_form | closure | 4266.6 | 4248.9 | 0.4161% | 185278 | 0 | 2 | 0.0017% | 1.08 |

Admitted 6/6; worst mesh volume error
**0.5387 %** (JBC); worst Richardson-extrapolated
error 0.1713 %; observed orders 1.07, 1.06, 1.38, 1.06, 1.06, 1.08
(first-order, as expected for linear chords of a concave section shape). Two independent geometry laws
agree on KCS to 0.023 %.

Negative controls (must all be rejected):

| control | rejected | failed gates |
|---|---|---|
| placeholder_304_triangles | True | "" |
| byte_identical_pair | True | [["G1_not_placeholder", "G2_no_duplicate_geometry", "G3_watertight", "G4_manifold_oriented", "G5_single_closed_shell", "G7_volume_vs_reference", "G8_principal_dimensions"], ["G1_not_placeholder", "G2_no_duplicate_geometry", "G3_watertight", "G4_manifold_oriented", "G5_single_closed_shell", "G7_volume_vs_reference", "G8_principal_dimensions"]] |
| clipped_cb | True | ["G6_particulars_consistent", "G7_volume_vs_reference"] |
| incoherent_particulars | True | ["G6_particulars_consistent"] |

## 4. Figure

`runs/m17335_fleet/fleet_summary.png`: left — mesh-vs-reference volume error per ship on a log axis
against the 1 % mission tolerance (green = literature reference, blue = closure reference); right — the
negative-control verdicts.

## 5. Caveats (what this is not)

* **Official geometry is WITHHELD.** No NMRI T2015 / SIMMAN / Gothenburg IGES or offset file is in the
  repository and none could be downloaded. Five of six ships are built by *surrogate* laws that carry the
  published Cb/Cm/LCB and nothing else (no bulb, no transom, no local fairing); Series 60 uses this
  repository's Series-60 law. The gate is therefore proven on surrogates; admitting a real file only needs
  `geometry_route: "imported_stl"` or `"offset_table"` plus its sha256.
* **Reference basis**: KCS, KVLCC2 and ONRT compare against recalled published displacement volumes
  (`literature_value_recalled`, `verified=false`), JBC/DTC/Series 60 against the closure Cb·L·B·T
  (`closure`) — a consistency check, *not* an external validation. The manifest and the CSV carry that
  label per ship; `ships_with_external_reference` is empty by construction, and printing it as passed
  would be dishonest.
* Literature digits were entered offline from recollection because no source could be fetched; they are
  flagged `retrieval=offline_recollection, verified=false`. Zero-handcopy discipline is respected in the
  sense that every number in this report and in the artefacts is injected from `evidence.json`/the
  manifest, but the *literature* numbers themselves could not be re-read from the source in this
  environment.
* LCB is carried by the `parametric_sections` law and is reported, but not gated (no tolerant reference
  value was available offline); the Series-60 route does not control LCB at all.
* Placeholder replacement: the plan is the gate itself — `G1_not_placeholder` (304-triangle signature,
  zero volume, 8-corner box) and `G2_no_duplicate_geometry` (sha256 collision across the fleet) name the
  defect and refuse the file, so `/generate` can route every STL through `tensorcad.fleet` before the
  geometry reaches a solver.

## 6. Error judgement (explicit)

| quantity | worst case | criterion | verdict |
|---|---|---|---|
| mesh volume vs reference | 0.5387 % (JBC) | < 1 % (mission), < 3 % | **PASS**, an order below the 3 % bar |
| Richardson-extrapolated volume | 0.1713 % | < 1 % | **PASS** |
| watertight / manifold / Euler | 0 boundary edges, 0 non-manifold edges, χ = 2 on all six | exact | **PASS** |
| L/B/T span | exact to 1e-4 relative on beam, 1e-9 on length/draft | < 1 % | **PASS** |
| negative controls rejected | 4/4 | all | **PASS** |

Is the key quantity within 3 %? Yes — 0.54 % worst case, and within 1 % as the anchor demands. What that
does **not** prove is that the surrogate meshes are the benchmark ships; that claim needs the official
files and is reported as WITHHELD.
