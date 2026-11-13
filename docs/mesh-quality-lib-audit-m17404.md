# Board audit — m17404 / `GAP-CAD-D1-MESH-QUALITY-LIB` (mesh-quality metric library)

*Auditor:* m17437 (board-audit, 交差核验单 — not a redo). *Repo:* TensorCAD.
*Audited commit:* `9b66af8` (`geomquality: canonical mesh-quality metric library + differential
cross-check + CI gate`), tree `/root/tp-agent/work/m17404`.
*Machine-readable artefact of this audit:* `artifacts/mesh_quality_audit/dialect.json`.

## 1. Scope and what was checked

Ten measures (`area`, `min/max_angle_deg`, `edge_ratio`, `aspect_ratio`, `aspect_frobenius`,
`gamma`, `radius_ratio`, `scaled_jacobian`, `condition_number_inverse`), the aggregate report,
the criteria/gate path, and the four claimed deliverables of the mission: artefact, benchmark
numbers, `repro_cmd`, six-section report.

Executed by the auditor in the audited tree (not read from the report):

| check | command | result |
|---|---|---|
| reproduction | `python tools/mesh_quality_crosscheck.py` | exit 0, `PASS`, all six meshes recomputed, artefact rewritten |
| unit/invariant suite | `python -m pytest tests/test_mesh_quality_lib.py -q` | 10 passed |
| no regression in neighbouring gates | `... tests/test_gate.py tests/test_mesh_conditioning.py tests/test_geometry_audit.py -q` | 77 passed, 1 skipped |
| anchors vs ideal values | `ideal_values` in the artefact | worst deviation 7.10543e-15 |
| route agreement | `differential.*.metrics.*.difference` | worst 2.40452e-08 (acos-conditioned angle on `suboff_coarse`) |
| production reachability | `grep -rn geomquality --include=*.py tensorcad` | only inside `tensorcad/geomquality/*`; no solver/gate path imports `geomquality.quality` |

Independence limitations: the auditor used the same host, no gmsh/VTK/meshmagick
(gmsh: absent, meshmagick: absent, vtk: absent), and no external reference values are available for the same meshes.

## 2. Anchor criterion, item by item

Anchor text: *与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异<1e-6 或差异解释）；可作 CI 门用*.

| item | criterion | measured | verdict |
|---|---|---|---|
| ≥3 meshes | ≥3 | 6 (`box_unit`, `series60`, `wigley`, `suboff_coarse`, `suboff_taubin30`, `sliver_patch`) | PASS |
| differential < 1e-6 | vs upstream | **not run** — all three upstream packages absent, `pip index` on this host (loong64) has no `gmsh` wheel | NOT RUN |
| "或差异解释" | explanation recorded | `upstream.*.reason` records the verbatim import error for each; the doc re-words acceptance as `|library − independent route| < 1e-6` | PASS (explained) |
| usable as CI gate | gate + thresholds + sensitivity | `quality_gate` + `DEFAULT_CRITERIA` + sweep in `gate.json.sensitivity`; sliver/raw-Suboff FAIL, smoothed PASS | PASS at library level, **not wired into CI** (see F2) |
| `tensor_benchmark_cases` row for the anchor | benchmark judge | query returns 0 rows for `anchor_slug='GAP-CAD-D1-MESH-QUALITY-LIB'` | ABSENT |

Because the benchmark table holds no row for this anchor, no external reference `got` exists to
compare against: the only `ref` available is analytic (ideal values and closed forms), and the
only cross-check is an internally written independent route. That is a real verification, but it
is *verification*, not *benchmark validation*.

## 3. Numbers reproduced by the auditor

Ideal (closed-form) values vs measured:

| measure | expected | got | difference |
|---|---|---|---|
| aspect_frobenius | 1 | 1 | 2.22e-16 |
| aspect_ratio | 1 | 1 | 2.22e-16 |
| condition_number_inverse | 0.57735 | 0.57735 | 1.11e-16 |
| edge_ratio | 1 | 1 | 0 |
| gamma | 1 | 1 | 2.22e-16 |
| max_angle_deg | 60 | 60 | 7.11e-15 |
| min_angle_deg | 60 | 60 | 7.11e-15 |
| radius_ratio | 1 | 1 | 1.11e-16 |
| scaled_jacobian | 0.866025 | 0.866025 | 0 |

Worst |library − independent route| per mesh:

| mesh | max abs deviation |
|---|---|
| box_unit | 2.132e-14 |
| series60 | 3.908e-14 |
| sliver_patch | 2.555e-11 |
| suboff_coarse | 2.405e-08 |
| suboff_taubin30 | 3.553e-14 |
| wigley | 4.263e-14 |

## 4. Findings (severity-ranked)

### F1 (medium, blocking for the "same threshold" claim) — the two aspect-ratio dialects are not the same quantity, so the shared limit 25 is not transferable

`mesh_conditioning`/`validation` ship `aspect = l_max²/(2A)`; the new library ships
`aspect = l_max·p/(4√3·A)`. The ratio is `p/(2√3·l_max) ∈ (1/√3, √3/2] = (0.577, 0.866]`, i.e.
bounded away from 1 and shape-dependent, so the library's number is systematically *smaller*:

| mesh | validation l_max²/2A | canonical aspect | ratio | gate/validation |
|---|---|---|---|---|
| box_unit | 2 | 1.3938 | 0.69692 | 2 |
| series60 | 6.1679 | 3.8244 | 0.62005 | 2 |
| sliver_patch | 400 | 231.23 | 0.57807 | 2 |
| suboff_coarse | 14826 | 8559.9 | 0.57736 | 2 |
| suboff_taubin30 | 8.4281 | 4.9337 | 0.58539 | 2 |
| wigley | 5.0966 | 3.2068 | 0.6292 | 2 |

Measured consequence, per mesh, of applying both dialects at the same numeric limit 25:

| mesh | triangles | elements whose verdict flips | of which pass every other canonical check |
|---|---|---|---|
| series60 | 1998 | 0 | 0 |
| suboff_coarse | 9088 | 33 | 8 |
| suboff_taubin30 | 9088 | 0 | 0 |
| wigley | 1998 | 0 | 0 |

Single-element counterexample (no mesh, no threshold tuning, scale-invariant to 2e-16):

```
triangle [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0741, 0.0]]
legacy = 26.9906  -> FAIL at 25
canonical = 15.6044  -> PASS at 25
min_angle = 4.23787 deg (>2), gamma = 0.085407 (>0.05),
edge_ratio = 1.99453 (<100)
```

The code comment in `quality.py` (`DEFAULT_CRITERIA`) says the limits are "the thresholds this
repository already ships elsewhere … both re-expressed here in the canonical definitions", and
`docs/mesh-quality-lib.md` §3 says "the repository's existing decimation guard (25), canonically
expressed". The first half is true (2° and 25 are the shipped numbers); "re-expressed" is not:
adopting the library's gate relaxes the aspect guard by a factor of 1.15–1.73. The library's own
artefact already *prints* both numbers side by side (`dialect_audit`), so this is a wording and
contract defect, not a hidden arithmetic error — but a reader who replaces the old gate with the
new one silently changes behaviour.

**Remediation:** either (a) state the conversion `canonical = legacy·p/(2√3·l_max)` next to the
limit and re-derive the canonical limit from 25 (i.e. 14.4–21.7, or simply drop the equivalence
claim), or (b) keep 25 and say plainly that the new gate is the more permissive one.

### F2 (medium) — the "CI gate" is not in CI, and the library is not reachable from production

The gate exists as `quality_gate`/`gate_report` returning a verdict dict; `tools/mesh_quality_crosscheck.py`
exits 1 on violation, so it *can* be a door. But: no production module imports
`tensorcad.geomquality.quality` (only `tools/` and `tests/`), no workflow in `.github/` invokes the
driver, and `mesh_conditioning`/`gate` keep their own dialects. So the mission's "可作 CI 门用" is
readiness, not integration.

### F3 (low) — definition provenance of `gamma` is transcribed, not verified

`gamma` is labelled *Gmsh 'gamma'* in `METRIC_DEFINITIONS` and the doc, but gmsh is absent and the
host's package index has no `gmsh` wheel, so the transcription was never diffed. The formula
satisfies the only locally checkable constraint (ideal = 1 at the equilateral element). Gmsh's own
normalisation of Γ must be re-checked on a host where `import gmsh` works; until then the label
carries the same status as the rest of the upstream column: transcribed. (Recorded as an
unverified claim, not as an error — I could not exhibit a counterexample.)

### F4 (low) — duplicate definitions left in place

`tensorcad/mesh_conditioning._aspect_ratio` and `tensorcad/validation` still define aspect ratio
their own way; the new module adds a fourth definition without deprecating or aliasing the old
ones. The doc names the four dialects (good) but nothing prevents a fifth from appearing.

## 5. Report quality (six-section requirement)

`tensor_reports` id 16393 for m17404: the six sections are present (case/method/parameters/
comparison/sensitivity/failure modes), 1 figure, `qa_status=pending`. The comparison table states
`ref` as ideal/closed forms and `got` as measured values — consistent with §2 above, i.e. it is
honest about *not* having an upstream reference. Two wording items overstate: the "same measure"
claim (F1) and "可作 CI 门用" (F2).

## 6. Error judgement against a baseline (explicit)

| key quantity | reference | measured | relative error | < 3 % ? |
|---|---|---|---|---|
| `aspect_ratio`, `gamma`, `edge_ratio`, `radius_ratio`, `scaled_jacobian`, `frobenius`, `condition_number_inverse` on the equilateral element | closed form (ideal 1 / 0.866 / 0.577) | see §3 | <= 2.3e-16 | yes (by 14 orders of magnitude) |
| 10 measures on the six hull meshes | auditor-independent route (Heron, law of cosines, SVD) | see §3 | worst 2.4e-08 absolute on an angle of ~0.0039 deg, i.e. ~6e-6 relative; all other measures <= 4.3e-14 absolute | yes |
| the same on rational meshes | exact `fractions.Fraction` arithmetic | `box_unit` 2.13e-14, `sliver_patch` 2.54e-12 | <= 1e-12 | yes |
| **vs gmsh / meshmagick / VTK** | **external reference** | **not computed** | **n/a** | **cannot be answered — no baseline was run (all three packages absent, no `gmsh` wheel on the loong64 index). This is the honest answer, not a pass.** |
| aspect dialect vs the shipped gate | `mesh_conditioning` limit 25 | library value is 0.577–0.866x the shipped one | systematic, not round-off | **no — the two numbers are not the same quantity (F1)** |

## 7. Disposition

**Approve with notes** for the library itself: the implementation is real, the arithmetic is
independently verified to 2.40452e-08 and exact to
1e-12 on rational meshes, the artefacts regenerate, the tests assert invariants
(ideal values, exact arithmetic, rigid-motion/scale invariance, degenerate→worst-not-NaN, gate
direction, aggregate consistency) rather than merely executing code.

**Maturity level: L2** — independent functional implementation. Not L3: no production caller.
Not L4: the anchor's own benchmark (gmsh/meshmagick differential) was not executed and no
`tensor_benchmark_cases` row exists, so there is no external reference and no error against a
baseline; internal verification against closed forms is L1-grade evidence carried along, not a
benchmark.

Requested changes (F1, F2) are documentation/contract level and do not invalidate the deliverable.
