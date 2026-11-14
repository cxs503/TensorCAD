# Board audit — m20384 / GAP-CAD-T1-SOURCE-HUNT (TensorCAD)

Audit mission: **m20425** · auditor: `board-audit` · subject: mission m20384, branch
`agent/m20384`, commit `460859e0927cddfa5ab6db8a1ea5a2b0450f643a`
("TensorCAD: source-hunt harness for the baseline fleet geometry provenance (DTC/ONRT/JBC)").
This is a **verdict/verification** round, not a redo: nothing in the subject tree was modified.

Deliverable under audit (four files, +979 lines):

| file | role |
|---|---|
| `tools/baseline_fleet/hunt_sources.py` (425 L) | probe → digest → parse → record tool, CLI entry |
| `artifacts/baseline_fleet/source_hunt.json` (13 408 B) | committed record: 11 endpoints, 4 ship blocks, cross-check, open questions |
| `tests/test_source_hunt.py` (124 L) | 9 offline tests (record invariants + parser on unseen markup) |
| `docs/audit_m20384_GAP-CAD-T1-SOURCE-HUNT.md` (80 L) | prose report |

Record SHA-256 `7b5474ff2a64609f8ac18fc69c3fef0a1b7037050a8d5be13f4ddd3473bf8264`.

## 1. Verdict

**Level 2 (L2 — independent functional implementation).** Real code with a real
input→compute→output path (HTTP probe → SHA-256 → table parse with matched row text →
machine-readable record), independently re-executable by the auditor today, covered by 9
assertion-bearing tests. It is **not** L3 (nothing in the platform or solver calls it; its
only in-repo callers are its own tests and its doc) and **not** L4 (see F1: the benchmark
referee has no case for this anchor, and the report's "0 %" rows are transcriptions of the
same bytes, not accuracy validation). The falsified DTC lead and the withheld JBC scale
factor do not lower the grade — the tool works and the negatives are booked honestly — but
they do mean the anchor's second acceptance item is unmet.

## 2. Checks performed (all re-run by the auditor, not taken from the report)

| # | check | command / method | result |
|---|---|---|---|
| C1 | digest binding to bytes | recomputed SHA-256 of the 6 cached pages in the subject cache | **6/6 match** record (bytes and digests) |
| C2 | numbers really in the page bytes | independent regex pass over `simman_summary.bin` / `t2015_onrt_gc.bin` (no subject parser) | **11/11 canonical tokens present**, record values equal them |
| C3 | scale rows self-consistent | publish scale row vs Lpp/Lwl pair in the same table | rel err 5.3e-7 / 6.3e-6 / 1.2e-7 (KCS/KVLCC2/ONRT) |
| C4 | ONRT cross-publisher scale | auditor parses T2015 `L WL (m) 3.147 154.0` itself | 154.0/3.147 = 48.935494 vs SIMMAN 48.9355 → **1.2e-7** |
| C5 | subject's own gate | `python tools/baseline_fleet/hunt_sources.py --verify-only` | rc=0, "0 problem(s), 11 sources" |
| C6 | subject's tests | `pytest tests/test_source_hunt.py -q` | **9 passed** (0.12 s) |
| C7 | clean re-run of the tool | `hunt_sources.py --out <tmp>` from a blank dir | rc=0, **10/10 network endpoints reproduce the record digests** (geometry archive only because `--fetch-geometry` was not passed) |
| C8 | live re-probe of every endpoint | `audit_independent.py --live` | **11/11 endpoints match** recorded status + byte count + full SHA-256 |
| C9 | geometry archive, full digest | independent `curl` of `onrt_wo_gap.zip`, 752 918 B | `bd747eff…c5c695685` — **exact match**; archive holds `onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs` (3 374 136 B) |
| C10 | benchmark referee | `SELECT … FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T1-SOURCE-HUNT'` | **0 rows** (the neighbouring case #52 belongs to `GAP-CAD-T1-BASELINE-FLEET`) |
| C11 | CI wiring claim | `grep` `.github/` | no dedicated `--verify-only` step, but `tests.yml:25 pytest --cov` runs `test_verify_only_is_green`, i.e. the claim holds indirectly |

Not run / not verifiable here: the `web.archive.org` route (the record reports three
30–90 s timeouts; the auditor did not re-attempt it — marked **not_run**, not "still dead");
whether SIMMAN's `onrt_wo_gap.zip` is *still* byte-identical to the T2015 copy (the record
compares against a 16-hex-digit digest stored in `fetch_sources.py`; the auditor confirmed
the T2015 side only).

## 3. Findings, severity-ranked

**F1 — MEDIUM (blocks L4): the benchmark referee is absent, and 4 of the 5 "0 %" rows are
identities, not validations.** `tensor_benchmark_cases` holds no row for anchor
`GAP-CAD-T1-SOURCE-HUNT` (C10), so no got-vs-tolerance comparison can be machine-adjudicated.
The report's §4 table is honest about this for two rows ("identity by construction") but the
column still reads "0 % … yes", which a reader can mistake for accuracy evidence. What the
round actually proves is *provenance and transcription* (the digits in the record are the
digits in the publisher bytes), which is a weaker and different claim than the L4 criterion.
Remediation: register the anchor's machine-checkable case(s) — e.g. `KCS Lpp`/`KVLCC2 ∇`
read vs canonical, `ONRT scale` vs the T2015 LWL pair (C4 gives 1.2e-7, well inside 3 %) —
and label identity rows as `transcription` rather than `0 %`.

**F2 — LOW/MEDIUM (reproducibility today, not durable): the committed record cannot
re-verify its own digests.** The page and archive bytes live in
`artifacts/baseline_fleet/raw/hunt_pages/` which `.gitignore:14` excludes, so a fresh clone
can check digest *syntax* only — which is exactly what `--verify-only` and
`test_verify_only_is_green` do. The gap is closed **empirically** by this audit (C7/C8/C9
reproduce every digest from the live publishers, full 64-hex), so the record is currently
true, but the evidence chain breaks the day a host changes bytes and archive.org is
unreachable. Remediation: commit the small page set (60 kB total, the 752 kB archive can stay
out) or a signed digest bundle, so `--verify-only` can hash bytes offline.

**F3 — LOW: `network_calls: false` means the record is a replay.** The committed record was
written from the cached bytes (`generated_utc` 12:41:02, cache mtime 12:39–12:40,
`fetched_this_run: false`), so it does not, by itself, prove the endpoints were contacted in
that run. Honest but easy to misread; the auditor's live re-probe (C8) is what turns it into
evidence. Remediation: store the fetch `checked_utc` per source (it does) **and** the
`fetched_this_run` flag in the report prose, or write the record on the fetching run.

**F4 — LOW: mass/volume convention under the ONRT "0 %" agreement.** The record stores
T2015's `Delta 72.6 kg / 8507 ton` as `displacement_m3`, and SIMMAN's column header is
`Displacement (m3)`. The two numbers agree (0 %) **only** under ρ = 1.000 t/m³; at seawater
ρ = 1.025 the same pair is **2.44 %** apart — which is precisely the "earlier 2.43 % gap"
the report attributes to "which column is read". Nothing is numerically wrong; the density
convention that makes the statement true is simply never stated. Remediation: rename the
field `displacement` with an explicit `unit`/`rho_basis` tag.

**F5 — NITPICK: `geometry_cross_check.identical` compares 16 hex digits.** The code compares
`fetch_sources.TARGETS["ONRT"][3]`, a 64-bit prefix, while the prose says "byte-identical".
The auditor upgraded this to a real 64-hex statement (C9), so the conclusion stands; the code
should compare full digests or say "digest-prefix match".

**F6 — POSITIVE: the tests assert invariants, not just execution.** Digest format per
reachable source, `page_sha256` binding for every extracted value, matched row text present,
literature-canonical values (230/320/52030/312622/8507/48.9355), the parser on *unseen*
markup (so a pass cannot come from re-reading the cache), and `main([...--verify-only]) == 0`.
Weakness: `test_published_principal_particulars` encodes the same source page's values, so it
is a transcription regression guard, not an independent accuracy check.

**F7 — POSITIVE: acceptance bookkeeping is honest.** 2 of 4 declared items MET (scratch-built
manifest; SIMMAN tables scraped), 2 NOT MET and flagged `WITHHELD` in both the report and
`open_questions` (DTC geometry; JBC scale). The DTC lead of the anchor is *falsified* with
digested negative evidence (T2015 `dtc.html` 404, geometry dir 403, uni-due 404/404 with an
identical generic-error digest, archive.org unreachable). No inflation found in either the
report or the record.

**F8 — NITPICK: no caller outside tests/docs.** `grep` finds no platform or solver side that
consumes `source_hunt.json`; the provenance is therefore not *used* by the fleet build yet
(`benchmarks/verified/baseline_fleet/geometry_diagnostics.json` still points at another work
tree, `/root/tp-agent/work/m20163/...`, as the report's §5 admits).

## 4. Overstated claims found

* "CI runs `--verify-only`" — true only through `tests.yml`'s pytest run (C11). Nitpick.
* "0 % error … `<3 %`: yes" for KCS Lpp, KVLCC2 ∇, ONRT ∇, ONRT ∇-cross — these are
  transcription identities; the table's own parenthetical admits it for one row only.
* "byte-identical" stated at a stronger level than the code checks (16 hex digits, F5).

No instance of a fabricated number, a missing artifact, or a claimed-but-absent run was found.

## 5. Disposition

**Approve with notes — L2 confirmed, no rework demanded of the subject; the anchor itself is
partially resolved.** Blocking for a higher grade (not for this verdict): F1 (register the
benchmark cases / relabel identity rows) and F2 (durable offline-verifiable digest bundle).
Both are small follow-ups that a new mission can pick up; neither invalidates the delivered
provenance record, which this audit reproduced end to end from the live publishers.

## 6. Machine artifacts

* `runs/m20425/audit_independent.py` — the auditor's independent verifier (own parser, own
  hashing, live re-probe), `--live` flag.
* `runs/m20425/audit_evidence.json` — full evidence object (digests, numbers, scale rows,
  units, subject test/gate output, clean re-run diff, live probe, checks matrix).
* `runs/m20425/fig_m20425_audit.png` — panel (a) live re-probe map, panel (b) identity rows
  vs the one real 2.44 % unit-convention sensitivity.