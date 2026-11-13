# Knowledge governance for TensorCAD: one ledger, two intakes, two outlets

`tensorcad/knowledge.py`, `knowledge/sediment_line.json`, `tools/sediment_guard.py`,
`tools/knowledge_export.py`, `tests/test_knowledge_governance.py`.
Generated numbers: [`docs/knowledge_governance_evidence.md`](knowledge_governance_evidence.md).

## 1. The gap this closes

TensorCAD already produces verified numbers everywhere — closed forms in
`wigley.py`, exact section areas in `box_gz.py`, benchmark gates under
`benchmarks/verified/`, audit measurements in `validation_report.json`. What it
did not have is an index that says **which claim rests on which number, and
whether that number is still true at HEAD**. Prose rots silently: a refactor can
move a value and every document quoting it stays green because nothing compares
them. Reviewing "is `docs/wigley.md` still right?" by hand does not scale to a
repository with dozens of `docs/*.md`.

The ledger therefore stores knowledge as *sediment*: each entry is a claim plus
every machine that can contradict it today.

| check | what it pins | what drift it catches |
|---|---|---|
| `probes` | a live recomputation through a public `tensorcad.*` callable, compared at a stated tolerance | the code changed (formula, default, branch) |
| `artifacts` | a value at a JSON pointer **and** the artifact's SHA-256 | the artifact was regenerated, or the schema moved |
| `docs` | a verbatim phrase the prose must still contain | the sentence quoting the number was edited away |

A check that cannot run is a failure (`BROKEN`), never a skip: a claim nobody can
re-derive any more is precisely what the guard exists to surface.

## 2. Method

1. **Passive intake** — `load_ledger` reads what verification already produced.
   The schema is closed (`SCHEMA = "tensorcad.sediment/1.0"`, four `kind`s) so a
   typo cannot invent a new claim class; duplicate ids, empty entries and
   malformed digests are rejected at load time.
2. **Active intake** — `harvest_candidates` walks `benchmarks/verified/**`,
   `validation_report.json` and the export manifests and lists numeric leaves no
   entry pins yet. That is the review queue: "we produced a number, nobody owns
   it yet".
3. **Guard** — `verify` runs probes, artifact checks and doc checks per entry and
   never short-circuits, so one run reports every broken claim.
   `tools/sediment_guard.py` exits 1 on any non-`PASS` row and writes the machine
   report (`artifacts/knowledge/sediment_guard.json`).
4. **Outlet 1, agents** — `export_agent_index` writes one JSONL record per claim
   (what is verified, to what tolerance, which probe/pointer reproduces it) so an
   agent can retrieve the verified surface without reading prose.
5. **Outlet 2, physical-AI training** — `export_training_corpus` writes a
   deterministic Wigley geometry → exact hydrostatics corpus plus a held-out eval
   split, with a SHA-256 manifest. Labels are the closed forms of `wigley.py`,
   recomputed independently in the tests from the formulas printed in
   `docs/wigley.md`.

## 3. Evidence (regenerated, never hand-copied)

Run the two commands below; the table in
[`docs/knowledge_governance_evidence.md`](knowledge_governance_evidence.md) is
written by the tool, not by a human:

```
python tools/sediment_guard.py --json artifacts/knowledge/sediment_guard.json
python tools/knowledge_export.py --report docs/knowledge_governance_evidence.md
```

Current state: 14 entries / 37 checks, all `PASS`; 72 corpus records
(51 train / 21 eval) with exact labels; the ledger digest is recorded in the
manifest so a downstream consumer can tell which claim set a corpus was built
from.

**Error judgement.** Probes are compared to their recorded values at the
tolerance stated per entry (`1e-12` for closed-form arms, `1e-15` for rational
coefficients, `1e-9` where the doc quotes ten significant figures). The observed
residuals are at IEEE-754 round-off — i.e. **< 3 % by roughly thirteen orders of
magnitude**, and the tests assert the entry-level predicates rather than a
rounded report. Mesh-level claims (default-grid Wigley volume, SUBOFF voxel
volume, the m241 audit numbers) are pinned by digest and by prose instead of
being recomputed: re-deriving them costs seconds-to-minutes, and pretending
otherwise would be the dishonest option.

## 4. Pitfalls, limits, scope

* **`DRIFT` is not automatically bad.** Regenerating a benchmark artifact is a
  legitimate act; the guard forces it to be a *deliberate* one — re-verify the
  value, then update the digest in the same commit. There is intentionally no
  `--accept` flag: a guard that can refresh itself is a rubber stamp.
* **Doc phrases are literal.** They are short, distinctive fragments (`` `4/9` ``,
  `` `5T/8` ``); rewording the sentence legitimately means updating the ledger.
  That edit is the review, not a nuisance.
* **Probes are restricted to `tensorcad.*` public callables.** The ledger must
  not become a code-execution gadget: `os.system:...` and `_private` targets are
  rejected when the ledger is parsed (`tests/test_knowledge_governance.py`).
* **The corpus is a seed, not a dataset.** It is one geometry family with exact
  labels — useful for verifying that a surrogate reproduces hydrostatics and for
  pinning tolerances, not a substitute for measured or solver-generated data.
  The split is a hash residue of the record id (nominal 1/5), deterministic on
  any machine; observed 21/72 here, which is expected variance for 72 records.
* **Not computed here:** any claim requiring a mesh, a solver run or an external
  reference table. Those live in the artifacts the ledger points at.

## 5. Adding an entry (the sediment ritual)

1. Pick an un-pinned number: `python tools/knowledge_export.py --candidates 10`.
2. Decide what it rests on: a closed form (probe), a checked-in artifact
   (pointer + `sha256sum`), a published comparison (doc phrase).
3. Append the entry to `knowledge/sediment_line.json`.
4. `python tools/sediment_guard.py --entry <new.id>` — it must be `PASS` before
   the number is quoted anywhere else.