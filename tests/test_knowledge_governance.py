"""Tests for the sediment ledger, the guard, and the two knowledge outlets.

The guard is the deliverable, so most of these tests are about *failure*: a
stale artifact, an edited sentence, a pointer that moved, a probe outside the
package.  A green run therefore means the failure modes were exercised, not just
that the happy path was executed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.knowledge import (
    BROKEN,
    DRIFT,
    MISMATCH,
    MISSING,
    PASS,
    KnowledgeError,
    corpus_records,
    export_agent_index,
    export_training_corpus,
    harvest_candidates,
    load_ledger,
    parse_ledger,
    read_jsonl,
    resolve_pointer,
    sha256_file,
    verify,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = REPO_ROOT / "knowledge" / "sediment_line.json"


@pytest.fixture(scope="module")
def ledger():
    return load_ledger(LEDGER_PATH)


def _document() -> dict:
    return json.loads(LEDGER_PATH.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# ledger loading and validation
# --------------------------------------------------------------------------- #
def test_ledger_loads_with_unique_ids_and_checkable_entries(ledger):
    assert ledger.sha256 == sha256_file(LEDGER_PATH)
    assert len(ledger.entries) >= 10
    assert len(set(ledger.ids())) == len(ledger.entries)
    for entry in ledger.entries:
        assert entry.checks >= 1, f"{entry.id} has no check"
        assert entry.acceptance.strip(), f"{entry.id} has no acceptance line"


def test_every_entry_reproduces_today(ledger):
    report = verify(ledger, REPO_ROOT)
    failed = report.failed()
    assert not failed, json.dumps(failed, indent=2)[:4000]
    assert report.counts(PASS) == len(ledger.entries)
    assert report.as_dict()["counts"][DRIFT] == 0


def test_probe_outside_package_is_rejected():
    document = _document()
    document["entries"] = [
        {
            "id": "bad.namespace",
            "title": "t",
            "claim": "c",
            "kind": "analytic_identity",
            "acceptance": "a",
            "probes": [{"call": "os.system:echo", "expected": 0}],
        }
    ]
    with pytest.raises(KnowledgeError):
        parse_ledger(document)


def test_private_probe_target_is_rejected():
    document = _document()
    document["entries"] = [
        {
            "id": "bad.private",
            "title": "t",
            "claim": "c",
            "kind": "analytic_identity",
            "acceptance": "a",
            "probes": [{"call": "tensorcad.wigley:_resolve", "expected": 0}],
        }
    ]
    with pytest.raises(KnowledgeError):
        parse_ledger(document)


def test_entry_without_any_check_is_rejected():
    document = _document()
    document["entries"] = [
        {
            "id": "empty.claim",
            "title": "t",
            "claim": "c",
            "kind": "analytic_identity",
            "acceptance": "a",
        }
    ]
    with pytest.raises(KnowledgeError):
        parse_ledger(document)


def test_duplicate_ids_are_rejected(ledger):
    document = _document()
    document["entries"] = [document["entries"][0], document["entries"][0]]
    with pytest.raises(KnowledgeError):
        parse_ledger(document)


def test_bad_schema_is_rejected():
    with pytest.raises(KnowledgeError):
        parse_ledger({"schema": "something/9", "entries": []})


# --------------------------------------------------------------------------- #
# the four failure modes the guard must catch
# --------------------------------------------------------------------------- #
def _single_entry_ledger(entry_id: str):
    return type(load_ledger(LEDGER_PATH))(
        entries=tuple(e for e in load_ledger(LEDGER_PATH).entries if e.id == entry_id),
        sha256="synthetic",
    )


def test_artifact_digest_drift_is_caught(tmp_path):
    """A regenerated artifact must not be able to change a number unnoticed."""
    source = REPO_ROOT / "validation_report.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["checks"]["hull_volume_default"]["measured_volume_m3"] = 7.5  # value moved
    (tmp_path / "validation_report.json").write_text(json.dumps(document), encoding="utf-8")
    (tmp_path / "README.md").write_text(
        (REPO_ROOT / "README.md").read_text(encoding="utf-8"), encoding="utf-8"
    )
    report = verify(_single_entry_ledger("hull.analytic_volume_default"), tmp_path)
    statuses = {row["status"] for row in report.rows[0]["checks"]}
    assert DRIFT in statuses
    assert not report.ok


def test_artifact_value_mismatch_with_matching_digest_is_caught(tmp_path):
    """Even if somebody refreshes the digest, the value itself must still match."""
    source = REPO_ROOT / "validation_report.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["checks"]["hull_volume_default"]["measured_volume_m3"] = 7.5
    target = tmp_path / "validation_report.json"
    target.write_text(json.dumps(document), encoding="utf-8")
    (tmp_path / "README.md").write_text("nothing relevant here", encoding="utf-8")

    ledger = load_ledger(LEDGER_PATH)
    entry = ledger.by_id("hull.analytic_volume_default")
    patched_artifacts = tuple(
        type(ref)(
            path=ref.path,
            json_pointer=ref.json_pointer,
            expected=ref.expected,
            sha256=sha256_file(target) if ref.path.endswith("validation_report.json") else ref.sha256,
            rel_tol=ref.rel_tol,
            abs_tol=ref.abs_tol,
            note=ref.note,
        )
        for ref in entry.artifacts
    )
    patched = type(entry)(
        id=entry.id,
        title=entry.title,
        claim=entry.claim,
        kind=entry.kind,
        acceptance=entry.acceptance,
        probes=entry.probes,
        artifacts=patched_artifacts,
        docs=entry.docs,
        tags=entry.tags,
    )
    from tensorcad.knowledge import verify_entry

    row = verify_entry(patched, tmp_path)
    assert row["status"] == MISMATCH
    assert any("7.5" in check["detail"] for check in row["checks"] if check["status"] == MISMATCH)


def test_removed_doc_sentence_is_caught(tmp_path):
    """Prose drift: the sentence quoting the number is evidence, not decoration."""
    (tmp_path / "docs").mkdir()
    text = (REPO_ROOT / "docs" / "wigley.md").read_text(encoding="utf-8")
    assert "`4/9`" in text
    (tmp_path / "docs" / "wigley.md").write_text(text.replace("`4/9`", "`0.444`"), encoding="utf-8")
    report = verify(_single_entry_ledger("wigley.cb_rational"), tmp_path)
    row = report.rows[0]
    # probes still pass (the code is fine); the doc check is what fails
    assert row["status"] == MISSING
    assert any(check["check"] == "doc:docs/wigley.md" for check in row["checks"])


def test_moved_artifact_pointer_is_broken_not_skipped(tmp_path):
    document = {"checks": {"hull_volume_default": {}}}
    (tmp_path / "validation_report.json").write_text(json.dumps(document), encoding="utf-8")
    (tmp_path / "README.md").write_text("x", encoding="utf-8")
    report = verify(_single_entry_ledger("hull.analytic_volume_default"), tmp_path)
    statuses = {row["status"] for row in report.rows[0]["checks"]}
    assert BROKEN in statuses
    assert not report.ok


def test_raising_probe_is_a_failure(ledger):
    entry = ledger.by_id("wigley.volume_ittc_case")
    probe = type(entry.probes[0])(
        call=entry.probes[0].call,
        expected=entry.probes[0].expected,
        args=("not-a-number",),
        rel_tol=entry.probes[0].rel_tol,
    )
    patched = type(entry)(
        id=entry.id,
        title=entry.title,
        claim=entry.claim,
        kind=entry.kind,
        acceptance=entry.acceptance,
        probes=(probe,),
        artifacts=(),
        docs=(),
        tags=entry.tags,
    )
    from tensorcad.knowledge import verify_entry

    assert verify_entry(patched, REPO_ROOT)["status"] == BROKEN


def test_resolve_pointer_follows_lists_and_rejects_scalars():
    document = {"a": {"b": [10, 20]}}
    assert resolve_pointer(document, "/a/b/1") == 20
    assert resolve_pointer(document, "") is document
    with pytest.raises(KnowledgeError):
        resolve_pointer(document, "/a/b/9")
    with pytest.raises(KnowledgeError):
        resolve_pointer(document, "/a/b/0/x")


# --------------------------------------------------------------------------- #
# outlet 1: agent index
# --------------------------------------------------------------------------- #
def test_agent_index_covers_every_entry(ledger, tmp_path):
    path = tmp_path / "agent_index.jsonl"
    assert export_agent_index(ledger, path) == len(ledger.entries)
    records = read_jsonl(path)
    assert [record["id"] for record in records] == sorted(ledger.ids())
    for record in records:
        reproduce = record["reproduce"]
        assert reproduce["probes"] or reproduce["artifacts"] or reproduce["docs"]
        assert record["acceptance"]


# --------------------------------------------------------------------------- #
# outlet 2: training corpus + held-out eval split
# --------------------------------------------------------------------------- #
def test_corpus_targets_match_independent_closed_forms(ledger):
    """Recompute in the test from the formulas written in docs/wigley.md."""
    for record in corpus_records(ledger):
        inputs = record["inputs"]
        length, beam, draft = inputs["length_m"], inputs["beam_m"], inputs["draft_m"]
        n, m = inputs["longitudinal_exponent"], inputs["vertical_exponent"]
        cb = n * m / ((n + 1) * (m + 1))
        expected = {
            "block_coefficient": cb,
            "prismatic_coefficient": n / (n + 1),
            "midship_coefficient": m / (m + 1),
            "waterplane_coefficient": n / (n + 1),
            "displacement_volume_m3": cb * length * beam * draft,
            "waterplane_area_m2": length * beam * n / (n + 1),
            "vertical_centre_buoyancy_m": draft * (1.0 - (m + 1) / (2.0 * (m + 2))),
        }
        assert set(record["targets"]) == set(expected)
        for name, value in expected.items():
            got = record["targets"][name]["value"]
            assert math.isclose(got, value, rel_tol=1e-14, abs_tol=0.0), (record["id"], name, got, value)
            assert record["targets"][name]["rel_tol"] == 1e-9


def test_classic_wigley_record_is_present_and_exact(ledger):
    records = {record["id"]: record for record in corpus_records(ledger)}
    classic = records["wigley-L10-B1-T0.625-n2-m2"]
    assert math.isclose(classic["targets"]["block_coefficient"]["value"], 4 / 9, rel_tol=1e-15)
    assert math.isclose(classic["targets"]["displacement_volume_m3"]["value"], 25 / 9, rel_tol=1e-15)


def test_export_is_byte_deterministic(ledger, tmp_path):
    first = export_training_corpus(tmp_path / "a", ledger)
    second = export_training_corpus(tmp_path / "b", ledger)
    assert first == second
    for name in first["files"]:
        left = (tmp_path / "a" / name).read_bytes()
        right = (tmp_path / "b" / name).read_bytes()
        assert left == right, f"{name} is not deterministic"


def test_manifest_digests_match_the_files(ledger, tmp_path):
    manifest = export_training_corpus(tmp_path, ledger)
    for name, info in manifest["files"].items():
        path = tmp_path / name
        assert sha256_file(path) == info["sha256"], name
        assert path.stat().st_size == info["bytes"], name
    assert manifest["ledger"]["sha256"] == ledger.sha256
    assert json.loads((tmp_path / "manifest.json").read_text()) == manifest


def test_split_is_partition_deterministic_and_nonempty(ledger, tmp_path):
    manifest = export_training_corpus(tmp_path, ledger)
    train = read_jsonl(tmp_path / "training_corpus.jsonl")
    held_out = read_jsonl(tmp_path / "eval_set.jsonl")
    assert len(train) + len(held_out) == manifest["counts"]["records_total"]
    assert held_out, "the eval split must not be empty"
    assert not ({record["id"] for record in train} & {record["id"] for record in held_out})
    for record in train:
        assert record["split"] == "train"
    for record in held_out:
        assert record["split"] == "eval"
    # re-derive the residue rule independently
    import hashlib

    for record in held_out:
        residue = int(hashlib.sha256(record["id"].encode()).hexdigest()[:8], 16) % 5
        assert residue == 0


# --------------------------------------------------------------------------- #
# active intake
# --------------------------------------------------------------------------- #
def test_harvest_lists_unpinned_numbers_and_skips_pinned_ones(ledger):
    candidates = harvest_candidates(REPO_ROOT, ledger, limit=200)
    assert candidates, "the benchmark artifacts should hold at least one un-pinned number"
    pointed = {(item["path"], item["json_pointer"]) for item in candidates}
    for entry in ledger.entries:
        for ref in entry.artifacts:
            assert (ref.path, ref.json_pointer) not in pointed
    for path, pointer in pointed:
        assert (REPO_ROOT / path).is_file()
        # every harvested value must be reachable again
        document = json.loads((REPO_ROOT / path).read_text())
        resolve_pointer(document, pointer)