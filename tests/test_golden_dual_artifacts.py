"""The committed dual-arm artefacts must be fresh, judged, and machine-free.

``artifacts/gap-cad-d2-golden/`` vendors the outputs of the dual-derived golden
suite (the TensorPapers ``papers/gap-cad-d2-golden-regression`` reference side)
so that the gate's own tests in ``test_geometry_golden_gate.py`` run in this
repository instead of depending on a sibling checkout that only CI has.

Three ways a vendored artefact set can lie, all guarded here:

* **Unjudged.** The files exist but no longer satisfy the acceptance reading
  (``tools/golden/gap_d2_assert_report.py``): a case stopped being measured, a
  detection rate fell, a WITHHELD leaked a verdict.
* **Stale.** The geometry code moved after the artefacts were generated.  The
  provenance sidecar records the sha256 of every ``tensorcad/*.py`` at
  generation time; this test recomputes it and fails on any drift, so a stale
  artefact set cannot keep the gate green.
* **Machine-bound.** An absolute worktree path of the generating machine is
  baked into the JSON, the same defect class
  ``test_golden_evidence_freshness.py`` polices in ``docs/golden/``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts" / "gap-cad-d2-golden"
SUITE = ART / "golden_suite.json"
MUTANTS = ART / "mutation_report.json"
PROV = ART / "provenance.json"

MACHINE_MARKERS = ("/root/", "/home/", "/tmp/", "/nfs/", "tp-agent")


def _judge() -> importlib.util.module_from_spec:
    path = ROOT / "tools" / "golden" / "gap_d2_assert_report.py"
    spec = importlib.util.spec_from_file_location("gap_d2_assert_report", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _package_sha256() -> str:
    """sha256 over (relative path, content) of every .py under tensorcad/."""
    digest = hashlib.sha256()
    for dirpath, _dirnames, filenames in os.walk(ROOT / "tensorcad"):
        for name in sorted(filenames):
            if not name.endswith(".py"):
                continue
            full = Path(dirpath) / name
            rel = full.relative_to(ROOT).as_posix()
            digest.update(rel.encode("utf-8"))
            digest.update(b"\0")
            digest.update(full.read_bytes())
            digest.update(b"\0")
    return digest.hexdigest()


def test_artefacts_are_committed() -> None:
    for path in (SUITE, MUTANTS, PROV):
        assert path.exists(), f"missing committed artefact: {path.relative_to(ROOT)}"


def test_committed_artefacts_pass_the_acceptance_judge() -> None:
    judge = _judge()
    suite = json.loads(SUITE.read_text(encoding="utf-8"))
    mutants = json.loads(MUTANTS.read_text(encoding="utf-8"))
    assert "8 cases PASS" in judge.check_suite(suite)
    assert "detection power 1.00" in judge.check_mutants(mutants)


def test_provenance_matches_the_current_tree() -> None:
    prov = json.loads(PROV.read_text(encoding="utf-8"))
    recorded = prov["tensorcad_sha256"]
    current = _package_sha256()
    assert recorded == current, (
        "the committed dual-arm artefacts were generated from a different "
        f"tensorcad tree (recorded {recorded[:12]}..., current {current[:12]}...): "
        "re-run the dual suite and re-vendor the artefacts"
    )


def test_artefacts_are_not_machine_bound() -> None:
    for path in (SUITE, MUTANTS, PROV):
        text = path.read_text(encoding="utf-8")
        for marker in MACHINE_MARKERS:
            assert marker not in text, f"{path.name} is machine-bound: {marker!r}"
