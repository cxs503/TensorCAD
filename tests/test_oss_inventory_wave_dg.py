"""Gates for the open-source DG inventory wave (FEAT-OSS-INVENTORY-DG).

These tests never re-run the network survey.  They verify that the *committed*
artifacts are complete, internally consistent and honest, so the inventory cannot
degrade back into a hand-written wish list:

* the census covers >=20 upstream repositories with a valid, unique tier;
* the deliberately unreachable negative-control entry is recorded as failed
  (failures are recorded, not silently dropped);
* every ``GAP-DG-*`` anchor cites >=2 reference implementations drawn from the
  census, carries a numeric acceptance criterion, and an external anchor that is
  not the target repository itself;
* every current-state claim is a grep probe satisfied at the pinned commit, and
  ``present`` probes really point at a file:line sample;
* the licence verdicts follow the declared policy (no strongly-copyleft repo may
  be advertised as code-portable);
* the markdown views are byte-identical to a fresh render of ``registry.json``
  (no hand editing), and the README quotes counts that exist in the registry;
* the wave engine's own ``verify`` reports no problem.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WAVE_PY = ROOT / "tools" / "oss_inventory_wave" / "wave.py"
DOMAIN = ROOT / "tools" / "oss_inventory_wave" / "domains" / "dg"
DOCS = ROOT / "docs" / "oss_inventory_dg"
REGISTRY = DOCS / "registry.json"
EVIDENCE = DOCS / "evidence" / "target_grep.json"


def _engine():
    spec = importlib.util.spec_from_file_location("oss_inventory_wave", WAVE_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def engine():
    return _engine()


@pytest.fixture(scope="module")
def cfg(engine):
    return engine.load_domain("dg")


@pytest.fixture(scope="module")
def registry():
    assert REGISTRY.exists(), f"missing {REGISTRY}"
    return json.loads(REGISTRY.read_text())


def _upstream(registry):
    return [r for r in registry["repos"] if r["role"] == "upstream"]


# --------------------------------------------------------------------------- #
# census completeness
# --------------------------------------------------------------------------- #
def test_census_covers_at_least_twenty_upstream_repos(registry):
    up = _upstream(registry)
    assert len(up) >= 20, f"only {len(up)} upstream repositories"


def test_every_repo_has_a_valid_tier_and_counts_agree(registry):
    tiers = registry["tiers"]
    for rec in registry["repos"]:
        assert rec["tier"] in tiers, f"{rec['id']}: invalid tier {rec['tier']}"
    counted = {}
    for rec in _upstream(registry):
        counted[rec["tier"]] = counted.get(rec["tier"], 0) + 1
    assert counted == registry["counts"]["by_tier"]
    assert sum(registry["counts"]["by_tier"].values()) == len(_upstream(registry))


def test_survey_status_is_recorded_for_every_pool_entry(registry):
    for rec in registry["repos"]:
        assert rec["status"] in {"ok", "failed"}, f"{rec['id']}: unknown status"
        if rec["status"] == "ok":
            assert (rec["metadata"] or {}).get("full_name"), rec["id"]
        else:
            assert rec["error"], f"{rec['id']}: failure without an error string"


def test_negative_control_failure_is_recorded(registry):
    ctl = [r for r in registry["repos"] if r["role"] == "negative_control"]
    assert ctl, "no negative control in the pool"
    for rec in ctl:
        assert rec["status"] == "failed", "the unreachable control must be recorded as failed"
        assert "404" in (rec["error"] or "") or "not_on_github" in (rec["error"] or "")
    assert registry["counts"]["negative_controls_recorded_failed"] == len(ctl)


def test_off_github_upstreams_are_marked_not_guessed(registry):
    for rec in registry["repos"]:
        if rec["status"] == "failed" and rec["role"] == "upstream":
            assert "not_on_github" in (rec["error"] or "") or "HTTP" in (rec["error"] or "")
            assert rec["metadata"] is None, f"{rec['id']}: failed entry must not carry metadata"


# --------------------------------------------------------------------------- #
# measured evidence, not memory
# --------------------------------------------------------------------------- #
def test_feature_coverage_is_measured_for_surveyed_repos(registry):
    surveyed = [r for r in _upstream(registry) if r["status"] == "ok"]
    assert surveyed, "no repo was surveyed"
    covered = [r for r in surveyed if sum(h["n"] for h in r["feature_paths"].values()) > 0]
    assert len(covered) >= 15, f"only {len(covered)} surveyed repos produced feature hits"
    feats = {f["id"] for f in registry["features"]}
    for rec in surveyed:
        assert set(rec["feature_paths"]) == feats, f"{rec['id']}: feature set mismatch"


def test_licence_verdicts_follow_the_policy(registry, cfg):
    policy = cfg["licence_policy"]
    for rec in _upstream(registry):
        verdict = rec["portability"]
        spdx = (rec["metadata"] or {}).get("license_spdx")
        if verdict == "port_ok":
            assert spdx in policy["port_ok"], f"{rec['id']}: {spdx} is not permissive"
        if spdx in policy["read_only_copyleft"]:
            assert verdict == "read_only_reimplement", f"{rec['id']}: copyleft marked {verdict}"
    assert sum(registry["counts"]["by_portability"].values()) == len(_upstream(registry))


# --------------------------------------------------------------------------- #
# gap anchors
# --------------------------------------------------------------------------- #
def test_gap_anchors_meet_the_wave_acceptance(registry):
    gaps = registry["gaps"]
    assert len(gaps) >= 4, f"only {len(gaps)} GAP-DG anchors"
    ids = [g["id"] for g in gaps]
    assert len(ids) == len(set(ids))
    for gap in gaps:
        assert re.fullmatch(r"GAP-DG-\d{3}", gap["id"]), gap["id"]
        assert gap["priority"] in {"P0", "P1", "P2"}
        assert gap["status"] in {"absent", "partial"}
        assert len(gap["refs"]) >= 2, f"{gap['id']}: needs >=2 reference implementations"
        assert gap["feature_points"], f"{gap['id']}: no function point"
        assert re.search(r"\d", gap["acceptance"]), f"{gap['id']}: acceptance has no number"
        assert re.search(r"(<=|>=|<|>|低于|不低于|误差|增长)", gap["acceptance"]), \
            f"{gap['id']}: acceptance has no comparator"
        assert registry["target_repo"].lower() not in gap["external_anchor"].lower(), \
            f"{gap['id']}: external anchor cites the target repo"
        assert len(gap["external_anchor"]) > 30, f"{gap['id']}: anchor is empty"
        assert gap["target_evidence"], f"{gap['id']}: no current-state probes"


def test_gap_reference_implementations_are_in_the_census(registry):
    pool = {r["id"]: r for r in registry["repos"]}
    for gap in registry["gaps"]:
        for ref in gap["refs"]:
            assert ref["id"] in pool, f"{gap['id']}: ref {ref['id']} not in the census"
        assert any(ref["tier"] == "T1" for ref in gap["refs"]), \
            f"{gap['id']}: no T1 reference implementation"


def test_current_state_claims_are_pinned_and_satisfied(registry):
    assert registry["counts"]["gaps_evidence_satisfied"] == len(registry["gaps"])
    for gap in registry["gaps"]:
        assert gap["evidence_satisfied"], f"{gap['id']}: probes unsatisfied"
        for probe in gap["target_evidence"]:
            assert probe["pattern"], gap["id"]
            assert probe["expect"] in {"absent", "present"}
            if probe["expect"] == "absent":
                assert probe["hits"] == 0, f"{gap['id']}: {probe['pattern']} unexpectedly present"
            else:
                assert probe["hits"] > 0, f"{gap['id']}: {probe['pattern']} not found"
                assert probe["samples"], f"{gap['id']}: present probe has no file:line sample"
                assert re.match(r"^\S+:\d+: ", probe["samples"][0]), gap["id"]


def test_committed_evidence_matches_the_registry_commit(registry):
    assert EVIDENCE.exists()
    blob = json.loads(EVIDENCE.read_text())
    assert blob["commit"] == registry["target_commit"]
    assert blob["target_repo"] == registry["target_repo"]
    assert {g["id"] for g in blob["gaps"]} == {g["id"] for g in registry["gaps"]}


# --------------------------------------------------------------------------- #
# views are generated, not hand written
# --------------------------------------------------------------------------- #
def test_markdown_views_are_reproducible_from_the_json(engine, registry):
    assert (DOCS / "registry.md").read_text() == engine.render_registry_md(registry)
    assert (DOCS / "gap-proposals.md").read_text() == engine.render_gaps_md(registry)


def test_readme_quotes_only_numbers_that_exist_in_the_registry(registry):
    readme = (DOCS / "README.md").read_text()
    counts = registry["counts"]
    assert f"{counts['upstream_repos']} 条上游" in readme
    assert f"`{registry['target_commit']}`" in readme
    for gap in registry["gaps"]:
        assert gap["id"] in readme, f"{gap['id']} missing from the README table"
    # no unresolved template placeholder may survive rendering
    assert "{{" not in readme and "}}" not in readme


def test_figures_exist_and_are_not_empty(registry):
    for name in ("coverage_heatmap.png", "tier_survey.png"):
        path = DOCS / "figures" / name
        assert path.exists(), name
        assert path.stat().st_size > 10_000, f"{name} looks empty"
        assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", name


# --------------------------------------------------------------------------- #
# the engine's own verifier
# --------------------------------------------------------------------------- #
def test_wave_verifier_reports_no_problem(engine, cfg):
    result = engine.verify(cfg)
    assert result["problems"] == [], result["problems"]