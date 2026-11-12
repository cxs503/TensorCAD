"""Invariants for the repo-patrol (upstream diff-set + commit-theme evolution) tooling.

Two layers of assertions:

* **synthetic repos** built in ``tmp_path`` -- classification rules, risk rules,
  landing-mode detection, upstream sync/behind detection, unmerged-branch
  inventory and determinism are exercised on history we control, so the expected
  values are known independently of this repository's real history;
* **real artefacts** in ``artifacts/repo-patrol/upstream-diffset-v1`` -- the
  published report must not drift from ``summary.json``/``patrol.json`` and the
  manifest must pin every file by sha256.

Run: ``pytest -q tests/test_repo_patrol.py`` (offline; the live ``ls-remote``
probe is skipped everywhere here).
"""
# ruff: noqa: E501  (long lines are mostly Chinese report prose)
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools" / "repo_patrol"
ART = ROOT / "artifacts" / "repo-patrol" / "upstream-diffset-v1"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"repo_patrol_{name}", TOOLS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


patrol = _load("patrol")
build_report = _load("build_report")
CFG = patrol.load_config(TOOLS / "themes.json")

GIT_ENV = {
    "GIT_AUTHOR_NAME": "Patrol Test",
    "GIT_AUTHOR_EMAIL": "patrol@example.invalid",
    "GIT_COMMITTER_NAME": "Patrol Test",
    "GIT_COMMITTER_EMAIL": "patrol@example.invalid",
    "GIT_AUTHOR_DATE": "2026-11-01T00:00:00+0000",
    "GIT_COMMITTER_DATE": "2026-11-01T00:00:00+0000",
}


def git(repo: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, **GIT_ENV},
    )
    assert proc.returncode == 0, f"git {' '.join(args)}: {proc.stderr}"
    return proc.stdout


def write(repo: Path, rel: str, body: str) -> None:
    path = repo / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)


def init_repo(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    git(path, "init", "-q", "-b", "main")
    git(path, "config", "user.name", "Patrol Test")
    git(path, "config", "user.email", "patrol@example.invalid")
    return path


def commit(repo: Path, message: str, files: dict[str, str] | None = None, allow_empty: bool = False) -> str:
    for rel, body in (files or {}).items():
        write(repo, rel, body)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message, *(["--allow-empty"] if allow_empty else []))
    return git(repo, "rev-parse", "HEAD").strip()


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def synthetic(tmp_path_factory) -> dict:
    """A repo whose landing history deliberately contains every shape we classify."""
    base = tmp_path_factory.mktemp("patrol")
    origin = base / "origin.git"
    repo = init_repo(base / "work")

    c_first = commit(repo, "m100: nurbs degree elevation", {"tensorcad/nurbs.py": "x = 1\n",
                                                            "tests/test_nurbs.py": "def test_x(): pass\n"})
    # core change without any test change -> core_without_tests
    c_core = commit(repo, "m101: hull waterline offset", {"tensorcad/hull.py": "y = 2\n"})
    # squash landing: merge-looking subject, single parent
    c_squash = commit(repo, "Merge #2 (m200) via merge-master", {"docs/note.md": "note\n"})
    # empty squash landing: merge-looking subject, no tree change at all
    c_empty = commit(repo, "Merge #3 (m201) via merge-master", allow_empty=True)
    # large churn
    c_big = commit(repo, "m102: bulk tables", {"docs/big.md": "line\n" * 1200})
    # true merge
    git(repo, "checkout", "-q", "-b", "side")
    c_side = commit(repo, "m300: side branch work", {"tensorcad/side.py": "z = 3\n",
                                                     "tests/test_side.py": "def test_z(): pass\n"})
    git(repo, "checkout", "-q", "main")
    git(repo, "merge", "-q", "--no-ff", "-m", "Merge #4 (m300) via push-merge", "side")
    c_merge = git(repo, "rev-parse", "HEAD").strip()
    # low-signal subject carrying a gate-filtered declaration
    c_gate = commit(repo, "m301: matured mission deliverables（6 个文件被门禁滤除）",
                    {"docs/gate.md": "g\n"})

    subprocess.run(["git", "clone", "-q", "--bare", str(repo), str(origin)], check=True)
    git(repo, "remote", "add", "origin", str(origin))
    git(repo, "fetch", "-q", "origin")
    git(repo, "branch", "-q", "--set-upstream-to=origin/main", "main")

    # a merged-away branch and an unmerged one, both living on the remote
    git(repo, "checkout", "-q", "-b", "feature/open")
    c_open = commit(repo, "m400: open feature experiment", {"docs/open.md": "o\n"})
    git(repo, "push", "-q", "origin", "feature/open")
    git(repo, "checkout", "-q", "main")

    return {
        "origin": origin,
        "repo": repo,
        "commits": {
            "first": c_first,
            "core": c_core,
            "squash": c_squash,
            "empty": c_empty,
            "big": c_big,
            "side": c_side,
            "merge": c_merge,
            "gate": c_gate,
            "open": c_open,
        },
    }


@pytest.fixture(scope="module")
def payload(synthetic) -> dict:
    return patrol.collect(synthetic["repo"], "origin/main", window=10, network=False)


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
def test_area_classification_uses_config_prefixes():
    assert patrol.classify_area("tensorcad/nurbs.py", CFG) == "core"
    assert patrol.classify_area("tests/test_nurbs.py", CFG) == "tests"
    assert patrol.classify_area("docs/x.md", CFG) == "docs"
    assert patrol.classify_area("CHANGELOG.md", CFG) == "repo-meta"
    assert patrol.classify_area("unheard/of.txt", CFG) == "other"


def test_topics_come_from_subject_and_file_stems():
    topics = patrol.classify_topics("m1: nurbs degree elevation", ["tensorcad/nurbs.py"], CFG)
    assert topics[0] == "geometry-kernel"
    assert "unclassified" not in topics


def test_bookkeeping_files_do_not_tag_every_commit():
    """CHANGELOG is part of nearly every landing; it must not set the topic."""
    topics = patrol.classify_topics("m9: matured deliverables", ["CHANGELOG.md"], CFG)
    assert topics == ["unclassified"]


def test_test_file_stems_lose_the_test_prefix():
    """``tests/test_nurbs.py`` must contribute ``nurbs``, not ``test``."""
    topics = patrol.classify_topics("m2: wire it up", ["tests/test_nurbs.py"], CFG)
    assert "geometry-kernel" in topics


def test_gate_filtered_count_is_parsed_or_absent():
    assert patrol.parse_gate_filtered("m301: matured mission deliverables（6 个文件被门禁滤除）", CFG) == 6
    assert patrol.parse_gate_filtered("m301: nothing dropped", CFG) is None


# --------------------------------------------------------------------------- #
# landing modes / risks on the synthetic history
# --------------------------------------------------------------------------- #
def _by_subject(payload: dict, needle: str) -> dict:
    hits = [c for c in payload["window"]["commits"] if needle in c["subject"]]
    assert len(hits) == 1, f"{needle!r}: {len(hits)} matches"
    return hits[0]


def test_true_merge_vs_squash_is_decided_by_parents_not_subject(payload):
    merge = _by_subject(payload, "Merge #4")
    squash = _by_subject(payload, "Merge #2")
    assert merge["is_merge"] and merge["landing_mode"] == "true-merge"
    assert not squash["is_merge"] and squash["landing_mode"] == "squash"
    assert merge["packaged_count"] >= 1
    assert "squash_landing_no_ancestry" in squash["risks"]
    assert "squash_landing_no_ancestry" not in merge["risks"]


def test_empty_landing_is_flagged(payload):
    empty = _by_subject(payload, "Merge #3")
    assert empty["files_changed"] == 0 and empty["insertions"] == 0
    assert {"empty_landing", "empty_delta"} <= set(empty["risks"])
    assert payload["window"]["totals"]["empty_landings"] == 1


def test_core_change_without_tests_is_flagged(payload):
    core = _by_subject(payload, "m101")
    assert "core_without_tests" in core["risks"]
    ok = _by_subject(payload, "m100")
    assert "core_without_tests" not in ok["risks"]


def test_large_churn_threshold_comes_from_config(payload):
    big = _by_subject(payload, "m102")
    assert big["insertions"] >= CFG["risk_rules"]["large_churn_lines"]
    assert "large_churn" in big["risks"]


def test_gate_filtered_and_low_signal_risks(payload):
    gate = _by_subject(payload, "m301")
    assert gate["gate_filtered_files"] == 6
    assert "gate_filtered_delivery" in gate["risks"]
    assert "low_signal_subject" in gate["risks"]


# --------------------------------------------------------------------------- #
# upstream lens
# --------------------------------------------------------------------------- #
def test_upstream_in_sync_when_local_head_is_the_fetched_ref(payload):
    up = payload["upstream"]
    assert up["verdict"] == "IN_SYNC"
    assert (up["ahead"], up["behind"]) == (0, 0)
    assert up["diffset_files"] == []
    assert up["live_probe"] == "skipped"


def test_offline_collection_marks_sync_verdict_withheld(payload):
    assert payload["evidence_grades"]["upstream_sync_verdict"] == "WITHHELD"


def test_behind_upstream_is_detected(synthetic):
    origin, repo = synthetic["origin"], synthetic["repo"]
    other = origin.parent / "peer"
    subprocess.run(["git", "clone", "-q", str(origin), str(other)], check=True)
    git(other, "config", "user.name", "Peer")
    git(other, "config", "user.email", "peer@example.invalid")
    commit(other, "peer: land something", {"docs/peer.md": "p\n"})
    git(other, "push", "-q", "origin", "main")
    git(repo, "fetch", "-q", "origin")
    try:
        up = patrol.collect(repo, "origin/main", window=3, network=False)
        assert up["upstream"]["behind"] == 1
        assert up["upstream"]["ahead"] == 0
        assert up["upstream"]["verdict"] == "BEHIND"
        assert up["upstream"]["diffset_files"] == ["docs/peer.md"]
        assert up["upstream"]["diffset_stat"]["insertions"] == 1
    finally:
        # restore the fixture's world for later tests
        git(other, "reset", "-q", "--hard", "HEAD~1")
        git(other, "push", "-q", "--force", "origin", "main")
        git(repo, "fetch", "-q", "origin")


def test_unmerged_branches_are_inventoried(payload):
    branches = {b["branch"]: b for b in payload["unmerged_branches"]}
    assert "origin/feature/open" in branches
    b = branches["origin/feature/open"]
    assert b["commits_ahead"] == 1
    assert b["missions"] == ["m400"]
    assert b["files_changed"] == 1


# --------------------------------------------------------------------------- #
# self-checks / determinism
# --------------------------------------------------------------------------- #
def test_all_crosschecks_agree(payload):
    checks = payload["crosschecks"]
    assert len(checks) >= 3
    for c in checks:
        assert c["match"], c
        assert c["delta"] == 0


def test_payload_is_deterministic(synthetic, payload):
    again = patrol.collect(synthetic["repo"], "origin/main", window=10, network=False)
    assert patrol.canonical(again) == patrol.canonical(payload)
    assert patrol.summary_of(again)["data_sha256"] == patrol.summary_of(payload)["data_sha256"]


def test_churn_span_is_reported_as_amplification_not_equality(payload):
    span = payload["churn_span"]
    assert span["per_commit_lines"] >= span["span_lines"]
    assert span["amplification"] is not None and span["amplification"] >= 1.0


# --------------------------------------------------------------------------- #
# artefacts + report of the real run
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not ART.is_dir(), reason="patrol artefacts not generated")
def test_manifest_pins_every_artifact():
    manifest = json.loads((ART / "manifest.json").read_text())
    assert manifest["files"], "manifest is empty"
    for name, sha in manifest["files"].items():
        assert hashlib.sha256((ART / name).read_bytes()).hexdigest() == sha, f"{name}: hash drift"
    assert len(manifest["data_sha256"]) == 64


@pytest.mark.skipif(not ART.is_dir(), reason="patrol artefacts not generated")
def test_tampered_artifacts_are_rejected(tmp_path):
    for name in ("patrol.json", "summary.json", "manifest.json"):
        (tmp_path / name).write_bytes((ART / name).read_bytes())
    data = json.loads((tmp_path / "patrol.json").read_text())
    data["window"]["totals"]["commits"] = 999
    (tmp_path / "patrol.json").write_text(json.dumps(data))
    with pytest.raises(SystemExit):
        build_report.load_bundle(tmp_path)


@pytest.mark.skipif(not ART.is_dir(), reason="patrol artefacts not generated")
def test_real_run_numbers_are_consistent_between_artifacts():
    payload_real = json.loads((ART / "patrol.json").read_text())
    summary = json.loads((ART / "summary.json").read_text())
    assert summary["data_sha256"] == patrol.summary_of(payload_real)["data_sha256"]
    assert summary["commits"] == len(payload_real["window"]["commits"])
    assert summary["merge_commits"] == sum(1 for c in payload_real["window"]["commits"] if c["is_merge"])
    assert summary["squash_landings"] + summary["merge_commits"] == summary["commits"]
    assert all(c["match"] for c in payload_real["crosschecks"])
    # the upstream verdict must be a declared value, never silently missing
    assert payload_real["upstream"]["verdict"] in {"NO_UPSTREAM", "IN_SYNC", "BEHIND", "AHEAD_ONLY", "DIVERGED"}


@pytest.mark.skipif(not ART.is_dir(), reason="patrol artefacts not generated")
def test_report_cannot_drift_from_artifacts():
    payload_real, summary, manifest = build_report.load_bundle(ART)
    body = build_report.render(payload_real, summary, manifest, "docs/repo_patrol", payload_real["repo"]["head_date"])
    report = ROOT / "docs/repo_patrol" / f"REPORT-{payload_real['repo']['head_date']}.md"
    assert report.is_file(), report
    published = report.read_text()
    for token in (
        summary["upstream_verdict"],
        f"真合并 {summary['merge_commits']} 个",
        f"squash 落地 {summary['squash_landings']} 个",
        f"{sum(1 for c in payload_real['crosschecks'] if c['match'])}/{len(payload_real['crosschecks'])}",
        "是否 <3%",
    ):
        assert token in published, f"report does not carry {token!r}"
    for commit in payload_real["window"]["commits"]:
        assert commit["short"] in published, f"commit {commit['short']} missing from report table"
    # re-rendering from the same artifacts must reproduce the published report
    assert body == published
