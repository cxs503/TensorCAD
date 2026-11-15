"""The two arms of the GAP-CAD-D2 geometry gate must watch the same paths.

The gate installed here has two independent carriers: the CI workflow
``.github/workflows/geometry-golden.yml`` (authoritative, fails closed) and the
husky-style pre-commit step ``.githooks/pre-commit-geometry-golden.sh`` (local, ~35 s).
Both encode "which edits can move a golden number" as a hard-coded list, and the lists
are maintained in different files by different mechanisms -- a path filter in YAML and a
regex in shell.  Drift between them is silent in both directions: a path only CI watches
gives no local warning, a path only the hook watches is unenforced for anyone who did
not install it.  These tests make the drift loud, and pin the two acceptance properties
that are cheap to check without running the 35 s suite: documentation-only edits do not
trigger the gate, and the workflow asserts the suite's own report rather than only its
exit code.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKFLOW = os.path.join(REPO_ROOT, ".github", "workflows", "geometry-golden.yml")
HOOK = os.path.join(REPO_ROOT, ".githooks", "pre-commit-geometry-golden.sh")
PRE_COMMIT = os.path.join(REPO_ROOT, ".githooks", "pre-commit")

# Edits that cannot move a geometry number.  Documented in
# docs/gap-cad-d2-geometry-golden-gate.md as the "does not trigger" column.
DOC_ONLY = [
    "docs/baseline_fleet.md",
    "docs/gap-cad-d2-geometry-golden-gate.md",
    "README.md",
    "CHANGELOG.md",
    "docs/golden/mutation_audit.md",
]


def _workflow_paths() -> dict[str, list[str]]:
    """Return {"push": [...], "pull_request": [...]} from the workflow's on: block."""
    with open(WORKFLOW, encoding="utf-8") as handle:
        lines = handle.read().splitlines()

    sections: dict[str, list[str]] = {"push": [], "pull_request": []}
    current = None
    in_on_block = False
    for line in lines:
        if line.startswith("on:"):
            in_on_block = True
            continue
        if not in_on_block:
            continue
        if line and not line.startswith(" "):
            break  # left the on: block (e.g. "jobs:")
        match = re.match(r"^  (push|pull_request):$", line)
        if match:
            current = match.group(1)
            continue
        match = re.match(r"^      - \"(.+)\"$", line)
        if match and current:
            sections[current].append(match.group(1))
            continue
        if re.match(r"^  \w+:", line):
            current = None  # a different trigger (workflow_dispatch:, schedule:, ...)
    return sections


def _hook_watched() -> set[str]:
    """Expand the shell hook's WATCHED regex into the file list it covers."""
    with open(HOOK, encoding="utf-8") as handle:
        source = handle.read()
    match = re.search(r"^WATCHED='\^\((?P<body>.+)\)'$", source, re.MULTILINE)
    assert match, "WATCHED pattern not found in the pre-commit hook"
    body = match.group("body")
    prefix, suffix = "tensorcad/(", r")\.py"
    assert body.startswith(prefix) and body.endswith(suffix), (
        f"unrecognised WATCHED form: {body!r}"
    )
    alternatives = body[len(prefix):-len(suffix)].split("|")
    assert all(re.fullmatch(r"[A-Za-z0-9_]+", name) for name in alternatives), alternatives
    return {f"tensorcad/{name}.py" for name in alternatives}


@pytest.fixture(scope="module")
def workflow_paths():
    return _workflow_paths()


@pytest.fixture(scope="module")
def hook_watched():
    return _hook_watched()


def test_hook_and_workflow_watch_the_same_geometry_modules(workflow_paths, hook_watched):
    watched_by_workflow = {
        path for path in workflow_paths["push"] if path.startswith("tensorcad/")
    }
    assert watched_by_workflow == hook_watched, (
        "CI and the pre-commit hook disagree about which modules can move a golden "
        f"number; CI-only: {sorted(watched_by_workflow - hook_watched)}, "
        f"hook-only: {sorted(hook_watched - watched_by_workflow)}"
    )
    assert len(hook_watched) >= 11


def test_push_and_pull_request_filters_agree(workflow_paths):
    assert workflow_paths["push"], "push path filter is empty -- the gate would never run"
    assert workflow_paths["pull_request"], "pull_request path filter is empty"
    assert set(workflow_paths["push"]) == set(workflow_paths["pull_request"])


def test_gate_self_tests_its_own_workflow_file(workflow_paths):
    for trigger in ("push", "pull_request"):
        assert ".github/workflows/geometry-golden.yml" in workflow_paths[trigger]


def test_every_watched_path_exists(hook_watched):
    missing = sorted(path for path in hook_watched if not os.path.isfile(
        os.path.join(REPO_ROOT, path)))
    assert not missing, f"the gate watches paths that do not exist: {missing}"


def test_documentation_only_edits_are_not_watched(workflow_paths, hook_watched):
    watched = set(workflow_paths["push"]) | set(workflow_paths["pull_request"]) | hook_watched
    triggered = sorted(path for path in DOC_ONLY if path in watched)
    assert not triggered, f"documentation-only paths would trigger the gate: {triggered}"


def test_workflow_asserts_the_acceptance_criteria_not_just_the_exit_code():
    """The workflow must run the suite, run the sweep, and then *judge* the reports."""
    with open(WORKFLOW, encoding="utf-8") as handle:
        source = handle.read()
    assert "python suite.py" in source and "python mutants.py" in source
    assert "gap_d2_assert_report.py" in source
    # the criteria themselves live in the shared judge, not in YAML prose
    with open(os.path.join(REPO_ROOT, "tools", "golden", "gap_d2_assert_report.py"),
              encoding="utf-8") as handle:
        judge = handle.read()
    for needle in (
        "control_arm_all_pass",
        "strong_detection_rate",
        "gate_detection_rate",
        "withheld_refused",
        "unexpected_misses",
        "false_alarms",
    ):
        assert needle in judge, f"the shared judge never checks {needle}"


def test_pre_commit_chains_the_geometry_golden_step():
    with open(PRE_COMMIT, encoding="utf-8") as handle:
        source = handle.read()
    assert "pre-commit-geometry-golden.sh" in source
    assert os.access(HOOK, os.X_OK), "the chained hook must be executable"


def _run_hook(staged: list[str], golden_dir: str | None = None) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO_ROOT
    env["GOLDEN_STAGED_FILE"] = "/tmp/golden-gate-test-staged.txt"
    if golden_dir is not None:
        env["GOLDEN_DIR"] = golden_dir
    else:
        env.pop("GOLDEN_DIR", None)
    with open(env["GOLDEN_STAGED_FILE"], "w", encoding="utf-8") as handle:
        handle.write("\n".join(staged) + "\n")
    return subprocess.run(
        [HOOK], cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=180
    )


def test_hook_lets_a_documentation_only_commit_through_without_running_the_suite():
    result = _run_hook(DOC_ONLY, golden_dir="/nonexistent-suite")
    assert result.returncode == 0, result.stderr
    assert "geometry-affecting" not in result.stdout
    assert "suite.py" not in result.stdout + result.stderr


def test_hook_detects_a_geometry_edit():
    # A misconfigured GOLDEN_DIR exits non-zero *after* the path filter has fired, which
    # is how this test observes the trigger without paying for the 35 s suite.
    result = _run_hook(["tensorcad/mesh.py"], golden_dir="/nonexistent-suite")
    assert "geometry-affecting change detected" in result.stdout
    assert "tensorcad/mesh.py" in result.stdout
    assert result.returncode != 0


def test_hook_refuses_a_named_but_missing_suite_instead_of_downgrading():
    result = _run_hook(["tensorcad/hydrostatics.py"], golden_dir="/nonexistent-suite")
    assert result.returncode == 1
    assert "refusing to guess" in result.stderr


def _assert_report_module():
    path = os.path.join(REPO_ROOT, "tools", "golden", "gap_d2_assert_report.py")
    spec = importlib.util.spec_from_file_location("gap_d2_assert_report", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _artefact(name: str) -> dict:
    with open(os.path.join(REPO_ROOT, "artifacts", "gap-cad-d2-golden", name),
              encoding="utf-8") as handle:
        return json.load(handle)


def test_assert_report_accepts_the_committed_run():
    module = _assert_report_module()
    suite_line = module.check_suite(_artefact("golden_suite.json"))
    mutant_line = module.check_mutants(_artefact("mutation_report.json"))
    assert "8 cases PASS" in suite_line
    assert "detection power 1.00" in mutant_line


@pytest.mark.parametrize(
    "mutation, reason",
    [
        (lambda suite: suite["results"][3].__setitem__("verdict", "SKIP"),
         "not measured"),
        (lambda suite: suite["results"].pop(), "expected 8 golden cases"),
        (lambda suite: suite["results"][0]["findings"][0].__setitem__("verdict", "FAIL"),
         "unexpected finding verdicts"),
    ],
)
def test_assert_report_rejects_a_case_that_stopped_being_measured(mutation, reason):
    """The +3% beam regression exits suite.py 0 with two SKIP cases -- it must still fail."""
    module = _assert_report_module()
    suite = json.loads(json.dumps(_artefact("golden_suite.json")))
    mutation(suite)
    with pytest.raises(module.Rejected) as rejection:
        module.check_suite(suite)
    assert reason in str(rejection.value)


@pytest.mark.parametrize(
    "mutation, reason",
    [
        (lambda summary: summary.__setitem__("strong_detection_rate", 0.99), "must be 1.0"),
        (lambda summary: summary.__setitem__("control_arm_all_pass", False), "control_arm_all_pass"),
        (lambda summary: summary.__setitem__("withheld_refused", 15), "not refused"),
        (lambda summary: summary.__setitem__("false_alarms", 2), "false alarms"),
        (lambda summary: summary["control_arm"].__setitem__("wigley_coarse", "SKIP"), "control arm"),
    ],
)
def test_assert_report_rejects_lost_detection_power(mutation, reason):
    module = _assert_report_module()
    report = json.loads(json.dumps(_artefact("mutation_report.json")))
    mutation(report["summary"])
    with pytest.raises(module.Rejected) as rejection:
        module.check_mutants(report)
    assert reason in str(rejection.value)


def test_assert_report_runs_as_the_gate_invokes_it(tmp_path):
    """The CLI reads <suite-dir>/outputs/*.json -- the layout the run leaves behind."""
    suite_dir = tmp_path / "gap-cad-d2-golden-regression"
    outputs = suite_dir / "outputs"
    outputs.mkdir(parents=True)
    (outputs / "golden_suite.json").write_text(
        json.dumps(_artefact("golden_suite.json")), encoding="utf-8")
    (outputs / "mutation_report.json").write_text(
        json.dumps(_artefact("mutation_report.json")), encoding="utf-8")
    script = os.path.join(REPO_ROOT, "tools", "golden", "gap_d2_assert_report.py")
    accepted = subprocess.run([sys.executable, script, str(suite_dir)],
                              capture_output=True, text=True, timeout=120)
    assert accepted.returncode == 0, accepted.stderr
    assert "accepted:" in accepted.stdout

    (outputs / "golden_suite.json").write_text("{}", encoding="utf-8")
    rejected = subprocess.run([sys.executable, script, str(suite_dir)],
                              capture_output=True, text=True, timeout=120)
    assert rejected.returncode == 1
    assert rejected.stderr.startswith("REJECTED:")


def test_both_carriers_use_the_shared_report_judge():
    """CI and the hook must judge with one script, or the two arms can drift apart."""
    for path in (WORKFLOW, HOOK):
        with open(path, encoding="utf-8") as handle:
            source = handle.read()
        assert "gap_d2_assert_report.py" in source, path


def test_gate_document_is_generated_from_the_committed_artefacts():
    """docs/gap-cad-d2-geometry-golden-gate.md must match its own generator.

    The document's tolerance and detection tables are injected from
    artifacts/gap-cad-d2-golden/*.json, so a refreshed artefact set without a regenerated
    document would publish numbers no run produced.  ``--check`` is the same command a
    reviewer runs.
    """
    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "tools", "golden", "gap_d2_gate_docs.py"),
         "--check"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, (result.stdout + result.stderr).strip()


if __name__ == "__main__":  # pragma: no cover - manual probe
    sys.exit(pytest.main([__file__, "-q"]))