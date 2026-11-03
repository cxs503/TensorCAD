"""Reader-side format contract: hostile manifests and the stdlib-only integrity path.

Two audit items were closed by code review only in the v0.4.0 audit
(``expand_audit_tensorcad.md`` §3.2 items 4 and 7):

* #4 "``verify`` needs torch" -- the claim is that the *integrity* checks
  (manifest, revision, weights hash, history, derived artifacts) are
  standard-library only and therefore keep working when the optional ``neural``
  extra is absent; only the live field query degrades to ``skipped``.
* #7 "``parse_manifest`` does not validate ``derived[]``" -- the reader must not
  be weaker than the writers it accepts files from, so a hand-edited manifest
  cannot smuggle ``../../etc/passwd`` past the parser.

Both are *measured* here instead of argued:

* the CLI cases run in a subprocess whose import system refuses ``torch``,
  ``numpy`` and ``safetensors``, which is exactly the "core install" environment
  the docstring in :mod:`tensorcad.cli_neural` promises to support;
* the parser cases are a hostile corpus driven straight at ``parse_manifest``.

Nothing in this module imports torch, numpy or safetensors, so it runs on the
minimal installation it is testing.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad.neural_geometry import (
    COMMITTED,
    Lineage,
    NeuralGeometryFormatError,
    ObjectStore,
    PayloadRef,
    TransactionRecord,
    build_manifest,
    derived_record,
    dumps_canonical,
    parse_manifest,
    sha256_bytes,
)

REPO_ROOT = Path(__file__).resolve().parent.parent

PAYLOAD = b"weights-v1"
MESH = b"o mesh"

#: Modules the core install does not have. ``safetensors`` is the one that is
#: actually missing from a core install in the wild; ``torch`` and ``numpy`` are
#: blocked as well so the subprocess cannot silently satisfy the claim through an
#: extra that happens to be present on a developer machine.
BANNED_MODULES = ("torch", "numpy", "safetensors")

BLOCKER_PREFIX = r"""
import sys


class _Blocker:
    def __init__(self, banned):
        self.banned = set(banned)

    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in self.banned:
            raise ImportError(f"BLOCKED: {name} is not installed in this environment")
        return None


sys.meta_path.insert(0, _Blocker(BANNED))
"""

BLOCKER_SCRIPT = BLOCKER_PREFIX + """
from tensorcad.cli_neural import main

raise SystemExit(main(sys.argv[1:]))
"""


def blocked_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    return env


# --------------------------------------------------------------------------- #
# fixtures: a valid bundle built with the standard library only
# --------------------------------------------------------------------------- #


def make_manifest(*, revision: int = 1, payload: bytes = PAYLOAD) -> dict[str, object]:
    return build_manifest(
        object_id="urn:tensorcad:test:sphere",
        revision=revision,
        weights=PayloadRef("weights.safetensors", sha256_bytes(payload), len(payload)),
        architecture_config={"coordinate_dim": 3, "output_dim": 1, "state_dim": 1,
                             "hidden_width": 8, "hidden_layers": 2},
        field_spec={"quantity": "truncated_signed_distance", "inside_sign": "negative",
                    "surface_isovalue": 0.0, "length_unit": "m", "truncation_distance": 0.15},
        coordinates={"dimension": 3, "frame": "object_local",
                     "bounds": {"min": [-0.6, -0.6, -0.6], "max": [0.6, 0.6, 0.6]}},
        state={"schema_version": "0.2.0", "field_state": {"radius": 0.5}, "controls": {}},
        training={"schema_version": "0.2.0", "steps": 1},
        validation={"schema_version": "0.2.0", "status": "passed", "revision": revision},
        lineage=Lineage(parent_revision=None, root_revision=1, operation="fit"),
        provenance={"producer": "test"},
    )


def bundle_files(payload: bytes = PAYLOAD) -> dict[str, bytes | str]:
    return {
        "weights.safetensors": payload,
        "state.json": dumps_canonical({"schema_version": "0.2.0", "field_state": {"radius": 0.5}}),
        "training.json": dumps_canonical({"schema_version": "0.2.0", "steps": 1}),
        "validation.json": dumps_canonical({"schema_version": "0.2.0", "status": "passed"}),
    }


def record(revision: int = 1) -> TransactionRecord:
    return TransactionRecord(
        revision=revision,
        parent_revision=None,
        status=COMMITTED,
        operation="fit",
        weights_sha256=sha256_bytes(PAYLOAD),
        recorded_utc="2026-11-02T00:00:00+00:00",
    )


def publish(tmp_path: Path, *, derived: list[dict[str, object]] | None = None) -> Path:
    """Publish revision 1 and return the bundle root."""
    root = tmp_path / "sphere.ngo"
    store = ObjectStore(root)
    files = bundle_files()
    manifest = make_manifest()
    if derived is not None:
        files["derived/surface.obj"] = MESH
        manifest["derived"] = derived
    store.publish(files=files, manifest=manifest, expected_parent=None, record=record())
    return root


def mesh_record(payload: bytes = MESH, path: str = "derived/surface.obj") -> dict[str, object]:
    return derived_record(
        kind="surface_mesh",
        fmt="obj",
        path=path,
        digest=sha256_bytes(payload),
        size_bytes=len(payload),
        metrics={"triangles": 4, "vertices": 6},
    )


def run_cli(*arguments: str) -> subprocess.CompletedProcess[str]:
    """Run ``tensorcad-neural`` with torch/numpy/safetensors imports blocked."""
    script = f"BANNED = {BANNED_MODULES!r}\n" + BLOCKER_SCRIPT
    return subprocess.run(
        [sys.executable, "-c", script, *arguments],
        capture_output=True,
        text=True,
        env=blocked_env(),
        cwd=str(REPO_ROOT),
        timeout=120,
    )


# --------------------------------------------------------------------------- #
# audit item #4: the integrity path must not need the neural extra
# --------------------------------------------------------------------------- #


def test_the_core_install_really_cannot_import_the_neural_extra(tmp_path: Path) -> None:
    """Guard the guard: if the blocker stopped working the next tests are void."""
    for name in BANNED_MODULES:
        probe = subprocess.run(
            [sys.executable, "-c", f"BANNED = {BANNED_MODULES!r}\n" + BLOCKER_PREFIX + f"\nimport {name}"],
            capture_output=True, text=True, env=blocked_env(), cwd=str(REPO_ROOT), timeout=120,
        )
        assert probe.returncode != 0, f"{name} was importable -- the blocked environment is a lie"
        assert "BLOCKED" in probe.stderr, probe.stderr
    # and the CLI itself imports cleanly inside that environment
    probe = subprocess.run(
        [sys.executable, "-c", f"BANNED = {BANNED_MODULES!r}\n" + BLOCKER_PREFIX
         + "\nimport tensorcad, tensorcad.cli_neural, tensorcad.neural_geometry\nprint('imported')"],
        capture_output=True, text=True, env=blocked_env(), cwd=str(REPO_ROOT), timeout=120,
    )
    assert probe.returncode == 0, probe.stderr
    assert "imported" in probe.stdout


def test_info_reports_a_bundle_without_the_neural_extra(tmp_path: Path) -> None:
    bundle = publish(tmp_path, derived=[mesh_record()])
    result = run_cli("info", "--object", str(bundle))
    assert result.returncode == 0, result.stdout + result.stderr
    assert "urn:tensorcad:test:sphere" in result.stdout
    assert "derived       surface_mesh" in result.stdout
    assert "Traceback" not in result.stderr


def test_verify_passes_integrity_checks_without_the_neural_extra(tmp_path: Path) -> None:
    bundle = publish(tmp_path, derived=[mesh_record()])
    result = run_cli("verify", "--json", "--object", str(bundle))
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "verified"
    checks = {check["name"]: check for check in payload["checks"]}
    assert checks["manifest_parses"]["passed"] is True
    assert checks["weights_hash_matches"]["passed"] is True
    assert checks["history_matches_revision"]["passed"] is True
    assert checks["derived_hash:derived/surface.obj"]["passed"] is True
    assert checks["derived_size:derived/surface.obj"]["passed"] is True
    # the only thing the core install cannot do: query the live field
    assert checks["field_query"]["severity"] == "skipped"
    assert checks["field_query"]["passed"] is True
    # nothing is downgraded to a warning: every integrity check genuinely ran
    assert payload["warnings"] == []


def test_verify_still_fails_on_a_bad_hash_without_the_neural_extra(tmp_path: Path) -> None:
    """A skipped field query must not launder a real integrity failure into a pass."""
    bundle = publish(tmp_path, derived=[mesh_record()])
    weights = bundle / "weights.safetensors"
    weights.write_bytes(PAYLOAD + b"tampered")
    result = run_cli("verify", "--json", "--object", str(bundle))
    assert result.returncode == 1, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload["status"] == "failed"
    failed = {check["name"] for check in payload["checks"] if not check["passed"]}
    assert "weights_hash_matches" in failed


def test_verify_reports_a_tampered_derived_artifact_without_the_neural_extra(
    tmp_path: Path,
) -> None:
    bundle = publish(tmp_path, derived=[mesh_record()])
    (bundle / "derived" / "surface.obj").write_bytes(b"o tampered mesh")
    result = run_cli("verify", "--json", "--object", str(bundle))
    assert result.returncode == 1, result.stdout + result.stderr
    checks = {check["name"]: check for check in json.loads(result.stdout)["checks"]}
    assert checks["derived_hash:derived/surface.obj"]["passed"] is False


def test_a_derived_symlink_leaving_the_bundle_is_not_hashed(tmp_path: Path) -> None:
    """A string check cannot see a symlink; the reader must re-check after resolving.

    The bundle ships ``derived/surface.obj`` as a symlink to a file *outside* the
    bundle, and the manifest records that outside file's real digest. A reader
    that resolves only at hash time would report the artifact as intact, which is
    exactly the "reader weaker than writer" failure the containment rules exist
    to prevent.
    """
    outside = tmp_path / "outside.obj"
    outside.write_bytes(b"o secret mesh")
    bundle = publish(tmp_path, derived=[mesh_record(outside.read_bytes())])
    os.remove(bundle / "derived" / "surface.obj")
    os.symlink(outside, bundle / "derived" / "surface.obj")
    # the digest in the manifest matches the file the symlink points at ...
    assert sha256_bytes(outside.read_bytes()) == mesh_record(outside.read_bytes())["sha256"]
    sha256_bytes(outside.read_bytes())  # ... so only a containment check catches it
    result = run_cli("verify", "--json", "--object", str(bundle))
    assert result.returncode == 1, result.stdout + result.stderr
    checks = {check["name"]: check for check in json.loads(result.stdout)["checks"]}
    assert checks["derived_containment:derived/surface.obj"]["passed"] is False
    assert "escapes bundle root" in checks["derived_containment:derived/surface.obj"]["detail"]
    # the escaped payload is never opened, so no hash verdict is published for it
    assert "derived_hash:derived/surface.obj" not in checks


# --------------------------------------------------------------------------- #
# audit item #7: the reader must not be weaker than the writer
# --------------------------------------------------------------------------- #


def mutating(**overrides: object) -> dict[str, object]:
    """A manifest whose ``derived`` list is the single record in ``overrides``."""
    manifest = make_manifest()
    manifest["derived"] = [dict(mesh_record(), **overrides)]
    return manifest


HOSTILE_DERIVED: tuple[tuple[str, object, str], ...] = (
    ("derived_is_not_a_list", {"kind": "surface_mesh"}, "derived must be a list"),
    ("record_is_not_an_object", ["surface_mesh"], "derived[0] must be an object"),
    ("kind_missing", {"kind": None}, "derived[0].kind must be a non-empty string"),
    ("kind_is_numeric", {"kind": 7}, "derived[0].kind must be a non-empty string"),
    ("format_missing", {"format": ""}, "derived[0].format must be a non-empty string"),
    ("path_traversal", {"path": "derived/../../etc/passwd"}, "unsafe bundle path in derived[0].path"),
    ("path_is_absolute", {"path": "/etc/passwd"}, "unsafe bundle path in derived[0].path"),
    ("path_is_dot_relative", {"path": "./derived/surface.obj"}, "unsafe bundle path in derived[0].path"),
    ("path_is_just_dot", {"path": "."}, "invalid bundle path in derived[0].path"),
    ("path_is_empty", {"path": "  "}, "must be a non-empty relative path"),
    ("path_is_numeric", {"path": 3}, "must be a non-empty relative path"),
    ("sha256_missing", {"sha256": None}, "derived[0].sha256 must be a lowercase SHA-256 hex digest"),
    ("sha256_is_uppercase", {"sha256": "A" * 64}, "lowercase SHA-256 hex digest"),
    ("sha256_is_short", {"sha256": "0" * 63}, "lowercase SHA-256 hex digest"),
    ("sha256_is_numeric", {"sha256": 0}, "lowercase SHA-256 hex digest"),
    ("bytes_negative", {"bytes": -1}, "derived[0].bytes must be a non-negative integer"),
    ("bytes_is_bool", {"bytes": True}, "derived[0].bytes must be a non-negative integer"),
    ("bytes_is_float", {"bytes": 1.5}, "derived[0].bytes must be a non-negative integer"),
    ("bytes_missing", {"bytes": None}, "derived[0].bytes must be a non-negative integer"),
    ("metrics_is_a_list", {"metrics": []}, "derived[0].metrics must be an object"),
)


@pytest.mark.parametrize(
    ("case", "overrides", "expected"), HOSTILE_DERIVED, ids=[c[0] for c in HOSTILE_DERIVED]
)
def test_hostile_derived_records_are_rejected(case: str, overrides: object, expected: str) -> None:
    manifest = make_manifest()
    if case == "derived_is_not_a_list" or case == "record_is_not_an_object":
        manifest["derived"] = overrides
    else:
        manifest["derived"] = [dict(mesh_record(), **overrides)]  # type: ignore[arg-type]
    with pytest.raises(NeuralGeometryFormatError, match=expected.replace("[", r"\[").replace("]", r"\]")):
        parse_manifest(manifest)


def test_a_valid_derived_record_parses_identically() -> None:
    """Writer and reader agree bit for bit: the reader may not reshape the record."""
    written = mesh_record()
    parsed = parse_manifest(mutating()).derived[0]
    assert parsed == written
    assert parsed["path"] == "derived/surface.obj"


def test_absent_and_null_derived_are_both_empty() -> None:
    """Forward compatibility: an older writer that has no derived artifacts still loads."""
    no_key = make_manifest()
    null_key = make_manifest()
    null_key["derived"] = None
    assert parse_manifest(no_key).derived == ()
    assert parse_manifest(null_key).derived == ()


def test_the_writer_rejects_escapes_before_they_reach_the_manifest() -> None:
    with pytest.raises(NeuralGeometryFormatError, match="unsafe bundle path"):
        derived_record(kind="surface_mesh", fmt="obj", path="../escape.obj",
                       digest=sha256_bytes(MESH), size_bytes=len(MESH), metrics={})


def test_duplicate_derived_paths_are_carried_verbatim() -> None:
    """The parser is a format gate, not a policy engine: duplicates surface downstream."""
    first = mesh_record()
    second = dict(mesh_record(b"other"), path="derived/surface.obj")
    manifest = make_manifest()
    manifest["derived"] = [first, second]
    paths = [item["path"] for item in parse_manifest(manifest).derived]
    assert paths == ["derived/surface.obj", "derived/surface.obj"]


# --------------------------------------------------------------------------- #
# robustness: no mutation may raise anything but the documented format error
# --------------------------------------------------------------------------- #


#: Number of deterministic mutations in the robustness corpus (evidence in the audit).
FUZZ_CASES = 400


def _mutations(manifest: dict[str, object], count: int, seed: int = 20261102) -> list[dict[str, object]]:
    """Deterministic hostile mutations: retype, delete and inject, one edit each."""
    import copy
    import random

    rng = random.Random(seed)
    replacements = [
        None, True, False, 0, -1, 1.5, "", "x", [], {}, ["x"], {"x": 1}, float("nan"), float("inf")
    ]
    mutations: list[dict[str, object]] = []
    keys = [*manifest, "derived"]
    while len(mutations) < count:
        candidate = copy.deepcopy(manifest)
        key = rng.choice(keys)
        action = rng.choice(("retype", "delete", "retype", "inject_derived"))
        if action == "delete":
            candidate.pop(key, None)
        elif action == "inject_derived":
            candidate["derived"] = [dict(mesh_record(), **{rng.choice(
                ("kind", "format", "path", "sha256", "bytes", "metrics")
            ): rng.choice(replacements)})]
        else:
            candidate[key] = rng.choice(replacements)
        mutations.append(candidate)
    return mutations


def test_no_manifest_mutation_raises_an_undocumented_exception() -> None:
    """A parser that raises KeyError/TypeError on hostile input is a crash, not a rejection."""
    manifest = make_manifest()
    manifest["derived"] = [mesh_record()]
    accepted = rejected = 0
    for candidate in _mutations(manifest, FUZZ_CASES):
        try:
            parse_manifest(candidate)
        except NeuralGeometryFormatError:
            rejected += 1
        else:
            accepted += 1
    # the property under test is "no undocumented exception"; the split between
    # accepted and rejected depends on the corpus and is reported, not asserted
    assert accepted + rejected == FUZZ_CASES
    assert rejected > 0 and accepted > 0, "the corpus is degenerate: one verdict only"


def test_the_fuzz_corpus_is_stable_across_runs() -> None:
    """The seed is fixed, so the evidence quoted in the audit is reproducible."""
    manifest = make_manifest()
    manifest["derived"] = [mesh_record()]
    first = _mutations(manifest, 5)
    second = _mutations(manifest, 5)
    # ``repr`` is used because the corpus deliberately contains NaN, and NaN != NaN
    assert [repr(item) for item in first] == [repr(item) for item in second]


def test_canonical_dumps_survive_a_round_trip() -> None:
    """Manifests written by the store parse back to the same object."""
    manifest = parse_manifest(mutating())
    reloaded = json.loads(dumps_canonical(mutating()))
    assert parse_manifest(reloaded).derived == manifest.derived
