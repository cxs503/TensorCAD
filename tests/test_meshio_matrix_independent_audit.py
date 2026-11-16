"""Tests for ``tools/meshio_matrix_independent_audit.py`` (GAP-CAD-B2).

The audit tool is the independent re-check of the 6-format x 3-geometry meshio
exchange matrix.  Two properties are pinned here:

* **all-green** -- every audited cell (36 = 18 combos x 2 directions) passes the
  volume / topology / orientation / closedness / spec gates, i.e. the green
  table in ``docs/meshio_format_matrix_green_table_m21249.md`` is reproducible.
* **portability** -- the tool must not assume its artefacts land inside the
  checkout.  ``--out`` / ``--fig`` may point anywhere (e.g. a scratch dir under
  ``/tmp``); the JSON/manifest must record a usable path instead of raising
  ``ValueError`` from ``Path.relative_to``.  This is the completed form of the
  portability fix that shipped as PR#737 (m21249).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
TOOL = REPO_ROOT / "tools" / "meshio_matrix_independent_audit.py"


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "meshio_matrix_independent_audit", TOOL
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_rel_helper_inside_and_outside_repo():
    """``_rel`` is repo-relative inside the checkout, absolute outside it."""
    mod = _load_tool()
    root = mod.ROOT.resolve()

    inside = root / "tools" / "meshio_matrix_independent_audit.py"
    assert mod._rel(inside) == "tools/meshio_matrix_independent_audit.py"

    # A path that is not under the checkout must not raise; it degrades to the
    # absolute path (the case that used to crash with ValueError).
    outside = Path("/definitely/not/under/the/repo/out.stl")
    assert mod._rel(outside) == str(outside.resolve())


def test_audit_all_green_and_portable(tmp_path: Path):
    """Full 36-cell audit: all green, and portable to an out-of-repo --out."""
    mod = _load_tool()

    # tmp_path lives outside the checkout on this host; that is exactly the
    # case the portability fix targets.  Guard the assumption explicitly.
    assert not tmp_path.resolve().is_relative_to(REPO_ROOT.resolve())

    out_dir = tmp_path / "matrix18"
    rc = mod.main(["--out", str(out_dir), "--fig", str(out_dir / "figures")])
    assert rc == 0, "audit reported failing cells"

    payload = json.loads((out_dir / "matrix18.json").read_text())
    rows = payload["rows"]
    assert len(rows) == 36, "expected 18 combos x 2 directions"
    assert all(r["green"] for r in rows), [
        (r["direction"], r["format"], r["geometry"])
        for r in rows if not r["green"]
    ]

    # Artefacts the green table indexes must exist.
    for name in ("matrix18.json", "matrix.md", "manifest.json"):
        assert (out_dir / name).is_file(), name
    assert (out_dir / "files").is_dir()

    # Manifest sha256 must match the bytes on disk (re-checkable without rerun).
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest, "empty manifest"
    for rel, digest in manifest.items():
        target = Path(rel)
        if not target.is_absolute():
            target = REPO_ROOT / target
        assert target.is_file(), f"manifest entry missing on disk: {rel}"
        assert _sha256(target) == digest, f"sha256 mismatch for {rel}"

    # Every recorded artefact path must be resolvable (the portability
    # guarantee): inside the repo it is repo-relative, outside it is absolute.
    for rel in manifest:
        target = Path(rel)
        if not target.is_absolute():
            target = REPO_ROOT / target
        assert target.is_file(), f"unresolvable manifest path: {rel}"


def test_audit_exit_code_reflects_failures(monkeypatch, tmp_path):
    """A non-green cell must surface as a non-zero exit code."""
    mod = _load_tool()
    real_audit_out = mod.audit_out

    def failing_audit_out(fmt, geom_name, files_dir):
        record = real_audit_out(fmt, geom_name, files_dir)
        record["green"] = False  # force a red cell
        return record

    monkeypatch.setattr(mod, "audit_out", failing_audit_out)
    rc = mod.main(["--out", str(tmp_path / "forced_fail"),
                   "--fig", str(tmp_path / "forced_fail_fig")])
    assert rc == 1
