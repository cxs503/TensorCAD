"""The Gmsh *source* oracle build script must fail honestly, and must parse.

`docs/audit_m22616_gmsh_source_oracle.md` records that the script shipped with a `sh`
syntax error, so the oracle it is supposed to build never existed. Nothing on this host can
build the oracle (Gmsh's element factory closure is absent), so the two properties that
*are* testable are pinned here:

* the script parses under `sh -n` -- a parse error is invisible to every reader who never
  runs it, which is exactly how the original defect survived;
* with no Gmsh source tree it exits 2 (not 0, not 3) and says so, so an absent oracle can
  never be mistaken for a built one.

Plus the audit's own central fact: the harness never touches `MElementFactory::create`,
which is what makes the link failure a build-isolation problem rather than a missing piece
of the measurement.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ORACLES = REPO / "tools" / "mesh_quality_oracles"
SCRIPT = ORACLES / "build_gmsh_source_oracle.sh"


def test_build_script_parses_under_posix_sh():
    sh = shutil.which("sh")
    if sh is None:
        pytest.skip("no /bin/sh on this host")
    completed = subprocess.run([sh, "-n", str(SCRIPT)], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


def test_absent_source_tree_exits_2(tmp_path):
    sh = shutil.which("sh")
    if sh is None:
        pytest.skip("no /bin/sh on this host")
    completed = subprocess.run(
        [sh, str(SCRIPT), str(tmp_path / "not-a-gmsh-tree"), str(tmp_path / "out")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 2, (completed.returncode, completed.stdout, completed.stderr)
    assert "gmsh source tree not found" in completed.stderr
    assert not (tmp_path / "out").exists()


def test_harness_does_not_reach_the_element_factory():
    for name in ("gmsh_source_oracle.cpp", "gmsh_source_shims.cpp"):
        text = (ORACLES / name).read_text()
        assert "MElementFactory" not in text.replace("MElementFactoryNotOnQualityPath", ""), (
            f"{name} references the element factory; the link failure would then be part of "
            "the measurement path, not a build-isolation problem"
        )