"""Every example script must be runnable from a clean checkout without installation.

README promises that the examples are executable evidence, not prose. A checkout that
has not run ``pip install -e .`` puts the repository root on ``sys.path`` for nothing:
``python examples/<name>/run_*.py`` starts with ``sys.path[0] == examples/<name>``, so a
top-level ``from tensorcad import ...`` raises ``ModuleNotFoundError`` before the script
computes anything. ``audit_report.md`` (m14491) caught three such scripts; a sweep of all
48 entry points at m16106 found 21.

Two checks:

* **static** -- for every ``examples/*/run_*.py`` the ``sys.path`` bootstrap appears
  *before* the first ``tensorcad`` import. Cheap, runs on every ``pytest`` invocation.
* **dynamic** -- the module body (imports + definitions, ``__main__`` guard not executed)
  is executed with ``PYTHONPATH`` removed and the working directory outside the
  checkout, in a subprocess. This is the check the old examples failed, and it is the
  only one that cannot be satisfied by a comment.
"""
from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted(REPO_ROOT.glob("examples/*/run_*.py"))

#: entry points that raised ``ModuleNotFoundError`` in a clean checkout before m16106
PREVIOUSLY_BROKEN = [
    "analytic_hydrostatics/run_analytic_check.py",
    "asymmetric_hull/run_example.py",
    "bonjean/run_example.py",
    "csg_boolean/run_example.py",
    "damage/run_example.py",
    "fairness/run_example.py",
    "freeboard_seam/run_example.py",
    "geometry_core/run_voxel_export.py",
    "gz/run_example.py",
    "hydrostatics_verification/run_verification.py",
    "intact_stability/run_example.py",
    "nurbs_degree_elevation/run_example.py",
    "nurbs_edge_alignment/run_example.py",
    "nurbs_piecewise_elevation/run_example.py",
    "nurbs_surface_degree_elevation/run_example.py",
    "planing/run_example.py",
    "righting_arm/run_example.py",
    "series60_verification/run_verification.py",
    "stability/run_example.py",
    "transom_stern/run_example.py",
    "wigley_stability/run_example.py",
]

PROBE = (
    "import runpy, sys\n"
    "runpy.run_path(sys.argv[1], run_name='tc_example_probe')\n"
)


def _first_tensorcad_import(text: str) -> int | None:
    for node in ast.parse(text).body:
        if isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] == "tensorcad":
            return node.lineno
        if isinstance(node, ast.Import):
            if any(alias.name.split(".")[0] == "tensorcad" for alias in node.names):
                return node.lineno
    return None


@pytest.mark.parametrize("script", EXAMPLES, ids=lambda p: p.parent.name + "/" + p.name)
def test_the_bootstrap_precedes_the_first_tensorcad_import(script: Path) -> None:
    text = script.read_text(encoding="utf-8")
    first_import = _first_tensorcad_import(text)
    assert first_import is not None, f"{script} no longer imports tensorcad; drop it from this gate"
    preamble = "\n".join(text.splitlines()[: first_import - 1])
    assert "sys.path" in preamble, (
        f"{script.relative_to(REPO_ROOT)} imports tensorcad before arranging sys.path; "
        "running it from a checkout without an installed package raises ModuleNotFoundError"
    )


@pytest.mark.parametrize("relative", PREVIOUSLY_BROKEN)
def test_the_module_body_imports_without_the_package_being_installed(
    relative: str, tmp_path: Path
) -> None:
    script = REPO_ROOT / "examples" / relative
    assert script.is_file(), f"missing example {relative}"
    env = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    completed = subprocess.run(
        [sys.executable, "-c", PROBE, str(script)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert "ModuleNotFoundError" not in completed.stderr, completed.stderr[-2000:]
    assert completed.returncode == 0, completed.stderr[-2000:]
