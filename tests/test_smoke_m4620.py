"""Pipeline E2E smoke test (task m4620).

Checks that the code path (a module at the repository root) and the test path
(``tests/``) of the TensorCAD checkout are both reachable and importable, i.e.
that an agent-generated marker module can be imported by a test without any
installation step. The marker module lives outside the ``tensorcad`` package, so
the repository root is prepended to ``sys.path`` explicitly: that makes the test
independent of whether it is launched as ``pytest`` or ``python -m pytest``.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import smoke_marker_m4620  # noqa: E402  (import after the sys.path bootstrap)


def test_marker_constant_matches_task_number() -> None:
    assert smoke_marker_m4620.PIPELINE_SMOKE_M4620 == "4620"


def test_trivial_arithmetic() -> None:
    assert 1 + 1 == 2
