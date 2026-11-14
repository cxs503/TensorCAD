"""Tests for the waterline scale calibration (:mod:`tensorcad.fleet.scale_calibration`).

What is asserted, and why it is an invariant rather than an execution:

* **Recovery** — on a synthetic hull whose scale is *planted*, the two waterline
  conditions must return that scale.  A solver that executes but returns 742 (the
  first version of this module did, by locking onto a skeg edge at z = 22 mm) fails
  here.
* **Refusal** — a target ratio that no immersion satisfies must come back
  ``bracketed = False`` with ``root_m is None`` rather than a confident wrong root.
* **Patch accounting** — a patch that collapses to a point is reported as degenerate,
  never folded into the statistics.
* **The JBC record on disk** — the committed calibration report must state a scale,
  the two independent estimates, and the volume gate that failed; and the digest it
  pins must equal the digest of the file when the (git-ignored) artifact is present.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from tensorcad.fleet.scale_calibration import (
    SCHEMA,
    PatchGrids,
    solve_immersion,
    waterline_metrics,
    waterline_points,
)

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "tensorcad" / "fleet" / "data" / "scale_calibration.json"


def _synthetic_grids(
    scale_true: float,
    draft: float,
    *,
    waterline_length: float,
    waterline_beam: float,
    n_u: int = 41,
    n_v: int = 161,
    z_max: float = 0.62,
) -> PatchGrids:
    """A prismatic surface patch whose full-scale waterline is known by construction.

    The model frame is the full-scale hull divided by ``scale_true``: at immersion
    ``draft / scale_true`` the surface is ``waterline_length`` long and
    ``waterline_beam`` wide.  Below that immersion the length tapers so a cut near the
    keel is *not* a valid solution, which is the trap the ranking rule must survive.
    """
    immersion = draft / scale_true
    x = np.linspace(-0.5 * waterline_length, 0.5 * waterline_length, n_u)
    z = np.linspace(0.0, z_max, n_v)
    block = np.empty((n_u, n_v, 3))
    for j, zj in enumerate(z):
        taper = min(1.0, (zj / immersion) ** 0.5) if immersion > 0 else 0.0
        half = 0.5 * waterline_beam * taper
        for i, xi in enumerate(x):
            block[i, j] = (xi * taper, np.sign(xi) * half if xi else 0.0, zj)
    return PatchGrids(
        grids=[block],
        degenerate=[],
        nonzero=[0],
        factor=1.0,
        bbox=(
            (float(x.min()), -0.5 * waterline_beam, 0.0),
            (float(x.max()), 0.5 * waterline_beam, z_max),
        ),
    )


def test_solver_recovers_planted_scale_from_each_waterline_condition():
    scale_true, draft = 40.0, 16.5
    length_wl, beam_wl = 286.0, 45.0
    grids = _synthetic_grids(
        scale_true, draft, waterline_length=length_wl / scale_true, waterline_beam=beam_wl / scale_true
    )

    by_length = solve_immersion(grids, "length_m", length_wl / draft)
    by_beam = solve_immersion(grids, "beam_m", beam_wl / draft)

    assert by_length["bracketed"] and by_length["converged"]
    assert by_beam["bracketed"] and by_beam["converged"]
    scale_length = draft / by_length["root_m"]
    scale_beam = draft / by_beam["root_m"]
    assert scale_length == pytest.approx(scale_true, rel=0.01)
    assert scale_beam == pytest.approx(scale_true, rel=0.01)
    # Both conditions must pick the *design* waterline, not the tapered cut near the keel.
    assert by_length["root_m"] > 0.5 * (draft / scale_true)


def test_solver_refuses_an_unreachable_target():
    grids = _synthetic_grids(
        40.0, 16.5, waterline_length=286.0 / 40.0, waterline_beam=45.0 / 40.0
    )
    impossible = solve_immersion(grids, "length_m", 1.0e6)
    assert impossible["bracketed"] is False
    assert impossible["root_m"] is None


def test_waterline_metrics_read_the_known_waterplane():
    length_wl, beam_wl, scale = 286.0, 45.0, 40.0
    draft = 16.5
    grids = _synthetic_grids(
        scale, draft, waterline_length=length_wl / scale, waterline_beam=beam_wl / scale
    )
    metrics = waterline_metrics(grids, draft / scale)
    assert metrics["length_m"] == pytest.approx(length_wl / scale, rel=0.01)
    assert metrics["beam_m"] == pytest.approx(beam_wl / scale, rel=0.01)
    assert metrics["area_m2"] > 0.0
    points = waterline_points(grids, draft / scale)
    assert len(points) > 10


def test_record_on_disk_is_self_consistent():
    assert RECORD.exists(), "the JBC calibration record must be committed"
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    assert record["schema"] == SCHEMA
    calibration = record["calibration"]
    assert calibration["scale"] is not None
    # The two independent conditions must be reported, and their spread quoted.
    assert calibration["scale_by_beam"] is not None
    assert calibration["scale_by_length"] is not None
    spread = abs(calibration["scale_by_beam"] - calibration["scale_by_length"]) / calibration[
        "scale_by_beam"
    ]
    assert spread < 0.05, f"the two waterline conditions disagree by {100 * spread:.1f} %"
    # Beam closure is the sharpest single condition on this file and lands on 1:40.
    assert calibration["scale_by_beam"] == pytest.approx(40.0, abs=0.05)
    # The volume gate must stay open: the evaluated surface is less full than the
    # publisher's displacement, and hiding that would be the failure mode this whole
    # exercise exists to prevent.
    closure = calibration["volume_closure"]
    assert closure["within_tolerance"] is False
    assert closure["volume_error_percent"] < -3.0
    assert calibration["gate"] == "WITHHELD"
    assert calibration["patch_quality"]["patches_degenerate"] > 0


def test_record_digest_matches_the_artifact_when_present():
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    iges = ROOT / record["artifact"]["iges_path_in_worktree"]
    if not iges.exists():
        pytest.skip("JBC IGES artifact is git-ignored and absent in this checkout")
    digest = hashlib.sha256(iges.read_bytes()).hexdigest()
    assert digest == record["artifact"]["iges_sha256"]