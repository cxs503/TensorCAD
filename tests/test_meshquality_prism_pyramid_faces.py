"""Interior-face non-orthogonality for prisms and pyramids, and the two face tables.

``docs/mesh_quality_prism_pyramid.md`` §6 used to record "non-orthogonality is not
extended to prisms/pyramids" as a known gap, and ``door.admit`` therefore skipped the
measure for a prism-only or pyramid-only mesh -- a *silent* pass, which the same document
calls "worse than no gate". This file pins the closure of that gap:

* the face tables themselves (:data:`CELL_FACES`), by two properties that a mistyped
  winding fails while still producing plausible numbers -- every face normal points away
  from the cell centroid, and the divergence-theorem volume of the table is positive and
  equal to the metric module's independent tet-split volume;
* the closed form on two conforming stacks, where the exact answer is ``atan(shear)`` for
  a mesh sheared in x proportional to z (the same closed form the hexahedral stack uses),
  checked against a second route that computes the angle with ``atan2(|n x d|, n . d)``
  and therefore shares no normalisation step with the library's ``acos(n_hat . d_hat)``;
* the gate and the two callers: ``door.admit`` now emits ``non_orthogonality_prism`` /
  ``non_orthogonality_pyramid``, the CLI's ``--max-non-orth`` accepts the two kinds, and a
  mesh that is worse than the limit fails with the kind named.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad.meshquality import NODES_PER_CELL, non_orthogonality, prism_metrics, pyramid_metrics
from tensorcad.meshquality.door import CELL_KINDS_WITH_INTERNAL_FACES, admit
from tensorcad.meshquality.reference import reference_non_orthogonality
from tensorcad.meshquality.topology import (
    CELL_FACES,
    CELL_NODES,
    _centroid,
    _newell_normal,
    _sub,
    gate_non_orthogonality,
)

#: A right prism translated off the origin: a table wound the wrong way would still give a
#: positive volume if the cell sat at the origin, so the nodes are deliberately elsewhere.
PRISM = ((1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (1.0, 2.0, 1.0), (1.0, 1.0, 3.0), (2.0, 1.0, 3.0), (1.0, 2.0, 3.0))
#: A right square pyramid, likewise translated: base 1, height 2, volume 2/3.
PYRAMID = ((1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (2.0, 2.0, 1.0), (1.0, 2.0, 1.0), (1.5, 1.5, 3.0))

VOLUME_KINDS = ("prism", "pyramid")


def _divergence_volume(points, faces) -> float:
    """Signed volume from the divergence theorem, fan-triangulating each face."""
    total = 0.0
    for face in faces:
        a = points[face[0]]
        for i in range(1, len(face) - 1):
            b, c = points[face[i]], points[face[i + 1]]
            total += (
                a[0] * (b[1] * c[2] - b[2] * c[1])
                - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
            )
    return total / 6.0


def prism_stack(layers: int = 3, shear: float = 0.0, height: float = 1.0):
    """A conforming column of right prisms, sheared in x proportional to z."""
    points: list[tuple[float, float, float]] = []
    for level in range(layers + 1):
        z = level * height
        x = shear * z
        points += [(x, 0.0, z), (x + 1.0, 0.0, z), (x, 1.0, z)]
    cells = [
        (3 * k, 3 * k + 1, 3 * k + 2, 3 * k + 3, 3 * k + 4, 3 * k + 5) for k in range(layers)
    ]
    return points, cells


def pyramid_pair(shear: float = 0.0):
    """Two pyramids sharing one base quad -- apex up and apex down, a hex/tet transition."""
    base = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    apex_up, apex_down = (0.5, 0.5, 1.0), (0.5, 0.5, -1.0)
    points = [(x + shear * z, y, z) for x, y, z in base + [apex_up, apex_down]]
    return points, [(0, 1, 2, 3, 4), (0, 1, 2, 3, 5)]


# --------------------------------------------------------------------------------------
# the tables
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("kind", VOLUME_KINDS)
def test_face_table_winds_every_face_outward(kind: str) -> None:
    points = PRISM if kind == "prism" else PYRAMID
    centre = _centroid(points, range(len(points)))
    for face in CELL_FACES[kind]:
        assert len(set(face)) == len(face), f"{kind} face {face} repeats a node"
        normal = _newell_normal([points[i] for i in face])
        outward = _sub(_centroid(points, face), centre)
        assert sum(normal[i] * outward[i] for i in range(3)) > 0.0, (
            f"{kind} face {face} is wound inward"
        )


@pytest.mark.parametrize("kind", VOLUME_KINDS)
def test_face_table_volume_matches_metrics_module(kind: str) -> None:
    """The divergence-theorem volume of the table equals the metric's tet-split volume."""
    points = PRISM if kind == "prism" else PYRAMID
    from_table = _divergence_volume(points, CELL_FACES[kind])
    metrics = prism_metrics if kind == "prism" else pyramid_metrics
    cell = tuple(range(CELL_NODES[kind]))
    from_metrics = metrics(points, [cell])["volume"][0]
    # prism: triangle of area 1/2 x height 2; pyramid: unit square base x height 2 / 3
    assert from_table == pytest.approx(1.0 if kind == "prism" else 2.0 / 3.0, abs=1e-12)
    assert from_table == pytest.approx(from_metrics, abs=1e-12)
    assert from_table > 0.0


@pytest.mark.parametrize("kind", VOLUME_KINDS)
def test_node_table_agrees_with_metrics_and_face_count(kind: str) -> None:
    assert CELL_NODES[kind] == NODES_PER_CELL[kind]
    expected_faces = 5 if kind == "pyramid" else 5  # 1 quad + 4 tri, or 2 tri + 3 quad
    assert len(CELL_FACES[kind]) == expected_faces
    assert all(max(face) < CELL_NODES[kind] for face in CELL_FACES[kind])
    # prism: 2 triangles + 3 quads; pyramid: 1 quad + 4 triangles
    assert sorted({len(face) for face in CELL_FACES[kind]}) == [3, 4]


def test_every_gated_kind_has_a_face_table() -> None:
    """A kind the door measures must have a table, otherwise the door would raise."""
    assert set(CELL_KINDS_WITH_INTERNAL_FACES) <= set(CELL_FACES)
    assert set(CELL_KINDS_WITH_INTERNAL_FACES) == {"tet", "hex", "prism", "pyramid"}


# --------------------------------------------------------------------------------------
# closed forms on conforming stacks
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("shear", [0.0, 0.25, 0.5, 1.0, 2.0])
def test_prism_stack_reproduces_closed_form(shear: float) -> None:
    points, cells = prism_stack(layers=3, shear=shear)
    faces = non_orthogonality(points, cells, "prism")
    exact = math.degrees(math.atan(shear))
    assert faces.cells == 3
    assert faces.internal_faces == 2
    assert faces.boundary_faces == 11  # 3 columns x 5 faces - 2 shared
    assert faces.degenerate_faces == 0
    assert faces.maximum == pytest.approx(exact, abs=1e-12)
    assert faces.mean == pytest.approx(exact, abs=1e-12)


@pytest.mark.parametrize("shear", [0.0, 0.5, 1.0, 3.0])
def test_pyramid_pair_reproduces_closed_form(shear: float) -> None:
    points, cells = pyramid_pair(shear=shear)
    faces = non_orthogonality(points, cells, "pyramid")
    exact = math.degrees(math.atan(shear))
    assert faces.internal_faces == 1
    assert faces.boundary_faces == 8
    assert faces.maximum == pytest.approx(exact, abs=1e-12)


@pytest.mark.parametrize("kind", VOLUME_KINDS)
@pytest.mark.parametrize("shear", [0.0, 0.25, 0.5, 1.0, 3.0])
def test_library_angle_agrees_with_reference_route(kind: str, shear: float) -> None:
    """Second route: face normal by a triangle fan, angle by ``atan2``, no normalisation.

    :func:`tensorcad.meshquality.reference.reference_non_orthogonality` walks the same
    :data:`CELL_FACES` table but builds each normal as a fan of triangles around the face
    centroid and takes ``atan2(|n x d|, n . d)``. It was already generic over the table, so
    extending the table is what makes it a check for these two kinds at all -- and it is
    the route the audit tool reports, not a copy of the library's arithmetic.
    """
    points, cells = prism_stack(3, shear) if kind == "prism" else pyramid_pair(shear)
    library = non_orthogonality(points, cells, kind)
    second = reference_non_orthogonality(points, cells, kind)
    assert len(second) == library.internal_faces
    assert library.maximum == pytest.approx(max(second), abs=1e-12)
    assert library.maximum == pytest.approx(math.degrees(math.atan(shear)), abs=1e-12)


# --------------------------------------------------------------------------------------
# boundary faces, the gate, and the conforming-mesh guard
# --------------------------------------------------------------------------------------
def test_single_prism_has_no_internal_face_and_the_gate_says_so() -> None:
    points, cells = prism_stack(layers=1)
    faces = non_orthogonality(points, cells, "prism")
    assert faces.internal_faces == 0
    assert faces.boundary_faces == 5
    assert faces.maximum != faces.maximum  # NaN, not 0: nothing was measured
    (check,) = gate_non_orthogonality(faces)
    assert not check.passed and "no internal faces" in check.detail


def test_gate_fails_a_prism_stack_over_the_limit_and_names_the_measure() -> None:
    points, cells = prism_stack(layers=2, shear=3.0)  # atan(3) = 71.565 deg > 70
    faces = non_orthogonality(points, cells, "prism")
    (check,) = gate_non_orthogonality(faces)
    assert not check.passed
    assert check.name == "volume.max_non_orthogonality"
    assert check.value == pytest.approx(math.degrees(math.atan(3.0)), abs=1e-12)
    assert check.limit == "<= 70 deg"
    assert "1 internal face" in check.detail


def test_pyramid_gate_passes_a_moderate_shear() -> None:
    points, cells = pyramid_pair(shear=0.5)  # 26.57 deg
    (check,) = gate_non_orthogonality(non_orthogonality(points, cells, "pyramid"))
    assert check.passed


def test_three_cells_sharing_a_face_still_raises_for_prisms() -> None:
    points = PRISM
    cells = [(0, 1, 2, 3, 4, 5), (0, 1, 2, 3, 4, 5), (0, 1, 2, 3, 4, 5)]
    with pytest.raises(ValueError, match="non-conforming"):
        non_orthogonality(points, cells, "prism")


def test_wrong_node_count_is_rejected() -> None:
    with pytest.raises(ValueError, match="prism cells need 6 nodes"):
        non_orthogonality(PRISM, [(0, 1, 2, 3)], "prism")
    with pytest.raises(ValueError, match="pyramid cells need 5 nodes"):
        non_orthogonality(PYRAMID, [(0, 1, 2, 3, 4, 5)], "pyramid")


def test_unknown_kind_still_raises() -> None:
    with pytest.raises(ValueError, match="non_orthogonality supports"):
        non_orthogonality(PRISM, [(0, 1, 2, 3, 4, 5)], "wedge")


# --------------------------------------------------------------------------------------
# the two callers: door.admit and the CLI
# --------------------------------------------------------------------------------------
def test_door_measures_prism_non_orthogonality_instead_of_skipping_it() -> None:
    points, cells = prism_stack(layers=2, shear=3.0)
    verdict = admit(points, cells={"prism": cells})
    payload = verdict.payload
    assert "non_orthogonality_prism" in payload
    record = payload["non_orthogonality_prism"]
    assert record["internal_faces"] == 1
    assert record["max_deg"] == pytest.approx(math.degrees(math.atan(3.0)), abs=1e-12)
    failed = {check.name for check in verdict.checks if not check.passed}
    assert "volume.max_non_orthogonality" in failed
    assert verdict.passed is False


def test_door_measures_pyramid_non_orthogonality() -> None:
    points, cells = pyramid_pair(shear=0.25)
    verdict = admit(points, cells={"pyramid": cells})
    assert verdict.payload["non_orthogonality_pyramid"]["internal_faces"] == 1
    non_orth = [
        c for c in verdict.checks if c.name == "volume.max_non_orthogonality"
    ]
    assert len(non_orth) == 1 and non_orth[0].passed


def _run_cli(spec: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tensorcad.meshquality.cli", "--input", str(spec), *args],
        capture_output=True,
        text=True,
    )


def test_cli_max_non_orth_accepts_prism_cells(tmp_path: Path) -> None:
    points, cells = prism_stack(layers=2, shear=0.5)
    spec = tmp_path / "prisms.json"
    spec.write_text(
        json.dumps(
            {"points": [list(p) for p in points], "elements": {"prism": [list(c) for c in cells]}}
        ),
        encoding="utf-8",
    )
    out = tmp_path / "report.json"
    result = _run_cli(spec, "--max-non-orth", "70", "--json", str(out), "--quiet")
    assert result.returncode == 0, result.stderr
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["non_orthogonality"]["kind"] == "prism"
    assert payload["non_orthogonality"]["max_deg"] == pytest.approx(26.56505117707799, abs=1e-9)


def test_cli_max_non_orth_fails_a_steep_prism_stack(tmp_path: Path) -> None:
    points, cells = prism_stack(layers=2, shear=3.0)
    spec = tmp_path / "prisms.json"
    spec.write_text(
        json.dumps(
            {"points": [list(p) for p in points], "elements": {"prism": [list(c) for c in cells]}}
        ),
        encoding="utf-8",
    )
    out = tmp_path / "report.json"
    result = _run_cli(spec, "--max-non-orth", "70", "--json", str(out), "--quiet")
    assert result.returncode == 1
    gate = json.loads(out.read_text(encoding="utf-8"))["gate"]
    assert gate["passed"] is False
    assert "volume.max_non_orthogonality" in gate["failures"]