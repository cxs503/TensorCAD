"""Acceptance budget for a trim loop that does not invert onto its surface exactly.

GAP-CAD-T1-IGES-PIPE: the KP505 propeller (a third-party CATIA export) carries one hub
face whose 3-D trim curve misses its revolved surface by 2.2e-4 mm -- 1.5e-6 of the
250 mm part.  The old gate refused that face at ``1e-6`` of the model extent, i.e. it
dropped 1/8 of the model for a deviation 500x below the chord tolerance the mesh is
built to.  These tests pin both sides of the new contract:

* a face inside the **transfer budget** (``1e-4`` of the extent, a tenth of the chord
  budget) is tessellated and *labelled* -- the report never calls it exact;
* a face outside it is still refused (a wrong surface or a mis-paired curve misses by
  order 1, four decades above the budget);
* the mesh does not inherit the deviation: boundary vertices are emitted from the
  inverted parameters, so they sit on the surface to round-off;
* a refused face cannot leave a clean ``max_inversion_residual`` behind.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geomquality.iges_trim import Surface, TrimmedFace, TrimmedModel
from tensorcad.geomquality.trim_mesh import tessellate_trimmed_model

ROOT_EXTENT = 1.0  # largest |coordinate| of the trim loop, i.e. the model extent


class Bump(Surface):
    """A non-planar patch: ``S(u, v) = (u, v, h exp(-(u^2 + v^2)))`` over (-2, 2)^2."""

    kind = "bump"

    def __init__(self, height: float = 0.1) -> None:
        self.height = height

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:
        return ((-2.0, 2.0), (-2.0, 2.0))

    def evaluate(self, u: float, v: float) -> tuple[float, float, float]:
        return (u, v, self.height * math.exp(-(u * u + v * v)))


def _loop(offset: float, samples: int = 48) -> tuple[tuple[float, float, float], ...]:
    """The unit circle in parameter space, raised off the surface by ``offset`` in z."""
    surface = Bump()
    points = []
    for index in range(samples):
        angle = 2.0 * math.pi * index / samples
        u, v = math.cos(angle), math.sin(angle)
        x, y, z = surface.evaluate(u, v)
        points.append((x, y, z + offset))
    points.append(points[0])
    return tuple(points)


def _model(*offsets: float) -> TrimmedModel:
    faces = tuple(
        TrimmedFace(index=index, surface=Bump(), loops=(_loop(offset),),
                    name=f"bump{index}", source_entity=100 + index, surface_kind="bump")
        for index, offset in enumerate(offsets)
    )
    return TrimmedModel(faces=faces, units_name="MM", metres_per_unit=0.001,
                        source="synthetic", census={})


def _on_surface_deviation(vertices) -> float:
    surface = Bump()
    worst = 0.0
    for x, y, z in vertices:
        worst = max(worst, abs(z - surface.evaluate(x, y)[2]))
    return worst


def test_budgets_are_ordered_and_tied_to_the_model_extent():
    result = tessellate_trimmed_model(_model(0.0))
    extent = ROOT_EXTENT
    assert result.residual_tolerance == pytest.approx(1e-6 * extent)
    assert result.transfer_tolerance == pytest.approx(1e-4 * extent)
    assert result.transfer_tolerance >= result.residual_tolerance
    assert result.chord_tolerance == pytest.approx(1e-3 * extent)
    # The acceptance budget must stay an order of magnitude tighter than the accuracy the
    # mesh is built to, otherwise "inside the budget" would mean nothing for the mesh.
    assert result.transfer_tolerance <= 0.1 * result.chord_tolerance


def test_exact_trim_is_classed_exact_and_never_banded():
    result = tessellate_trimmed_model(_model(0.0))
    assert result.refused == {}
    assert result.faces_in_transfer_band == {}
    assert result.face_stats[0]["residual_class"] == "exact"
    assert result.max_residual <= result.residual_tolerance
    assert len(result.faces) > 100


def test_face_inside_the_transfer_budget_is_tessellated_and_labelled():
    offset = 1e-5  # between the exact budget (1e-6) and the transfer budget (1e-4)
    result = tessellate_trimmed_model(_model(offset))
    assert result.refused == {}
    assert set(result.faces_in_transfer_band) == {0}
    assert result.faces_in_transfer_band[0] == pytest.approx(offset, rel=0.05)
    assert result.face_stats[0]["residual_class"] == "transfer_budget"
    assert len(result.faces) > 100
    # No mesh vertex inherits the deviation: boundary vertices come from (u, v).
    assert _on_surface_deviation(result.vertices) <= 1e-12


def test_face_beyond_the_transfer_budget_is_still_refused():
    offset = 1e-2  # 1e0 of the extent, the scale a mis-paired curve misses by
    result = tessellate_trimmed_model(_model(offset))
    assert result.faces == []
    assert set(result.refused) == {0}
    assert "transfer budget" in result.refused[0]
    # The refused face's residual is still measured, so G2 cannot report a clean 0.0.
    assert result.max_residual == pytest.approx(offset, rel=0.05)
    assert result.max_residual > result.transfer_tolerance


def test_a_caller_can_tighten_the_budget_but_never_excuse_a_face_through_it():
    result = tessellate_trimmed_model(_model(1e-5), transfer_tolerance=1e-6)
    assert set(result.refused) == {0}
    assert result.transfer_tolerance == pytest.approx(1e-6)


def test_one_refused_face_does_not_hide_behind_one_good_face():
    result = tessellate_trimmed_model(_model(0.0, 1e-2))
    assert set(result.refused) == {1}
    assert set(result.faces_in_transfer_band) == set()
    assert result.max_residual == pytest.approx(1e-2, rel=0.05)
    assert len(result.faces) > 100


def test_kp505_hub_face_is_recovered_inside_the_transfer_budget():
    """The third-party specimen the budget was calibrated on: 7/8 -> 8/8 faces."""
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1]
              / "benchmarks" / "anchors" / "sources" / "KP505_modelScale_250mm.igs")
    if not source.exists():
        pytest.skip("KP505 propeller specimen not present in this checkout")
    from tensorcad.geomquality.trim_mesh import import_iges

    result = import_iges(str(source), repair=False)
    report = result.report
    imported = report["import"]
    assert imported["faces_in_file"] == 8
    assert imported["faces_refused"] == 0, imported["refused"]
    assert imported["faces_tessellated"] == 8
    assert set(imported["faces_in_transfer_band"]) == {"1"}
    residual = imported["faces_in_transfer_band"]["1"]
    assert imported["residual_tolerance"] < residual <= imported["transfer_tolerance"]
    assert report["gates"]["G1_all_faces_read"]["passed"] is True
    assert report["gates"]["G2_trim_on_surface"]["passed"] is True
    # The recovery is declared, not disguised as an exact inversion.
    assert report["gates"]["G2_trim_on_surface"]["faces_in_transfer_band"] == 1
    assert report["gates"]["G2_trim_on_surface"]["exact_limit"] == imported["residual_tolerance"]