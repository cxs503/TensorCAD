"""Entity 142 (curve on a parametric surface): parametric image vs 3-D curve.

A 142 record normally carries the same trim twice -- once as a 3-D curve and once as its
image in the surface's parameter space.  Reading the parameter-space image as a 3-D curve
is geometrically *plausible* (it is a closed loop somewhere inside the model bounding
box), which is how the 250 mm KP505 propeller ended up with eight faces, four of them
carrying a second, nonsense boundary loop and refused at a 1.26e2 model-unit inversion
residual (8.3e5 times the gate) instead of being meshed.

These tests pin the measured classification and the parametric trim route:

* helpers that decide "is this a parameter-space image" and how far it is from the 3-D
  curve (agreement is measured to the nearest *segment*, not the nearest vertex);
* the mesh stage consuming a parametric loop without inverting anything, and still
  gating a 3-D loop that is off its surface;
* the real file, when the asset is present (it is not redistributable, so the test skips).
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tensorcad.geomquality.iges_trim import (
    PlaneSurface,
    TrimmedFace,
    TrimmedModel,
    _close_uv,
    _uv_deviation,
    _uv_in_domain,
    read_iges_trimmed,
)
from tensorcad.geomquality.trim_mesh import tessellate_trimmed_model

KP505 = Path(__file__).resolve().parents[1] / "benchmarks/anchors/sources/KP505_modelScale_250mm.igs"
KP505_SHA256 = "0bb54ae7618e8d8b1b7c23de8d384e04b0465071e34fe32e1d606ddbb5df462f"


def _plane(margin: float = 4.0) -> PlaneSurface:
    return PlaneSurface(
        origin=(0.0, 0.0, 0.0),
        u_axis=(1.0, 0.0, 0.0),
        v_axis=(0.0, 1.0, 0.0),
        normal=(0.0, 0.0, 1.0),
        margin=margin,
    )


def _frame(uv_loop, loop3d, cross_check=0.0) -> TrimmedModel:
    face = TrimmedFace(
        index=0,
        surface=_plane(),
        loops=(tuple(loop3d),),
        name="face_000",
        source_entity=1,
        surface_kind="plane",
        uv_loops=(tuple(uv_loop),),
        cross_checks=(cross_check,),
    )
    return TrimmedModel(faces=(face,), units_name="MM", metres_per_unit=0.001,
                        source="synthetic", census={})


# --------------------------------------------------------------------------- helpers
def test_uv_in_domain_accepts_only_parameter_range() -> None:
    surface = _plane(margin=1.0)  # domain is +-4
    assert _uv_in_domain(surface, [(0.0, 0.0, 0.0), (3.9, -3.9, 0.0)])
    assert not _uv_in_domain(surface, [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0)])


def test_uv_deviation_is_measured_to_segments_not_vertices() -> None:
    """A coarse 3-D reference must not be reported as a disagreement.

    The polyline below is the same straight line as the parameter image, sampled at two
    vertices only.  A vertex metric reads that sampling gap (0.5) as a disagreement; the
    segment metric reads the truth (0).
    """
    surface = _plane(margin=1.0)
    surface = PlaneSurface(
        origin=(0.0, 0.0, 0.0), u_axis=(1.0, 0.0, 0.0), v_axis=(0.0, 1.0, 0.0),
        normal=(0.0, 0.0, 1.0), margin=1.0,
    )
    uv = [(0.25, 0.0), (0.5, 0.0), (0.75, 0.0)]
    polyline = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    assert _uv_deviation(surface, uv, polyline) < 1e-12
    vertex_metric = max(
        min(math.dist(surface.evaluate(u, v), point) for point in polyline) for u, v in uv
    )
    assert vertex_metric > 0.4  # what the naive metric would have claimed


def test_uv_deviation_detects_a_real_mismatch() -> None:
    surface = _plane(margin=1.0)
    uv = [(0.25, 0.5), (0.75, 0.5)]
    polyline = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    assert _uv_deviation(surface, uv, polyline) == pytest.approx(0.5, rel=1e-12)


def test_close_uv_closes_an_open_loop_and_reports_the_gap() -> None:
    surface = _plane(margin=1.0)
    closed, gap = _close_uv(surface, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.2, 0.0)])
    assert closed[0] == closed[-1]
    assert gap == pytest.approx(math.hypot(0.2, 0.0), rel=1e-12)
    already, gap_zero = _close_uv(surface, [(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)])
    assert already[0] == already[-1] and len(already) == 3 and gap_zero == 0.0


# ------------------------------------------------------- the mesh stage honours it
def test_parametric_loop_is_used_without_inversion() -> None:
    """A parametric loop meshes even when its 3-D polyline is nowhere near the surface.

    That is the point of the route: on a 126 x 54 control-point patch the inversion walk
    over a 1 825-point loop costs minutes, and the parameter-space image is exact.
    """
    uv = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]
    off_surface = [(x, y, 5.0) for x, y, _ in [(u, v, 0.0) for u, v in uv]]
    tessellation = tessellate_trimmed_model(_frame(uv, off_surface))
    assert tessellation.refused == {}
    assert tessellation.max_residual == 0.0
    stats = tessellation.face_stats[0]
    assert stats["trim_routes"] == {"parametric": 1, "inverted": 0}
    assert stats["parametric_cross_check"] == 0.0
    assert stats["triangles"] >= 2


def test_three_d_loop_off_its_surface_is_still_refused() -> None:
    """Without a parametric image the inversion gate is unchanged: no silent pass."""
    face = TrimmedFace(
        index=0,
        surface=_plane(),
        loops=(((0.0, 0.0, 5.0), (2.0, 0.0, 5.0), (2.0, 2.0, 5.0), (0.0, 0.0, 5.0)),),
        name="face_000",
        source_entity=1,
        surface_kind="plane",
    )
    model = TrimmedModel(faces=(face,), units_name="MM", metres_per_unit=0.001,
                         source="synthetic", census={})
    tessellation = tessellate_trimmed_model(model)
    assert tessellation.faces == []
    assert 0 in tessellation.refused
    assert "does not lie on its surface" in tessellation.refused[0]


# ------------------------------------------------------------- real third-party file
@pytest.mark.skipif(not KP505.exists(), reason="KP505 IGES asset is not redistributable")
def test_kp505_entity_142_parametric_image_is_classified() -> None:
    """The record that used to be read as a second 3-D loop now yields one exact trim."""
    import hashlib

    digest = hashlib.sha256(KP505.read_bytes()).hexdigest()
    assert digest == KP505_SHA256
    model = read_iges_trimmed(str(KP505))
    assert len(model.faces) == 8
    face = next(item for item in model.faces if item.source_entity == 123)
    assert face.trim_routes() == ["parametric"]
    cross_check = face.cross_checks[0]
    assert cross_check is not None and cross_check < 0.05  # model units (mm)
    # The parametric image replaces the 3-D curve of the same 142 record, and the
    # surplus reference is recorded as a refusal instead of becoming a boundary loop.
    assert any("same trim as its parameter-space image" in text
               for text in model.refused.values())