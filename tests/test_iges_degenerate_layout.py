"""Regression tests: type-128 control-tail tie-break on degenerate (2x2 / 2xn) nets and DE pointer numbering.

Both defects were found on the authoritative workshop bytes (ONRT Rhino export,
JBC) while working the GAP-CAD-T1-BASELINE-FLEET anchor; see
``docs/baseline_fleet_iges_layout_fix.md`` for the measured numbers.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.iges import (
    ENTITY_BSPLINE_SURFACE,
    BSplineSurface,
    DirectoryEntry,
    IgesFile,
    _split_control,
)


def _record_2x2_surface() -> list[float]:
    """The parameter list of ONRT DE 27 (P records 14-23), verbatim token order.

    A 2x2 bilinear patch written by Rhino: header, two knot vectors, then a
    constant weight block ``1,1,1,1``, then the four control points, then the
    writer's trailing parameter bounds ``0,0.0338,0,0.0338`` and ``0,1,25``.
    """

    knots = [0.0, 0.0, 0.03375527426160339, 0.03375527426160339]
    weights = [1.0, 1.0, 1.0, 1.0]
    points = [
        0.9224234297863333, 0.04348254925436184, -0.01833581102116458,
        0.9194814637864123, 0.04348254925436184, -0.05196263627320664,
        0.9224234297863333, 0.009727274992758455, -0.01833581102116458,
        0.9194814637864123, 0.009727274992758455, -0.05196263627320664,
    ]
    return [ENTITY_BSPLINE_SURFACE, 1, 1, 1, 1, 0, 0, 1, 0, 0] + knots + knots + weights + points + [
        0.0, 0.03375527426160339, 0.0, 0.03375527426160339, 0, 1, 25,
    ]


def test_degenerate_patch_prefers_the_weight_block_present_in_the_record():
    """A 2x2 net has no interior control point, so every reading scores 0.0.

    The tie must not fall back to insertion order: the interleaved reading of this
    record turns the weight block into a control point at the origin of the unit
    cube and leaves a negative weight (-0.052) behind, which is not representable.
    """

    surface = BSplineSurface.from_parameters(_record_2x2_surface())
    assert surface.layout == "weight_block"
    assert surface.layout_ambiguous is False
    assert all(w == 1.0 for w in surface.weights)

    (x_lo, y_lo, z_lo), (x_hi, y_hi, z_hi) = surface.bounding_box()
    # The true patch is 3 mm x 34 mm x 34 mm (model units, 1:48.9 scale); a
    # mis-read one spans the unit cube and reports max|y| = 1.0 m on a hull whose
    # half-beam is 0.19 m.
    assert x_hi - x_lo < 0.01
    assert y_hi - y_lo < 0.05
    assert z_hi - z_lo < 0.05
    assert max(abs(y_lo), abs(y_hi)) < 0.1


def test_non_positive_weight_candidate_is_rejected():
    """A candidate whose weights contain 0 or a negative value is not a reading."""

    body = [1.0, 2.0, 3.0, -0.5, 0.1, 0.2, 0.3, 1.0, 0.4, 0.5, 0.6, 1.0]
    control, weights, layout, _ambiguous = _split_control(
        body, points=3, rational=True, count_u=3, count_v=1
    )
    assert layout == "weight_block"
    assert all(w > 0.0 for w in weights)
    assert len(control) == 3


@pytest.mark.parametrize("pointer, expected", [(27, 28), (28, 28), (1, 2), (2, 2)])
def test_pointer_resolution_accepts_first_and_second_record_numbering(pointer, expected):
    """ONRT/Rhino numbers the first record of a DE, this reader keys the second."""

    entries = {
        2: DirectoryEntry(sequence=2, entity_type=314, parameter_pointer=1, parameter_lines=1,
                          form=0, label="COLOR", level=0, status="00000000"),
        28: DirectoryEntry(sequence=28, entity_type=ENTITY_BSPLINE_SURFACE, parameter_pointer=14,
                           parameter_lines=10, form=0, label="Shell", level=7, status="00000000"),
        30: DirectoryEntry(sequence=30, entity_type=141, parameter_pointer=24, parameter_lines=1,
                           form=1, label="", level=1, status="00000000"),
    }
    model = IgesFile(path="memory", global_fields=(), entries=entries, raw_parameters={},
                     units="m", unit_flag=6)
    resolved = model.resolve(pointer)
    assert resolved is not None and resolved.sequence == expected


def test_pointer_resolution_never_walks_past_the_next_entry():
    """An even pointer with no exact entry must not silently borrow a neighbour."""

    entries = {
        4: DirectoryEntry(sequence=4, entity_type=126, parameter_pointer=1, parameter_lines=1,
                          form=0, label="", level=0, status="00000000"),
    }
    model = IgesFile(path="memory", global_fields=(), entries=entries, raw_parameters={},
                     units="m", unit_flag=6)
    assert model.resolve(3) is not None      # odd pointer: its second record is key 4
    assert model.resolve(6) is None          # even pointer: no neighbour lookup
