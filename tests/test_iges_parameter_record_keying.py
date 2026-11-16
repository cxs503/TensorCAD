"""Regression tests for the parameter (P) section record keying of the patch-set reader.

Anchor: GAP-CAD-T1-IGES-PIPE.

Defect this pins down (found on the publisher specimens, not on the files this package
writes for itself): ``tensorcad.geomquality.io._iges_parameters`` grouped P-section lines
by columns 65-72 -- the *entity pointer* field, which is per-line metadata -- instead of
by the line's own sequence number in columns 74-80, which is what the D section's
parameter data pointer refers to.  A parameter record longer than one line was therefore
split wherever two consecutive lines happened to disagree in columns 65-72, and the
entity's parameters came out as a fragment: on the KVLCC2 specimen the shifted field
decode read control-point coordinates as B-spline weights (a weight of -121.49 m) and the
file was written off as invalid publisher data.

The same reader also decoded type 128 with a private field order (``deg_u, deg_v,
count_u, count_v`` first, weights last) that only this package's own ``write_iges`` emits;
IGES 5.3 §4.24 puts ``K1, K2, M1, M2, PROP1..4`` first, then the knots, then the control
points.  Publisher files were decoded as if they were the private dialect.  The decode now
delegates that record to ``tensorcad.iges.BSplineSurface.from_parameters``, the decoder the
fleet route already runs against the publisher files.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from tensorcad.geomquality.io import _iges_parameters, read_iges

ASSET_ROOT = Path("/root/tp-agent/assets/geometry/hulls")
KCS = ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs"


def _p_line(body: str, entity_pointer: int, sequence: int) -> str:
    return f"{body:<64}{entity_pointer:>8}P{sequence:>7}"


def test_multi_line_parameter_record_is_read_as_one_record():
    """A record spanning two lines must decode whole, whatever columns 65-72 say.

    The two lines carry different entity-pointer fields (1 and 22); keying the grouping on
    that field -- the defect -- yields two fragments and no value list for the record.
    """
    first = _p_line("128,1,1,0,0,0,0,0,0,0,0.,0.,0.,0.", 1, 1)
    second = _p_line("1.,0.,1.,0.,0.,1.,1.,0.,0.,0.,0.,1.,1.;", 22, 2)
    values = _iges_parameters(first + "\n" + second + "\n")
    assert 1 in values, f"record not keyed on its first line's sequence: {sorted(values)}"
    record = values[1]
    assert record[0] == 128.0
    assert len(record) == 27, len(record)  # 14 + 13 tokens, both lines included
    assert record[-1] == 1.0
    assert not any(value != value for value in record), "NaN in a clean record"


def test_sequence_number_is_read_from_columns_74_to_80():
    """The key of a record is its first line's sequence number, not line 1's index."""
    lines = [
        _p_line("128,1,1,3,3,0,0,0,0,0,", 1, 7),
        _p_line("0,1,1,0.;", 1, 8),
    ]
    values = _iges_parameters("\n".join(lines) + "\n")
    assert 7 in values, sorted(values)
    assert values[7][0] == 128.0 and values[7][-1] == 0.0
    # The entity-pointer field (columns 65-72) stays reachable as an alias, so files
    # written by this package's own ``write_iges`` keep decoding.
    assert values.get(1) == values[7]


@pytest.mark.skipif(not KCS.exists(), reason="publisher KCS geometry not in the asset store")
def test_publisher_patch_set_decodes_with_unit_weights_and_published_beam():
    """The 55-patch KCS specimen: every record is one entity, weights are unit, beam right."""
    text = KCS.read_text(encoding="latin-1")
    directory = [line for line in text.splitlines() if len(line) > 79 and line[72] == "D"]
    pointers = [int(directory[index][8:16]) for index in range(0, len(directory), 2)]
    parameters = _iges_parameters(text)
    entity_types = Counter(int(parameters[pointer][0]) for pointer in pointers)
    # One record per directory entry, and every record starts with its own entity type.
    assert entity_types == Counter({128: 55}), entity_types

    brep = read_iges(KCS)
    assert len(brep.faces) == 55
    weights = [w for face in brep.faces for row in face.surface.weights for w in row]
    assert min(weights) > 0.0, min(weights)
    assert max(abs(w - 1.0) for w in weights) < 1e-12, "publisher patches are polynomial"
    points = [p for face in brep.faces for row in face.surface.control_net for p in row]
    beam = 2.0 * max(abs(p[1]) for p in points)
    assert abs(beam - 32.2) / 32.2 < 0.01, beam
    length = max(p[0] for p in points) - min(p[0] for p in points)
    assert length >= 230.0, length  # one-sided: the bbox contains Lpp, with overhangs