"""The surface census: which entity types a hull's geometry actually arrives in.

Measured fact this suite pins (asset store, 2026-11-14):

| hull   | type 128 (decoded) | type 126 (not decoded) |
|--------|--------------------|------------------------|
| KCS    | 55                 | 0                      |
| JBC    | 45                 | 172                    |
| KVLCC2 | 19                 | 304                    |

``tensorcad.iges`` decodes type 128 only.  KCS is therefore read whole and gates
VERIFIED (+0.30 %); JBC and KVLCC2 arrive in residue and were recorded WITHHELD
with a *downstream* message ("field produced fewer than two usable waterline
slices"), which reads like a shape verdict.  These tests hold the reader to
reporting the census, and the fleet producer to attributing such a refusal to the
reader instead of to the hull.

Synthetic files only: no publisher bytes are needed, so this runs offline.
"""

from __future__ import annotations

import pytest

from tensorcad.fleet_baseline import ShipSpec, build_ship
from tensorcad.iges import read_iges


def _iges(entities: list[tuple[int, str]]) -> str:
    """A minimal but format-faithful IGES file carrying ``entities``.

    ``entities`` is ``[(entity_type, parameter_payload), ...]``.  Directory and
    parameter sections are emitted in the 80-column fixed layout the reader
    parses: directory pointers in columns 1-8, parameter-line counts in columns
    25-32 of the second directory record, parameter sequence numbers in columns
    74-80.
    """

    def _rec(text: str, letter: str, seq: int) -> str:
        return text[:72].ljust(72) + letter + f"{seq:7d}"

    directory = []
    parameters = []
    para_seq = 1
    for index, (entity_type, payload) in enumerate(entities):
        chunks = [payload[i:i + 64] for i in range(0, len(payload), 64)] or [""]
        pointer = para_seq
        for chunk in chunks:
            parameters.append(_rec(chunk, "P", para_seq))
            para_seq += 1
        de_seq = 2 * index + 1
        first = (f"{entity_type:8d}{pointer:8d}{0:8d}{0:8d}{0:8d}"
                 f"{0:8d}{0:8d}{0:8d}{'00000000':>8s}")
        second = (f"{entity_type:8d}{0:8d}{0:8d}{0:8d}{len(chunks):8d}"
                  f"{0:8d}{0:8d}{0:8d}{'00000000':>8s}")
        directory.append(_rec(first, "D", de_seq))
        directory.append(_rec(second, "D", de_seq + 1))

    start = _rec("S      1", "S", 1)
    global_section = _rec(
        "1H,,1H;,4HTEST,4HT.igs,8Htensorca,8Htensorca,32,38,6,38,15,"
        "4HTEST,1.0,2,2HMM,1000,1.0,13H20261114.1200,1.0E-07,",
        "G",
        1,
    )
    return "\n".join([start, global_section, *directory, *parameters,
                      _rec("S      1", "T", 1)]) + "\n"


def _write(tmp_path, entities, name="synthetic.igs"):
    path = tmp_path / name
    path.write_text(_iges(entities))
    return path


#: One 128 (rational B-spline surface) and one 126 (B-spline surface) entity.
_P128 = ("128,1,1,1,1,0,0,0,0,0,"                              # type, upper u/v, degree u/v, props
         "0,0,1,1,0,0,1,1,"                                      # knots u, knots v
         "0,0,0,1, 1,0,0,1, 0,1,0,1, 1,1,0,1;")                  # 4 points x (x,y,z,w)
_P126 = "126,1,1,0,0,1,0,0.0,0.0,2.0,1.0,1.0,0.0,0.0,0.0,2.0,0.0,0.0,1.0,2.0;"


def test_census_separates_decoded_from_dropped_surfaces(tmp_path):
    path = _write(tmp_path, [(128, _P128), (126, _P126)])
    model = read_iges(str(path))
    census = model.surface_census()
    assert census["readable_type128"] == 1
    assert census["unreadable_type126"] == 1
    assert census["by_type"][126] == 1 and census["by_type"][128] == 1
    # The decoded geometry is the 128 residue only -- this is the mechanism by
    # which a hull can be read "successfully" in residuum.
    assert len(model.surfaces()) == census["readable_type128"]


def test_census_counts_every_type_in_the_file(tmp_path):
    path = _write(tmp_path, [(126, _P126), (126, _P126), (126, _P126)])
    census = read_iges(str(path)).surface_census()
    assert census["readable_type128"] == 0
    assert census["unreadable_type126"] == 3
    assert census["by_type"] == {126: 3}


def _spec(path) -> ShipSpec:
    return ShipSpec(
        name="SYNTH",
        draft=5.0,
        length_pp=100.0,
        beam=20.0,
        published_volume_m3=None,
        geometry_scale=1.0,
        iges_path=str(path),
    )


def test_refusal_names_the_dropped_entity_type_not_the_hull(tmp_path):
    """A 126-only file must fail as a reader gap, with the count in the message."""
    path = _write(tmp_path, [(126, _P126)])
    with pytest.raises(ValueError) as excinfo:
        build_ship(_spec(path))
    message = str(excinfo.value)
    assert "type-126" in message
    assert "1 type-126" in message
    # and it must not read as a hydrostatic/shape verdict
    assert "waterline slices" not in message


def test_refusal_stays_plain_when_no_surfaces_were_dropped(tmp_path):
    """An empty file is a plain absence -- no 126 excuse is invented for it."""
    path = _write(tmp_path, [(102, "102,0,0,0,0,0;")])
    with pytest.raises(ValueError) as excinfo:
        build_ship(_spec(path))
    assert "type-126" not in str(excinfo.value)