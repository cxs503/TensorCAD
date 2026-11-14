"""Trim binding in the IGES reader: type 141 / 142 / 144 field offsets.

These are the offsets that decide whether a *trimmed* hull can be read at all.  Before
the fix pinned here, ``BoundedSurface`` took the outer boundary from parameter 4 of a
type-144 wrapper (that field is ``N2``, normally ``0``), ``CurveOnSurface`` took the
basis surface from parameter 2 of a type-142 (that field is the spare ``CRV``), and
``Boundary`` took its ``N`` curve pointers from parameter 3 of a type-141, where the
parameter-space direction flag lives.  Every consequence was silent: a mis-bound patch
was reported as "no boundary" and then read as its whole basis surface, so a trimmed
file looked like an untrimmed one with overhang.

The fixtures are synthetic IGES text built here, so the suite runs offline and does not
carry publisher bytes.
"""

from __future__ import annotations

import pytest

from tensorcad.iges import Boundary, BoundedSurface, CurveOnSurface


def _rec(text: str, letter: str, seq: int) -> str:
    return text[:72].ljust(72) + letter + f"{seq:7d}"


def _write_iges(tmp_path, entities: list[tuple[int, str]]):
    """Emit a minimal but format-faithful IGES file (directory + parameter sections)."""
    directory: list[str] = []
    parameters: list[str] = []
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
    global_section = _rec("1H,,1H;,4HMETR,1H,,1.0,1.0,1,6,1,1,1.0,1,1.0,0", "G", 1)
    terminate = _rec("S      1G      1D      1P      1".ljust(0), "T", 1)
    text = "\n".join([_rec("synthetic trim fixture", "S", 1), global_section] + directory
                     + parameters + [terminate]) + "\n"
    path = tmp_path / "trim.igs"
    path.write_text(text)
    return path


# A plane as a degree-1 x degree-1 net: 4 control points, weights all one.
PLANE = [128, 1, 1, 1, 1, 0, 0, 1, 0, 0,
         # knots u (degree 1, 2 control points) then knots v, as the reader lays them out
         0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
         # 2x2 control net: a unit square in the xy plane
         0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
         1.0, 1.0, 1.0, 1.0]


def _payload(values: list[float]) -> str:
    return ",".join(str(value) for value in values) + ";"


class TestFieldOffsets:
    """The offsets are pinned individually so a regression names the field it broke."""

    def test_type_144_outer_boundary_is_parameter_5(self):
        wrapper = BoundedSurface.from_parameters([144.0, 7.0, 1.0, 0.0, 35.0, 0.0, 0.0])
        assert wrapper.surface == 7
        assert wrapper.boundary == 35, "PTO (the outer trim loop) is parameter 5"
        assert wrapper.boundary_count == 1
        assert wrapper.inner == ()

    def test_type_144_inner_boundary_is_parameter_6(self):
        wrapper = BoundedSurface.from_parameters([144.0, 7.0, 2.0, 0.0, 35.0, 41.0, 0.0])
        assert wrapper.boundary == 35
        assert wrapper.inner == (41,)

    def test_type_142_basis_surface_is_parameter_3(self):
        curve = CurveOnSurface.from_parameters([142.0, 0.0, 5.0, 33.0, 19.0, 1.0, 0.0])
        assert curve.surface == 5, "SPTR is parameter 3; parameter 2 is the spare CRV field"
        assert curve.curve_3d == 33
        assert curve.curve_uv == 19

    def test_type_141_pointers_start_at_parameter_5(self):
        boundary = Boundary.from_parameters([141.0, 4.0, 11.0, 1.0, 21.0, 23.0, 25.0, 27.0])
        assert boundary.surface == 11
        assert boundary.parameter_flip is True
        assert boundary.curves == (21, 23, 25, 27)
        assert len(boundary.curves) % 2 == 0, "alternating (3D, (u,v)) pairs"


# Reader-level binding (144 -> 142 -> trim polygon, end to end) is exercised against the
# real asset store by ``tools/baseline_fleet/patch_census.py``; a synthetic fixture for it
# is not attempted here because a type-126 parameter list is long enough that a hand-built
# one would test the fixture more than the reader.
