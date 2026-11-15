"""Hidden-surface removal on an overlapping patch soup (``wetted_area``).

The invariant under test is *duplication invariance*: describing the same skin
twice must not change the wetted area.  It is the property a publisher IGES file
violates (``KCS_hull_SVA.igs`` raw triangle sum = +19.4 % against the published
9 424 m²), and it is checked here on geometry with an exact answer, so the
de-duplication is not merely "closer to the reference".
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.geomquality.hull_surface import (
    PARTICULARS,
    hull_characteristics,
    wetted_area,
)

ASSET_ROOT = Path("/root/tp-agent/assets/geometry/hulls")
KCS = ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs"

LENGTH, HALF_BEAM, DRAFT = 10.0, 2.5, 1.0


def _box_half_hull(length=LENGTH, half_beam=HALF_BEAM, draft=DRAFT):
    """Half-breadth box: bottom ``z = 0``, side ``y = half_beam``, ends at ``x = 0, L``.

    The ``y = 0`` face is the centreplane cut, not wetted surface; there is no deck
    (a face lying in ``z = T`` would be at the waterline, not under it).
    """
    index = {}
    vertices = []
    for x in (0.0, length):
        for y in (0.0, half_beam):
            for z in (0.0, draft):
                index[(x, y, z)] = len(vertices)
                vertices.append((x, y, z))

    def quad(*corners):
        ids = [index[c] for c in corners]
        return [(ids[0], ids[1], ids[2]), (ids[0], ids[2], ids[3])]

    L, B, T = length, half_beam, draft
    faces = []
    faces += quad((0, 0, 0), (L, 0, 0), (L, B, 0), (0, B, 0))     # bottom
    faces += quad((0, B, 0), (L, B, 0), (L, B, T), (0, B, T))     # starboard side
    faces += quad((0, 0, 0), (0, B, 0), (0, B, T), (0, 0, T))     # aft end
    faces += quad((L, 0, 0), (L, B, 0), (L, B, T), (L, 0, T))     # fore end
    faces += quad((0, 0, 0), (L, 0, 0), (L, 0, T), (0, 0, T))     # centreplane cut
    return vertices, faces


#: side + bottom + both ends, doubled for symmetry; the centreplane adds no area
BOX_AREA = 2.0 * (LENGTH * DRAFT + LENGTH * HALF_BEAM + 2.0 * HALF_BEAM * DRAFT)


def test_box_wetted_area_is_exact_and_unchanged_by_dedup():
    vertices, faces = _box_half_hull()
    area, discarded = wetted_area(vertices, faces, waterline=DRAFT)
    assert area == pytest.approx(BOX_AREA, rel=1e-12)
    assert discarded == pytest.approx(0.0, abs=1e-12)


def test_duplicated_skin_does_not_double_the_area():
    """A second, exactly coincident copy of the starboard side adds no area."""
    vertices, faces = _box_half_hull()
    duplicate = faces[2:4]  # the starboard side quad
    offset = len(vertices)
    both = list(vertices) + list(vertices)
    doubled = list(faces) + [tuple(i + offset for i in face) for face in duplicate]

    area, discarded = wetted_area(both, doubled, waterline=DRAFT)
    assert area == pytest.approx(BOX_AREA, rel=1e-9), "duplicated skin changed the area"
    assert discarded == pytest.approx(2.0 * LENGTH * DRAFT, rel=1e-6)

    raw, _ = wetted_area(both, doubled, waterline=DRAFT, dedup=False)
    assert raw == pytest.approx(BOX_AREA + 2.0 * LENGTH * DRAFT, rel=1e-12)
    assert raw > area  # the defect the de-duplication removes


def test_sheet_hidden_under_the_outer_skin_is_dropped():
    """A patch one centimetre inside the hull skin is not wetted surface."""
    vertices, faces = _box_half_hull()
    inner = 0.01
    offset = len(vertices)
    patch = []
    xs, zs = (2.0, 3.0, 4.0), (0.1, 0.5, 0.9)
    grid = {}
    for x in xs:
        for z in zs:
            grid[(x, z)] = len(vertices)
            vertices.append((x, HALF_BEAM - inner, z))
    for xi in range(len(xs) - 1):
        for zi in range(len(zs) - 1):
            a = grid[(xs[xi], zs[zi])]
            b = grid[(xs[xi + 1], zs[zi])]
            c = grid[(xs[xi + 1], zs[zi + 1])]
            d = grid[(xs[xi], zs[zi + 1])]
            patch += [(a, b, c), (a, c, d)]
    area, discarded = wetted_area(vertices, list(faces) + patch, waterline=DRAFT)
    assert area == pytest.approx(BOX_AREA, rel=1e-9)
    assert discarded == pytest.approx(2.0 * (len(xs) - 1) * (len(zs) - 1) * 1.0 * 0.4,
                                      rel=1e-6)


def test_full_hull_model_is_not_doubled():
    """A full (both-sides) model must not be multiplied by the symmetry factor."""
    vertices, faces = _box_half_hull()
    wetted_faces = faces[:-2]  # drop the centreplane cut from both sides
    mirrored = [(-x, -y, z) for (x, y, z) in vertices]
    offset = len(vertices)
    hull = list(vertices) + mirrored
    full_faces = wetted_faces + [tuple(i + offset for i in face) for face in wetted_faces]
    area, discarded = wetted_area(hull, full_faces, waterline=DRAFT)
    # BOX_AREA is already the full-hull (both sides) wetted area
    assert area == pytest.approx(BOX_AREA, rel=1e-12)
    assert discarded == pytest.approx(0.0, abs=1e-12)


@pytest.mark.skipif(not KCS.exists(), reason=f"workshop asset absent: {KCS}")
def test_kcs_dedup_flips_the_pinned_wetted_area_defect():
    """The published KCS wetted area is reached once duplicate sheets are removed."""
    from tensorcad.geomquality.trim_mesh import import_cad

    result = import_cad(KCS, weld_tolerance=1e-3)
    mesh = result.mesh
    particulars = PARTICULARS["KCS"]
    got = hull_characteristics(list(mesh.vertices), list(mesh.faces),
                               waterline=particulars["draft"], stations=61, levels=48)
    # the raw soup sum stays on record as the measured defect (+19.4 %)
    assert got.s_wetted_naive > 1.15 * particulars["s_wetted"]
    assert got.s_wetted_discarded > 0.0
    assert got.s_wetted == pytest.approx(got.s_wetted_naive - got.s_wetted_discarded, rel=1e-9)
    error = 100.0 * (got.s_wetted / particulars["s_wetted"] - 1.0)
    assert abs(error) <= 3.0, f"de-duplicated wetted area error {error:+.2f} %"
    # the de-duplication must not move the volume/Cb, which never summed the soup
    assert got.volume == pytest.approx(51870.7, rel=2e-3)