"""Why an exact divergence integral over an open hull skin is *short*, in closed form.

Anchor ``GAP-CAD-T1-IGES-PIPE``.  A published hull IGES is a **skin**, not a solid:
it carries the wetted surface and usually nothing else.  The centreplane cut and the
design-waterline cap are not part of the file, so ``∮ n dA`` over the file is *not*
zero and ``V = (1/3) ∮ r·n dA`` over the file is *not* the displacement.

That is exactly what happened to KVLCC2: an exact divergence-theorem integral over the
parsed patch set returned 256 657 m³ vs the published 312 622 m³ (−17.9 %) with
``|∮ n dA| = 11 272 m²``, and the shortfall was read as a defective patch set.  These
tests pin the arithmetic instead: on a skin whose displacement is known exactly, the
divergence integral is short by **exactly** ``(1/3) · T · A_wl`` — the missing waterline
cap — and adding that one term recovers the volume to machine precision.  The
centreplane cap of a half-breadth model contributes ``r·n = y = 0`` and therefore
nothing at all.

The second consequence, which is what makes the module's own route honest: the
envelope-breadth integral in :func:`hull_surface.hull_characteristics` needs no closure
and lands on the exact volume on the very same open skin.  Geometry and route are
therefore separable, and a short divergence integral says nothing about the shape.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geomquality.hull_surface import hull_characteristics

LENGTH, HALF_BEAM, DRAFT = 10.0, 2.5, 1.0


def _box_skin(half_beam: float = HALF_BEAM, *, half: bool = True,
              cap: bool = False, half_cap: bool = False):
    """A box-shaped hull carried as an *open* skin (no centreplane, no waterline cap).

    Vertex index = ``4*xi + 2*yi + zi`` over ``x in (0, L)``, ``y`` in the breadth pair,
    ``z in (0, T)``.  Outward normals throughout: ``-x`` at ``x=0``, ``+x`` at ``x=L``,
    ``-z`` on the bottom, ``±y`` on the side shells.  ``cap=True`` adds the design
    waterline (``z = T``, outward ``+z``), which is the only cap a hull skin is missing
    apart from the centreplane cut.
    """
    ys = (0.0, half_beam) if half else (-half_beam, half_beam)
    vertices = [(x, y, z) for x in (0.0, LENGTH) for y in ys for z in (0.0, DRAFT)]
    far = 2 * len(ys)
    faces = [
        (0, 3, 2), (0, 1, 3),                                  # x = 0, outward -x
        (far, far + 2, far + 3), (far, far + 3, far + 1),      # x = L, outward +x
        (0, far + 2, far), (0, 2, far + 2),                    # z = 0, outward -z
        (2, far + 3, far + 2), (2, 3, far + 3),                # y = +b, outward +y
    ]
    if not half:
        faces += [(0, far, far + 1), (0, far + 1, 1)]          # y = -b, outward -y
    if cap:
        faces += [(1, far + 1, far + 3), (1, far + 3, 3)]      # z = T, outward +z
    if half_cap:
        faces += [(0, far, far + 1), (0, far + 1, 1)]          # y = 0, outward -y
    return vertices, faces


def _divergence(vertices, faces):
    """``(∮ n dA, (1/3) ∮ r·n dA)`` for the faces exactly as given."""
    nx = ny = nz = volume = 0.0
    for face in faces:
        a, b, c = (vertices[i] for i in face)
        u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        cross = (u[1] * v[2] - u[2] * v[1],
                 u[2] * v[0] - u[0] * v[2],
                 u[0] * v[1] - u[1] * v[0])
        nx += 0.5 * cross[0]
        ny += 0.5 * cross[1]
        nz += 0.5 * cross[2]
        volume += (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]) / 6.0
    return (nx, ny, nz), volume


def test_open_half_skin_divergence_integral_is_short_by_the_waterline_cap():
    """``V_skin = V - (1/3)·T·A_wl`` exactly; the centreplane cap adds nothing."""
    vertices, faces = _box_skin()
    closure, volume = _divergence(vertices, faces)
    half_volume = LENGTH * HALF_BEAM * DRAFT                       # 25 m³
    waterplane = LENGTH * HALF_BEAM                                # 25 m² (one side)
    cap_term = DRAFT * waterplane / 3.0                            # 8.333… m³

    assert volume == pytest.approx(half_volume - cap_term, rel=1e-12)
    assert volume + cap_term == pytest.approx(half_volume, rel=1e-12)
    # The open skin is short by exactly the cap term: 33 1/3 % of this body.
    assert volume / half_volume == pytest.approx(2.0 / 3.0, rel=1e-12)

    # The closure residual is made of cap areas, not of a defect: the side shell leaves
    # ``+A_side`` in y and the missing waterline leaves ``-A_wl`` in z.  The centreplane
    # cut would remove the y term but contributes ``r·n = y = 0`` to the volume.
    assert closure[0] == pytest.approx(0.0, abs=1e-12)
    assert closure[1] == pytest.approx(LENGTH * DRAFT, rel=1e-12)
    assert closure[2] == pytest.approx(-waterplane, rel=1e-12)


def test_waterline_cap_closes_the_divergence_integral():
    """Adding the one missing cap takes the closure residual to zero and the volume exact."""
    skin_vertices, skin_faces = _box_skin()
    solid_vertices, solid_faces = _box_skin(cap=True, half_cap=True)
    closure, volume = _divergence(solid_vertices, solid_faces)
    assert max(abs(component) for component in closure) < 1e-12
    assert volume == pytest.approx(LENGTH * HALF_BEAM * DRAFT, rel=1e-12)
    # The waterline cap is the term that carries volume; the centreplane cap closes the
    # y residual but contributes nothing (``r·n = y = 0`` on that plane).
    waterline_only_vertices, waterline_only_faces = _box_skin(cap=True)
    closure_wl, volume_wl = _divergence(waterline_only_vertices, waterline_only_faces)
    assert volume_wl == pytest.approx(volume, rel=1e-12)
    assert closure_wl[1] == pytest.approx(LENGTH * DRAFT, rel=1e-12)
    assert len(solid_faces) == len(skin_faces) + 4
    assert skin_vertices == solid_vertices


def test_both_sides_skin_is_short_by_the_same_single_cap():
    """A both-sides skin misses only the waterline cap: the same formula, full breadth."""
    vertices, faces = _box_skin(half=False)
    _, volume = _divergence(vertices, faces)
    full_volume = LENGTH * 2.0 * HALF_BEAM * DRAFT                 # 50 m³
    waterplane = LENGTH * 2.0 * HALF_BEAM                          # 50 m²
    assert volume == pytest.approx(full_volume - DRAFT * waterplane / 3.0, rel=1e-12)


def test_envelope_route_needs_no_closure_on_the_very_same_skin():
    """``hull_characteristics`` is exact on the open skin: geometry ≠ route.

    This is the statement that lets KVLCC2 be judged on its shape even though its
    ``∮ n dA`` is non-zero (as it must be for any hull-surface file).
    """
    vertices, faces = _box_skin()
    got = hull_characteristics(vertices, faces, waterline=DRAFT, stations=40, levels=8)
    # reported for the *full* hull: the half-breadth skin is doubled by symmetry
    assert got.volume == pytest.approx(2.0 * LENGTH * HALF_BEAM * DRAFT, rel=1e-12)
    assert got.symmetry_factor == 2
    assert got.cb == pytest.approx(1.0, rel=1e-12)


def test_shortfall_is_a_cap_term_not_noise():
    """Scaling the waterplane area scales the shortfall one-for-one — it *is* the cap."""
    shortfalls = []
    for beam in (1.0, 2.0, 4.0):
        vertices, faces = _box_skin(beam)
        _, volume = _divergence(vertices, faces)
        shortfalls.append(LENGTH * beam * DRAFT - volume)
    ratios = [shortfalls[i + 1] / shortfalls[i] for i in range(len(shortfalls) - 1)]
    assert all(math.isclose(ratio, 2.0, rel_tol=1e-12) for ratio in ratios)