"""Uniform 1-to-4 triangle subdivision: a *single connected* large-mesh source.

Why this module exists
----------------------
The parametric hull generator tops out at ``nx <= 4001`` stations and
``nz <= 513`` levels, i.e. ``4*(nx-1)*(nz-1) = 8,199,998`` triangles for one
body.  The 10^7-triangle band therefore had to be reached with
:func:`tensorcad.perf_benchmark.build_fleet_mesh`, which assembles *k* disjoint
copies -- that measures a fleet, not "one 10^7-triangle body".  Subdividing a
coarser hull closes that gap with no new geometry authoring: the result is one
connected closed surface, and because the four children of a planar parent tile
it exactly, the enclosed volume (and the surface area) is preserved to
floating-point, which makes the integral a *check* rather than an estimate.

Algorithm
---------
Each triangle ``(a, b, c)`` becomes four (winding preserved)::

    (a, m_ab, m_ca)  (m_ab, b, m_bc)  (m_ca, m_bc, c)  (m_ab, m_bc, m_ca)

where ``m_xy`` is the arithmetic midpoint of the edge ``{x, y}``.  Midpoints are
keyed on the *undirected index pair* ``(min, max)``, so two faces sharing an
edge get the **same** new vertex index by construction -- welding uses index
identity, never a tolerance.  No spatial search, no float comparison, no
dependency: the module is stdlib-only like the rest of the core.

Topology bookkeeping (a closed mesh with ``E = 3F/2``)::

    V' = V + E      F' = 4F      E' = 2E + 3F      V' - E' + F' = 2  (unchanged)

so watertightness, edge manifoldness and the Euler characteristic survive
subdivision; the parent/child orientation stays outward because the child
windings are the parent's winding restricted to each sub-triangle.

Cost is O(F) time and memory; the 10^6-10^7 rungs are dominated by Python
object construction (``Mesh`` holds tuples), not by the arithmetic.
"""
from __future__ import annotations

from .geometry import Face, Mesh, Vec3

__all__ = ["subdivide_mesh", "expected_face_count", "subdivision_edges"]


def expected_face_count(face_count: int, levels: int = 1) -> int:
    """Faces after ``levels`` uniform subdivisions (``4**levels`` growth)."""
    if levels < 0:
        raise ValueError(f"levels must be >= 0, got {levels!r}")
    return int(face_count) * 4 ** levels


def subdivision_edges(mesh: Mesh) -> dict[tuple[int, int], int]:
    """Map undirected edge ``(min, max)`` -> new vertex index (after ``vertices``)."""
    base = len(mesh.vertices)
    edges: dict[tuple[int, int], int] = {}
    for a, b, c in mesh.faces:
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            if key not in edges:
                edges[key] = base + len(edges)
    return edges


def subdivide_mesh(mesh: Mesh, levels: int = 1) -> Mesh:
    """Return ``mesh`` refined ``levels`` times with 1-to-4 midpoint subdivision.

    ``levels=0`` returns an equal-but-copied mesh (identity in everything but
    object identity).  The input is never mutated.
    """
    if levels < 0:
        raise ValueError(f"levels must be >= 0, got {levels!r}")
    if not mesh.faces:
        raise ValueError("cannot subdivide a mesh with no faces")

    vertices: list[Vec3] = list(mesh.vertices)
    faces: list[Face] = list(mesh.faces)

    for _ in range(int(levels)):
        base = len(vertices)
        edges: dict[tuple[int, int], int] = {}
        # Midpoints are created once per undirected edge; adjacent faces share
        # the very same tuple object, so the weld is exact by identity.
        midpoints: list[Vec3] = []
        for a, b, c in faces:
            for u, v in ((a, b), (b, c), (c, a)):
                key = (u, v) if u < v else (v, u)
                if key not in edges:
                    pu, pv = vertices[u], vertices[v]
                    edges[key] = base + len(midpoints)
                    midpoints.append((
                        0.5 * (pu[0] + pv[0]),
                        0.5 * (pu[1] + pv[1]),
                        0.5 * (pu[2] + pv[2]),
                    ))
        vertices = vertices + midpoints

        new_faces: list[Face] = []
        append = new_faces.append
        for a, b, c in faces:
            mab = edges[(a, b) if a < b else (b, a)]
            mbc = edges[(b, c) if b < c else (c, b)]
            mca = edges[(c, a) if c < a else (a, c)]
            append((a, mab, mca))
            append((mab, b, mbc))
            append((mca, mbc, c))
            append((mab, mbc, mca))
        faces = new_faces

    return Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )