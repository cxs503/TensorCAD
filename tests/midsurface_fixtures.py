"""Fixtures for the mid-surface tests: closed solids built from their own triangles.

Every mesh here is constructed rather than imported, so a test never depends on the module it
is judging. Winding is built outward and vertices are shared, so the meshes are closed and
edge-manifold by construction (the module refuses anything else, which is tested separately).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.section_gz import prism_mesh


def _key(point: Sequence[float]) -> tuple[float, float, float]:
    return (round(point[0], 12), round(point[1], 12), round(point[2], 12))


def box(
    low: Sequence[float] = (0.0, 0.0, 0.0),
    high: Sequence[float] = (1.0, 1.0, 1.0),
    *,
    divisions: int = 1,
) -> Mesh:
    """Axis-aligned box as a closed triangulated surface, each face an ``n x n`` grid."""
    builder = MeshBuilder()
    index_of: dict[tuple[float, float, float], int] = {}

    def vertex(point: Sequence[float]) -> int:
        key = _key(point)
        if key not in index_of:
            index_of[key] = builder.add_vertex((point[0], point[1], point[2]))
        return index_of[key]

    dx = high[0] - low[0]
    dy = high[1] - low[1]
    dz = high[2] - low[2]
    faces = (
        ((high[0], low[1], low[2]), (0.0, dy, 0.0), (0.0, 0.0, dz)),  # +x
        ((low[0], low[1], low[2]), (0.0, 0.0, dz), (0.0, dy, 0.0)),  # -x
        ((low[0], high[1], low[2]), (0.0, 0.0, dz), (dx, 0.0, 0.0)),  # +y
        ((low[0], low[1], low[2]), (dx, 0.0, 0.0), (0.0, 0.0, dz)),  # -y
        ((low[0], low[1], high[2]), (dx, 0.0, 0.0), (0.0, dy, 0.0)),  # +z
        ((low[0], low[1], low[2]), (0.0, dy, 0.0), (dx, 0.0, 0.0)),  # -z
    )
    steps = max(1, int(divisions))
    for origin, du, dv in faces:
        grid = []
        for i in range(steps + 1):
            row = []
            for j in range(steps + 1):
                point = (
                    origin[0] + du[0] * i / steps + dv[0] * j / steps,
                    origin[1] + du[1] * i / steps + dv[1] * j / steps,
                    origin[2] + du[2] * i / steps + dv[2] * j / steps,
                )
                row.append(vertex(point))
            grid.append(row)
        for i in range(steps):
            for j in range(steps):
                a, b, c, d = grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]
                builder.add_triangle(a, b, c)
                builder.add_triangle(a, c, d)
    return builder.build()


def cube_solid(length: float = 1.0, *, divisions: int = 1) -> Mesh:
    """A solid cube: its faces meet faces whose normals point the same way, so no wall pairs."""
    return box((0.0, 0.0, 0.0), (length, length, length), divisions=divisions)


def plate(length: float = 1.0, breadth: float = 1.0, thickness: float = 0.06, *, divisions: int = 1) -> Mesh:
    """A solid plate: two large faces pair at ``thickness``, its edge strips do not."""
    return box((0.0, 0.0, 0.0), (length, breadth, thickness), divisions=divisions)


def regular_section(inradius: float, sides: int) -> tuple[tuple[float, float], ...]:
    """Regular ``sides``-gon in the (y, z) plane, ``inradius`` = centre-to-facet distance."""
    angle = math.pi / sides
    circumradius = inradius / math.cos(angle)
    return tuple(
        (
            circumradius * math.cos(2.0 * angle * index + angle),
            circumradius * math.sin(2.0 * angle * index + angle),
        )
        for index in range(sides)
    )


def polygon_prism(section: Sequence[Sequence[float]], length: float, *, x0: float = 0.0) -> Mesh:
    """Convex section extruded along +x (the package's own prism builder)."""
    return prism_mesh(tuple((point[0], point[1]) for point in section), length, x0=x0)


def regular_prism(inradius: float, sides: int, length: float) -> Mesh:
    return polygon_prism(regular_section(inradius, sides), length)


def taper_prism(width: float = 1.0, thin: float = 0.05, thick: float = 0.3, length: float = 1.0) -> Mesh:
    """A quadrilateral prism whose wall thins from ``thick`` to ``thin`` along +y."""
    section = ((0.0, 0.0), (width, 0.0), (width, thick), (0.0, thin))
    return polygon_prism(section, length)


def l_prism(leg: float = 1.0, width: float = 0.3, length: float = 1.0) -> Mesh:
    """An L-shaped section extruded along +x: convex facets, non-convex caps on purpose.

    The caps are triangulated as two rectangles (the classic two-leg decomposition), so the
    mesh is closed and outward with a genuinely non-convex cap polygon -- which is exactly the
    panel geometry the mid-surface extractor has to refuse instead of over-covering.
    """
    section = (
        (0.0, 0.0),
        (leg, 0.0),
        (leg, width),
        (width, width),
        (width, leg),
        (0.0, leg),
    )
    builder = MeshBuilder()
    near = [builder.add_vertex((0.0, point[0], point[1])) for point in section]
    far = [builder.add_vertex((length, point[0], point[1])) for point in section]
    count = len(section)
    for index in range(count):
        nxt = (index + 1) % count
        a, b = near[index], near[nxt]
        c, d = far[nxt], far[index]
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    # Cap triangulation: the two legs as rectangles, wound outward (+x on the far cap).
    rectangles = ((0, 1, 2, 3), (0, 3, 4, 5))
    for i, j, k, m in rectangles:
        builder.add_triangle(far[i], far[j], far[k])
        builder.add_triangle(far[i], far[k], far[m])
        builder.add_triangle(near[i], near[k], near[j])
        builder.add_triangle(near[i], near[m], near[k])
    return builder.build()


def merge(*meshes: Mesh) -> Mesh:
    """Concatenate meshes into one surface (their components stay disjoint)."""
    builder = MeshBuilder()
    for mesh in meshes:
        index_map = [builder.add_vertex(point) for point in mesh.vertices]
        for face in mesh.faces:
            builder.add_triangle(index_map[face[0]], index_map[face[1]], index_map[face[2]])
    return builder.build()


def quadrant_plate(length: float = 1.0, thickness: float = 0.06, gap: float = 1e-6) -> Mesh:
    """Four plate quadrants, each its own solid, separated by ``gap``: a CAD import seam."""
    half = 0.5 * length
    quad = half - 0.5 * gap
    boxes = []
    for i in (0, 1):
        for j in (0, 1):
            x0 = 0.0 if i == 0 else half + 0.5 * gap
            y0 = 0.0 if j == 0 else half + 0.5 * gap
            boxes.append(box((x0, y0, 0.0), (x0 + quad, y0 + quad, thickness)))
    return merge(*boxes)