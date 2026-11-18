"""Regression for the O(components x edges) closed-component scan in audit_mesh.

The old implementation recomputed, for every component, the set of edges touching
it by scanning *all* edges and testing set membership per edge use.  On a
fragmented mesh that is O(C x E): the ONRT publisher file (132 faces, 1794
components, ~240k edges) spent the bulk of its 414 s import wall clock there
(cluster job cj5190 profile: 947 s of 1068 s profiled time in audit_mesh).

The replacement is a single pass over the edges: a component is closed iff no
edge touching it has a use count other than two.  These tests pin both the
semantics (closed_component_count on a mesh with a known answer) and the cost
(the old O(C x E) scan must not come back).
"""
from __future__ import annotations

import time

from tensorcad.geomquality.trim_mesh import audit_mesh


def _subdivided_box(n: int, origin: tuple[float, float, float] = (0.0, 0.0, 0.0)):
    """A closed box [0,1]^3: (n+1)^3 shared grid vertices, boundary triangles only.

    Every box edge is shared by exactly two triangles, so the box is one closed
    component regardless of the (outward) orientation chosen per face.
    """
    ox, oy, oz = origin
    vertices = [(ox + i / n, oy + j / n, oz + k / n) for i in range(n + 1)
                for j in range(n + 1) for k in range(n + 1)]

    def idx(i: int, j: int, k: int) -> int:
        return (i * (n + 1) + j) * (n + 1) + k

    faces: list[tuple[int, int, int]] = []
    # (fixed axis, fixed value, u axis, v axis, outward sign)
    for axis, fixed, u, v, sign in (
        (2, 0, 0, 1, -1), (2, 1, 0, 1, 1),
        (0, 0, 1, 2, -1), (0, 1, 1, 2, 1),
        (1, 0, 0, 2, -1), (1, 1, 0, 2, 1),
    ):
        for iu in range(n):
            for iv in range(n):
                def pt(a: int, b: int) -> int:
                    c = [0, 0, 0]
                    c[axis] = fixed
                    c[u] = iu + a
                    c[v] = iv + b
                    return idx(*c)
                a, b, c, d = pt(0, 0), pt(1, 0), pt(1, 1), pt(0, 1)
                if sign < 0:
                    faces.append((a, b, c))
                    faces.append((a, c, d))
                else:
                    faces.append((a, c, b))
                    faces.append((a, d, c))
    return vertices, faces


def test_closed_component_count_on_mixed_mesh():
    """One closed box + three open triangles: exactly one closed component."""
    vertices, faces = _subdivided_box(4)
    base = len(vertices)
    for k in range(3):
        x = 10.0 + 5.0 * k
        vertices += [(x, 0.0, 0.0), (x + 1.0, 0.0, 0.0), (x, 1.0, 0.0)]
        faces += ((base + 3 * k, base + 3 * k + 1, base + 3 * k + 2),)
    audit = audit_mesh(vertices, faces)
    assert audit.component_count == 4
    assert audit.closed_component_count == 1
    assert audit.free_edge_count == 9  # 3 open triangles x 3 free edges
    assert audit.watertight is False


def test_closed_component_count_many_components_is_fast():
    """The old O(C x E) scan must not come back.

    200 disjoint closed boxes (200 components) plus one open triangle: the old
    code paid 201 x E edge scans with a Python-level membership test per edge
    use; the single-pass replacement is O(E).  A 10 s budget is generous by
    orders of magnitude for the new code and unreachable for the old one.
    """
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for k in range(200):
        box_v, box_f = _subdivided_box(2, origin=(0.0, 0.0, 100.0 * k))
        offset = len(vertices)
        vertices.extend(box_v)
        faces.extend(tuple(i + offset for i in tri) for tri in box_f)
    vertices += [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    faces.append((len(vertices) - 3, len(vertices) - 2, len(vertices) - 1))

    started = time.time()
    audit = audit_mesh(vertices, faces)
    elapsed = time.time() - started
    assert audit.component_count == 201
    assert audit.closed_component_count == 200
    assert audit.free_edge_count == 3
    assert elapsed < 10.0, f"audit_mesh took {elapsed:.1f} s: O(C x E) scan regression?"
