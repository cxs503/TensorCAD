"""Independent exact reference for the shell tests: 2-D half-plane erosion of a convex section.

The module under test offsets a *surface*: per vertex, it fits the incident facet planes in 3-D
and solves a 3x3 least-squares system. This oracle never looks at a vertex, a plane system or a
least-squares solve. For a prism it uses the other definition of the same object -- the eroded
solid is the intersection of the offset half-planes -- and computes it with Sutherland-Hodgman
clipping in the section plane, then multiplies by the eroded extrusion length. Two independent
routes to one number, which is what makes the agreement in ``test_shell.py`` worth something.

Only convex sections are admissible (clipping is exact for them); a concave section would need
the offset surfaces to be trimmed against each other, which is the very collision problem the
module reports rather than solves.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

Point2 = tuple[float, float]


def _area(polygon: Sequence[Point2]) -> float:
    total = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        total += y0 * z1 - y1 * z0
    return 0.5 * total


def _clip(polygon: Sequence[Point2], normal: Point2, bound: float) -> list[Point2]:
    """Keep the part of the (convex) polygon with ``normal . p <= bound``."""
    result: list[Point2] = []
    count = len(polygon)
    for index in range(count):
        current = polygon[index]
        following = polygon[(index + 1) % count]
        inside_current = normal[0] * current[0] + normal[1] * current[1] - bound
        inside_following = normal[0] * following[0] + normal[1] * following[1] - bound
        if inside_current <= 0.0:
            result.append(current)
        if (inside_current < 0.0 < inside_following) or (inside_following < 0.0 < inside_current):
            ratio = inside_current / (inside_current - inside_following)
            result.append(
                (
                    current[0] + ratio * (following[0] - current[0]),
                    current[1] + ratio * (following[1] - current[1]),
                )
            )
    return result


def eroded_section(section: Sequence[Point2], thickness: float) -> list[Point2]:
    """Section clipped by every edge's inward-offset half-plane, in the section's own order."""
    polygon = list(section)
    count = len(section)
    for index in range(count):
        start = section[index]
        end = section[(index + 1) % count]
        edge_y, edge_z = end[0] - start[0], end[1] - start[1]
        length = math.hypot(edge_y, edge_z)
        # A counter-clockwise section has (ey, -ex) as the outward normal of this edge.
        normal = (edge_z / length, -edge_y / length)
        bound = normal[0] * start[0] + normal[1] * start[1] - thickness
        polygon = _clip(polygon, normal, bound)
        if len(polygon) < 3:
            raise ValueError(f"the section is eroded away at thickness {thickness!r}")
    return polygon


def eroded_prism_volume(section: Sequence[Point2], length: float, thickness: float) -> float:
    """Volume of ``section x [0, length]`` eroded by ``thickness`` on all faces.

    Exact for a convex section: the eroded solid is the product of the eroded section and the
    shortened extrusion ``length - 2 * thickness``.
    """
    if length - 2.0 * thickness <= 0.0:
        raise ValueError("the extrusion is shorter than twice the thickness")
    return _area(eroded_section(section, thickness)) * (length - 2.0 * thickness)


def section_area(section: Sequence[Point2]) -> float:
    return abs(_area(section))
