"""Exact-arithmetic oracle for triangle--triangle intersection (verification only).

Independent of `tensorcad.intersect`: the library decides intersection by testing
whether a closed edge of one triangle pierces the other, while this oracle clips
each triangle against the *other triangle's plane* and intersects the resulting
parameter intervals. Coplanar pairs are clipped as 2-D convex polygons. All
coordinates are `Fraction`, so the oracle is exact and cannot share a rounding
failure with the floating-point implementation under test.

Conventions that mirror the library contract:
  * a triangle with zero area carries no intersection (degenerate input);
  * triangles are closed sets, so touching at a point or along a segment counts.
"""
from __future__ import annotations

from fractions import Fraction

Tri = tuple[tuple[Fraction, Fraction, Fraction], ...]


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _is_degenerate(tri) -> bool:
    n = _cross(_sub(tri[1], tri[0]), _sub(tri[2], tri[0]))
    return n == (0, 0, 0)


def _canonical(normal, scale):
    """Unit-ish normal with the dominant component positive (sign convention)."""
    unit = tuple(Fraction(c) / scale for c in normal)
    dominant = max(range(3), key=lambda i: abs(unit[i]))
    if unit[dominant] < 0:
        unit = (-unit[0], -unit[1], -unit[2])
    return unit


def _same_plane(normal, offset, tri) -> bool:
    """True when every vertex of `tri` lies on {normal . p == offset} exactly."""
    return all(_dot(normal, point) == offset for point in tri)


def _clip_interval(lo, hi, coeff, rhs):
    """Intersect [lo, hi] with {t : coeff * t >= rhs}; None when empty."""
    if coeff == 0:
        return (lo, hi) if Fraction(0) >= rhs else None
    bound = Fraction(rhs, 1) / coeff
    if coeff > 0:
        return (max(lo, bound), hi) if max(lo, bound) <= hi else None
    return (lo, min(hi, bound)) if lo <= min(hi, bound) else None


def _line_in_triangle(point, direction, tri):
    """Parameter interval of {point + t*direction} inside `tri` (exact)."""
    normal = _cross(_sub(tri[1], tri[0]), _sub(tri[2], tri[0]))
    lo, hi = None, None
    for i in range(3):
        a, b = tri[i], tri[(i + 1) % 3]
        # Edge a->b: inside half-plane is cross(b-a, p-a) . normal >= 0.
        edge = _sub(b, a)
        w = _cross(edge, _sub(point, a))
        coeff = _dot(_cross(edge, direction), normal)
        rhs = -_dot(w, normal)
        if coeff == 0 and rhs > 0:
            return None
        if coeff != 0:
            bound = Fraction(rhs, 1) / coeff
            if coeff > 0:
                lo = bound if lo is None else max(lo, bound)
            else:
                hi = bound if hi is None else min(hi, bound)
    if lo is None:
        lo = Fraction(-1)
    if hi is None:
        hi = Fraction(1)
    if lo > hi:
        return None
    return (lo, hi)


def _project(tri, axis):
    u, v = [(1, 2), (2, 0), (0, 1)][axis]
    return [(p[u], p[v]) for p in tri]


def _clip_polygon(subject, clipper):
    """Clip convex polygon `subject` by convex polygon `clipper` (2-D, exact)."""
    output = list(subject)
    n = len(clipper)
    # Orientation of the clipper decides which side is interior.
    area = Fraction(0)
    for i in range(n):
        (x0, y0), (x1, y1) = clipper[i], clipper[(i + 1) % n]
        area += x0 * y1 - x1 * y0
    sign = 1 if area > 0 else -1
    for i in range(n):
        if not output:
            return []
        a, b = clipper[i], clipper[(i + 1) % n]
        edge = (b[0] - a[0], b[1] - a[1])

        def side(p):
            return sign * (edge[0] * (p[1] - a[1]) - edge[1] * (p[0] - a[0]))

        new = []
        m = len(output)
        for j in range(m):
            cur = output[j]
            nxt = output[(j + 1) % m]
            sc, sn = side(cur), side(nxt)
            if sc >= 0:
                new.append(cur)
            if (sc > 0 and sn < 0) or (sc < 0 and sn > 0):
                t = Fraction(sc, 1) / (sc - sn)
                new.append((cur[0] + t * (nxt[0] - cur[0]), cur[1] + t * (nxt[1] - cur[1])))
        output = new
    return output


def triangles_intersect_exact(first, second) -> bool:
    """Exact answer for 'do these two closed triangles share a point'."""
    t1 = tuple(tuple(Fraction(c) for c in p) for p in first)
    t2 = tuple(tuple(Fraction(c) for c in p) for p in second)
    if _is_degenerate(t1) or _is_degenerate(t2):
        return False
    n1 = _cross(_sub(t1[1], t1[0]), _sub(t1[2], t1[0]))
    n2 = _cross(_sub(t2[1], t2[0]), _sub(t2[2], t2[0]))
    # Plane offsets, sign-canonicalised so coplanar planes compare equal.
    d1 = _dot(n1, t1[0])
    d2 = _dot(n2, t2[0])
    scale1 = max(abs(c) for c in n1)
    scale2 = max(abs(c) for c in n2)
    n1c = _canonical(n1, scale1)
    n2c = _canonical(n2, scale2)
    parallel = _cross(n1c, n2c) == (0, 0, 0)
    if parallel and _cross(n1c, n2c) == (0, 0, 0) and _same_plane(n1, _dot(n1, t1[0]), t2):
        axis = max(range(3), key=lambda i: abs(n1c[i]))
        p1 = _project(t1, axis)
        p2 = _project(t2, axis)
        return bool(_clip_polygon(p1, p2))
    if parallel:
        return False
    # Direction of the plane-plane intersection line.
    direction = _cross(n1, n2)
    # A point on both planes: solve for the point closest to the origin in the
    # rational sense -- pick a coordinate to eliminate and solve the 2x2 system.
    axis = max(range(3), key=lambda i: abs(direction[i]))
    keep = [i for i in range(3) if i != axis]
    u, v = keep
    det = n1[u] * n2[v] - n1[v] * n2[u]
    if det == 0:  # pragma: no cover - parallel already excluded
        return False
    point = [Fraction(0), Fraction(0), Fraction(0)]
    point[u] = Fraction(d1 * n2[v] - d2 * n1[v], det)
    point[v] = Fraction(n1[u] * d2 - n2[u] * d1, det)
    point = tuple(point)
    span1 = _line_in_triangle(point, direction, t1)
    if span1 is None:
        return False
    span2 = _line_in_triangle(point, direction, t2)
    if span2 is None:
        return False
    lo = max(span1[0], span2[0])
    hi = min(span1[1], span2[1])
    return lo <= hi

# ---------------------------------------------------------------------------
# Second exact route: 3-D half-space clipping (Sutherland--Hodgman).
#
# The two routes above and this one partition the problem differently: route 1
# clips the plane/plane line against each triangle, route 2 clips each triangle
# as a polygon against the other's three half-spaces. They agree on every
# exactly-representable (dyadic) pair tested, which is what makes them usable as
# a reference for the floating-point implementation. On inputs that are *not*
# exactly representable (thirds) they can disagree on near-degenerate
# configurations; the library is compared against this route there, because
# clipping both triangles against both sets of half-spaces is the more direct
# statement of the question.
# ---------------------------------------------------------------------------


def _clip_by_halfspace(poly, normal, offset):
    out = []
    count = len(poly)
    if count == 0:
        return out
    distances = [_dot(normal, point) - offset for point in poly]
    for i in range(count):
        cur, nxt = poly[i], poly[(i + 1) % count]
        d_cur, d_next = distances[i], distances[(i + 1) % count]
        if d_cur <= 0:
            out.append(cur)
        if (d_cur < 0 < d_next) or (d_cur > 0 > d_next):
            t = d_cur / (d_cur - d_next)
            out.append(tuple(cur[k] + t * (nxt[k] - cur[k]) for k in range(3)))
    return out


def _clip_triangle(subject, clipper):
    normal = _cross(_sub(clipper[1], clipper[0]), _sub(clipper[2], clipper[0]))
    poly = list(subject)
    for i in range(3):
        a, b = clipper[i], clipper[(i + 1) % 3]
        inward = _cross(normal, _sub(b, a))
        poly = _clip_by_halfspace(poly, inward, _dot(inward, a))
        if not poly:
            return []
    return poly


def triangles_intersect_clip3d(first, second) -> bool:
    """Exact answer by clipping each triangle against the other's half-spaces."""
    t1 = tuple(tuple(Fraction(c) for c in p) for p in first)
    t2 = tuple(tuple(Fraction(c) for c in p) for p in second)
    if _is_degenerate(t1) or _is_degenerate(t2):
        return False
    return bool(_clip_triangle(t1, t2)) or bool(_clip_triangle(t2, t1))
