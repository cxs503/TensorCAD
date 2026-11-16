"""One-shot patch: add exact nearest-surface queries to TriangleGrid."""
from pathlib import Path

p = Path("tensorcad/spatial_index.py")
src = p.read_text()

# --- 1. dataclass -----------------------------------------------------------
anchor = '''@dataclass(frozen=True)
class RayHit:'''
surface_hit = '''@dataclass(frozen=True)
class SurfaceHit:
    """Nearest point on the mesh surface to a query point."""

    distance: float
    triangle: int
    point: Vec3

    def to_dict(self) -> dict[str, object]:
        return {
            "distance": self.distance,
            "triangle": self.triangle,
            "point": list(self.point),
        }


@dataclass(frozen=True)
class RayHit:'''
assert src.count(anchor) == 1
src = src.replace(anchor, surface_hit)

# --- 2. method --------------------------------------------------------------
anchor2 = '''    def cast_all(self, origin: Vec3, direction: Vec3) -> tuple[RayHit, ...]:'''
method = '''    def nearest_surface(self, point: Vec3) -> SurfaceHit | None:
        """Exact Euclidean distance from ``point`` to the mesh surface.

        Shell search over the same uniform grid the ray queries use: the cell
        holding the point is tested first, then successive Chebyshev shells, and
        the search stops as soon as the best distance found is no greater than
        the distance from the query point to the outside of the searched block.
        Every triangle is stored in *all* cells its bounding box touches, so an
        untested triangle lies wholly outside the searched block -- which is what
        makes that bound valid and the answer exact (not approximate, unlike a
        nearest-vertex or k-ring heuristic).

        Returns ``None`` only for an empty mesh, which the constructor rejects.
        """
        px, py, pz = point
        coords = self._coords
        starts = self._starts
        flat = self._flat
        nx, ny, nz = self._nx, self._ny, self._nz
        lo, cell, inv = self._lo, self._cell, self._inv_cell

        best = math.inf
        best_tri = -1
        best_q = (px, py, pz)

        def consider(k: int) -> None:
            nonlocal best, best_tri, best_q
            d2, qx, qy, qz = _closest_point_triangle(coords, k, px, py, pz)
            if d2 < best:
                best = d2
                best_tri = k
                best_q = (qx, qy, qz)

        # Oversized triangles live in no cell; they must be tested unconditionally.
        for k in self._oversized:
            consider(k)

        def cell_of(value: float, axis: int) -> int:
            c = int(math.floor((value - lo[axis]) * inv[axis]))
            if c < 0:
                return 0
            n = self._resolution[axis]
            return n - 1 if c >= n else c

        bx = cell_of(px, 0)
        by = cell_of(py, 1)
        bz = cell_of(pz, 2)

        r = 0
        while True:
            ix0 = bx - r if bx - r > 0 else 0
            ix1 = bx + r if bx + r < nx - 1 else nx - 1
            iy0 = by - r if by - r > 0 else 0
            iy1 = by + r if by + r < ny - 1 else ny - 1
            iz0 = bz - r if bz - r > 0 else 0
            iz1 = bz + r if bz + r < nz - 1 else nz - 1
            complete = (
                ix0 == 0 and ix1 == nx - 1
                and iy0 == 0 and iy1 == ny - 1
                and iz0 == 0 and iz1 == nz - 1
            )

            for iz in range(iz0, iz1 + 1):
                z_edge = iz == iz0 or iz == iz1
                zbase = nx * ny * iz
                for iy in range(iy0, iy1 + 1):
                    cellbase = zbase + nx * iy
                    if z_edge or iy == iy0 or iy == iy1:
                        xs = range(ix0, ix1 + 1)
                    elif ix1 > ix0:
                        xs = (ix0, ix1)
                    else:
                        xs = (ix0,)
                    for ix in xs:
                        c = cellbase + ix
                        for i in range(starts[c], starts[c + 1]):
                            consider(flat[i])

            if complete:
                break

            # Distance from the query point to the *outside* of the searched
            # block: no unexamined triangle can be closer than this.
            bound = math.inf
            if ix0 > 0:
                d = px - (lo[0] + ix0 * cell[0])
                bound = d * d if d > 0.0 else 0.0
            if ix1 < nx - 1:
                d = lo[0] + (ix1 + 1) * cell[0] - px
                dd = d * d if d > 0.0 else 0.0
                if dd < bound:
                    bound = dd
            if iy0 > 0:
                d = py - (lo[1] + iy0 * cell[1])
                dd = d * d if d > 0.0 else 0.0
                if dd < bound:
                    bound = dd
            if iy1 < ny - 1:
                d = lo[1] + (iy1 + 1) * cell[1] - py
                dd = d * d if d > 0.0 else 0.0
                if dd < bound:
                    bound = dd
            if iz0 > 0:
                d = pz - (lo[2] + iz0 * cell[2])
                dd = d * d if d > 0.0 else 0.0
                if dd < bound:
                    bound = dd
            if iz1 < nz - 1:
                d = lo[2] + (iz1 + 1) * cell[2] - pz
                dd = d * d if d > 0.0 else 0.0
                if dd < bound:
                    bound = dd
            if best <= bound:
                break
            r += 1

        if best_tri < 0:
            return None
        return SurfaceHit(distance=math.sqrt(best), triangle=best_tri, point=best_q)

    def cast_all(self, origin: Vec3, direction: Vec3) -> tuple[RayHit, ...]:'''
assert src.count(anchor2) == 1
src = src.replace(anchor2, method)

# --- 3. helpers -------------------------------------------------------------
anchor3 = '''def _intersect(
    coords: "array",'''
helpers = '''def _closest_point_triangle(
    coords: "array", k: int, px: float, py: float, pz: float
) -> tuple[float, float, float, float]:
    """Squared distance and closest point from ``p`` to triangle ``k``.

    Ericson's barycentric region test (*Real-Time Collision Detection* §5.1.5):
    one branch per Voronoi region of the triangle, no iteration, exact for
    non-degenerate triangles.  Degenerate (zero-area) triangles fall back to the
    closest of their three corner points.
    """
    base = 9 * k
    ax = coords[base]; ay = coords[base + 1]; az = coords[base + 2]
    bx = coords[base + 3]; by = coords[base + 4]; bz = coords[base + 5]
    cx = coords[base + 6]; cy = coords[base + 7]; cz = coords[base + 8]

    abx = bx - ax; aby = by - ay; abz = bz - az
    acx = cx - ax; acy = cy - ay; acz = cz - az
    apx = px - ax; apy = py - ay; apz = pz - az

    d1 = abx * apx + aby * apy + abz * apz
    d2 = acx * apx + acy * apy + acz * apz
    if d1 <= 0.0 and d2 <= 0.0:
        qx, qy, qz = ax, ay, az
        dx, dy, dz = apx, apy, apz
        return dx * dx + dy * dy + dz * dz, qx, qy, qz

    bpx = px - bx; bpy = py - by; bpz = pz - bz
    d3 = abx * bpx + aby * bpy + abz * bpz
    d4 = acx * bpx + acy * bpy + acz * bpz
    if d3 >= 0.0 and d4 <= d3:
        dx, dy, dz = bpx, bpy, bpz
        return dx * dx + dy * dy + dz * dz, bx, by, bz

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        denom = d1 - d3
        if denom != 0.0:
            v = d1 / denom
            qx = ax + abx * v; qy = ay + aby * v; qz = az + abz * v
            dx = px - qx; dy = py - qy; dz = pz - qz
            return dx * dx + dy * dy + dz * dz, qx, qy, qz

    cpx = px - cx; cpy = py - cy; cpz = pz - cz
    d5 = abx * cpx + aby * cpy + abz * cpz
    d6 = acx * cpx + acy * cpy + acz * cpz
    if d6 >= 0.0 and d5 <= d6:
        dx, dy, dz = cpx, cpy, cpz
        return dx * dx + dy * dy + dz * dz, cx, cy, cz

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        denom = d2 - d6
        if denom != 0.0:
            w = d2 / denom
            qx = ax + acx * w; qy = ay + acy * w; qz = az + acz * w
            dx = px - qx; dy = py - qy; dz = pz - qz
            return dx * dx + dy * dy + dz * dz, qx, qy, qz

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        denom = (d4 - d3) + (d5 - d6)
        if denom != 0.0:
            w = (d4 - d3) / denom
            qx = bx + (cx - bx) * w; qy = by + (cy - by) * w; qz = bz + (cz - bz) * w
            dx = px - qx; dy = py - qy; dz = pz - qz
            return dx * dx + dy * dy + dz * dz, qx, qy, qz

    total = va + vb + vc
    if total > 0.0:
        denom = 1.0 / total
        v = vb * denom
        w = vc * denom
        qx = ax + abx * v + acx * w
        qy = ay + aby * v + acy * w
        qz = az + abz * v + acz * w
        dx = px - qx; dy = py - qy; dz = pz - qz
        return dx * dx + dy * dy + dz * dz, qx, qy, qz

    # Degenerate triangle: closest corner wins.
    best_d2 = apx * apx + apy * apy + apz * apz
    best_q = (ax, ay, az)
    d2b = bpx * bpx + bpy * bpy + bpz * bpz
    if d2b < best_d2:
        best_d2 = d2b
        best_q = (bx, by, bz)
    d2c = cpx * cpx + cpy * cpy + cpz * cpz
    if d2c < best_d2:
        best_d2 = d2c
        best_q = (cx, cy, cz)
    return best_d2, best_q[0], best_q[1], best_q[2]


def brute_force_surface_distance(mesh: Mesh, point: Vec3) -> float:
    """Reference O(N) distance from ``point`` to the mesh surface.

    The cross-check counterpart of :meth:`TriangleGrid.nearest_surface`, kept in
    the library so the acceleration structure can be validated on any mesh.
    """
    px, py, pz = point
    coords = array("d")
    for face in mesh.faces:
        for index in face:
            v = mesh.vertices[index]
            coords.append(v[0]); coords.append(v[1]); coords.append(v[2])
    best = math.inf
    for k in range(len(mesh.faces)):
        d2 = _closest_point_triangle(coords, k, px, py, pz)[0]
        if d2 < best:
            best = d2
    return math.sqrt(best)


def _intersect(
    coords: "array",'''
assert src.count(anchor3) == 1
src = src.replace(anchor3, helpers)

p.write_text(src)
print("patched", p, len(src.splitlines()), "lines")