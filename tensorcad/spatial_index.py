"""Uniform-grid spatial index over a triangle mesh: ray casting and point location.

Why this exists
---------------
:class:`tensorcad.benchmark._RayCaster` casts rays by testing *every* triangle,
which is fine for a 2k-triangle smoke mesh and unusable at 10^6--10^7 triangles
(GAP-CAD-D4-PERF-SCALE measured 13 queries/s at 1.0e5 triangles, i.e. a
ten-thousand-ray workload would take a quarter of an hour).  A uniform grid with
3-D DDA traversal reduces a query to ``O(cells crossed)`` triangle tests instead
of ``O(N)``, which is what makes point location at 10^6--10^7 triangles
practical on one machine.

Design
------
* **CSR buckets, no Python object lists.**  Triangle references live in a flat
  :class:`array.array` of 32-bit ints addressed by per-cell prefix sums (two-pass
  counting sort).  At 10^7 triangles that is hundreds of MB instead of the
  multi-GB a ``dict[cell] -> list`` layout would need.
* **Near-cubic cells.**  The cell count along each axis is proportional to the
  bounding-box extent, so cells stay near-cubic for an elongated hull box
  instead of degenerating into slivers.
* **Oversized triangles go to a fallback list.**  A triangle whose bbox would
  touch more than ``oversized_limit`` cells is *not* copied into every cell
  (that is how one huge triangle turns a grid into a memory bomb); it is tested
  for every ray instead, which keeps the answer exact.
* **Same geometry as a brute-force test.**  ``cast_ray`` is Möller--Trumbore in
  double precision and returns the nearest positive intersection, so it agrees
  with a brute-force loop to rounding error.  ``tests/test_spatial_index.py``
  pins that equivalence, because an acceleration structure that changes answers
  is worthless.
* **Standard library only.**  No numpy, consistent with the rest of the core.

Conventions
-----------
* Ray directions need not be normalised; distances are in mesh units either way.
* A ray starting inside a closed mesh hits its far side first; use
  :meth:`TriangleGrid.contains` for inside/outside tests.
* ``contains`` is undefined (may return either) for points exactly on the
  surface, as for any parity test in floating point.

Usage
-----
::

    from tensorcad.spatial_index import TriangleGrid
    grid = TriangleGrid(mesh)                       # build; build_seconds recorded
    hit = grid.cast_ray((0.0, 0.0, 0.0), direction)  # RayHit | None
    inside = grid.contains((0.5, 0.0, -0.1))         # bool
"""
from __future__ import annotations

import math
from array import array
from dataclasses import dataclass
from typing import Iterator, Sequence

from .geometry import Mesh, Vec3

#: Grid resolutions are clamped to this range so a sliver mesh cannot allocate
#: an absurd number of cells.
MIN_CELLS_PER_AXIS = 1
MAX_CELLS_PER_AXIS = 2048

#: Default cap on total grid cells (32M cells = 128 MB of int32 counters).
DEFAULT_MAX_CELLS = 32_000_000

#: A triangle whose bbox covers more cells than this stays in the fallback list
#: instead of being spread over the grid.
DEFAULT_OVERSIZED_LIMIT = 4096

#: Geometric tolerances.  Deliberately tiny: the point of the acceleration
#: structure is to return the *same* answer as the brute-force test, so no slack
#: is introduced beyond double-precision rounding.
_EDGE_EPS = 1e-12
_PARALLEL_EPS = 1e-14
_MIN_T = 1e-12

#: Probe direction used by :meth:`TriangleGrid.contains`.  Deliberately
#: irrational-ish so a query point never lines up with a mesh vertex or edge.
_CONTAIN_PROBE = (1.0, 0.0013, 0.0007)


@dataclass(frozen=True)
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
class RayHit:
    """Nearest ray--triangle intersection."""

    distance: float
    triangle: int
    point: Vec3

    def to_dict(self) -> dict[str, object]:
        return {"distance": self.distance, "triangle": self.triangle, "point": list(self.point)}


class TriangleGrid:
    """Uniform grid over the triangles of ``mesh`` with DDA ray traversal.

    Parameters
    ----------
    mesh:
        Any :class:`~tensorcad.geometry.Mesh`; vertices may be shared.
    resolution:
        Optional ``(nx, ny, nz)`` cell counts.  When omitted the resolution is
        chosen so an *average* triangle bbox touches about
        ``target_triangles_per_cell`` cells, then clamped by ``max_cells``.
    target_triangles_per_cell:
        Occupancy target for the automatic resolution.
    max_cells:
        Upper bound on ``nx * ny * nz``.
    oversized_limit:
        Triangles covering more cells than this are tested by every ray.
    """

    def __init__(
        self,
        mesh: Mesh,
        resolution: tuple[int, int, int] | None = None,
        target_triangles_per_cell: float = 4.0,
        max_cells: int = DEFAULT_MAX_CELLS,
        oversized_limit: int = DEFAULT_OVERSIZED_LIMIT,
    ) -> None:
        n_tris = len(mesh.faces)
        if n_tris == 0:
            raise ValueError("cannot index a mesh with no triangles")
        self.mesh = mesh
        self.n_triangles = n_tris

        xs = [v[0] for v in mesh.vertices]
        ys = [v[1] for v in mesh.vertices]
        zs = [v[2] for v in mesh.vertices]
        lo = (min(xs), min(ys), min(zs))
        hi = (max(xs), max(ys), max(zs))
        # Pad so surface triangles never sit exactly on the outer planes.
        ext = [max(hi[i] - lo[i], 1e-12) for i in range(3)]
        pad = [1e-9 * e for e in ext]
        self._lo = tuple(lo[i] - pad[i] for i in range(3))
        self._ext = tuple(ext[i] + 2.0 * pad[i] for i in range(3))
        self._hi = tuple(self._lo[i] + self._ext[i] for i in range(3))

        if resolution is None:
            self._resolution = self._auto_resolution(n_tris, target_triangles_per_cell, max_cells)
        else:
            self._resolution = self._clamp_resolution(resolution)
        nx, ny, nz = self._resolution
        self._nx, self._ny, self._nz = nx, ny, nz
        self._cell = tuple(self._ext[i] / self._resolution[i] for i in range(3))
        self._inv_cell = tuple(1.0 / c for c in self._cell)
        self._n_cells = nx * ny * nz

        # Flattened vertex triples, 9 doubles per triangle: the query inner loop
        # touches this array directly, which is cheaper than chasing
        # mesh.vertices[...] tuples per intersection test.
        self._coords = array("d")
        append = self._coords.append
        vertex_lookup = mesh.vertices
        for face in mesh.faces:
            a = vertex_lookup[face[0]]
            b = vertex_lookup[face[1]]
            c = vertex_lookup[face[2]]
            append(a[0]); append(a[1]); append(a[2])
            append(b[0]); append(b[1]); append(b[2])
            append(c[0]); append(c[1]); append(c[2])

        self._oversized: list[int] = []
        self._starts: array = array("i")
        self._flat: array = array("i")
        self._counts: array = array("i", bytes(4 * self._n_cells))
        self.build_seconds: float = 0.0
        self._build(oversized_limit)

    # -- construction -------------------------------------------------------

    @staticmethod
    def _clamp_resolution(resolution: Sequence[int]) -> tuple[int, int, int]:
        return (
            max(MIN_CELLS_PER_AXIS, min(MAX_CELLS_PER_AXIS, int(resolution[0]))),
            max(MIN_CELLS_PER_AXIS, min(MAX_CELLS_PER_AXIS, int(resolution[1]))),
            max(MIN_CELLS_PER_AXIS, min(MAX_CELLS_PER_AXIS, int(resolution[2]))),
        )

    def _auto_resolution(
        self, n_tris: int, target: float, max_cells: int
    ) -> tuple[int, int, int]:
        """Cell counts proportional to the bbox extents, sized by occupancy."""
        target = max(1.0, float(target))
        volume_budget = max(1.0, n_tris / target)
        ex, ey, ez = self._ext
        h = (ex * ey * ez / volume_budget) ** (1.0 / 3.0)
        res = self._clamp_resolution(
            tuple(max(1, int(math.ceil(e / h))) for e in (ex, ey, ez))
        )
        while res[0] * res[1] * res[2] > max_cells:
            res = self._clamp_resolution(tuple(max(1, int(round(r * 0.8))) for r in res))
            if res == (1, 1, 1):
                break
        return res

    def _cell_index_range(self, axis: int, vmin: float, vmax: float) -> tuple[int, int]:
        lo = self._lo[axis]
        inv = self._inv_cell[axis]
        n = self._resolution[axis]
        i0 = int(math.floor((vmin - lo) * inv))
        i1 = int(math.floor((vmax - lo) * inv))
        if i0 < 0:
            i0 = 0
        if i1 >= n:
            i1 = n - 1
        if i1 < i0:
            i1 = i0
        return i0, i1

    def _build(self, oversized_limit: int) -> None:
        import time

        t0 = time.perf_counter()
        coords = self._coords
        ranges = array("i")
        counts = self._counts
        nx, ny = self._nx, self._ny
        n_tris = self.n_triangles

        for k in range(n_tris):
            base = 9 * k
            x0 = coords[base]; x1 = coords[base + 3]; x2 = coords[base + 6]
            y0 = coords[base + 1]; y1 = coords[base + 4]; y2 = coords[base + 7]
            z0 = coords[base + 2]; z1 = coords[base + 5]; z2 = coords[base + 8]
            ix0, ix1 = self._cell_index_range(0, min(x0, x1, x2), max(x0, x1, x2))
            iy0, iy1 = self._cell_index_range(1, min(y0, y1, y2), max(y0, y1, y2))
            iz0, iz1 = self._cell_index_range(2, min(z0, z1, z2), max(z0, z1, z2))
            ranges.extend((ix0, ix1, iy0, iy1, iz0, iz1))
            if (ix1 - ix0 + 1) * (iy1 - iy0 + 1) * (iz1 - iz0 + 1) > oversized_limit:
                self._oversized.append(k)
                continue
            for iz in range(iz0, iz1 + 1):
                zbase = nx * ny * iz
                for iy in range(iy0, iy1 + 1):
                    cell = zbase + nx * iy
                    for ix in range(ix0, ix1 + 1):
                        counts[cell + ix] += 1

        starts = array("i", bytes(4 * (self._n_cells + 1)))
        total = 0
        for cell in range(self._n_cells):
            starts[cell] = total
            total += counts[cell]
        starts[self._n_cells] = total

        flat = array("i", bytes(4 * total))
        cursor = array("i", starts)
        oversized = set(self._oversized)
        r = ranges
        for k in range(n_tris):
            if k in oversized:
                continue
            base = 6 * k
            ix0 = r[base]; ix1 = r[base + 1]
            iy0 = r[base + 2]; iy1 = r[base + 3]
            iz0 = r[base + 4]; iz1 = r[base + 5]
            for iz in range(iz0, iz1 + 1):
                zbase = nx * ny * iz
                for iy in range(iy0, iy1 + 1):
                    cell = zbase + nx * iy
                    for ix in range(ix0, ix1 + 1):
                        c = cell + ix
                        flat[cursor[c]] = k
                        cursor[c] += 1

        self._starts = starts
        self._flat = flat
        self.build_seconds = time.perf_counter() - t0

    # -- introspection ------------------------------------------------------

    @property
    def resolution(self) -> tuple[int, int, int]:
        return self._resolution

    @property
    def n_cells(self) -> int:
        return self._n_cells

    @property
    def n_entries(self) -> int:
        """Total triangle references stored in cells (duplicates included)."""
        return len(self._flat)

    @property
    def n_oversized(self) -> int:
        return len(self._oversized)

    def stats(self) -> dict[str, object]:
        occupied = 0
        max_per_cell = 0
        for c in self._counts:
            if c:
                occupied += 1
                if c > max_per_cell:
                    max_per_cell = c
        return {
            "triangles": self.n_triangles,
            "resolution": list(self._resolution),
            "cells": self._n_cells,
            "entries": self.n_entries,
            "entries_per_triangle": (
                self.n_entries / self.n_triangles if self.n_triangles else 0.0
            ),
            "occupied_cells": occupied,
            "mean_triangles_per_occupied_cell": (
                self.n_entries / occupied if occupied else 0.0
            ),
            "max_triangles_per_cell": max_per_cell,
            "oversized_triangles": self.n_oversized,
            "index_bytes": (
                4 * self._n_cells + 4 * (self._n_cells + 1)
                + 4 * self.n_entries + 8 * 9 * self.n_triangles
            ),
            "build_seconds": self.build_seconds,
        }

    # -- queries ------------------------------------------------------------

    def cast_ray(
        self,
        origin: Vec3,
        direction: Vec3,
        max_distance: float | None = None,
    ) -> RayHit | None:
        """Nearest intersection along ``origin + t * direction`` (``t > 0``)."""
        ox, oy, oz = origin
        dx, dy, dz = direction
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        if norm == 0.0:
            raise ValueError("ray direction must be non-zero")
        dx /= norm
        dy /= norm
        dz /= norm

        clip = self._clip_to_bounds(ox, oy, oz, dx, dy, dz)
        if clip is None:
            return None
        t_enter, t_exit = clip
        limit = math.inf if max_distance is None else max_distance
        if t_enter > limit:
            return None
        if t_exit > limit:
            t_exit = limit

        best_t = min(limit, math.inf)
        best_tri = -1
        coords = self._coords
        starts = self._starts
        flat = self._flat

        for k in self._oversized:
            t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, best_t)
            if t is not None:
                best_t = t
                best_tri = k

        for cell, t_cell_exit in self._cells_crossed(ox, oy, oz, dx, dy, dz, t_enter, t_exit):
            for i in range(starts[cell], starts[cell + 1]):
                k = flat[i]
                t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, best_t)
                if t is not None:
                    best_t = t
                    best_tri = k
            if best_t <= t_cell_exit:
                break

        if best_tri < 0:
            return None
        return RayHit(
            distance=best_t,
            triangle=best_tri,
            point=(ox + dx * best_t, oy + dy * best_t, oz + dz * best_t),
        )

    def nearest_distance(self, origin: Vec3, direction: Vec3) -> float | None:
        """First hit distance, or ``None`` when the ray misses."""
        hit = self.cast_ray(origin, direction)
        return None if hit is None else hit.distance

    def nearest_surface(self, point: Vec3, *, band: float | None = None) -> SurfaceHit | None:
        """Exact Euclidean distance from ``point`` to the mesh surface.

        Shell search over the same uniform grid the ray queries use: the cell
        holding the point is tested first, then successive Chebyshev shells, and
        the search stops as soon as the best distance found is no greater than
        the distance from the query point to the outside of the searched block.
        Every triangle is stored in *all* cells its bounding box touches, so an
        untested triangle lies wholly outside the searched block -- which is what
        makes that bound valid and the answer exact (not approximate, unlike a
        nearest-vertex or k-ring heuristic).

        ``band`` turns the query into a *narrow-band* one: when no triangle lies
        within ``band`` of the point the search stops and returns
        ``SurfaceHit(distance=band, triangle=-1)`` instead of the (much larger)
        true distance.  That is the only way this stays affordable on a 10^7
        triangle mesh, where an interior query with a 0.5 m answer would
        otherwise have to clear a 0.5 m ball of triangles; the truncation is
        sign-agnostic and the caller keeps the sign test.  ``triangle == -1``
        marks a truncated result, so a caller can tell "no triangle within the
        band" from "the nearest triangle is at distance band".

        Returns ``None`` only for an empty mesh, which the constructor rejects.
        """
        band2 = None if band is None else float(band) ** 2
        if band2 is not None and band2 <= 0.0:
            raise ValueError(f"band must be positive, got {band!r}")
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
            if band2 is not None and bound > band2:
                # Nothing unexamined can be within the band.
                if best <= band2:
                    break
                return SurfaceHit(distance=math.sqrt(band2), triangle=-1, point=(px, py, pz))
            r += 1

        if best_tri < 0:
            return None
        return SurfaceHit(distance=math.sqrt(best), triangle=best_tri, point=best_q)

    def cast_all(self, origin: Vec3, direction: Vec3) -> tuple[RayHit, ...]:
        """All intersections along the ray, sorted by distance.

        Duplicate hits at shared edges/vertices (the same geometric crossing
        reported by two adjacent triangles) are collapsed.
        """
        ox, oy, oz = origin
        dx, dy, dz = direction
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        if norm == 0.0:
            raise ValueError("ray direction must be non-zero")
        dx /= norm
        dy /= norm
        dz /= norm
        clip = self._clip_to_bounds(ox, oy, oz, dx, dy, dz)
        if clip is None:
            return ()
        t_enter, t_exit = clip
        starts = self._starts
        flat = self._flat
        coords = self._coords

        found: list[tuple[float, int]] = []
        for k in self._oversized:
            t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, math.inf)
            if t is not None:
                found.append((t, k))
        for cell, _t_exit in self._cells_crossed(ox, oy, oz, dx, dy, dz, t_enter, t_exit):
            for i in range(starts[cell], starts[cell + 1]):
                k = flat[i]
                t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, math.inf)
                if t is not None:
                    found.append((t, k))

        found.sort()
        deduped: list[tuple[float, int]] = []
        for t, k in found:
            if deduped and abs(t - deduped[-1][0]) <= 1e-9 * max(1.0, abs(t)):
                continue
            deduped.append((t, k))
        return tuple(
            RayHit(distance=t, triangle=k, point=(ox + dx * t, oy + dy * t, oz + dz * t))
            for t, k in deduped
        )

    def contains(self, point: Vec3) -> bool:
        """Inside/outside test by ray-crossing parity.

        Undefined (may return either) for points lying exactly on the surface.
        """
        return len(self.cast_all(point, _CONTAIN_PROBE)) % 2 == 1

    def query_many(
        self,
        origin: Vec3,
        directions: Sequence[Vec3],
        max_distance: float | None = None,
    ) -> list[RayHit | None]:
        """Batch wrapper: one query per direction, in order."""
        return [self.cast_ray(origin, d, max_distance) for d in directions]

    # -- traversal helpers --------------------------------------------------

    def _clip_to_bounds(
        self, ox: float, oy: float, oz: float,
        dx: float, dy: float, dz: float,
    ) -> tuple[float, float] | None:
        """Slab test; returns the parametric span inside the grid bounds."""
        t_enter = 0.0
        t_exit = math.inf
        for axis, (o, d) in enumerate(((ox, dx), (oy, dy), (oz, dz))):
            lo = self._lo[axis]
            hi = self._hi[axis]
            if d == 0.0:
                if o < lo or o > hi:
                    return None
                continue
            t1 = (lo - o) / d
            t2 = (hi - o) / d
            if t1 > t2:
                t1, t2 = t2, t1
            if t1 > t_enter:
                t_enter = t1
            if t2 < t_exit:
                t_exit = t2
            if t_enter > t_exit:
                return None
        return t_enter, t_exit

    def _cells_crossed(
        self, ox: float, oy: float, oz: float,
        dx: float, dy: float, dz: float,
        t_enter: float, t_exit: float,
    ) -> Iterator[tuple[int, float]]:
        """Yield ``(cell_index, t_leaving_cell)`` in DDA order (Amanatides-Woo)."""
        nx, ny, nz = self._nx, self._ny, self._nz
        lo = self._lo
        inv = self._inv_cell
        res = self._resolution
        cell_size = self._cell

        t_start = t_enter
        px = ox + dx * t_start
        py = oy + dy * t_start
        pz = oz + dz * t_start

        def cell_of(p: float, axis: int) -> int:
            i = int(math.floor((p - lo[axis]) * inv[axis]))
            if i < 0:
                return 0
            if i >= res[axis]:
                return res[axis] - 1
            return i

        ix = cell_of(px, 0)
        iy = cell_of(py, 1)
        iz = cell_of(pz, 2)

        step_x = 0 if dx == 0.0 else (1 if dx > 0.0 else -1)
        step_y = 0 if dy == 0.0 else (1 if dy > 0.0 else -1)
        step_z = 0 if dz == 0.0 else (1 if dz > 0.0 else -1)

        if step_x == 0:
            t_max_x = math.inf
            t_delta_x = math.inf
        else:
            boundary = lo[0] + (ix + 1) * cell_size[0] if step_x > 0 else lo[0] + ix * cell_size[0]
            t_max_x = (boundary - ox) / dx
            t_delta_x = cell_size[0] / abs(dx)
        if step_y == 0:
            t_max_y = math.inf
            t_delta_y = math.inf
        else:
            boundary = lo[1] + (iy + 1) * cell_size[1] if step_y > 0 else lo[1] + iy * cell_size[1]
            t_max_y = (boundary - oy) / dy
            t_delta_y = cell_size[1] / abs(dy)
        if step_z == 0:
            t_max_z = math.inf
            t_delta_z = math.inf
        else:
            boundary = lo[2] + (iz + 1) * cell_size[2] if step_z > 0 else lo[2] + iz * cell_size[2]
            t_max_z = (boundary - oz) / dz
            t_delta_z = cell_size[2] / abs(dz)

        # Rounding can put the start cell one step behind the entry point.
        while step_x and t_max_x < t_start:
            ix += step_x
            t_max_x += t_delta_x
        while step_y and t_max_y < t_start:
            iy += step_y
            t_max_y += t_delta_y
        while step_z and t_max_z < t_start:
            iz += step_z
            t_max_z += t_delta_z
        if not (0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz):  # pragma: no cover
            return

        nxny = nx * ny
        while True:
            cell = ix + nx * iy + nxny * iz
            if t_max_x <= t_max_y and t_max_x <= t_max_z:
                yield cell, t_max_x
                if t_max_x > t_exit:
                    return
                ix += step_x
                if ix < 0 or ix >= nx:
                    return
                t_max_x += t_delta_x
            elif t_max_y <= t_max_z:
                yield cell, t_max_y
                if t_max_y > t_exit:
                    return
                iy += step_y
                if iy < 0 or iy >= ny:
                    return
                t_max_y += t_delta_y
            else:
                yield cell, t_max_z
                if t_max_z > t_exit:
                    return
                iz += step_z
                if iz < 0 or iz >= nz:
                    return
                t_max_z += t_delta_z


def _closest_point_triangle(
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
    coords: "array",
    k: int,
    ox: float, oy: float, oz: float,
    dx: float, dy: float, dz: float,
    best: float,
) -> float | None:
    """Möller--Trumbore for triangle ``k``; ``None`` unless it is a closer hit."""
    base = 9 * k
    ax = coords[base]; ay = coords[base + 1]; az = coords[base + 2]
    e1x = coords[base + 3] - ax
    e1y = coords[base + 4] - ay
    e1z = coords[base + 5] - az
    e2x = coords[base + 6] - ax
    e2y = coords[base + 7] - ay
    e2z = coords[base + 8] - az

    px = dy * e2z - dz * e2y
    py = dz * e2x - dx * e2z
    pz = dx * e2y - dy * e2x
    det = e1x * px + e1y * py + e1z * pz
    if -_PARALLEL_EPS < det < _PARALLEL_EPS:
        return None
    inv_det = 1.0 / det

    tx = ox - ax
    ty = oy - ay
    tz = oz - az
    u = (tx * px + ty * py + tz * pz) * inv_det
    if u < -_EDGE_EPS or u > 1.0 + _EDGE_EPS:
        return None

    qx = ty * e1z - tz * e1y
    qy = tz * e1x - tx * e1z
    qz = tx * e1y - ty * e1x
    v = (dx * qx + dy * qy + dz * qz) * inv_det
    if v < -_EDGE_EPS or u + v > 1.0 + _EDGE_EPS:
        return None

    t = (e2x * qx + e2y * qy + e2z * qz) * inv_det
    if t <= _MIN_T or t >= best:
        return None
    return t


def brute_force_nearest(mesh: Mesh, origin: Vec3, direction: Vec3) -> float | None:
    """Reference O(N) nearest-hit distance, for cross-checking the grid.

    Kept in the library (not only in tests) so users can validate the
    acceleration structure on their own meshes.  Same convention as
    :meth:`TriangleGrid.cast_ray`: ``t > 0``, direction normalised internally.
    """
    dx, dy, dz = direction
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm == 0.0:
        raise ValueError("ray direction must be non-zero")
    dx /= norm
    dy /= norm
    dz /= norm

    coords = array("d")
    for face in mesh.faces:
        for index in face:
            v = mesh.vertices[index]
            coords.append(v[0]); coords.append(v[1]); coords.append(v[2])

    best = math.inf
    ox, oy, oz = origin
    for k in range(len(mesh.faces)):
        t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, best)
        if t is not None:
            best = t
    return None if best == math.inf else best


def brute_force_all(mesh: Mesh, origin: Vec3, direction: Vec3) -> tuple[float, ...]:
    """All brute-force intersections along the ray, sorted (reference)."""
    dx, dy, dz = direction
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm == 0.0:
        raise ValueError("ray direction must be non-zero")
    dx /= norm
    dy /= norm
    dz /= norm

    coords = array("d")
    for face in mesh.faces:
        for index in face:
            v = mesh.vertices[index]
            coords.append(v[0]); coords.append(v[1]); coords.append(v[2])

    ox, oy, oz = origin
    hits = []
    for k in range(len(mesh.faces)):
        t = _intersect(coords, k, ox, oy, oz, dx, dy, dz, math.inf)
        if t is not None:
            hits.append(t)
    hits.sort()
    return tuple(hits)