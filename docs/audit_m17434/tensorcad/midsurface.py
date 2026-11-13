"""Mid-surface extraction from a thin-walled solid, and planar defeaturing.

Structural idealisation starts by throwing away the wall: a stiffened panel, a hull
side shell, a bracket are all *plates*, and what a shell/plate solver wants is the
surface half-way between the two faces of the wall, plus the local thickness as an
element property. ``GAP-CAD-013`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md``
records that TensorCAD could audit a geometry but could not perform this
idealisation; this module closes it.

**The method, and why this one.** A thin wall is a pair of surfaces a distance ``t``
apart with antiparallel normals. So for every triangle take its outward normal
``n`` and its centroid ``c``, shoot the ray ``c - s n`` (``s > 0``) into the solid,
and keep the *nearest* triangle the ray strikes. If that triangle's normal is
antiparallel to ``n`` within an angle tolerance and the distance lies inside the
admissible band, the two triangles are the two faces of one wall, the measured
distance **is** the local thickness there, and the mid-surface triangle is the
originating triangle displaced by ``-n t/2``. Nothing is fitted and nothing is
smoothed: the mid-surface is exactly the set of wall mid-planes, which is why its
area can be compared with a closed form and *is* compared in
``examples/midsurface``.

The alternative -- skinning the zero level set of ``|f| - t/2`` -- needs a field
for the solid, resamples, and loses the one-to-one correspondence between a
thickness and a face. The ray cast keeps that correspondence, so the thickness is
published per pair instead of as a scalar average, and an unpaired face is a
reported fact rather than a smoothed-over hole.

**The thickness gate is not optional.** A solid cube's opposite faces *are*
antiparallel, so the ray always finds a partner; without a band the whole cube
"idealises" into a mid-plane through its centre. :func:`extract_midsurface`
therefore takes ``max_thickness``, defaulting to
:data:`DEFAULT_MAX_THICKNESS_FRACTION` of the bounding-box diagonal (5 %): thick
enough for a genuine thin wall, too thin for a solid body. Faces whose nearest
antiparallel partner is outside the band are listed in
:attr:`MidSurfaceResult.unpaired` -- the caller sees the solid body it handed in,
instead of a silently half-hollowed one.

**Walls that a cap does not simplify.** Mid-surface extraction is a local,
per-face operation, so a *capped* tube's end disk is paired with the other end disk
and offset as a whole disk, rather than replaced by its annulus. The honest
statement is therefore: extract mid-surfaces from a **skin** (two wall faces, no
caps), or from a body whose caps are thinner than the band; do not expect caps to
become annuli. ``examples/midsurface`` builds skins for exactly this reason, and
the hull case uses :func:`tensorcad.shell.offset_inward`, whose face ``i``
corresponds to face ``i`` of the outer surface.

**Defeaturing.** Idealisation also wants fewer faces: a flat panel meshed with 800
triangles is one *face* to a solver. :func:`defeature_planar_regions` grows regions
over triangles that are coplanar within an angle and a plane-offset tolerance,
takes each region's boundary loop, and re-triangulates that loop (ear clipping in
the region's own plane). The re-triangulation is accepted only if the region's area
is preserved to round-off and the loop is a single simple polygon; otherwise the
region keeps its original triangles and is counted in
:attr:`DefeatureResult.skipped_regions`. Holes are *not* stitched (a region with an
inner boundary is skipped), and curved regions are never merged -- a cylinder keeps
as many triangles as it came in with, which the example reports rather than hides.

Everything here is standard-library only, like the rest of the package.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

from . import mesh as meshlib
from .geometry import Face, Mesh, Vec3

#: Default admissible wall thickness, as a fraction of the bounding-box diagonal.
DEFAULT_MAX_THICKNESS_FRACTION = 0.05
#: Default tolerance on "the two faces of a wall are antiparallel", in degrees.
DEFAULT_ANGLE_TOLERANCE_DEG = 20.0
#: Default tolerance on "these triangles are coplanar", in degrees.
DEFAULT_COPLANAR_ANGLE_DEG = 1.0
#: Plane-offset tolerance of the defeaturing region grow, relative to the model size.
DEFAULT_PLANE_TOLERANCE_FRACTION = 1e-6


class MidSurfaceError(ValueError):
    """Raised for input that cannot describe a mid-surface problem."""


# --------------------------------------------------------------------------- #
# small vector helpers (the package stays standard-library only)
# --------------------------------------------------------------------------- #
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _mul(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec3) -> Vec3:
    length = _norm(a)
    if length == 0.0:
        raise MidSurfaceError("cannot normalise a zero-length vector")
    return _mul(a, 1.0 / length)


def _positive(value: float, what: str) -> float:
    if not (value > 0.0) or not math.isfinite(value):
        raise MidSurfaceError(f"{what} must be positive and finite: {value!r}")
    return value


def _scale_of(mesh: Mesh) -> float:
    box = mesh.bounding_box()
    diagonal = _norm(_sub(box.maximum, box.minimum))
    if not (diagonal > 0.0):
        raise MidSurfaceError("the mesh has zero extent; a mid-surface needs a size")
    return diagonal


# --------------------------------------------------------------------------- #
# ray casting against the mesh, on a uniform bucket grid
# --------------------------------------------------------------------------- #
class _FaceGrid:
    """Uniform bucket grid over triangle bounding boxes, for nearest-hit queries."""

    def __init__(self, mesh: Mesh, cell: float) -> None:
        self.mesh = mesh
        self.cell = _positive(cell, "the grid cell size")
        self.buckets: dict[tuple[int, int, int], list[int]] = {}
        self.normals = [meshlib.face_normal(mesh, face) for face in mesh.faces]
        self.centroids = [meshlib.face_centroid(mesh, face) for face in mesh.faces]
        for index, face in enumerate(mesh.faces):
            minimum, maximum = self._face_box(face)
            for key in self._keys(minimum, maximum):
                self.buckets.setdefault(key, []).append(index)

    def _face_box(self, face: Face) -> tuple[Vec3, Vec3]:
        points = [self.mesh.vertices[i] for i in face]
        minimum = tuple(min(p[axis] for p in points) for axis in range(3))
        maximum = tuple(max(p[axis] for p in points) for axis in range(3))
        return minimum, maximum  # type: ignore[return-value]

    def _keys(self, minimum: Vec3, maximum: Vec3) -> Iterable[tuple[int, int, int]]:
        lo = [int(math.floor(minimum[axis] / self.cell)) for axis in range(3)]
        hi = [int(math.floor(maximum[axis] / self.cell)) for axis in range(3)]
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                for k in range(lo[2], hi[2] + 1):
                    yield (i, j, k)

    def query(self, origin: Vec3, direction: Vec3, length: float) -> list[int]:
        """Candidate faces whose bucket the ray segment passes through."""
        candidates: set[int] = set()
        steps = max(1, int(math.ceil(length / self.cell)))
        for step in range(steps + 1):
            point = _add(origin, _mul(direction, length * step / steps))
            key = tuple(int(math.floor(point[axis] / self.cell)) for axis in range(3))
            for neighbour in self._neighbours(key):
                candidates.update(self.buckets.get(neighbour, ()))
        return sorted(candidates)

    @staticmethod
    def _neighbours(key: tuple[int, int, int]) -> Iterable[tuple[int, int, int]]:
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    yield (key[0] + di, key[1] + dj, key[2] + dk)


def _ray_triangle(origin: Vec3, direction: Vec3, a: Vec3, b: Vec3, c: Vec3) -> float | None:
    """Moller-Trumbore: distance along ``direction`` to the triangle, or ``None``."""
    edge1 = _sub(b, a)
    edge2 = _sub(c, a)
    pvec = _cross(direction, edge2)
    det = _dot(edge1, pvec)
    if abs(det) < 1e-300:
        return None
    inv_det = 1.0 / det
    tvec = _sub(origin, a)
    u = _dot(tvec, pvec) * inv_det
    if u < -1e-12 or u > 1.0 + 1e-12:
        return None
    qvec = _cross(tvec, edge1)
    v = _dot(direction, qvec) * inv_det
    if v < -1e-12 or u + v > 1.0 + 1e-12:
        return None
    return _dot(edge2, qvec) * inv_det


# --------------------------------------------------------------------------- #
# mid-surface
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ThicknessSample:
    """One wall measurement: ``face`` paired with ``partner`` at ``thickness``."""

    face: int
    partner: int
    thickness: float
    normal_dot: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "face": self.face,
            "partner": self.partner,
            "thickness": self.thickness,
            "normal_dot": self.normal_dot,
        }


@dataclass(frozen=True)
class MidSurfaceResult:
    """The extracted mid-surface plus the evidence for it."""

    mesh: Mesh
    samples: tuple[ThicknessSample, ...]
    unpaired: tuple[int, ...]
    max_thickness: float
    min_thickness: float
    angle_tolerance_deg: float
    source_face_count: int
    source_area: float

    @property
    def paired_face_count(self) -> int:
        return len(self.samples)

    @property
    def mean_thickness(self) -> float:
        if not self.samples:
            return 0.0
        return sum(sample.thickness for sample in self.samples) / len(self.samples)

    @property
    def min_measured(self) -> float:
        return min((sample.thickness for sample in self.samples), default=0.0)

    @property
    def max_measured(self) -> float:
        return max((sample.thickness for sample in self.samples), default=0.0)

    @property
    def area(self) -> float:
        return meshlib.surface_area(self.mesh)

    @property
    def paired_fraction(self) -> float:
        if self.source_face_count == 0:
            return 0.0
        return len(self.samples) / self.source_face_count

    @property
    def area_change_fraction(self) -> float:
        if not self.source_area:
            return 0.0
        return (self.area - self.source_area) / self.source_area

    def to_dict(self) -> dict[str, object]:
        return {
            "source_face_count": self.source_face_count,
            "source_area": self.source_area,
            "mid_face_count": len(self.mesh.faces),
            "mid_vertex_count": len(self.mesh.vertices),
            "mid_area": self.area,
            "paired_face_count": self.paired_face_count,
            "unpaired_face_count": len(self.unpaired),
            "paired_fraction": self.paired_fraction,
            "mean_thickness": self.mean_thickness,
            "min_measured_thickness": self.min_measured,
            "max_measured_thickness": self.max_measured,
            "max_thickness_band": self.max_thickness,
            "min_thickness_band": self.min_thickness,
            "angle_tolerance_deg": self.angle_tolerance_deg,
            "area_change_fraction": self.area_change_fraction,
        }


def _emit_mid_triangles(
    mesh: Mesh, grid: _FaceGrid, samples: Sequence[ThicknessSample]
) -> tuple[list[Vec3], list[Face]]:
    """One mid triangle per measured wall, copied from a face chosen by the geometry.

    A mid triangle is a *translated copy* of one of the two facings of a wall, and the two
    facings of a curved wall do not have the same area, so the side the copy comes from
    decides the mid-surface area by up to the wall's own ``t/R``. Taking "whichever face
    the file lists first" makes that answer a property of the labelling: relabelling the
    faces of the same geometry moved the mid area of a t/R = 2 % tube by 1.2 % and of a
    sphere by 2.1 % (pinned in ``tests/test_midsurface_oracle.py``). The side is therefore
    chosen per wall, from the faces alone:

    * the facing whose outward normal points away from the model's bounding-box centre
      owns the wall (its normal against the centroid-minus-centre vector) -- on a skin
      that is the outer facing, which is also the side the ``t/(2R)`` surplus is stated
      for;
    * a tie -- a flat plate's two facings, which are congruent up to the translation -- is
      broken by the normal itself, so one winding covers the whole sheet instead of a
      checkerboard of them;
    * each face emits at most one mid triangle, so a wall that both facings measured is
      not turned into two overlapping triangles, and the walls are visited in an order
      fixed by the chosen face's centroid so that the greedy match does not depend on
      face numbering either.
    """
    box = mesh.bounding_box()
    middle = _mul(_add(box.minimum, box.maximum), 0.5)
    measured = {sample.face: sample for sample in samples}
    outward: dict[int, float] = {}
    normals: dict[int, Vec3] = {}
    for face in measured:
        raw = grid.normals[face]
        if _norm(raw) == 0.0:
            outward[face] = -math.inf
            normals[face] = (0.0, 0.0, 0.0)
            continue
        normal = _unit(raw)
        normals[face] = normal
        outward[face] = _dot(normal, _sub(grid.centroids[face], middle))

    # A wall is *reciprocal* when each facing measured the other. One-sided measurements
    # (face A's ray finds B, but B's ray finds something else) are ambiguous about which
    # facing owns the mid triangle, and on a hull plating they are exactly the faces at
    # the folded waterline rim, where the two facings disagree about the wall's size. They
    # stay in ``samples`` -- the thickness was really measured -- but they do not emit.
    partner_of = {sample.face: sample.partner for sample in samples}
    walls: dict[tuple[int, int], None] = {}
    for sample in samples:
        if partner_of.get(sample.partner) != sample.face:
            continue
        walls[(min(sample.face, sample.partner), max(sample.face, sample.partner))] = None

    def _side_key(face: int) -> tuple[float, float, float, float]:
        # prefer the outward facing, then one fixed normal direction, then the centroid
        normal = normals[face]
        return (-outward[face], -normal[0], -normal[1], -normal[2])

    def _wall_order(item: tuple[int, int]) -> tuple[float, float, float, float]:
        owner = min((face for face in item if face in measured), key=_side_key)
        centroid = grid.centroids[owner]
        return (round(-outward[owner], 9), centroid[0], centroid[1], centroid[2])

    mid_vertices: list[Vec3] = []
    mid_faces: list[Face] = []
    used: set[int] = set()
    for wall in sorted(walls, key=_wall_order):
        # Only the owning facing emits. Falling back to the other facing when the owner
        # is taken would make the answer depend on which of two walls that share a face
        # happens to be visited first -- and two walls with the same owner would emit
        # *the same* mid triangle anyway, so skipping the second loses nothing.
        owner = min((face for face in wall if face in measured), key=_side_key)
        if owner in used:
            continue
        used.add(owner)
        offset = _mul(normals[owner], -0.5 * measured[owner].thickness)
        base = len(mid_vertices)
        mid_vertices.extend(_add(mesh.vertices[i], offset) for i in mesh.faces[owner])
        mid_faces.append((base, base + 1, base + 2))
    return mid_vertices, mid_faces


def extract_midsurface(
    mesh: Mesh,
    *,
    max_thickness: float | None = None,
    min_thickness: float = 0.0,
    angle_tolerance_deg: float = DEFAULT_ANGLE_TOLERANCE_DEG,
    weld_tolerance: float | None = None,
) -> MidSurfaceResult:
    """Extract the mid-surface of the thin walls in ``mesh``.

    Every face is ray-cast along its inward normal; the nearest antiparallel partner
    within ``[min_thickness, max_thickness]`` makes it a wall face, and the pair
    contributes one mid-surface triangle displaced by ``-n t/2``. Faces with no such
    partner are reported in :attr:`MidSurfaceResult.unpaired` and contribute nothing.

    ``max_thickness`` defaults to :data:`DEFAULT_MAX_THICKNESS_FRACTION` of the
    bounding-box diagonal. Raising it lets thicker walls through but also starts
    idealising solid bodies, whose opposite faces are antiparallel too.
    """
    if not mesh.faces:
        raise MidSurfaceError("the mesh has no faces; there is nothing to idealise")
    scale = _scale_of(mesh)
    if max_thickness is None:
        max_thickness = DEFAULT_MAX_THICKNESS_FRACTION * scale
    _positive(max_thickness, "max_thickness")
    if min_thickness < 0.0 or not math.isfinite(min_thickness):
        raise MidSurfaceError(f"min_thickness must be finite and non-negative: {min_thickness!r}")
    if min_thickness >= max_thickness:
        raise MidSurfaceError(
            f"the thickness band is empty: min_thickness {min_thickness!r} >= "
            f"max_thickness {max_thickness!r}"
        )
    if not (0.0 <= angle_tolerance_deg < 90.0):
        raise MidSurfaceError(
            f"angle_tolerance_deg must be in [0, 90): {angle_tolerance_deg!r}"
        )
    cos_limit = math.cos(math.radians(angle_tolerance_deg))

    grid = _FaceGrid(mesh, cell=max(max_thickness, 1e-6 * scale))
    samples: list[ThicknessSample] = []
    unpaired: list[int] = []

    for index, face in enumerate(mesh.faces):
        normal = grid.normals[index]
        if _norm(normal) == 0.0:
            unpaired.append(index)
            continue
        normal = _unit(normal)
        inward = _mul(normal, -1.0)
        centroid = grid.centroids[index]
        best: tuple[tuple[float, float, float, float], float, int] | None = None
        for candidate in grid.query(centroid, inward, max_thickness):
            if candidate == index:
                continue
            hit = _ray_triangle(
                centroid,
                inward,
                *(mesh.vertices[i] for i in mesh.faces[candidate]),
            )
            if hit is None or hit <= 1e-12 * scale:
                continue
            # Equidistant candidates are the rule rather than the exception on a
            # symmetric mesh (the facets either side of a sphere's pole ring are
            # equidistant to round-off), so "the first one the grid listed" is not a
            # decision: that list is sorted by face index. Rank by distance at the weld
            # scale first, then by the candidate's own centroid, which is geometry.
            rank = (round(hit / scale, 12), *grid.centroids[candidate])
            if best is None or rank < best[0]:
                best = (rank, hit, candidate)
        if best is None:
            unpaired.append(index)
            continue
        _, thickness, partner = best
        opponent = grid.normals[partner]
        if _norm(opponent) == 0.0:
            unpaired.append(index)
            continue
        normal_dot = _dot(normal, _unit(opponent))
        if normal_dot > -cos_limit or not (min_thickness <= thickness <= max_thickness):
            unpaired.append(index)
            continue
        samples.append(
            ThicknessSample(
                face=index, partner=partner, thickness=thickness, normal_dot=normal_dot
            )
        )

    mid_vertices, mid_faces = _emit_mid_triangles(mesh, grid, samples)

    tolerance = weld_tolerance if weld_tolerance is not None else 1e-9 * scale
    skin = Mesh(
        vertices=tuple(mid_vertices),
        faces=tuple(mid_faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    welded = meshlib.weld(skin, tolerance=tolerance) if skin.faces else skin
    return MidSurfaceResult(
        mesh=welded,
        samples=tuple(samples),
        unpaired=tuple(unpaired),
        max_thickness=max_thickness,
        min_thickness=min_thickness,
        angle_tolerance_deg=angle_tolerance_deg,
        source_face_count=len(mesh.faces),
        source_area=meshlib.surface_area(mesh),
    )


# --------------------------------------------------------------------------- #
# planar defeaturing: merge coplanar triangles into one polygon each
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CoplanarRegion:
    """A maximal set of adjacent, coplanar triangles, and its simplified form."""

    faces: tuple[int, ...]
    normal: Vec3
    offset: float
    area: float
    boundary_loops: tuple[tuple[int, ...], ...]
    simplified: bool
    reason: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "face_count": len(self.faces),
            "normal": list(self.normal),
            "offset": self.offset,
            "area": self.area,
            "boundary_loop_count": len(self.boundary_loops),
            "boundary_vertex_counts": [len(loop) for loop in self.boundary_loops],
            "simplified": self.simplified,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class DefeatureResult:
    """The defeatured mesh plus the region ledger behind it."""

    mesh: Mesh
    regions: tuple[CoplanarRegion, ...]
    angle_tolerance_deg: float
    plane_tolerance: float

    @property
    def face_count_before(self) -> int:
        return sum(len(region.faces) for region in self.regions)

    @property
    def face_count_after(self) -> int:
        return len(self.mesh.faces)

    @property
    def face_reduction(self) -> float:
        before = self.face_count_before
        return 0.0 if before == 0 else 1.0 - self.face_count_after / before

    @property
    def skipped_regions(self) -> tuple[CoplanarRegion, ...]:
        return tuple(region for region in self.regions if not region.simplified)

    @property
    def area_before(self) -> float:
        return sum(region.area for region in self.regions)

    @property
    def area_after(self) -> float:
        return meshlib.surface_area(self.mesh)

    @property
    def area_relative_error(self) -> float:
        if not self.area_before:
            return 0.0
        return abs(self.area_after - self.area_before) / self.area_before

    def to_dict(self) -> dict[str, object]:
        return {
            "region_count": len(self.regions),
            "simplified_region_count": len(self.regions) - len(self.skipped_regions),
            "skipped_region_count": len(self.skipped_regions),
            "face_count_before": self.face_count_before,
            "face_count_after": self.face_count_after,
            "face_reduction": self.face_reduction,
            "area_before": self.area_before,
            "area_after": self.area_after,
            "area_relative_error": self.area_relative_error,
            "angle_tolerance_deg": self.angle_tolerance_deg,
            "plane_tolerance": self.plane_tolerance,
            "regions": [region.to_dict() for region in self.regions],
        }


def _edge_map(faces: Sequence[Face]) -> dict[tuple[int, int], list[tuple[int, tuple[int, int]]]]:
    """Undirected edge -> list of (local face index, directed edge)."""
    table: dict[tuple[int, int], list[tuple[int, tuple[int, int]]]] = {}
    for local, face in enumerate(faces):
        for a, b in meshlib.face_edges(face):
            key = meshlib.undirected_edge((a, b))
            table.setdefault(key, []).append((local, (a, b)))
    return table


def _loops_from_boundary(boundary: Sequence[tuple[int, int]]) -> tuple[tuple[int, ...], ...]:
    """Chain directed boundary edges of a region into closed loops."""
    outgoing: dict[int, list[int]] = {}
    for a, b in boundary:
        outgoing.setdefault(a, []).append(b)
    loops: list[tuple[int, ...]] = []
    while outgoing:
        start = next(iter(outgoing))
        loop = [start]
        current = start
        closed = False
        while True:
            targets = outgoing.get(current)
            if not targets:
                break
            nxt = targets.pop()
            if not targets:
                del outgoing[current]
            if nxt == start:
                closed = True
                break
            if nxt in loop:  # pragma: no cover - defensive: non-simple boundary
                break
            loop.append(nxt)
            current = nxt
        if closed and len(loop) >= 3:
            loops.append(tuple(loop))
    return tuple(loops)


def _plane_basis(normal: Vec3) -> tuple[Vec3, Vec3]:
    helper = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _unit(_cross(normal, helper))
    v = _cross(normal, u)
    return u, v


def _segments_cross(p1, p2, p3, p4) -> bool:  # type: ignore[no-untyped-def]
    def orient(ax, ay, bx, by, cx, cy):  # type: ignore[no-untyped-def]
        value = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        if abs(value) < 1e-15:
            return 0
        return 1 if value > 0 else -1

    o1 = orient(*p1, *p2, *p3)
    o2 = orient(*p1, *p2, *p4)
    o3 = orient(*p3, *p4, *p1)
    o4 = orient(*p3, *p4, *p2)
    return o1 * o2 < 0 and o3 * o4 < 0


def _polygon_is_simple(points: Sequence[tuple[float, float]]) -> bool:
    """O(n^2) simplicity test in the region plane (no crossing edges)."""
    count = len(points)
    for i in range(count):
        a = points[i]
        b = points[(i + 1) % count]
        for j in range(i + 1, count):
            if j == i or (j + 1) % count == i or j == (i + 1) % count:
                continue
            if _segments_cross(a, b, points[j], points[(j + 1) % count]):
                return False
    return True


def _point_in_triangle(point, a, b, c) -> bool:  # type: ignore[no-untyped-def]
    def sign(p1, p2, p3):  # type: ignore[no-untyped-def]
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1 = sign(point, a, b)
    d2 = sign(point, b, c)
    d3 = sign(point, c, a)
    has_negative = d1 < 0 or d2 < 0 or d3 < 0
    has_positive = d1 > 0 or d2 > 0 or d3 > 0
    return not (has_negative and has_positive)


def _drop_collinear(points: Sequence[tuple[float, float]], tolerance: float) -> list[int]:
    """Indices of a polygon loop with the in-line vertices removed.

    A coplanar region's boundary carries a vertex wherever the *triangulation* had one,
    so a plate meshed 20x10 bounds with 60 vertices that describe four straight edges.
    Dropping them is what turns those 60 into 4 -- and the caller re-checks the area, so
    a wrong drop cannot pass silently.
    """
    keep = list(range(len(points)))
    changed = True
    while changed and len(keep) > 3:
        changed = False
        for position in range(len(keep)):
            previous = points[keep[position - 1]]
            current = points[keep[position]]
            following = points[keep[(position + 1) % len(keep)]]
            ax, ay = current[0] - previous[0], current[1] - previous[1]
            bx, by = following[0] - current[0], following[1] - current[1]
            scale = math.hypot(ax, ay) * math.hypot(bx, by)
            if scale == 0.0 or abs(ax * by - ay * bx) <= tolerance * scale:
                keep.pop(position)
                changed = True
                break
    return keep


def _ear_clip(points: Sequence[tuple[float, float]]) -> list[tuple[int, int, int]]:
    """Triangulate a simple polygon given in the region plane; returns local indices."""
    count = len(points)
    if count < 3:
        return []
    signed = 0.0
    for i in range(count):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % count]
        signed += x1 * y2 - x2 * y1
    order = list(range(count))
    if signed < 0.0:
        order.reverse()
    triangles: list[tuple[int, int, int]] = []
    guard = 0
    while len(order) > 3 and guard < 4 * count * count:
        guard += 1
        clipped = False
        for position in range(len(order)):
            i0 = order[position - 1]
            i1 = order[position]
            i2 = order[(position + 1) % len(order)]
            a, b, c = points[i0], points[i1], points[i2]
            cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if cross <= 0.0:
                continue
            blocked = False
            for other in order:
                if other in (i0, i1, i2):
                    continue
                if _point_in_triangle(points[other], a, b, c):
                    blocked = True
                    break
            if blocked:
                continue
            triangles.append((i0, i1, i2))
            order.pop(position)
            clipped = True
            break
        if not clipped:  # pragma: no cover - defensive: degenerate polygon
            return []
    if len(order) == 3:
        triangles.append((order[0], order[1], order[2]))
    return triangles


def _triangle_area_3d(a: Vec3, b: Vec3, c: Vec3) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def defeature_planar_regions(
    mesh: Mesh,
    *,
    angle_tolerance_deg: float = DEFAULT_COPLANAR_ANGLE_DEG,
    plane_tolerance: float | None = None,
) -> DefeatureResult:
    """Merge coplanar adjacent triangles into single polygons and re-triangulate them.

    Two adjacent triangles join when their normals differ by at most
    ``angle_tolerance_deg`` and their planes lie within ``plane_tolerance`` of each
    other (default: :data:`DEFAULT_PLANE_TOLERANCE_FRACTION` of the bounding-box
    diagonal). Each resulting region is replaced by an ear-clipping triangulation of
    its outer boundary loop -- accepted only when the region has exactly one
    boundary loop, the loop is simple, and the area is preserved to round-off.
    Curved regions are untouched (no two facets are coplanar), so a curved shell
    keeps its face count; that is reported, not hidden.
    """
    if not mesh.faces:
        raise MidSurfaceError("the mesh has no faces; there is nothing to defeature")
    if not (0.0 <= angle_tolerance_deg < 90.0):
        raise MidSurfaceError(
            f"angle_tolerance_deg must be in [0, 90): {angle_tolerance_deg!r}"
        )
    scale = _scale_of(mesh)
    if plane_tolerance is None:
        plane_tolerance = DEFAULT_PLANE_TOLERANCE_FRACTION * scale
    cos_limit = math.cos(math.radians(angle_tolerance_deg))
    unit_normals: list[Vec3] = []
    for face in mesh.faces:
        normal = meshlib.face_normal(mesh, face)
        length = _norm(normal)
        unit_normals.append(_mul(normal, 1.0 / length) if length else (0.0, 0.0, 0.0))
    offsets = [
        _dot(unit_normals[index], meshlib.face_centroid(mesh, face))
        for index, face in enumerate(mesh.faces)
    ]
    areas = [meshlib.triangle_area(mesh, face) for face in mesh.faces]

    adjacency: dict[int, set[int]] = {index: set() for index in range(len(mesh.faces))}
    for entries in _edge_map(mesh.faces).values():
        if len(entries) != 2:
            continue
        first, second = entries[0][0], entries[1][0]
        adjacency[first].add(second)
        adjacency[second].add(first)

    seen: set[int] = set()
    regions: list[CoplanarRegion] = []
    out_vertices: list[Vec3] = []
    out_faces: list[Face] = []

    for seed in range(len(mesh.faces)):
        if seed in seen:
            continue
        normal = unit_normals[seed]
        offset = offsets[seed]
        stack = [seed]
        seen.add(seed)
        group: list[int] = []
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbour in adjacency[current]:
                if neighbour in seen:
                    continue
                if _dot(normal, unit_normals[neighbour]) < cos_limit:
                    continue
                if abs(offsets[neighbour] - offset) > plane_tolerance:
                    continue
                seen.add(neighbour)
                stack.append(neighbour)
        group.sort()
        region_area = sum(areas[index] for index in group)
        local_faces = [mesh.faces[index] for index in group]
        boundary: list[tuple[int, int]] = []
        for entries in _edge_map(local_faces).values():
            if len(entries) == 1:
                boundary.append(entries[0][1])
        loops = _loops_from_boundary(boundary)

        simplified = False
        reason = ""
        if len(loops) != 1:
            reason = f"region has {len(loops)} boundary loops; holes are not stitched"
        elif len(loops[0]) < 3:
            reason = "boundary loop is degenerate"
        else:
            loop = loops[0]
            origin = mesh.vertices[loop[0]]
            u_axis, v_axis = _plane_basis(normal)
            planar = [
                (
                    _dot(_sub(mesh.vertices[index], origin), u_axis),
                    _dot(_sub(mesh.vertices[index], origin), v_axis),
                )
                for index in loop
            ]
            kept = _drop_collinear(planar, 1e-12)
            if len(kept) < 3:
                reason = "boundary loop collapses to a line after dropping collinear vertices"
            elif not _polygon_is_simple([planar[index] for index in kept]):
                reason = "boundary loop is not simple in the region plane"
            else:
                loop = [loop[index] for index in kept]
                planar = [planar[index] for index in kept]
                local = _ear_clip(planar)
                if not local:
                    reason = "ear clipping produced no triangles"
                else:
                    candidate_area = sum(
                        _triangle_area_3d(
                            mesh.vertices[loop[i0]],
                            mesh.vertices[loop[i1]],
                            mesh.vertices[loop[i2]],
                        )
                        for i0, i1, i2 in local
                    )
                    tolerance = 1e-9 * region_area if region_area > 0.0 else 0.0
                    if region_area <= 0.0 or abs(candidate_area - region_area) > tolerance:
                        reason = (
                            "re-triangulation does not preserve the region area "
                            f"({candidate_area!r} vs {region_area!r})"
                        )
                    else:
                        base = len(out_vertices)
                        out_vertices.extend(mesh.vertices[index] for index in loop)
                        out_faces.extend(
                            (base + i0, base + i1, base + i2) for i0, i1, i2 in local
                        )
                        simplified = True
        if not simplified:
            base = len(out_vertices)
            for face in local_faces:
                out_vertices.extend(mesh.vertices[index] for index in face)
                out_faces.append((base, base + 1, base + 2))
                base += 3
        regions.append(
            CoplanarRegion(
                faces=tuple(group),
                normal=normal,
                offset=offset,
                area=region_area,
                boundary_loops=loops,
                simplified=simplified,
                reason=reason,
            )
        )

    merged = Mesh(
        vertices=tuple(out_vertices),
        faces=tuple(out_faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    merged = meshlib.weld(merged, tolerance=1e-9 * scale) if merged.faces else merged
    return DefeatureResult(
        mesh=merged,
        regions=tuple(regions),
        angle_tolerance_deg=angle_tolerance_deg,
        plane_tolerance=plane_tolerance,
    )


# --------------------------------------------------------------------------- #
# skin builders used by the verification example (and by callers who need one)
# --------------------------------------------------------------------------- #
def plate_solid(
    length: float, breadth: float, thickness: float, *, nx: int = 20, ny: int = 10
) -> Mesh:
    """A closed box ``length x breadth x thickness``, large faces meshed ``nx x ny``.

    The two large faces are subdivided so mid-surface extraction has something to
    simplify; each of the four rim strips is one quad, which is the correct
    resolution of a plate rim at this scale.
    """
    for value, name in ((length, "length"), (breadth, "breadth"), (thickness, "thickness")):
        _positive(value, name)
    if nx < 1 or ny < 1:
        raise MidSurfaceError("nx and ny must be at least 1")
    half_t = 0.5 * thickness
    vertices: list[Vec3] = []
    faces: list[Face] = []
    lookup: dict[tuple[float, float, float], int] = {}

    def vertex(x: float, y: float, z: float) -> int:
        key = (x, y, z)
        if key not in lookup:
            lookup[key] = len(vertices)
            vertices.append(key)
        return lookup[key]

    def quad(a: int, b: int, c: int, d: int) -> None:
        faces.append((a, b, c))
        faces.append((a, c, d))

    for top in (True, False):
        z = half_t if top else -half_t
        for i in range(nx):
            for j in range(ny):
                x0, x1 = length * i / nx, length * (i + 1) / nx
                y0, y1 = breadth * j / ny, breadth * (j + 1) / ny
                a = vertex(x0, y0, z)
                b = vertex(x1, y0, z)
                c = vertex(x1, y1, z)
                d = vertex(x0, y1, z)
                if top:
                    quad(a, b, c, d)
                else:
                    quad(a, d, c, b)
    rim = [
        ((0.0, 0.0), (length, 0.0)),
        ((length, 0.0), (length, breadth)),
        ((length, breadth), (0.0, breadth)),
        ((0.0, breadth), (0.0, 0.0)),
    ]
    for (x0, y0), (x1, y1) in rim:
        a = vertex(x0, y0, half_t)
        b = vertex(x1, y1, half_t)
        c = vertex(x1, y1, -half_t)
        d = vertex(x0, y0, -half_t)
        quad(a, b, c, d)
    return Mesh(vertices=tuple(vertices), faces=tuple(faces), units="m")


def cylinder_skin(
    radius: float,
    thickness: float,
    length: float,
    *,
    segments: int = 128,
    z0: float = 0.0,
) -> Mesh:
    """Two concentric **open** cylinders: the skin of a tube of wall ``thickness``.

    Open at both ends on purpose -- see the module docstring: a capped tube's end
    disk pairs with the other end disk instead of becoming an annulus.
    """
    for value, name in ((radius, "radius"), (thickness, "thickness"), (length, "length")):
        _positive(value, name)
    if not (0.0 < thickness < radius):
        raise MidSurfaceError("the wall must be thinner than the cylinder's radius")
    if segments < 8:
        raise MidSurfaceError("segments must be at least 8")
    inner = radius - thickness
    vertices: list[Vec3] = []
    faces: list[Face] = []
    for surface, inverted in ((radius, False), (inner, True)):
        base = len(vertices)
        for k in range(segments):
            angle = 2.0 * math.pi * k / segments
            x = surface * math.cos(angle)
            y = surface * math.sin(angle)
            vertices.append((x, y, z0))
            vertices.append((x, y, z0 + length))
        for k in range(segments):
            k2 = (k + 1) % segments
            a = base + 2 * k
            b = base + 2 * k2
            c = base + 2 * k2 + 1
            d = base + 2 * k + 1
            if not inverted:
                faces.append((a, b, c))
                faces.append((a, c, d))
            else:
                faces.append((a, c, b))
                faces.append((a, d, c))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces), units="m")


def sphere_skin(radius: float, thickness: float, *, rings: int = 64, segments: int = 128) -> Mesh:
    """Two concentric sphere meshes: the skin of a spherical shell."""
    for value, name in ((radius, "radius"), (thickness, "thickness")):
        _positive(value, name)
    if not (0.0 < thickness < radius):
        raise MidSurfaceError("the wall must be thinner than the sphere's radius")
    if rings < 4 or segments < 8:
        raise MidSurfaceError("rings >= 4 and segments >= 8 are required")
    inner = radius - thickness
    vertices: list[Vec3] = []
    faces: list[Face] = []
    for surface, inverted in ((radius, False), (inner, True)):
        base = len(vertices)
        for i in range(rings + 1):
            polar = math.pi * i / rings
            for k in range(segments):
                azimuth = 2.0 * math.pi * k / segments
                vertices.append(
                    (
                        surface * math.sin(polar) * math.cos(azimuth),
                        surface * math.sin(polar) * math.sin(azimuth),
                        surface * math.cos(polar),
                    )
                )
        for i in range(rings):
            for k in range(segments):
                k2 = (k + 1) % segments
                a = base + i * segments + k
                b = base + i * segments + k2
                c = base + (i + 1) * segments + k2
                d = base + (i + 1) * segments + k
                if not inverted:
                    # polar increases downward (+z to -z): this order is outward
                    faces.append((a, c, b))
                    faces.append((a, d, c))
                else:
                    faces.append((a, b, c))
                    faces.append((a, c, d))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces), units="m")


def mid_surface_between(outer: Mesh, inner: Mesh, *, weld_tolerance: float | None = None) -> Mesh:
    """Mid-surface of a wall given as two index-aligned surfaces.

    ``tensorcad.shell.offset_inward`` preserves face indices, so when both surfaces
    are available there is no ray casting to do: the mid-surface vertex is simply the
    average of the two corresponding vertices. This is the independent route the
    example uses to cross-check :func:`extract_midsurface` on the hull.
    """
    if len(outer.vertices) != len(inner.vertices) or len(outer.faces) != len(inner.faces):
        raise MidSurfaceError(
            "the two surfaces must share topology: "
            f"outer {len(outer.vertices)}v/{len(outer.faces)}f vs "
            f"inner {len(inner.vertices)}v/{len(inner.faces)}f"
        )
    vertices = tuple(
        (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]))
        for a, b in zip(outer.vertices, inner.vertices)
    )
    merged = Mesh(
        vertices=vertices, faces=outer.faces, units=outer.units,
        coordinate_frame=outer.coordinate_frame,
    )
    scale = _scale_of(merged)
    tolerance = weld_tolerance if weld_tolerance is not None else 1e-9 * scale
    return meshlib.weld(merged, tolerance=tolerance)