"""Multi-patch assembly: stitch a network of NURBS patches into one conforming mesh.

:mod:`tensorcad.edge_alignment` makes the *representation* of two patch edges
compatible (equal degree, equal knot vector, therefore identical control points when
the two edges describe the same curve) and names the step it leaves out:

    "Curves only (edges), not whole patches: aligning a patch network means aligning
    the edge curves it shares and using them as the patch data, which is exactly what
    a stitching operation needs; multi-patch face assembly itself stays out of scope,
    as it does in :mod:`tensorcad.nurbs`."

This module is that step. A :class:`PatchNetwork` is a list of patches plus an
explicit list of :class:`Seam` declarations -- which edge of which patch is
identified with which edge of which patch, and in which direction. :func:`stitch`
samples the network on a common grid and **welds by index, not by tolerance**: the
node arrays of the two sides of a seam are unified, so a seam that the patch data
already closes becomes a seam that the *topology* closes. The residual geometric gap
is measured and reported rather than assumed away, which is the difference between
"watertight because the faces share an index" and "watertight because the numbers
happened to be inside 1e-9".

Why index welding and not position welding
------------------------------------------
Welding on position (which is what :meth:`tensorcad.nurbs.NURBSSurface.tessellate`
does, keyed on the exact coordinate tuple) requires the two patches to evaluate their
common edge to the *same bits*, and independent patches generally do not: a
meridian-parametrised edge and an arc-parametrised edge hit the same point set at
different parameters, and even identical rational data evaluated through two
different tensor-product denominators can differ in the last decimal. Index welding
is what a real CAD kernel does, and it keeps the two concerns separate: the
**declaration** says which nodes are the same node, the **report** says how far apart
the two patches actually put them.

Named limitations
-----------------
* **Untrimmed patches only.** Every side of a patch is a full boundary curve; a
  trimmed patch (a bulbous bow cut out of a hull patch, a frame opening) has boundary
  curves that are not patch edges at all, and is out of scope here. This is the same
  reason :mod:`tensorcad.nurbs` has no trimmed surfaces.
* **A common (``nu``, ``nv``) grid.** Seams are welded node-for-node, so a seam that
  joins a ``u``-run side to a ``v``-run side needs ``nu == nv``; the check is explicit
  and raises rather than silently sampling one side more finely than the other.
* **C0 seams and closed shells, not smooth spline surfaces.** The result is a mesh,
  not a single G1 NURBS surface: nothing here changes the patch data to raise the
  smoothness across a seam. :meth:`Stitched.seam_report` measures the dihedral angle
  across every seam, so a seam that is only G0 (a box edge, 90 degrees) is reported
  as such instead of being mistaken for a defect.
* **Caps are not synthesised.** A network is closed only when every side is consumed
  by a seam or is degenerate (a pole). A hole is reported as ``open_sides``; closing
  a planar hole with a fan is :mod:`tensorcad.csg` / caller territory, a curved hole
  would need a trimmed patch.
* **Orientation is inherited.** Winding follows ``S_u x S_v`` of each patch, so the
  sign of the enclosed volume is the construction's, not the stitcher's. Use
  :func:`tensorcad.mesh.orient_outward` on the result if the patches disagree.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder, Vec3
from .nurbs import NURBSCurve, NURBSError, NURBSSurface, circular_arc

SIDES: tuple[str, ...] = ("u0", "u1", "v0", "v1")
"""The four boundary sides of a patch: the edge at ``u_min``, ``u_max``, ``v_min``,
``v_max``.  ``"v0"`` and ``"v1"`` are curves in ``u``; ``"u0"`` and ``"u1"`` are
curves in ``v``."""

_URUN = ("v0", "v1")


class PatchNetworkError(ValueError):
    """The declared network cannot be stitched as written."""


def _side_samples(side: str, nu: int, nv: int) -> int:
    """Number of grid nodes on a side (``nu`` for a ``u``-run, ``nv`` for a ``v``-run)."""
    return nu if side in _URUN else nv


def _check_side(side: str) -> str:
    if side not in SIDES:
        raise PatchNetworkError(f"side must be one of {SIDES}, got {side!r}")
    return side


@dataclass(frozen=True)
class Seam:
    """A declared identification of two patch edges.

    ``first``/``second`` are ``(patch_index, side)`` pairs.  ``reversed=True`` states
    that the second edge is traced in the opposite direction to the first (a
    neighbour whose parameter runs the other way); the two are then welded node
    ``k`` to node ``N-1-k``.  The direction is *declared*, never guessed -- the same
    rule :func:`tensorcad.edge_alignment.align_edges` follows.
    """

    first: tuple[int, str]
    second: tuple[int, str]
    reversed: bool = False
    label: str = ""

    def __post_init__(self) -> None:
        first = (int(self.first[0]), _check_side(self.first[1]))
        second = (int(self.second[0]), _check_side(self.second[1]))
        object.__setattr__(self, "first", first)
        object.__setattr__(self, "second", second)
        if first[0] < 0 or second[0] < 0:
            raise PatchNetworkError(f"patch indices must be non-negative: {self!r}")


@dataclass(frozen=True)
class PatchNetwork:
    """Patches plus the seams that identify their edges.

    Validation happens at construction, so a network either names every edge at most
    once or does not exist: a side claimed by two seams would need one node to be
    welded to two different ones, which is a contradiction rather than a tolerance
    question.
    """

    surfaces: tuple[NURBSSurface, ...]
    seams: tuple[Seam, ...] = ()

    def __post_init__(self) -> None:
        surfaces = tuple(self.surfaces)
        if not surfaces:
            raise PatchNetworkError("a patch network needs at least one patch")
        for index, surface in enumerate(surfaces):
            if not isinstance(surface, NURBSSurface):
                raise PatchNetworkError(
                    f"patch {index} is not a NURBSSurface: {type(surface).__name__}"
                )
        seams = tuple(self.seams)
        claimed: dict[tuple[int, str], int] = {}
        for position, seam in enumerate(seams):
            if not isinstance(seam, Seam):
                raise PatchNetworkError(f"seam {position} is not a Seam: {seam!r}")
            for side in (seam.first, seam.second):
                if side[0] >= len(surfaces):
                    raise PatchNetworkError(
                        f"seam {position} names patch {side[0]}, but the network has "
                        f"{len(surfaces)} patches"
                    )
                if side in claimed:
                    raise PatchNetworkError(
                        f"side {side} is claimed by seams {claimed[side]} and {position}; "
                        "a side can be welded to exactly one neighbour"
                    )
                claimed[side] = position
            if seam.first == seam.second:
                raise PatchNetworkError(f"seam {position} joins a side to itself")
        object.__setattr__(self, "surfaces", surfaces)
        object.__setattr__(self, "seams", seams)

    @property
    def claimed_sides(self) -> tuple[tuple[int, str], ...]:
        return tuple(side for seam in self.seams for side in (seam.first, seam.second))

    def unclaimed_sides(self) -> tuple[tuple[int, str], ...]:
        claimed = set(self.claimed_sides)
        return tuple(
            (index, side)
            for index in range(len(self.surfaces))
            for side in SIDES
            if (index, side) not in claimed
        )


@dataclass(frozen=True)
class SeamRecord:
    """The measured state of one seam.

    ``gap_max`` is the largest distance between the two patches' node positions at
    paired parameters, measured *before* welding (after welding the two are the same
    index by construction, so the number would be zero by definition and say nothing).
    ``angle_max_deg`` is the largest angle between the two patches' unit normals at
    those paired collar points: ``0`` for a seam that is tangent continuous (G1), the
    dihedral angle for a seam that is only positionally continuous (G0) -- a box edge
    is 90 degrees and correct.  ``inset_normals`` counts the nodes whose normal had to
    be taken at a neighbour because the exact point is singular (see
    :func:`_seam_normal`); a non-zero count is disclosed, not hidden.
    """

    label: str
    first: tuple[int, str]
    second: tuple[int, str]
    reversed: bool
    samples: int
    gap_max: float
    gap_mean: float
    angle_max_deg: float
    angle_mean_deg: float
    inset_normals: int = 0

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "first": list(self.first),
            "second": list(self.second),
            "reversed": self.reversed,
            "samples": self.samples,
            "gap_max_m": self.gap_max,
            "gap_mean_m": self.gap_mean,
            "angle_max_deg": self.angle_max_deg,
            "angle_mean_deg": self.angle_mean_deg,
            "inset_normals": self.inset_normals,
        }


@dataclass(frozen=True)
class Stitched:
    """The result of :func:`stitch`: one mesh plus everything measured about it."""

    mesh: Mesh
    seams: tuple[SeamRecord, ...]
    open_sides: tuple[tuple[int, str], ...]
    pole_sides: tuple[tuple[int, str], ...]
    nu: int
    nv: int
    welded_nodes: int
    patch_count: int

    def seam_report(self) -> dict[str, object]:
        gaps = [record.gap_max for record in self.seams]
        angles = [record.angle_max_deg for record in self.seams]
        return {
            "patches": self.patch_count,
            "samples": [self.nu, self.nv],
            "seams": [record.to_dict() for record in self.seams],
            "seam_count": len(self.seams),
            "gap_max_m": max(gaps) if gaps else 0.0,
            "angle_max_deg": max(angles) if angles else 0.0,
            "open_sides": [list(side) for side in self.open_sides],
            "pole_sides": [list(side) for side in self.pole_sides],
            "welded_nodes": self.welded_nodes,
        }

    def solid_report(
        self,
        *,
        reference_volume: float | None = None,
        reference_area: float | None = None,
    ) -> dict[str, object]:
        """Topology and integrals of the stitched shell, with optional references."""
        mesh = self.mesh
        counts = meshlib.edge_face_counts(mesh)
        volume = meshlib.total_signed_volume(mesh)
        area = meshlib.surface_area(mesh)
        report: dict[str, object] = {
            "vertices": len(mesh.vertices),
            "triangles": len(mesh.faces),
            "watertight": meshlib.is_watertight(mesh),
            "edge_manifold": meshlib.is_edge_manifold(mesh),
            "boundary_edges": len(meshlib.boundary_edges(mesh)),
            "non_manifold_edges": sum(1 for count in counts.values() if count > 2),
            "components": _component_count(mesh),
            "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
            "signed_volume_m3": volume,
            "surface_area_m2": area,
            "open_sides": len(self.open_sides),
        }
        if reference_volume is not None:
            report["reference_volume_m3"] = reference_volume
            report["volume_relative_error"] = (
                (volume - reference_volume) / reference_volume if reference_volume else math.inf
            )
        if reference_area is not None:
            report["reference_area_m2"] = reference_area
            report["area_relative_error"] = (
                (area - reference_area) / reference_area if reference_area else math.inf
            )
        return report


# ---------------------------------------------------------------------------
# grid helpers
# ---------------------------------------------------------------------------
def _grid_parameters(surface: NURBSSurface, nu: int, nv: int) -> tuple[list[float], list[float]]:
    (u_low, u_high), (v_low, v_high) = surface.domain
    u_values = [u_low + (u_high - u_low) * index / nu for index in range(nu + 1)]
    v_values = [v_low + (v_high - v_low) * index / nv for index in range(nv + 1)]
    u_values[0], u_values[-1] = u_low, u_high
    v_values[0], v_values[-1] = v_low, v_high
    return u_values, v_values


def _slot_ordering(side: str, nu: int, nv: int) -> list[tuple[int, int]]:
    """Slots ``(i, j)`` along a side, in increasing parameter order."""
    if side == "v0":
        return [(i, 0) for i in range(nu + 1)]
    if side == "v1":
        return [(i, nv) for i in range(nu + 1)]
    if side == "u0":
        return [(0, j) for j in range(nv + 1)]
    return [(nu, j) for j in range(nv + 1)]


def _collar_normal(surface: NURBSSurface, u: float, v: float) -> Vec3:
    """Unit normal a patch has at a boundary node (the same call as in the interior)."""
    return surface.normal(u, v)


def _seam_normal(
    surface: NURBSSurface,
    side: str,
    order: Sequence[tuple[int, int]],
    u_values: Sequence[float],
    v_values: Sequence[float],
    k: int,
) -> tuple[Vec3, bool]:
    """Normal at seam node ``k``, or at its nearest neighbour when the exact point is singular.

    A non-singular side can still end at a pole in the *other* parameter direction (a
    patch whose ``u0`` edge is an apex and whose ``v0`` edge is seamed), and there the
    Jacobian vanishes exactly.  The measurement then walks along the seam to the
    closest node with a defined normal and says so (``inset`` counts these), rather
    than reporting a fabricated angle or raising in the middle of a stitch.
    """
    for distance in range(len(order)):
        for index in (k - distance, k + distance):
            if not 0 <= index < len(order):
                continue
            i, j = order[index]
            try:
                return surface.normal(u_values[i], v_values[j]), index != k
            except NURBSError:
                continue
    raise PatchNetworkError(
        f"side {side} has no node with a defined normal; the patch is singular along it"
    )


class _Union:
    """Minimal union-find over node slots; the welder's only state."""

    def __init__(self, size: int) -> None:
        self._parent = list(range(size))

    def find(self, index: int) -> int:
        parent = self._parent
        root = index
        while parent[root] != root:
            root = parent[root]
        while parent[index] != root:
            parent[index], index = root, parent[index]
        return root

    def union(self, first: int, second: int) -> None:
        a, b = self.find(first), self.find(second)
        if a == b:
            return
        if a < b:
            self._parent[b] = a
        else:
            self._parent[a] = b


def stitch(network: PatchNetwork, *, nu: int = 16, nv: int = 16) -> Stitched:
    """Sample ``network`` on a common ``nu`` x ``nv`` grid and weld its seams.

    The mesh is built with the seam nodes already unified (index welding), and the
    distance that welding removed is reported per seam in :attr:`Stitched.seams`.
    """
    if nu < 1 or nv < 1:
        raise PatchNetworkError(f"stitching needs at least 1 x 1 samples, got {nu} x {nv}")
    surfaces = network.surfaces
    for position, seam in enumerate(network.seams):
        for patch, side in (seam.first, seam.second):
            axis = "u" if side.startswith("u") else "v"
            if surfaces[patch].singular_boundary_point(axis, upper=side.endswith("1")) is not None:
                raise PatchNetworkError(
                    f"seam {position} claims side {side} of patch {patch}, which is a pole: "
                    "every node of it is the same point, so leave it unclaimed and stitch "
                    "will collapse it (reported in pole_sides)"
                )
        first_count = _side_samples(seam.first[1], nu, nv)
        second_count = _side_samples(seam.second[1], nu, nv)
        if first_count != second_count:
            raise PatchNetworkError(
                f"seam {position} joins a side sampled with {first_count} nodes to one "
                f"sampled with {second_count}; a {seam.first[1]}/{seam.second[1]} seam "
                f"needs nu == nv (got {nu} x {nv})"
            )

    parameters: list[tuple[list[float], list[float]]] = []
    positions: list[list[list[Vec3]]] = []
    for surface in surfaces:
        u_values, v_values = _grid_parameters(surface, nu, nv)
        parameters.append((u_values, v_values))
        positions.append(
            [[surface.evaluate(u, v) for v in v_values] for u in u_values]
        )

    # slot (patch, i, j) <-> flat index
    width = (nu + 1) * (nv + 1)

    def slot(patch: int, i: int, j: int) -> int:
        return patch * width + i * (nv + 1) + j

    union = _Union(len(surfaces) * width)

    # --- weld the declared seams -------------------------------------------
    records: list[SeamRecord] = []
    for position, seam in enumerate(network.seams):
        first_patch, first_side = seam.first
        second_patch, second_side = seam.second
        first_order = _slot_ordering(first_side, nu, nv)
        second_order = _slot_ordering(second_side, nu, nv)
        count = len(first_order)
        gaps: list[float] = []
        angles: list[float] = []
        insets = 0
        for k in range(count):
            paired = count - 1 - k if seam.reversed else k
            (i_a, j_a) = first_order[k]
            (i_b, j_b) = second_order[paired]
            a = positions[first_patch][i_a][j_a]
            b = positions[second_patch][i_b][j_b]
            gaps.append(math.dist(a, b))
            normal_a, inset_a = _seam_normal(
                surfaces[first_patch], first_side, first_order,
                parameters[first_patch][0], parameters[first_patch][1], k,
            )
            normal_b, inset_b = _seam_normal(
                surfaces[second_patch], second_side, second_order,
                parameters[second_patch][0], parameters[second_patch][1], paired,
            )
            if inset_a or inset_b:
                insets += 1
            dot = sum(x * y for x, y in zip(normal_a, normal_b))
            angles.append(math.degrees(math.acos(max(-1.0, min(1.0, dot)))))
            union.union(slot(first_patch, i_a, j_a), slot(second_patch, i_b, j_b))
        records.append(
            SeamRecord(
                label=seam.label or f"seam{position}",
                first=seam.first,
                second=seam.second,
                reversed=seam.reversed,
                samples=count,
                gap_max=max(gaps),
                gap_mean=sum(gaps) / count,
                angle_max_deg=max(angles),
                angle_mean_deg=sum(angles) / count,
                inset_normals=insets,
            )
        )

    # --- collapse poles on the unclaimed sides ------------------------------
    pole_sides: list[tuple[int, str]] = []
    open_sides: list[tuple[int, str]] = []
    claimed = set(network.claimed_sides)
    for patch, surface in enumerate(surfaces):
        for side in SIDES:
            if (patch, side) in claimed:
                continue
            axis = "u" if side.startswith("u") else "v"
            upper = side.endswith("1")
            pole = surface.singular_boundary_point(axis, upper=upper)
            if pole is None:
                open_sides.append((patch, side))
                continue
            pole_sides.append((patch, side))
            order = _slot_ordering(side, nu, nv)
            root = slot(patch, *order[0])
            for i, j in order[1:]:
                union.union(root, slot(patch, i, j))

    # --- materialise vertices, one per weld class ---------------------------
    builder = MeshBuilder(units="m")
    vertex_of_root: dict[int, int] = {}
    nodes = [[[0] * (nv + 1) for _ in range(nu + 1)] for _ in surfaces]
    for patch in range(len(surfaces)):
        for i in range(nu + 1):
            for j in range(nv + 1):
                root = union.find(slot(patch, i, j))
                index = vertex_of_root.get(root)
                if index is None:
                    index = builder.add_vertex(positions[patch][i][j])
                    vertex_of_root[root] = index
                nodes[patch][i][j] = index

    for patch in range(len(surfaces)):
        for i in range(nu):
            for j in range(nv):
                builder.add_quad(
                    nodes[patch][i][j],
                    nodes[patch][i + 1][j],
                    nodes[patch][i + 1][j + 1],
                    nodes[patch][i][j + 1],
                )

    total_slots = len(surfaces) * width
    return Stitched(
        mesh=builder.build(),
        seams=tuple(records),
        open_sides=tuple(open_sides),
        pole_sides=tuple(pole_sides),
        nu=nu,
        nv=nv,
        welded_nodes=total_slots - len(vertex_of_root),
        patch_count=len(surfaces),
    )


def assemble(network: PatchNetwork, *, nu: int = 16, nv: int = 16) -> Stitched:
    """:func:`stitch` under the name a caller looks for."""
    return stitch(network, nu=nu, nv=nv)


def _component_count(mesh: Mesh) -> int:
    """Number of face-connected shells (union-find over shared vertices)."""
    parent: dict[int, int] = {}

    def find(index: int) -> int:
        parent.setdefault(index, index)
        root = index
        while parent[root] != root:
            root = parent[root]
        while parent[index] != root:
            parent[index], index = root, parent[index]
        return root

    for face in mesh.faces:
        a, b, c = face
        root = find(a)
        for other in (b, c):
            other_root = find(other)
            if other_root != root:
                parent[other_root] = root
    return len({find(index) for index in parent})


# ---------------------------------------------------------------------------
# constructors: the two shapes a caller needs to build a network by hand
# ---------------------------------------------------------------------------
def bilinear_patch(origin: Vec3, along_u: Vec3, along_v: Vec3) -> NURBSSurface:
    """A planar (bilinear) patch: the parallelogram spanned by two edge vectors.

    Its normal is ``along_u x along_v``, so the arguments also set the winding -- the
    sign of the solid's volume comes from these three vectors.  A bilinear patch is
    exact under tessellation, which is what makes it the anchor for the stitcher: a
    closed box of six bilinear patches has an analytic volume *and* an exact
    tessellation, so any volume error is a seam error.
    """
    ox, oy, oz = (float(origin[0]), float(origin[1]), float(origin[2]))
    ux, uy, uz = (float(along_u[0]), float(along_u[1]), float(along_u[2]))
    vx, vy, vz = (float(along_v[0]), float(along_v[1]), float(along_v[2]))
    net = tuple(
        tuple(
            (ox + i * ux + j * vx, oy + i * uy + j * vy, oz + i * uz + j * vz)
            for j in (0, 1)
        )
        for i in (0, 1)
    )
    return NURBSSurface(
        degree=(1, 1),
        knots=((0.0, 0.0, 1.0, 1.0), (0.0, 0.0, 1.0, 1.0)),
        control_net=net,
    )


def revolve_patch(profile: NURBSCurve, theta_start: float, theta_sweep: float) -> NURBSSurface:
    """Revolve a ``(radius, z)`` profile arc about the ``z`` axis through an exact arc.

    ``profile`` is a rational curve in the ``(rho, z)`` half-plane; its control points
    are read as ``(rho, z, unused)`` and only the first two components are used.  The
    patch is the tensor product of that profile with an **exact** rational quadratic
    arc of ``theta_sweep`` radians (``theta_sweep <= pi``, the same limit
    :func:`tensorcad.nurbs.circular_arc` documents), so ``u`` is the angle and ``v``
    is the profile parameter.

    Pairing a profile whose first or last control radius is zero (a cone apex, a pole)
    with the whole revolution gives a **singular** side, which :func:`stitch`
    collapses to one vertex and reports in ``pole_sides``.
    """
    if not 0.0 < theta_sweep <= math.pi * (1.0 + 1e-12):
        raise NURBSError(f"theta_sweep must lie in (0, pi], got {theta_sweep!r}")
    arc = circular_arc(1.0, theta_sweep, start_angle=theta_start)
    net = tuple(
        tuple(
            (radius * point[0], radius * point[1], height)
            for radius, height in ((profile.control_points[j][0], profile.control_points[j][1])
                                   for j in range(len(profile.control_points)))
        )
        for point in arc.control_points
    )
    weights = tuple(
        tuple(arc.weights[i] * profile.weights[j] for j in range(len(profile.control_points)))
        for i in range(len(arc.control_points))
    )
    return NURBSSurface(
        degree=(arc.degree, profile.degree),
        knots=(arc.knots, profile.knots),
        control_net=net,
        weights=weights,
    )


def _side_endpoints(surface: NURBSSurface, side: str) -> tuple[Vec3, Vec3]:
    """The two corners of a side, in increasing parameter order."""
    (u_low, u_high), (v_low, v_high) = surface.domain
    if side == "v0":
        return surface.evaluate(u_low, v_low), surface.evaluate(u_high, v_low)
    if side == "v1":
        return surface.evaluate(u_low, v_high), surface.evaluate(u_high, v_high)
    if side == "u0":
        return surface.evaluate(u_low, v_low), surface.evaluate(u_low, v_high)
    return surface.evaluate(u_high, v_low), surface.evaluate(u_high, v_high)


def _is_pole(surface: NURBSSurface, side: str) -> bool:
    axis = "u" if side.startswith("u") else "v"
    return surface.singular_boundary_point(axis, upper=side.endswith("1")) is not None


def derive_seams(
    surfaces: Sequence[NURBSSurface], *, tolerance: float = 1e-12
) -> tuple[Seam, ...]:
    """Pair the sides of ``surfaces`` whose corners coincide, and only those.

    The stitcher takes declarations, and a caller who already knows the topology
    should write them.  This helper exists for the case where the topology is a
    *consequence of the geometry* (six faces of a box, the quadrants of a revolved
    body) and writing thirty declarations by hand would be the error-prone part.

    Rules, all of them conservative:

    * pole sides are never candidates (they collapse, they do not seam);
    * a pair matches when both endpoints coincide within ``tolerance``, in either
      direction -- the direction becomes ``Seam.reversed``;
    * a side with **more than one** candidate is left unmatched rather than guessed,
      because a guess that welds the wrong two sides is worse than a report of an
      open side.  :func:`stitch` will then list it in ``open_sides``.
    """
    if tolerance <= 0.0:
        raise PatchNetworkError(f"tolerance must be positive, got {tolerance!r}")
    sides = [
        (patch, side)
        for patch in range(len(surfaces))
        for side in SIDES
        if not _is_pole(surfaces[patch], side)
    ]
    endpoints = {key: _side_endpoints(surfaces[key[0]], key[1]) for key in sides}
    used: set[tuple[int, str]] = set()
    seams: list[Seam] = []
    for position, first in enumerate(sides):
        if first in used:
            continue
        candidates: list[tuple[tuple[int, str], bool]] = []
        for second in sides[position + 1:]:
            if second in used:
                continue
            a = endpoints[first]
            b = endpoints[second]
            same = max(math.dist(a[0], b[0]), math.dist(a[1], b[1]))
            flipped = max(math.dist(a[0], b[1]), math.dist(a[1], b[0]))
            if min(same, flipped) <= tolerance:
                candidates.append((second, flipped < same))
        if len(candidates) != 1:
            continue
        second, reversed_ = candidates[0]
        seams.append(Seam(first, second, reversed=reversed_, label=f"derived{len(seams)}"))
        used.add(first)
        used.add(second)
    return tuple(seams)


def seam_network(
    surfaces: Iterable[NURBSSurface], seams: Iterable[Seam] = ()
) -> PatchNetwork:
    """Convenience constructor: ``PatchNetwork(tuple(surfaces), tuple(seams))``."""
    return PatchNetwork(tuple(surfaces), tuple(seams))
