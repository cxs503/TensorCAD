"""Body-fitted volume mesh: a hull surface turned into cells that fill the body.

The package could already describe the *outside* of a hull (:mod:`tensorcad.hull`) and
integrate the inside of it (:mod:`tensorcad.hydrostatics`), and it could rasterise an
implicit field into axis-aligned boxes (:mod:`tensorcad.voxel`). What it could not do is
hand a solver the volume: cells that fill the body, whose wall *is* the wetted surface,
and that are thin enough near that wall to resolve a boundary layer. A voxel grid is not
that -- its staircase wall is the wrong shape at exactly the scale the boundary layer
cares about.

The construction here is the one every hull modeller starts from: the hull families in
this package are all of the form *half breadth as a function of station and level*,
``b(x, z)``, so the body can be filled by a structured (index-space) grid whose nodes sit
at ``(x_i, y, z_j)``. The wall nodes are the hull surface's own nodes, and the *breadth*
direction is graded from each wall inward with a prescribed first-layer thickness and
growth ratio. The result is a conforming hexahedral mesh that can be split into
tetrahedra for solvers that want them, a VTK file any post-processor can read, and a
quality report that says out loud how good the cells are.

What it is not
--------------
* It is not a general B-rep tetrahedraliser. The input is a :class:`SectionGrid`: a
  half-breadth sampler on a rectilinear station/level grid, single-valued in ``y`` at
  every ``(x, z)``. Bodies that are not of that form are refused, not approximated.
* It does not smooth or optimise the cells. The quality numbers are *measured* so a
  caller can decide; there is no hidden projection pass that would make a bad grid look
  good.
* Collapsed columns (the keel line and the two end stations, where the hull genuinely has
  zero breadth) are kept as collapsed cells rather than nudged apart. A structured hull
  grid always contains those wedges; the summary counts them instead of hiding them.

Units are metres throughout, in the package's frame: ``x`` aft-to-fore, ``y`` to
starboard, ``z`` up.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Iterable, Sequence

if TYPE_CHECKING:  # pragma: no cover - avoids a runtime import cycle
    from .hull import HullParameters

Vec3 = tuple[float, float, float]

#: Cell quality below this scaled Jacobian is reported as distorted. A tetrahedron of a
#: cube (the Kuhn split) sits between roughly 0.1 and 0.5, a regular tetrahedron at
#: 1/sqrt(2) = 0.707, and a flat cell at 0.
DISTORTION_THRESHOLD = 0.1

#: The six ways to walk from one corner of a cell to the opposite one. Splitting each
#: hexahedron along all six paths tiles it exactly and is conforming across cell faces: on
#: a shared face every cell induces the diagonal that joins the two low corners of that
#: face, so neighbouring cells never disagree about how the face is cut.
_PERMUTATIONS: tuple[tuple[int, int, int], ...] = (
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0),
)


class VolumeMeshError(ValueError):
    """Raised when a grid or a mesh request cannot be built honestly."""


# ---------------------------------------------------------------------------------------
# Input: the surface the volume is built from
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SectionGrid:
    """Half breadth sampled on a rectilinear station/level grid.

    ``breadth(i, j)`` is the half breadth in metres at station ``xs[i]`` and level
    ``zs[j]``. It must be finite and non-negative; zero means the surface genuinely closes
    on the centreline there (keel line, stem, transom), which the builder treats as a
    collapsed column rather than as an error.

    The contract is deliberately narrow: a body whose half breadth is not single valued in
    ``y`` over its ``(x, z)`` box cannot be represented, and :meth:`validate` says so.
    """

    xs: tuple[float, ...]
    zs: tuple[float, ...]
    breadth: Callable[[int, int], float]

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if len(self.xs) < 2:
            raise VolumeMeshError(f"need at least two stations, got {len(self.xs)}")
        if len(self.zs) < 2:
            raise VolumeMeshError(f"need at least two levels, got {len(self.zs)}")
        for name, values in (("xs", self.xs), ("zs", self.zs)):
            for value in values:
                if not math.isfinite(value):
                    raise VolumeMeshError(f"{name} contains a non-finite coordinate {value!r}")
            steps = [b - a for a, b in zip(values, values[1:])]
            if any(step <= 0.0 for step in steps):
                raise VolumeMeshError(f"{name} must be strictly increasing")
        for i in range(len(self.xs)):
            for j in range(len(self.zs)):
                value = self.breadth(i, j)
                if not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise VolumeMeshError(
                        f"breadth at station {i} level {j} is not a finite number: {value!r}"
                    )
                if value < 0.0:
                    raise VolumeMeshError(
                        f"breadth at station {i} level {j} is negative ({value!r}); "
                        "a body that is single valued in y needs b(x, z) >= 0"
                    )

    @property
    def nx(self) -> int:
        return len(self.xs)

    @property
    def nz(self) -> int:
        return len(self.zs)

    @classmethod
    def from_hull(cls, parameters: "HullParameters") -> "SectionGrid":
        """The sampler of :func:`tensorcad.hull.generate_hull`, station for station.

        The hull's own section outline is walked keelward, so level ``j`` of the grid is
        ring point ``j`` of the section. Using the same sampler as the surface mesh is what
        makes the volume mesh's wall coincide with the hull mesh's shell nodes exactly,
        instead of to within a resampling error.
        """
        from .hull import _section_point  # same-package sampler, index for index

        parameters.validate()
        xs = tuple(parameters.length * i / (parameters.nx - 1) for i in range(parameters.nx))
        zs = tuple(
            -parameters.draft * (1.0 - j / (parameters.nz - 1)) for j in range(parameters.nz)
        )

        def breadth(i: int, j: int) -> float:
            return abs(_section_point(parameters, i, j)[1])

        return cls(xs=xs, zs=zs, breadth=breadth)


def box_section_grid(
    *, length: float, beam: float, draft: float, nx: int = 11, nz: int = 7
) -> SectionGrid:
    """A rectangular box as a :class:`SectionGrid` -- the analytic reference case.

    A box has constant half breadth, so its volume mesh is exact: the hexahedra tile the
    box with no discretisation error at all. That is what makes it a zero-error control.
    """
    for name, value in (("length", length), ("beam", beam), ("draft", draft)):
        if not math.isfinite(value) or value <= 0.0:
            raise VolumeMeshError(f"{name} must be a finite positive number, got {value!r}")
    if nx < 2 or nz < 2:
        raise VolumeMeshError("a box grid needs at least two stations and two levels")
    xs = tuple(length * i / (nx - 1) for i in range(nx))
    zs = tuple(-draft * (1.0 - j / (nz - 1)) for j in range(nz))
    return SectionGrid(xs=xs, zs=zs, breadth=lambda i, j: 0.5 * beam)


# ---------------------------------------------------------------------------------------
# The mesh
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class VolumeMesh:
    """A conforming hexahedral volume mesh, with an optional tetrahedral split.

    Node ``(i, j, k)`` is unique for every index triple, so the grid stays structured: the
    hexahedra of a structured grid share their corner nodes by construction, which is why
    :meth:`is_conforming` can pass without any welding step. Collapsed columns (zero
    breadth) leave *distinct* nodes at *equal* coordinates on purpose: welding them would
    turn the wedges into degenerate index arithmetic, and every consumer downstream would
    have to re-derive what is a wall and what is a fold.
    """

    vertices: tuple[Vec3, ...]
    hexes: tuple[tuple[int, int, int, int, int, int, int, int], ...]
    tets: tuple[tuple[int, int, int, int], ...]
    shape: tuple[int, int, int]
    columns: tuple[tuple[float, ...], ...]
    station_lines: tuple[float, ...]
    level_lines: tuple[float, ...]
    index_map: tuple[int, ...]
    grid: SectionGrid = field(repr=False, compare=False)
    layer_spec: "LayerSpec | None" = field(default=None, repr=False, compare=False)

    # -- geometry --------------------------------------------------------------------
    def _hex_tets(self, hex_nodes: Sequence[int]) -> tuple[tuple[int, int, int, int], ...]:
        corners = _corner_index_map(hex_nodes)
        out = []
        for permutation in _PERMUTATIONS:
            path = [corners[(0, 0, 0)]]
            steps = [0, 0, 0]
            for axis in permutation:
                steps[axis] = 1
                path.append(corners[(steps[0], steps[1], steps[2])])
            out.append((path[0], path[1], path[2], path[3]))
        return tuple(out)

    def hex_volume(self, hex_nodes: Sequence[int]) -> float:
        """Volume of one hexahedron, as the sum of its six Kuhn tetrahedra.

        The six paths from one corner to the opposite one tile the cell without overlap, so
        the volumes add. Each path's orientation depends on the permutation, which is why
        the absolute values are taken: what is being measured is the tiling, not a winding
        number.
        """
        return sum(
            abs(_tet_signed_volume(self.vertices, tet)) for tet in self._hex_tets(hex_nodes)
        )

    def volume(self) -> float:
        """Filled volume in cubic metres (absolute cell volumes, no cancellation)."""
        return sum(abs(self.hex_volume(hex_nodes)) for hex_nodes in self.hexes)

    def tet_volume(self) -> float:
        """Volume of the tetrahedral split; equals :meth:`volume` to round-off."""
        return sum(abs(_tet_signed_volume(self.vertices, tet)) for tet in self.tets)

    def node_id(self, i: int, j: int, k: int) -> int:
        """Node id of station ``i``, level ``j``, breadth ``k`` (0 = port wall)."""
        nx, nz, breadth_cells = self.shape
        if not (0 <= i < nx and 0 <= j < nz and 0 <= k <= breadth_cells):
            raise VolumeMeshError(f"node index out of range: {(i, j, k)} for shape {self.shape}")
        # Node ids follow first use, not index order (a cell adds the nodes of four
        # columns), so the map is explicit rather than arithmetic.
        return self.index_map[(i * nz + j) * (breadth_cells + 1) + k]

    def cell_nodes(self, i: int, j: int, k: int) -> tuple[int, int, int, int, int, int, int, int]:
        """The eight node ids of cell ``(i, j, k)``, in the mesh's VTK ordering."""
        if not (0 <= i < self.shape[0] - 1 and 0 <= j < self.shape[1] - 1 and 0 <= k < self.shape[2]):
            raise VolumeMeshError(f"cell index out of range: {(i, j, k)} for shape {self.shape}")
        return (
            self.node_id(i, j, k),
            self.node_id(i + 1, j, k),
            self.node_id(i + 1, j, k + 1),
            self.node_id(i, j, k + 1),
            self.node_id(i, j + 1, k),
            self.node_id(i + 1, j + 1, k),
            self.node_id(i + 1, j + 1, k + 1),
            self.node_id(i, j + 1, k + 1),
        )

    def cell_centroid(self, hex_nodes: Sequence[int]) -> Vec3:
        n = len(hex_nodes)
        return (
            sum(self.vertices[v][0] for v in hex_nodes) / n,
            sum(self.vertices[v][1] for v in hex_nodes) / n,
            sum(self.vertices[v][2] for v in hex_nodes) / n,
        )

    # -- topology --------------------------------------------------------------------
    def face_multiplicities(self) -> dict[tuple[int, int, int], int]:
        """How often each triangle occurs over the whole tetrahedral split.

        One is a boundary face, two is an interior face, three or more means a non-manifold
        fan. A caller can therefore check that the cells a solver sees form a closed tiling
        of the body without trusting the builder's word for it.
        """
        counts: dict[tuple[int, int, int], int] = {}
        for tet in self.tets:
            for a in range(4):
                for b in range(a + 1, 4):
                    for c in range(b + 1, 4):
                        face = tuple(sorted((tet[a], tet[b], tet[c])))
                        counts[face] = counts.get(face, 0) + 1
        return counts

    def is_conforming(self) -> bool:
        """True when every triangle is shared by at most two tetrahedra."""
        return all(count <= 2 for count in self.face_multiplicities().values())

    def boundary_faces(self) -> tuple[tuple[int, int, int], ...]:
        """Triangles that occur exactly once: the outer skin of the volume mesh."""
        return tuple(
            sorted(face for face, count in self.face_multiplicities().items() if count == 1)
        )

    def counts(self) -> dict[str, Any]:
        zero_volume = sum(
            1 for tet in self.tets if abs(_tet_signed_volume(self.vertices, tet)) == 0.0
        )
        return {
            "nodes": len(self.vertices),
            "hexes": len(self.hexes),
            "tets": len(self.tets),
            "zero_volume_tets": zero_volume,
            "shape_stations_levels_breadth": list(self.shape),
        }

    # -- quality ---------------------------------------------------------------------
    def quality(self, threshold: float = DISTORTION_THRESHOLD) -> dict[str, Any]:
        """Scaled-Jacobian statistics of the cells, split into wall cells and core cells.

        The scaled Jacobian of a cell is its volume over the volume of the regular cell
        with the same edges -- ``det(e1, e2, e3) / (|e1||e2||e3|)``, the best corner of a
        tetrahedron and the worst corner of a hexahedron. It is 1 for a regular
        hexahedron, ``1/sqrt(2) = 0.7071`` for a regular tetrahedron, and 0 for a flat
        cell. Negative means inverted, which is a defect and not a quality statement.

        Wall cells are the ones in the first breadth layer on either side; they are
        anisotropic by design, and reporting them together with the isotropic core would
        dress an aspect-ratio decision up as distortion.
        """
        hex_wall: list[float] = []
        hex_core: list[float] = []
        tet_wall: list[float] = []
        tet_core: list[float] = []
        collapsed = 0
        breadth_cells = self.shape[2]
        for index, hex_nodes in enumerate(self.hexes):
            k = index % breadth_cells
            value = min(_hex_corner_jacobians(self.vertices, hex_nodes))
            if value <= 0.0:
                # A wedge on the keel line or at a collapsed station: a real cell of a
                # structured hull grid, but not a cell a quality statistic can average in.
                collapsed += 1
                continue
            if k == 0 or k == breadth_cells - 1:
                hex_wall.append(value)
                bucket = tet_wall
            else:
                hex_core.append(value)
                bucket = tet_core
            bucket.extend(
                tet_scaled_jacobian(self.vertices, _positive_tet(self.vertices, tet))
                for tet in self._hex_tets(hex_nodes)
                if abs(_tet_signed_volume(self.vertices, tet)) > 0.0
            )
        return {
            "threshold": threshold,
            "hex": {
                "wall": _jacobian_stats(hex_wall, threshold),
                "core": _jacobian_stats(hex_core, threshold),
                "all": _jacobian_stats(hex_wall + hex_core, threshold),
            },
            "tet": {
                "wall": _jacobian_stats(tet_wall, threshold),
                "core": _jacobian_stats(tet_core, threshold),
                "all": _jacobian_stats(tet_wall + tet_core, threshold),
            },
            "collapsed_cells": collapsed,
        }

    # -- output ----------------------------------------------------------------------
    def to_vtk(self, path: str | Path, *, cells: str = "hex") -> Path:
        """Write a legacy VTK unstructured grid (ASCII), readable by ParaView et al."""
        target = Path(path)
        if cells not in ("hex", "tet"):
            raise VolumeMeshError(f"cells must be 'hex' or 'tet', got {cells!r}")
        if cells == "tet" and not self.tets:
            raise VolumeMeshError(
                "this mesh has no tetrahedral split (built with split_tets=False)"
            )
        if cells == "hex":
            connectivity: Iterable[Sequence[int]] = self.hexes
            vtk_type, per_cell = 12, 8  # VTK_HEXAHEDRON
        else:
            connectivity = self.tets
            vtk_type, per_cell = 10, 4  # VTK_TETRA
        cell_list = list(connectivity)
        lines = [
            "# vtk DataFile Version 3.0",
            "TensorCAD volume mesh",
            "ASCII",
            "DATASET UNSTRUCTURED_GRID",
            f"POINTS {len(self.vertices)} double",
        ]
        lines.extend(f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in self.vertices)
        lines.append(f"CELLS {len(cell_list)} {len(cell_list) * (per_cell + 1)}")
        lines.extend(f"{per_cell} " + " ".join(str(node) for node in cell) for cell in cell_list)
        lines.append(f"CELL_TYPES {len(cell_list)}")
        lines.extend(str(vtk_type) for _ in cell_list)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return target

    def summary(self) -> dict[str, Any]:
        """Everything a consumer or a document needs, computed rather than assumed."""
        volume = self.volume()
        counts = self.counts()
        collapsed = sum(1 for column in self.columns if not column or column[0] == 0.0)
        graded = sum(1 for column in self.columns if len(_graded_prefix(column)) >= 2)
        return {
            "units": "m",
            "counts": counts,
            "volume_m3": volume,
            "tet_volume_m3": self.tet_volume(),
            "volume_split_relative_difference": (
                abs(self.tet_volume() - volume) / volume if volume else None
            ),
            "conforming": self.is_conforming(),
            "boundary_triangles": len(self.boundary_faces()),
            "collapsed_columns": collapsed,
            "graded_columns": graded,
            "columns": len(self.columns),
            "quality": self.quality(),
            "first_layer_requested_m": (
                None if self.layer_spec is None else self.layer_spec.first_layer
            ),
            "growth_requested": None if self.layer_spec is None else self.layer_spec.growth,
            "layers_per_side": self.shape[2] // 2,
        }


@dataclass(frozen=True)
class LayerSpec:
    """The wall grading: first layer thickness, growth ratio, layers per side."""

    first_layer: float
    growth: float = 1.2
    layers: int = 4

    def __post_init__(self) -> None:
        if not math.isfinite(self.first_layer) or self.first_layer <= 0.0:
            raise VolumeMeshError(
                f"first_layer must be a finite positive length, got {self.first_layer!r}"
            )
        if not math.isfinite(self.growth) or self.growth < 1.0:
            raise VolumeMeshError(
                f"growth must be a finite ratio >= 1 (1 is uniform), got {self.growth!r}"
            )
        if self.layers < 1:
            raise VolumeMeshError(f"layers must be >= 1, got {self.layers!r}")


# ---------------------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------------------


def _corner_index_map(hex_nodes: Sequence[int]) -> dict[tuple[int, int, int], int]:
    """Map the eight ``(di, dj, dk)`` corners of a cell to node ids.

    The ordering is the VTK one: ``0, 1, 2, 3`` walk the low-**level** (low ``z``) face and
    ``4 ... 7`` the same corners one level up. From node 0 the three edges run along ``+x``
    (node 1), ``+y`` (node 3) and ``+z`` (node 4) in that order, which is right handed in
    the package frame ``(x fore, y starboard, z up)`` -- so the cells carry a positive
    Jacobian in any post-processor that honours the VTK convention.
    ``_corner_index_map`` maps ``(di, dj, dk)`` offsets -- station, level, breadth -- onto
    that sequence.
    """
    n = hex_nodes
    return {
        (0, 0, 0): n[0],
        (1, 0, 0): n[1],
        (1, 0, 1): n[2],
        (0, 0, 1): n[3],
        (0, 1, 0): n[4],
        (1, 1, 0): n[5],
        (1, 1, 1): n[6],
        (0, 1, 1): n[7],
    }


def _tet_signed_volume(vertices: Sequence[Vec3], tet: Sequence[int]) -> float:
    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = (vertices[i] for i in tet)
    ax, ay, az = x1 - x0, y1 - y0, z1 - z0
    bx, by, bz = x2 - x0, y2 - y0, z2 - z0
    cx, cy, cz = x3 - x0, y3 - y0, z3 - z0
    det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
    return det / 6.0


def _det3(a: Vec3, b: Vec3, c: Vec3) -> float:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _hex_corner_jacobians(
    vertices: Sequence[Vec3], hex_nodes: Sequence[int]
) -> tuple[float, ...]:
    """Scaled Jacobian at each of the eight corners of a hexahedron."""
    corners = _corner_index_map(hex_nodes)
    out = []
    for corner, node in corners.items():
        origin = vertices[node]
        vectors: list[Vec3] = []
        # The determinant is taken with the edges ordered (station, breadth, level) --
        # (x, y, z) -- because that triple is right handed in the package frame. The cell's
        # *node* order is a different thing and is left alone by this.
        for axis in (0, 2, 1):
            step = [corner[0], corner[1], corner[2]]
            step[axis] = 1 - step[axis]
            other = vertices[corners[(step[0], step[1], step[2])]]
            # Always orient the edge along the increasing cell axis: at a high-side corner
            # the only neighbour available is the low one, so the difference is negated.
            sign = 1.0 if corner[axis] == 0 else -1.0
            vectors.append(
                (
                    sign * (other[0] - origin[0]),
                    sign * (other[1] - origin[1]),
                    sign * (other[2] - origin[2]),
                )
            )
        denominator = _norm(vectors[0]) * _norm(vectors[1]) * _norm(vectors[2])
        out.append(_det3(vectors[0], vectors[1], vectors[2]) / denominator if denominator else 0.0)
    return tuple(out)


def tet_scaled_jacobian(vertices: Sequence[Vec3], tet: Sequence[int]) -> float:
    """Scaled Jacobian of a tetrahedron: the *worst* of its four corners.

    The corners are visited cyclically (``0 -> 1 -> 2 -> 3 -> 0``) with the handedness
    restored on every other one, so a correctly wound tetrahedron reports four positive
    corner Jacobians. The worst of them is the usual engineering convention, and it makes
    the number independent of where the vertex list happens to start.
    """
    worst = math.inf
    for corner in range(4):
        origin = vertices[tet[corner]]
        vectors = []
        for offset in (1, 2, 3):
            other = vertices[tet[(corner + offset) % 4]]
            vectors.append(
                (other[0] - origin[0], other[1] - origin[1], other[2] - origin[2])
            )
        denominator = _norm(vectors[0]) * _norm(vectors[1]) * _norm(vectors[2])
        value = _det3(vectors[0], vectors[1], vectors[2]) / denominator if denominator else 0.0
        # Walking the corners cyclically flips the handedness on every other one, so the
        # sign is restored before the worst corner is taken. Without this a correctly wound
        # tetrahedron would report two negative corners and be called inverted.
        worst = min(worst, value if corner % 2 == 0 else -value)
    return worst


def _positive_tet(vertices: Sequence[Vec3], tet: Sequence[int]) -> tuple[int, int, int, int]:
    """The same tet wound so that its volume is positive (a negative one is inverted)."""
    if _tet_signed_volume(vertices, tet) < 0.0:
        return (tet[0], tet[2], tet[1], tet[3])
    return (tet[0], tet[1], tet[2], tet[3])


def _jacobian_stats(values: Sequence[float], threshold: float) -> dict[str, Any]:
    if not values:
        return {
            "count": 0,
            "min": None,
            "p05": None,
            "median": None,
            "mean": None,
            "max": None,
            "below_threshold": 0,
            "below_threshold_fraction": None,
            "non_positive": 0,
        }
    ordered = sorted(values)
    count = len(ordered)
    below = sum(1 for value in ordered if value < threshold)
    return {
        "count": count,
        "min": ordered[0],
        "p05": ordered[max(0, int(0.05 * (count - 1)))],
        "median": ordered[count // 2],
        "mean": sum(ordered) / count,
        "max": ordered[-1],
        "below_threshold": below,
        "below_threshold_fraction": below / count,
        "non_positive": sum(1 for value in ordered if value <= 0.0),
    }


def _graded_prefix(thicknesses: Sequence[float]) -> list[float]:
    """The leading run of thicknesses that grows by more than 1 % -- the graded layers.

    The innermost cell of a column absorbs whatever breadth is left, so the *tail* of the
    list is not part of the geometric series. Returning the prefix lets a document state how
    many columns actually carry a full graded stack instead of implying that every column
    does.
    """
    out: list[float] = []
    for index, thickness in enumerate(thicknesses):
        if index == 0:
            out.append(thickness)
            continue
        if thickness > out[-1] * 1.01:
            out.append(thickness)
        else:
            break
    return out


def _layer_thicknesses(
    breadth: float, spec: LayerSpec, normal_scale: float
) -> tuple[float, ...]:
    """Thicknesses of the ``layers`` cells from one wall to the centreline.

    The first layer is requested as a *wall-normal* clearance, so it is stretched by
    ``normal_scale = sqrt(1 + (db/dx)^2 + (db/dz)^2)`` at that wall node: a step along
    ``y`` only clears the surface by ``step / normal_scale``. Every following layer grows by
    the requested ratio, and the innermost one absorbs the remainder so the column always
    reaches the centreline exactly. Columns too narrow for the requested stack keep a
    uniform fill, which :func:`_graded_prefix` reports rather than hides.
    """
    layers = spec.layers
    if breadth <= 0.0:
        return tuple(0.0 for _ in range(layers))
    first = spec.first_layer * normal_scale
    if first >= breadth * 0.5:
        return tuple(breadth / layers for _ in range(layers))
    graded: list[float] = []
    total = 0.0
    thickness = first
    for _ in range(layers - 1):
        if total + thickness > breadth * 0.9:
            break
        graded.append(thickness)
        total += thickness
        thickness *= spec.growth
    if not graded:
        return tuple(breadth / layers for _ in range(layers))
    remaining_slots = layers - len(graded)
    rest = (breadth - total) / remaining_slots
    if rest <= 0.0:  # pragma: no cover - guarded by the 0.9 test above
        return tuple(breadth / layers for _ in range(layers))
    graded.extend(rest for _ in range(remaining_slots))
    return tuple(graded)


def build_volume_mesh(
    grid: SectionGrid,
    *,
    first_layer: float,
    growth: float = 1.2,
    layers: int = 4,
    split_tets: bool = True,
) -> VolumeMesh:
    """Fill the body described by ``grid`` with graded hexahedra.

    Nodes are ``(station, level, breadth)``; the breadth direction runs port wall ->
    centreline -> starboard wall with ``layers`` cells per side, so a cross-section carries
    ``2 * layers`` cells across.
    """
    spec = LayerSpec(first_layer=first_layer, growth=growth, layers=layers)
    nx, nz, half = grid.nx, grid.nz, spec.layers
    breadth = [[grid.breadth(i, j) for j in range(nz)] for i in range(nx)]

    # Slope of the half breadth in each direction: central differences inside, one-sided at
    # the ends. It is used only to ask for a wall-normal first layer, so an O(dx^2) error
    # there moves the clearance by O(dx^2) -- measured in the example, not assumed.
    def slope(i: int, j: int, axis: int) -> float:
        if axis == 0:
            lo, hi = max(0, i - 1), min(nx - 1, i + 1)
            if lo == hi:
                return 0.0
            return (breadth[hi][j] - breadth[lo][j]) / (grid.xs[hi] - grid.xs[lo])
        lo, hi = max(0, j - 1), min(nz - 1, j + 1)
        if lo == hi:
            return 0.0
        return (breadth[i][hi] - breadth[i][lo]) / (grid.zs[hi] - grid.zs[lo])

    vertices: list[Vec3] = []
    node_id: dict[tuple[int, int, int], int] = {}
    columns: list[tuple[float, ...]] = []
    breadth_cells = 2 * half

    def add(i: int, j: int, k: int) -> int:
        key = (i, j, k)
        existing = node_id.get(key)
        if existing is not None:
            return existing
        y = column_positions[(i, j)][k]
        node_id[key] = len(vertices)
        vertices.append((grid.xs[i], y, grid.zs[j]))
        return node_id[key]

    column_positions: dict[tuple[int, int], list[float]] = {}
    for i in range(nx):
        for j in range(nz):
            b = breadth[i][j]
            normal_scale = math.sqrt(
                1.0 + slope(i, j, 0) ** 2 + slope(i, j, 1) ** 2
            )
            thicknesses = _layer_thicknesses(b, spec, normal_scale)
            columns.append(thicknesses)
            positions = [-b]
            running = -b
            for thickness in thicknesses:
                running += thickness
                positions.append(running)
            positions[half] = 0.0  # the centreline node is the centreline, exactly
            for k in range(half + 1, breadth_cells + 1):
                positions.append(-positions[breadth_cells - k])
            column_positions[(i, j)] = positions

    hexes: list[tuple[int, int, int, int, int, int, int, int]] = []
    for i in range(nx - 1):
        for j in range(nz - 1):
            for k in range(breadth_cells):
                # VTK ordering: nodes 0-3 walk the low level face, 4-7 the next level up.
                hexes.append(
                    (
                        add(i, j, k),
                        add(i + 1, j, k),
                        add(i + 1, j, k + 1),
                        add(i, j, k + 1),
                        add(i, j + 1, k),
                        add(i + 1, j + 1, k),
                        add(i + 1, j + 1, k + 1),
                        add(i, j + 1, k + 1),
                    )
                )

    mesh_stub = VolumeMesh(
        vertices=tuple(vertices),
        hexes=tuple(hexes),
        tets=(),
        shape=(nx, nz, breadth_cells),
        columns=tuple(columns),
        station_lines=grid.xs,
        level_lines=grid.zs,
        index_map=tuple(
            node_id[(i, j, k)]
            for i in range(nx)
            for j in range(nz)
            for k in range(breadth_cells + 1)
        ),
        grid=grid,
        layer_spec=spec,
    )

    # Every cell must be laid out with a consistent orientation (negative scaled Jacobian
    # at a corner means the corner order disagrees with the node layout, which is a builder
    # bug, not a property of the mesh). Zero is legal: that is a collapsed wedge.
    for hex_nodes in mesh_stub.hexes:
        if min(_hex_corner_jacobians(mesh_stub.vertices, hex_nodes)) < -1e-12:
            raise VolumeMeshError(
                "cell orientation is inconsistent with the node layout; this is a builder bug"
            )

    tets: list[tuple[int, int, int, int]] = []
    if split_tets:
        for hex_nodes in mesh_stub.hexes:
            for tet in mesh_stub._hex_tets(hex_nodes):
                tets.append(_positive_tet(mesh_stub.vertices, tet))

    return VolumeMesh(
        vertices=mesh_stub.vertices,
        hexes=mesh_stub.hexes,
        tets=tuple(tets),
        shape=mesh_stub.shape,
        columns=mesh_stub.columns,
        station_lines=grid.xs,
        level_lines=grid.zs,
        index_map=mesh_stub.index_map,
        grid=grid,
        layer_spec=spec,
    )


def volume_mesh_from_hull(
    parameters: "HullParameters",
    *,
    first_layer: float,
    growth: float = 1.2,
    layers: int = 4,
    split_tets: bool = True,
) -> VolumeMesh:
    """Convenience wrapper: the body of :func:`tensorcad.hull.generate_hull`'s hull."""
    return build_volume_mesh(
        SectionGrid.from_hull(parameters),
        first_layer=first_layer,
        growth=growth,
        layers=layers,
        split_tets=split_tets,
    )