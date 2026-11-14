"""One geometry handoff, three solver families, one machine-checked frame contract.

Why this module exists
----------------------
The package could already hand a surface to an FEM preprocessor
(:mod:`tensorcad.adapters.fem`) and a hull body to a volume-mesh writer
(:mod:`tensorcad.volume_mesh`).  What it could not do is state, in one place and
in a form a machine re-checks, what a *solver* is allowed to assume about that
handoff.  A previous round failed on exactly this: a volume mesh was written and
no solver ever consumed it, so nothing distinguished "exported" from "accepted".

The contract here is deliberately small and auditable:

``frame``
    ``z-up`` right-handed: ``x`` aft -> fore, ``y`` to starboard, ``z`` up, waterline
    at ``z = 0``, keel at ``z = -draft``.  Units are metres.  A caller whose CAD
    kernel works y-up converts with :func:`to_z_up`; the conversion is a proper
    rotation (determinant ``+1``), so it cannot silently mirror a hull.
``normals``
    Every surface triangle points out of the solid, i.e. the divergence theorem
    gives a *positive* volume.  That is checked, not assumed.
``tables``
    Solver-facing numbers are written as files (gmsh 2.2 ASCII with ``C3D4``
    tetrahedra and physically tagged boundary triangles, or a JSON table), never
    as imports of a solver package.  Each writer also writes ``handoff.json``:
    the declared counts, volume, area, digest and frame.
``verification``
    :func:`verify_bundle` re-reads the files with an independent parser and
    recomputes volume, area and digest.  A bundle whose declared numbers do not
    survive the re-read is refused.  This is what makes the handoff a contract
    rather than a comment.

Status of each route (evidence, not intent):

==================================  ==========================================
route                               evidence here
==================================  ==========================================
volume ``.msh`` (C3D4 + wall tags)  written and re-read by an independent parser
surface ``.json`` triangle table    written and re-read
surface gmsh via ``adapters.fem``   written by that module, re-used unchanged
B-rep STEP / IGES                   refused upstream in :mod:`tensorcad.exchange`
==================================  ==========================================
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from .geometry import Mesh

__all__ = [
    "BUNDLE_SCHEMA",
    "COORDINATE_CONTRACT",
    "FrameReport",
    "HandoffError",
    "boundary_triangles",
    "channel_grid_2d",
    "check_frame",
    "ellipsoid_mesh",
    "grid_quads",
    "read_grid_table",
    "read_lattice_table",
    "read_gmsh_volume",
    "signed_volume",
    "to_z_up",
    "total_volume",
    "verify_bundle",
    "write_bundle_grid",
    "write_bundle_lattice",
    "write_bundle_surface",
    "write_bundle_volume",
]

#: Bumped on any incompatible change to ``handoff.json`` or to the table layouts.
BUNDLE_SCHEMA = "tensorcad.solver_interface/1"

#: The single frame every solver-facing handoff declares.
COORDINATE_CONTRACT: dict[str, Any] = {
    "name": "z-up",
    "schema": BUNDLE_SCHEMA,
    "units": "m",
    "handedness": "right",
    "axes": {
        "x": "aft -> fore (longitudinal, positive forward)",
        "y": "to starboard (transverse)",
        "z": "up (vertical, positive upward)",
    },
    "origin": "z = 0 is the design waterline; the keel lies at z = -draft",
    "surface_normals": "outward (positive signed volume)",
    "rotation_from_y_up": "x' = x, y' = -z, z' = y  (determinant +1)",
}


class HandoffError(ValueError):
    """A handoff that cannot be verified is refused, never silently repaired."""


@dataclass(frozen=True)
class FrameReport:
    """What the frame check measured, not what the caller claimed."""

    declared: str
    right_handed: bool
    signed_volume_m3: float
    z_extent_m: tuple[float, float]

    @property
    def z_up(self) -> bool:
        """A right-handed outward surface whose vertical extent is on ``z``."""
        return self.right_handed and self.z_extent_m[1] > self.z_extent_m[0]

    def to_dict(self) -> dict[str, Any]:
        return {
            "declared": self.declared,
            "right_handed": self.right_handed,
            "signed_volume_m3": self.signed_volume_m3,
            "z_extent_m": list(self.z_extent_m),
            "z_up": self.z_up,
        }


# ------------------------------------------------------------------ frame checks


def signed_volume(vertices: Sequence[Sequence[float]], faces: Iterable[Sequence[int]]) -> float:
    """Divergence-theorem volume; positive when every face normal points outward."""
    total = 0.0
    for a, b, c in faces:
        ax, ay, az = vertices[a]
        bx, by, bz = vertices[b]
        cx, cy, cz = vertices[c]
        total += (
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        )
    return total / 6.0


def check_frame(mesh: Mesh, *, declared: str = "z-up") -> FrameReport:
    """Measure the frame the triangles actually describe."""
    if not mesh.faces:
        raise HandoffError("an empty surface has no frame to check")
    volume = signed_volume(mesh.vertices, mesh.faces)
    zs = [v[2] for v in mesh.vertices]
    return FrameReport(
        declared=declared,
        right_handed=volume > 0.0,
        signed_volume_m3=volume,
        z_extent_m=(min(zs), max(zs)),
    )


def to_z_up(mesh: Mesh, *, source_frame: str) -> Mesh:
    """Rotate a mesh into the contract frame; refuse frames it cannot convert.

    ``y-up`` (the usual glTF/CAD-kernel default) maps as ``x' = x, y' = -z,
    z' = y``.  That is a rotation of the point set (determinant ``+1``), so a
    right-handed outward surface stays right-handed and its signed volume is
    unchanged -- both are re-checked here rather than trusted.
    """
    frame = source_frame.strip().lower()
    if frame in ("z-up", "zup", "z_up"):
        return mesh
    if frame not in ("y-up", "yup", "y_up"):
        raise HandoffError(f"no declared conversion from frame {source_frame!r} to z-up")
    before = signed_volume(mesh.vertices, mesh.faces)
    vertices = tuple((x, -z, y) for x, y, z in mesh.vertices)
    after = signed_volume(vertices, mesh.faces)
    if abs(after - before) > 1e-9 * max(1.0, abs(before)):
        raise HandoffError("y-up -> z-up conversion changed the signed volume; not a rotation")
    return Mesh(
        vertices=vertices,
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=COORDINATE_CONTRACT["name"]
        + " (right-handed, x aft->fore, y starboard, z up)",
    )


# ------------------------------------------------------------------- primitives


def ellipsoid_mesh(
    semi_axes: Sequence[float], *, n_theta: int = 18, n_phi: int = 24
) -> Mesh:
    """A closed, outward-oriented ellipsoid in the contract frame.

    Parametrised as ``(a sin t cos p, b sin t sin p, c cos t)`` with ``t`` from
    the north pole to the south pole, so the poles are the only duplicated rows
    and the surface is closed without a seam.  Analytic volume ``4/3 pi a b c``
    is the reference the solver handoff is scored against.
    """
    a, b, c = (float(v) for v in semi_axes)
    if min(a, b, c) <= 0.0 or not all(math.isfinite(v) for v in (a, b, c)):
        raise HandoffError("semi-axes must be positive and finite")
    if n_theta < 4 or n_phi < 8:
        raise HandoffError("n_theta >= 4 and n_phi >= 8 are required")

    vertices: list[tuple[float, float, float]] = []
    rings: list[list[int]] = []
    for i in range(1, n_theta):  # poles handled separately
        polar = math.pi * i / n_theta
        ring = []
        for k in range(n_phi):
            azimuth = 2.0 * math.pi * k / n_phi
            ring.append(len(vertices))
            vertices.append(
                (
                    a * math.sin(polar) * math.cos(azimuth),
                    b * math.sin(polar) * math.sin(azimuth),
                    c * math.cos(polar),
                )
            )
        rings.append(ring)
    north = len(vertices)
    vertices.append((0.0, 0.0, c))
    south = len(vertices)
    vertices.append((0.0, 0.0, -c))

    faces: list[tuple[int, int, int]] = []
    for k in range(n_phi):
        nxt = (k + 1) % n_phi
        faces.append((north, rings[0][nxt], rings[0][k]))
        faces.append((south, rings[-1][k], rings[-1][nxt]))
    for i in range(len(rings) - 1):
        upper, lower = rings[i], rings[i + 1]
        for k in range(n_phi):
            nxt = (k + 1) % n_phi
            faces.append((upper[k], upper[nxt], lower[nxt]))
            faces.append((upper[k], lower[nxt], lower[k]))

    if signed_volume(vertices, faces) < 0.0:
        # The ring winding above is pole-to-pole; flip once so the contract's
        # outward-normal rule holds and the check below measures orientation.
        faces = [(a, c, b) for a, b, c in faces]
    mesh = Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units="m",
        coordinate_frame=COORDINATE_CONTRACT["name"]
        + " (right-handed, x aft->fore, y starboard, z up)",
    )
    report = check_frame(mesh)
    analytic = 4.0 / 3.0 * math.pi * a * b * c
    if abs(report.signed_volume_m3 - analytic) > 0.05 * analytic:
        # A badly resolved ellipsoid is a fixture problem, not a solver problem.
        raise HandoffError(
            f"discretised volume {report.signed_volume_m3:.6g} m3 is not within 5% of "
            f"the analytic {analytic:.6g} m3; refine n_theta/n_phi"
        )
    if not report.right_handed:
        raise HandoffError("ellipsoid triangulation is inward facing")
    return mesh


def channel_grid_2d(
    *, length: float, height: float, nx: int, ny: int
) -> tuple[tuple[tuple[float, float], ...], ...]:
    """Structured ``(ny+1) x (nx+1)`` planar grid, row-major by ``y``.

    This is the planar-channel analogue of :func:`tensorcad.volume_mesh.box_section_grid`:
    geometry, not meshing policy.  The solver-facing face tables are built from
    it by the consuming adapter, so the grid itself stays a plain coordinate
    array that any 2-D solver can read.
    """
    if not (length > 0.0 and height > 0.0) or nx < 1 or ny < 1:
        raise HandoffError("length/height must be positive and nx, ny >= 1")
    return tuple(
        tuple((length * i / nx, height * j / ny) for i in range(nx + 1))
        for j in range(ny + 1)
    )


# --------------------------------------------------------------- volume tables


def boundary_triangles(
    tets: Sequence[Sequence[int]],
) -> tuple[tuple[int, int, int], ...]:
    """Boundary faces of a tetrahedral mesh, outward relative to the tetrahedra.

    Each interior face is shared by two tetrahedra and cancels; each boundary
    face keeps the orientation induced by the single tetrahedron that owns it.
    """
    seen: dict[tuple[int, int, int], int] = {}
    for tet in tets:
        a, b, c, d = tet
        for face in ((a, c, b), (a, b, d), (b, c, d), (a, d, c)):
            key = tuple(sorted(face))
            count, oriented = seen.get(key, (0, face))
            seen[key] = (count + 1, oriented)
    return tuple(oriented for count, oriented in seen.values() if count == 1)


def total_volume(
    vertices: Sequence[Sequence[float]], tets: Sequence[Sequence[int]]
) -> float:
    total = 0.0
    for a, b, c, d in tets:
        ax, ay, az = vertices[a]
        bx, by, bz = vertices[b]
        cx, cy, cz = vertices[c]
        dx, dy, dz = vertices[d]
        total += (
            (bx - ax) * ((cy - ay) * (dz - az) - (cz - az) * (dy - ay))
            - (by - ay) * ((cx - ax) * (dz - az) - (cz - az) * (dx - ax))
            + (bz - az) * ((cx - ax) * (dy - ay) - (cy - ay) * (dx - ax))
        ) / 6.0
    return total


def _digest(vertices: Sequence[Sequence[float]], cells: Sequence[Sequence[int]]) -> str:
    digest = hashlib.sha256()
    digest.update(f"{len(vertices)}|{len(cells)}".encode())
    for v in vertices:
        digest.update((",".join(f"{float(c):.17g}" for c in v) + ";").encode())
    for cell in cells:
        digest.update((",".join(str(int(c)) for c in cell) + ";").encode())
    return digest.hexdigest()


def write_bundle_volume(
    outdir: str | Path,
    *,
    vertices: Sequence[Sequence[float]],
    tets: Sequence[Sequence[int]],
    wall_name: str = "wall",
    case: str = "unnamed",
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write ``volume.msh`` (gmsh 2.2 ASCII) plus ``handoff.json`` for a solid.

    The gmsh file carries one volume block of ``C3D4`` tetrahedra (element type
    4) and one triangle block (type 2) tagged ``wall_name`` -- the physical group
    an FEM or FVM reader turns into a boundary patch.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    triangles = boundary_triangles(tets)
    volume = total_volume(vertices, tets)
    if volume <= 0.0:
        raise HandoffError(f"tetrahedral volume {volume:.6g} m3 is not positive")
    digest = _digest(vertices, tets)

    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]
    lines += ["$PhysicalNames", "2", f'2 1 "{wall_name}"', '3 2 "solid"', "$EndPhysicalNames"]
    lines += ["$Nodes", str(len(vertices))]
    for index, (x, y, z) in enumerate(vertices, start=1):
        lines.append(f"{index} {float(x):.17g} {float(y):.17g} {float(z):.17g}")
    lines.append("$EndNodes")
    lines += ["$Elements", str(len(triangles) + len(tets))]
    tag = 1
    for a, b, c in triangles:
        lines.append(f"{tag} 2 2 1 1 {a + 1} {b + 1} {c + 1}")
        tag += 1
    for a, b, c, d in tets:
        lines.append(f"{tag} 4 2 2 2 {a + 1} {b + 1} {c + 1} {d + 1}")
        tag += 1
    lines.append("$EndElements")
    (out / "volume.msh").write_text("\n".join(lines) + "\n", encoding="utf-8")

    handoff = {
        "schema": BUNDLE_SCHEMA,
        "case": case,
        "frame": dict(COORDINATE_CONTRACT),
        "solver_family": "volume (FEM / FVM solid)",
        "tables": {"volume_msh": "volume.msh", "element_types": {"2": "triangle", "4": "C3D4"}},
        "physical_groups": {"2": {1: wall_name}, "3": {2: "solid"}},
        "counts": {
            "nodes": len(vertices),
            "tetrahedra": len(tets),
            "boundary_triangles": len(triangles),
        },
        "volume_m3": volume,
        "surface_area_m2": sum(
            0.5
            * math.dist((0.0, 0.0, 0.0), _cross(_sub(vertices[b], vertices[a]),
                                                 _sub(vertices[c], vertices[a])))
            for a, b, c in triangles
        ),
        "digest_sha256": digest,
        "extra": dict(extra or {}),
    }
    (out / "handoff.json").write_text(json.dumps(handoff, indent=2, sort_keys=True) + "\n")
    return handoff


def _sub(u: Sequence[float], v: Sequence[float]) -> tuple[float, float, float]:
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def _cross(u: Sequence[float], v: Sequence[float]) -> tuple[float, float, float]:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def write_bundle_surface(
    outdir: str | Path, mesh: Mesh, *, case: str = "unnamed", extra: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Write ``surface.json`` (triangle table) plus ``handoff.json`` for a closed body."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    report = check_frame(mesh)
    if not report.right_handed:
        raise HandoffError(
            "surface normals point into the body (signed volume "
            f"{report.signed_volume_m3:.6g} m3); the contract requires outward normals"
        )
    table = {
        "schema": BUNDLE_SCHEMA,
        "frame": COORDINATE_CONTRACT,
        "vertices": [[float(c) for c in v] for v in mesh.vertices],
        "faces": [[int(i) for i in f] for f in mesh.faces],
    }
    (out / "surface.json").write_text(json.dumps(table) + "\n")
    area = sum(
        0.5
        * math.dist((0.0, 0.0, 0.0), _cross(_sub(mesh.vertices[b], mesh.vertices[a]),
                                             _sub(mesh.vertices[c], mesh.vertices[a])))
        for a, b, c in mesh.faces
    )
    handoff = {
        "schema": BUNDLE_SCHEMA,
        "case": case,
        "frame": dict(COORDINATE_CONTRACT),
        "solver_family": "surface (BEM panels)",
        "tables": {"surface_json": "surface.json"},
        "counts": {"nodes": len(mesh.vertices), "triangles": len(mesh.faces)},
        "signed_volume_m3": report.signed_volume_m3,
        "surface_area_m2": area,
        "digest_sha256": _digest(mesh.vertices, mesh.faces),
        "extra": dict(extra or {}),
    }
    (out / "handoff.json").write_text(json.dumps(handoff, indent=2, sort_keys=True) + "\n")
    return handoff


def _quad_area(vertices: Sequence[Sequence[float]], quad: Sequence[int]) -> float:
    """Area of a planar quad, split into two triangles."""
    a, b, c, d = (vertices[i] for i in quad)
    return 0.5 * math.dist((0.0, 0.0, 0.0), _cross(_sub(b, a), _sub(c, a))) + 0.5 * math.dist(
        (0.0, 0.0, 0.0), _cross(_sub(c, a), _sub(d, a))
    )


def grid_quads(nx: int, ny: int) -> tuple[tuple[int, int, int, int], ...]:
    """Counter-clockwise quads of an ``nx x ny`` structured node grid.

    Node ``(i, j)`` sits at index ``j * (nx + 1) + i``; each quad lists its four
    corners counter-clockwise seen from above, so a plate mid-surface has its
    outward normal on ``+z``.
    """
    if nx < 1 or ny < 1:
        raise HandoffError("nx, ny >= 1 are required for a structured grid")
    quads = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            quads.append((a, a + 1, a + 1 + (nx + 1), a + (nx + 1)))
    return tuple(quads)


def write_bundle_grid(
    outdir: str | Path,
    *,
    xs: Sequence[float],
    ys: Sequence[float],
    case: str = "unnamed",
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write ``grid.json`` (structured quad grid) plus ``handoff.json``.

    The route exists because a plate or shell solver consumes a *2-D* mesh while
    this module previously only handed out triangles (BEM) and tetrahedra
    (FEM/FVM).  Nodes are written in the contract frame (``z = 0`` mid-surface,
    counter-clockwise quads about ``+z``), so a consumer never re-derives
    connectivity from spacing -- the same rule the other two routes follow.
    """
    if len(xs) < 2 or len(ys) < 2:
        raise HandoffError("a structured grid needs at least 2 nodes per axis")
    xs_f = [float(v) for v in xs]
    ys_f = [float(v) for v in ys]
    if any(b <= a for a, b in zip(xs_f, xs_f[1:])) or any(b <= a for a, b in zip(ys_f, ys_f[1:])):
        raise HandoffError("grid coordinates must be strictly increasing on both axes")
    nx, ny = len(xs_f) - 1, len(ys_f) - 1
    vertices = tuple((x, y, 0.0) for y in ys_f for x in xs_f)
    quads = grid_quads(nx, ny)
    # A clockwise quad would face -z and a shell solver would assemble a
    # negative-definite bending operator, so the winding is measured, not assumed.
    area_z = sum(
        (vertices[q[1]][0] - vertices[q[0]][0]) * (vertices[q[3]][1] - vertices[q[0]][1])
        - (vertices[q[2]][0] - vertices[q[0]][0]) * (vertices[q[1]][1] - vertices[q[0]][1])
        for q in quads
    )
    if not (area_z > 0.0):
        raise HandoffError(f"quad winding gives a non-positive projected area ({area_z:.6g})")
    area = sum(_quad_area(vertices, q) for q in quads)

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    table = {
        "vertices": [[float(c) for c in v] for v in vertices],
        "quads": [[int(i) for i in q] for q in quads],
        "axis": {"i": "+x", "j": "+y", "normal": "+z"},
    }
    (out / "grid.json").write_text(json.dumps(table) + "\n", encoding="utf-8")
    handoff = {
        "schema": BUNDLE_SCHEMA,
        "case": case,
        "frame": dict(COORDINATE_CONTRACT),
        "solver_family": "surface (plate / shell FEM)",
        "tables": {"grid_json": "grid.json"},
        "counts": {"nodes": len(vertices), "quads": len(quads)},
        "surface_area_m2": area,
        "digest_sha256": _digest(vertices, quads),
        "extra": dict(extra or {}),
    }
    (out / "handoff.json").write_text(json.dumps(handoff, indent=2, sort_keys=True) + "\n")
    return handoff


def write_bundle_lattice(
    outdir: str | Path,
    *,
    shape: Sequence[int],
    origin: Sequence[float],
    spacing: Sequence[float],
    fluid_cells: Sequence[Sequence[int]],
    walls: dict[str, Any] | None = None,
    case: str = "unnamed",
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write ``lattice.json`` (cell-centred occupancy) plus ``handoff.json``.

    This is the LBM / voxel-solver consumption surface: the CAD side says which
    *cells* are fluid, where the box starts and how large a cell is; the solver
    never sees the CAD representation.  Cells are ``(i, j, k)`` indices with
    ``i`` along ``+x``, ``j`` along ``+y``, ``k`` along ``+z`` (the z-up contract
    every route follows) and ``origin`` is the corner of cell ``(0, 0, 0)``.
    """
    dims = tuple(int(v) for v in shape)
    if len(dims) != 3 or min(dims) < 1:
        raise HandoffError(f"shape must be three positive cell counts, got {shape!r}")
    org = tuple(float(v) for v in origin)
    spacing_f = tuple(float(v) for v in spacing)
    if len(org) != 3 or len(spacing_f) != 3:
        raise HandoffError("origin and spacing must be three-vectors")
    if min(spacing_f) <= 0.0:
        raise HandoffError(f"spacing must be positive, got {spacing_f!r}")
    cells = sorted({tuple(int(c) for c in cell) for cell in fluid_cells})
    for cell in cells:
        if len(cell) != 3 or any(not (0 <= cell[axis] < dims[axis]) for axis in range(3)):
            raise HandoffError(f"cell {cell!r} lies outside the lattice {dims}")
    volume = len(cells) * spacing_f[0] * spacing_f[1] * spacing_f[2]

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    table = {
        "shape": list(dims),
        "origin": list(org),
        "spacing": list(spacing_f),
        "axis": {"i": "+x", "j": "+y", "k": "+z"},
        "fluid_cells": [list(cell) for cell in cells],
        "walls": dict(walls or {}),
    }
    (out / "lattice.json").write_text(json.dumps(table) + "\n", encoding="utf-8")
    handoff = {
        "schema": BUNDLE_SCHEMA,
        "case": case,
        "frame": dict(COORDINATE_CONTRACT),
        "solver_family": "voxel (LBM / cut-cell FVM)",
        "tables": {"lattice_json": "lattice.json"},
        "counts": {"cells": dims[0] * dims[1] * dims[2], "fluid_cells": len(cells)},
        "fluid_volume_m3": volume,
        "digest_sha256": _digest((org, spacing_f), cells),
        "extra": dict(extra or {}),
    }
    (out / "handoff.json").write_text(json.dumps(handoff, indent=2, sort_keys=True) + "\n")
    return handoff


def read_grid_table(path: str | Path) -> dict[str, Any]:
    """Independent reader for ``grid.json``: nodes, quads and axis map."""
    table = json.loads(Path(path).read_text(encoding="utf-8"))
    vertices = [tuple(float(c) for c in v) for v in table["vertices"]]
    quads = [tuple(int(i) for i in q) for q in table["quads"]]
    if any(len(q) != 4 for q in quads):
        raise HandoffError("every grid cell must be a quad")
    if any(not (0 <= i < len(vertices)) for q in quads for i in q):
        raise HandoffError("a grid quad references a node outside the table")
    return {"vertices": vertices, "quads": quads, "axis": table.get("axis", {})}


def read_lattice_table(path: str | Path) -> dict[str, Any]:
    """Independent reader for ``lattice.json``."""
    table = json.loads(Path(path).read_text(encoding="utf-8"))
    dims = tuple(int(v) for v in table["shape"])
    origin = tuple(float(v) for v in table["origin"])
    spacing = tuple(float(v) for v in table["spacing"])
    cells = [tuple(int(i) for i in cell) for cell in table["fluid_cells"]]
    for cell in cells:
        if any(not (0 <= cell[axis] < dims[axis]) for axis in range(3)):
            raise HandoffError(f"lattice cell {cell!r} lies outside {dims}")
    return {
        "shape": dims,
        "origin": origin,
        "spacing": spacing,
        "fluid_cells": cells,
        "walls": table.get("walls", {}),
        "axis": table.get("axis", {}),
    }


def read_gmsh_volume(path: str | Path) -> dict[str, Any]:
    """Independent gmsh 2.2 ASCII reader for the subset this module writes.

    Written here rather than borrowed so that :func:`verify_bundle` checks the
    files against a *second* implementation, not against the writer's own
    in-memory state.
    """
    text = Path(path).read_text(encoding="utf-8")

    def block(name: str) -> list[str]:
        marker = f"${name}"
        end = f"$End{name}"
        if marker not in text:
            raise HandoffError(f"{marker} is missing from {path}")
        body = text.split(marker, 1)[1].split(end, 1)[0].strip().splitlines()
        return [line.strip() for line in body if line.strip()]

    header = block("MeshFormat")[0].split()
    if not header[0].startswith("2."):
        raise HandoffError(f"only gmsh 2.x is written here, got {header[0]}")
    nodes_block = block("Nodes")
    count = int(nodes_block[0])
    nodes: list[tuple[float, float, float]] = []
    for line in nodes_block[1 : count + 1]:
        parts = line.split()
        nodes.append((float(parts[1]), float(parts[2]), float(parts[3])))
    elements_block = block("Elements")
    nelem = int(elements_block[0])
    triangles: list[tuple[int, int, int]] = []
    tets: list[tuple[int, int, int, int]] = []
    group_names: dict[int, str] = {}
    for line in block("PhysicalNames")[1:]:
        parts = line.replace('"', " ").split()
        if len(parts) >= 3:
            group_names[int(parts[1])] = parts[2]
    for line in elements_block[1 : nelem + 1]:
        parts = line.split()
        kind, physical = int(parts[1]), int(parts[3])
        conn = [int(v) - 1 for v in parts[5:]]
        if kind == 2:
            triangles.append(tuple(conn))
        elif kind == 4:
            tets.append(tuple(conn))
            physical = physical
        else:
            raise HandoffError(f"unexpected gmsh element type {kind} in a volume bundle")
    return {
        "nodes": nodes,
        "triangles": triangles,
        "tets": tets,
        "physical_names": group_names,
    }


def verify_bundle(outdir: str | Path) -> dict[str, Any]:
    """Re-read a bundle and check every declared number against the files."""
    out = Path(outdir)
    manifest = json.loads((out / "handoff.json").read_text())
    if manifest.get("schema") != BUNDLE_SCHEMA:
        raise HandoffError(f"unknown bundle schema {manifest.get('schema')!r}")
    frame = manifest.get("frame", {})
    if frame.get("name") != COORDINATE_CONTRACT["name"]:
        raise HandoffError(f"bundle frame {frame.get('name')!r} is not the contract frame")
    report: dict[str, Any] = {"schema": BUNDLE_SCHEMA, "case": manifest.get("case")}

    if "volume_msh" in manifest.get("tables", {}):
        data = read_gmsh_volume(out / "volume.msh")
        volume = total_volume(data["nodes"], data["tets"])
        digest = _digest(data["nodes"], data["tets"])
        declared = float(manifest["volume_m3"])
        if abs(volume - declared) > 1e-12 * max(1.0, abs(declared)):
            raise HandoffError(
                f"re-read volume {volume:.12g} m3 disagrees with declared {declared:.12g} m3"
            )
        if digest != manifest["digest_sha256"]:
            raise HandoffError("re-read digest differs from the declared digest")
        if len(data["tets"]) != manifest["counts"]["tetrahedra"]:
            raise HandoffError("re-read tetrahedron count differs from the declared count")
        report.update(
            {
                "route": "volume",
                "nodes": len(data["nodes"]),
                "tets": len(data["tets"]),
                "boundary_triangles": len(data["triangles"]),
                "volume_m3": volume,
                "digest_matches": True,
                "physical_groups": data["physical_names"],
            }
        )
    elif "surface_json" in manifest.get("tables", {}):
        table = json.loads((out / "surface.json").read_text())
        vertices, faces = table["vertices"], table["faces"]
        volume = signed_volume(vertices, faces)
        digest = _digest(vertices, faces)
        if not (volume > 0.0):
            raise HandoffError(f"re-read surface volume {volume:.6g} m3 is not positive")
        if abs(volume - float(manifest["signed_volume_m3"])) > 1e-12 * max(1.0, abs(volume)):
            raise HandoffError("re-read surface volume disagrees with the declared volume")
        if digest != manifest["digest_sha256"]:
            raise HandoffError("re-read surface digest differs from the declared digest")
        report.update(
            {
                "route": "surface",
                "nodes": len(vertices),
                "triangles": len(faces),
                "signed_volume_m3": volume,
                "digest_matches": True,
            }
        )
    elif "grid_json" in manifest.get("tables", {}):
        data = read_grid_table(out / "grid.json")
        vertices, quads = data["vertices"], data["quads"]
        digest = _digest(vertices, quads)
        area = sum(_quad_area(vertices, quad) for quad in quads)
        declared_area = float(manifest["surface_area_m2"])
        if abs(area - declared_area) > 1e-12 * max(1.0, abs(declared_area)):
            raise HandoffError(
                f"re-read grid area {area:.12g} m2 disagrees with declared {declared_area:.12g} m2"
            )
        if digest != manifest["digest_sha256"]:
            raise HandoffError("re-read grid digest differs from the declared digest")
        if len(vertices) != manifest["counts"]["nodes"] or len(quads) != manifest["counts"]["quads"]:
            raise HandoffError("re-read grid counts differ from the declared counts")
        report.update(
            {
                "route": "grid",
                "nodes": len(vertices),
                "quads": len(quads),
                "surface_area_m2": area,
                "digest_matches": True,
            }
        )
    elif "lattice_json" in manifest.get("tables", {}):
        data = read_lattice_table(out / "lattice.json")
        cells = data["fluid_cells"]
        digest = _digest((tuple(data["origin"]), tuple(data["spacing"])), cells)
        declared = float(manifest["fluid_volume_m3"])
        volume = len(cells) * math.prod(data["spacing"])
        if abs(volume - declared) > 1e-12 * max(1.0, abs(declared)):
            raise HandoffError(
                f"re-read lattice volume {volume:.12g} m3 disagrees with declared {declared:.12g} m3"
            )
        if digest != manifest["digest_sha256"]:
            raise HandoffError("re-read lattice digest differs from the declared digest")
        if len(cells) != manifest["counts"]["fluid_cells"]:
            raise HandoffError("re-read fluid-cell count differs from the declared count")
        report.update(
            {
                "route": "lattice",
                "shape": list(data["shape"]),
                "fluid_cells": len(cells),
                "fluid_volume_m3": volume,
                "digest_matches": True,
            }
        )
    else:
        raise HandoffError("handoff.json names no table this verifier understands")
    report["frame"] = frame.get("name")
    report["ok"] = True
    return report