"""Solver-neutral surface exchange for TensorFEM (and any other FEM code).

Where this sits: the core hands out a :class:`~tensorcad.geometry.Mesh` plus a
``tensorcad.geometry.handoff/1`` record; that record says what the surface *is*
(units, frame, digest, topology, measured volume, quality warnings) but it does not
carry the node and element tables a finite-element code needs. This module turns the
same surface into those tables and writes them in two formats:

* ``.msh`` (gmsh 2.2 ASCII) -- triangles with a *physical* tag per region, which is
  how an FEM preprocessor receives boundary patches (hull, sail, fins) as named
  groups rather than as one anonymous soup. Gmsh is a text format, so the writer and
  a reader live here in the standard library and a test can check the round trip
  without any extra install.
* ``.json`` -- the same tables plus the handoff manifest, byte-exact floats, for
  pipelines that would rather not parse a mesh format.

Nothing here imports TensorFEM, or OCP, or gmsh: the contract is the file, not a
private module import. That mirrors the ecosystem rule the other Tensor* projects
follow ("explicit coupling adapters, not private module imports"), and it keeps this
module usable in the standard-library-only core lane.

Conventions are declared rather than implied, because a surface handed to a solid
mechanics code is only meaningful with them:

``units``
    Metres. Both writers emit the numbers exactly as measured; gmsh has no unit
    concept, so the unit is recorded in the ``$Comments`` block and in the manifest.
``coordinate_frame``
    Copied from the mesh (the models here are x-aft / y-starboard / z-up).
``sign_convention``
    ``"outward"``: every triangle normal points out of the solid, so the outward
    pressure of a fluid solve is a positive normal traction. This is checked, not
    assumed -- :func:`fem_surface` measures the signed volume and refuses a surface
    whose orientation makes it negative.
``time_s``
    The instant the surface belongs to; static geometry is ``0.0``. A coupled run
    moves this number, nothing else.
``interpolation`` / ``integration``
    ``"P1 (linear) triangles"`` and ``"3-point symmetric Gauss per triangle"``, so a
    solver that needs a different quadrature knows what the exported normals and
    areas were computed with.
"""
from __future__ import annotations

import hashlib
import json
import math
import struct
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..geometry import Mesh
from ..geometry_field import GeometryField

__all__ = [
    "FEM_SCHEMA",
    "FemSurface",
    "fem_surface",
    "label_regions",
    "read_fem_json",
    "read_gmsh",
    "write_fem_json",
    "write_gmsh",
]

#: Schema tag written into every manifest; bump on incompatible change.
FEM_SCHEMA = "tensorcad.fem.surface/1"

#: gmsh element type for a 3-node triangle.
_GMSH_TRIANGLE = 2


@dataclass(frozen=True)
class RegionStats:
    """What one named region of the surface turned out to contain."""

    name: str
    physical_id: int
    triangles: int
    area_m2: float
    max_assignment_distance_m: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "physical_id": self.physical_id,
            "triangles": self.triangles,
            "area_m2": self.area_m2,
            "max_assignment_distance_m": self.max_assignment_distance_m,
        }


@dataclass(frozen=True)
class FemSurface:
    """Node/element tables for a boundary surface, with its declared conventions."""

    nodes: tuple[tuple[float, float, float], ...]
    triangles: tuple[tuple[int, int, int], ...]
    region_ids: tuple[int, ...]
    region_names: tuple[str, ...]
    regions: tuple[RegionStats, ...] = ()
    units: str = "m"
    coordinate_frame: str = "unspecified right-handed frame"
    time_s: float = 0.0
    sign_convention: str = "outward"
    interpolation: str = "P1 (linear) triangles"
    integration: str = "3-point symmetric Gauss per triangle"
    volume_m3: float = 0.0
    area_m2: float = 0.0
    watertight: bool = False
    edge_manifold: bool = False
    euler_characteristic: int = 0
    source_digest: str = ""
    quality: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "nodes", tuple(tuple(float(c) for c in n) for n in self.nodes))
        object.__setattr__(self, "triangles", tuple(tuple(int(i) for i in t) for t in self.triangles))
        object.__setattr__(self, "region_ids", tuple(int(r) for r in self.region_ids))
        object.__setattr__(self, "region_names", tuple(str(n) for n in self.region_names))
        object.__setattr__(self, "quality", dict(self.quality))
        if len(self.region_ids) != len(self.triangles):
            raise ValueError(
                f"{len(self.region_ids)} region tags for {len(self.triangles)} triangles"
            )
        for index, (a, b, c) in enumerate(self.triangles):
            for vertex in (a, b, c):
                if not 0 <= vertex < len(self.nodes):
                    raise ValueError(f"triangle {index} references out-of-range node {vertex}")
        if self.triangles:
            unknown = sorted(set(self.region_ids) - set(range(1, len(self.region_names) + 1)))
            if unknown:
                raise ValueError(f"region tags {unknown} have no name")

    @property
    def digest(self) -> str:
        """SHA-256 over nodes and elements, identical to the mesh digest convention."""
        digest = hashlib.sha256()
        digest.update(struct.pack("<II", len(self.nodes), len(self.triangles)))
        for node in self.nodes:
            digest.update(struct.pack("<3d", *node))
        for face in self.triangles:
            digest.update(struct.pack("<3I", *face))
        return digest.hexdigest()

    @property
    def region_id_of(self) -> dict[str, int]:
        return {name: index + 1 for index, name in enumerate(self.region_names)}

    def to_dict(self, *, include_tables: bool = False) -> dict[str, Any]:
        """The manifest: contract first, arrays only on request.

        The default is the small record a solver logs next to a run; ``include_tables``
        adds the coordinates and connectivity for JSON-only pipelines.
        """
        data: dict[str, Any] = {
            "schema": FEM_SCHEMA,
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "sign_convention": self.sign_convention,
            "time_s": self.time_s,
            "interpolation": self.interpolation,
            "integration": self.integration,
            "nodes": len(self.nodes),
            "triangles": len(self.triangles),
            "index_base": 0,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "euler_characteristic": self.euler_characteristic,
            "digest": self.digest,
            "mesh_sha256": self.source_digest,
            "signed_volume_m3": self.volume_m3,
            "surface_area_m2": self.area_m2,
            "regions": [stats.to_dict() for stats in self.regions],
            "quality": dict(self.quality),
        }
        if include_tables:
            data["node_table"] = [list(node) for node in self.nodes]
            data["element_table"] = [list(face) for face in self.triangles]
            data["region_table"] = list(self.region_ids)
        return data


def _triangle_areas_and_centroids(
    mesh: Mesh,
) -> tuple[list[float], list[tuple[float, float, float]]]:
    """Per-triangle area and centroid: the two quantities the handoff needs."""
    areas: list[float] = []
    centroids: list[tuple[float, float, float]] = []
    vertices = mesh.vertices
    for a, b, c in mesh.faces:
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        areas.append(0.5 * math.sqrt(cx * cx + cy * cy + cz * cz))
        centroids.append(
            ((pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0, (pa[2] + pb[2] + pc[2]) / 3.0)
        )
    return areas, centroids


def label_regions(
    mesh: Mesh,
    components: Mapping[str, GeometryField] | Sequence[GeometryField],
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[RegionStats, ...]]:
    """Assign every triangle to the analytic component whose surface it lies on.

    This is how a soup of triangles becomes boundary patches: for the SUBOFF full
    model the three fields (hull, sail, fins) are available analytically, so the
    triangle's region is the one whose signed distance is closest to zero at the
    triangle centroid. That is a geometric decision with a measurable error --
    :attr:`RegionStats.max_assignment_distance_m` reports the worst centroid distance
    to its assigned surface, so a reader can see whether the labelling is meaningful
    or whether the surface has drifted between two components.

    Requires NumPy (batch field evaluation) but nothing else. Triangle centroids are
    the reference point, not vertices or integration points: they are the interior
    point of the facet, so a triangle is never assigned because one of its corners
    happens to coincide with another component.
    """
    import numpy as np

    if isinstance(components, Mapping):
        named = list(components.items())
    else:
        named = [(f"region_{index + 1}", component) for index, component in enumerate(components)]
    if not named:
        raise ValueError("label_regions needs at least one component field")
    for name, component in named:
        if not isinstance(component, GeometryField):
            raise TypeError(f"component {name!r} is not a GeometryField")

    areas, centroids = _triangle_areas_and_centroids(mesh)
    if not centroids:
        return (), tuple(name for name, _ in named), ()
    points = np.asarray(centroids, dtype=float)
    distances = np.stack(
        [np.asarray(component.signed_distance(points), dtype=float) for _, component in named]
    )
    absolute = np.abs(distances)
    nearest = np.argmin(absolute, axis=0)
    ids = tuple(int(index) + 1 for index in nearest)
    names = tuple(name for name, _ in named)
    area_array = np.asarray(areas, dtype=float)
    stats: list[RegionStats] = []
    for index, name in enumerate(names):
        selected = nearest == index
        count = int(selected.sum())
        stats.append(
            RegionStats(
                name=name,
                physical_id=index + 1,
                triangles=count,
                area_m2=float(area_array[selected].sum()),
                max_assignment_distance_m=(
                    float(absolute[index][selected].max()) if count else 0.0
                ),
            )
        )
    return ids, names, tuple(stats)


def fem_surface(
    source: Any,
    *,
    components: Mapping[str, GeometryField] | Sequence[GeometryField] | None = None,
    units: str | None = None,
    coordinate_frame: str | None = None,
    time_s: float = 0.0,
) -> FemSurface:
    """Build the node/element tables for a sampled field or a bare mesh.

    ``source`` is a :class:`~tensorcad.geometry_field.SampledField` (preferred: the
    handoff record supplies the digest and quality warnings) or a bare :class:`Mesh`.
    ``components`` names the boundary patches; with ``None`` the whole surface is one
    patch called ``"wall"``, which is what a single-body fluid solve wants.

    The orientation is checked here rather than trusted: a surface whose signed volume
    came out negative is rejected with :class:`ValueError`, because an inward-facing
    surface turns every pressure load sign inside out, and that is exactly the class of
    error a solver cannot detect.
    """
    from ..mesh import edge_face_counts, is_edge_manifold, is_watertight, total_signed_volume

    if time_s < 0.0:
        raise ValueError(f"time_s must be non-negative, got {time_s!r}")
    mesh: Mesh
    digest = ""
    quality: dict[str, Any] = {}
    if isinstance(source, Mesh):
        mesh = source
    else:
        mesh = source.mesh
        handoff = source.handoff()
        digest = str(handoff["mesh_sha256"])
        quality = dict(handoff.get("quality", {}))
    if not mesh.faces:
        raise ValueError("cannot hand off an empty mesh")

    volume = float(total_signed_volume(mesh))
    if volume <= 0.0:
        raise ValueError(
            f"surface signed volume is {volume:.6g} m^3; the triangles face inward, "
            "so export would invert every load sign"
        )
    edges = edge_face_counts(mesh)
    areas, _ = _triangle_areas_and_centroids(mesh)
    region_ids: tuple[int, ...]
    region_names: tuple[str, ...]
    region_stats: tuple[RegionStats, ...]
    if components is None:
        region_ids = tuple(1 for _ in mesh.faces)
        region_names = ("wall",)
        region_stats = (RegionStats("wall", 1, len(mesh.faces), float(sum(areas)), 0.0),)
    else:
        region_ids, region_names, region_stats = label_regions(mesh, components)

    return FemSurface(
        nodes=tuple((float(x), float(y), float(z)) for x, y, z in mesh.vertices),
        triangles=mesh.faces,
        region_ids=region_ids,
        region_names=region_names,
        regions=region_stats,
        units=units if units is not None else mesh.units,
        coordinate_frame=coordinate_frame if coordinate_frame is not None else mesh.coordinate_frame,
        time_s=float(time_s),
        volume_m3=volume,
        area_m2=float(sum(areas)),
        # Measured on the exported tables, not copied from the producer: a consumer
        # that receives only this file must be able to trust these three numbers.
        watertight=is_watertight(mesh),
        edge_manifold=is_edge_manifold(mesh),
        euler_characteristic=len(mesh.vertices) - len(edges) + len(mesh.faces),
        source_digest=digest,
        quality=quality,
    )


def write_gmsh(surface: FemSurface, path: str | Path) -> Path:
    """Write ``surface`` as gmsh 2.2 ASCII with one physical tag per region.

    Floats use ``repr``, so the coordinates in the file are the doubles that were
    measured; a truncated ``%.6f`` would move nodes by up to a micrometre and, at
    these cell sizes, is the difference between two different meshes. The units and
    the digest go into a ``$Comments`` block because gmsh's node block has no place
    for them.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$Comments",
        f"generated by TensorCAD {FEM_SCHEMA}",
        f"units {surface.units}",
        f"coordinate_frame {surface.coordinate_frame}",
        f"sign_convention {surface.sign_convention}",
        f"time_s {surface.time_s!r}",
        f"digest {surface.digest}",
        "$EndComments",
        "$PhysicalNames",
        str(len(surface.region_names)),
    ]
    for index, name in enumerate(surface.region_names):
        # Dimension 2 (surface) with the physical id one-based; quotes are required.
        lines.append(f'{2} {index + 1} "{name}"')
    lines.append("$EndPhysicalNames")
    lines.append("$Nodes")
    lines.append(str(len(surface.nodes)))
    for node_id, (x, y, z) in enumerate(surface.nodes, start=1):
        lines.append(f"{node_id} {x!r} {y!r} {z!r}")
    lines.append("$EndNodes")
    lines.append("$Elements")
    lines.append(str(len(surface.triangles)))
    for element_id, (face, region) in enumerate(
        zip(surface.triangles, surface.region_ids, strict=True), start=1
    ):
        a, b, c = (index + 1 for index in face)
        lines.append(f"{element_id} {_GMSH_TRIANGLE} 2 {region} {region} {a} {b} {c}")
    lines.append("$EndElements")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def read_gmsh(path: str | Path) -> dict[str, Any]:
    """Parse a gmsh 2.2 file written by :func:`write_gmsh` and summarise it.

    Deliberately a *verifier*, not a general gmsh reader: it reads the node, element
    and physical-name blocks and returns counts, coordinates and per-physical-tag
    areas, which is exactly what a round-trip test needs to prove the file is
    readable and unchanged.
    """
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(f"no gmsh file at {source}")
    tokens = source.read_text(encoding="utf-8").split("\n")
    index = 0

    def seek(name: str) -> int:
        nonlocal index
        while index < len(tokens) and tokens[index].strip() != name:
            index += 1
        if index == len(tokens):
            raise ValueError(f"{source} has no {name} block")
        index += 1
        return index

    seek("$Nodes")
    node_count = int(tokens[index])
    index += 1
    nodes: dict[int, tuple[float, float, float]] = {}
    for _ in range(node_count):
        node_fields = tokens[index].split()
        nodes[int(node_fields[0])] = (
            float(node_fields[1]),
            float(node_fields[2]),
            float(node_fields[3]),
        )
        index += 1
    seek("$Elements")
    element_count = int(tokens[index])
    index += 1
    elements: dict[int, tuple[int, tuple[int, ...]]] = {}
    for _ in range(element_count):
        element_fields = [int(value) for value in tokens[index].split()]
        index += 1
        element_id, kind, tag_count = element_fields[0], element_fields[1], element_fields[2]
        tags = element_fields[3 : 3 + tag_count]
        element_corners = tuple(element_fields[3 + tag_count :])
        if kind != _GMSH_TRIANGLE:
            raise ValueError(f"element {element_id} is not a triangle (gmsh type {kind})")
        elements[element_id] = (tags[0], element_corners)
    names: dict[int, str] = {}
    index = 0
    seek("$PhysicalNames")
    for _ in range(int(tokens[index])):
        index += 1
        dimension, physical_id, label = tokens[index].split(maxsplit=2)
        names[int(physical_id)] = label.strip().strip('"')
        if int(dimension) != 2:
            raise ValueError(f"physical name {label.strip()} is not a surface")
    index += 1

    areas: dict[int, float] = {}
    for physical, corners in elements.values():
        a, b, c = (nodes[node] for node in corners)
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        areas[physical] = areas.get(physical, 0.0) + 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)
    return {
        "nodes": node_count,
        "triangles": element_count,
        "node_table": nodes,
        "physical_names": names,
        "area_by_region_m2": areas,
    }


def write_fem_json(surface: FemSurface, path: str | Path) -> Path:
    """Write the full manifest *and* the tables as JSON, floats at full precision."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = surface.to_dict(include_tables=True)
    target.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return target


def read_fem_json(path: str | Path) -> FemSurface:
    """Rebuild a :class:`FemSurface` from :func:`write_fem_json` output."""
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(f"no FEM JSON at {source}")
    payload = json.loads(source.read_text(encoding="utf-8"))
    if payload.get("schema") != FEM_SCHEMA:
        raise ValueError(f"{source} has schema {payload.get('schema')!r}, expected {FEM_SCHEMA!r}")
    nodes = tuple(tuple(node) for node in payload["node_table"])
    triangles = tuple(tuple(face) for face in payload["element_table"])
    entries = list(payload["regions"])
    stats = tuple(
        RegionStats(
            name=str(entry["name"]),
            physical_id=int(entry["physical_id"]),
            triangles=int(entry["triangles"]),
            area_m2=float(entry["area_m2"]),
            max_assignment_distance_m=float(entry["max_assignment_distance_m"]),
        )
        for entry in entries
    )
    return FemSurface(
        nodes=nodes,
        triangles=triangles,
        region_ids=tuple(payload["region_table"]),
        region_names=tuple(stats_entry.name for stats_entry in stats),
        regions=stats,
        units=payload["units"],
        coordinate_frame=payload["coordinate_frame"],
        time_s=float(payload["time_s"]),
        sign_convention=payload["sign_convention"],
        interpolation=payload["interpolation"],
        integration=payload["integration"],
        volume_m3=float(payload["signed_volume_m3"]),
        area_m2=float(payload["surface_area_m2"]),
        watertight=bool(payload["watertight"]),
        edge_manifold=bool(payload["edge_manifold"]),
        euler_characteristic=int(payload["euler_characteristic"]),
        source_digest=payload.get("mesh_sha256", ""),
        quality=payload.get("quality", {}),
    )
