"""OpenFOAM ``polyMesh``: read the mesh a solver actually gets, and gate it with the same library.

Why this module exists
---------------------
:mod:`tensorcad.meshquality.topology` already computes OpenFOAM's internal-face
non-orthogonality, but it takes *our* structured cell list. A solver consumer hands over
a ``constant/polyMesh`` directory instead: ``points``, ``faces``, ``owner``,
``neighbour``, ``boundary``. Until this module nothing in the tree could read that, so the
one quantity the documentation compares against ``checkMesh`` had no way to be measured on
a real case file. The gate could only judge meshes this repository had built itself.

What is implemented, and to whose definition
--------------------------------------------
* :func:`read_polymesh` / :func:`write_polymesh` -- ASCII ``polyMesh`` in and out, no
  third-party dependency. Binary files raise with the conversion command instead of
  silently producing garbage.
* :func:`polymesh_from_hexes` -- build a ``polyMesh`` from a hexahedral connectivity, so
  the repository's own meshes can be pushed through the solver-format route.
* :func:`face_geometry` -- face centre, area and area vector. Centre and area both come
  from the triangle fan about the first vertex, which is what ``Foam::face::centre`` and
  ``Foam::face::area`` do; a quad's centre is therefore its area centroid, not its vertex
  average.
* :func:`cell_geometry` -- volume and centre by OpenFOAM's pyramid decomposition
  (``primitiveMesh::calcCellCentresAndVols``): the cell is split into pyramids from an
  area-weighted estimate ``cEst`` of the centre, each pyramid's signed volume and centroid
  accumulated, and the centre taken as the volume-weighted mean. For a planar-faced
  closed polyhedron this is exact, not an approximation; it is *not* the same as the vertex
  average :func:`tensorcad.meshquality.topology.non_orthogonality` uses, and the benchmark
  reports how far apart the two routes are.
* :func:`check_polymesh` -- the three numbers ``checkMesh`` prints, as data:
  max non-orthogonality (deg), max skewness (dimensionless), cell volume range, plus a
  :class:`~tensorcad.meshquality.report.GateReport`.

Definitions, stated because they are where implementations disagree
-------------------------------------------------------------------
non-orthogonality
    ``acos( Sf . d / (|Sf| |d|) )`` in degrees, ``Sf`` the owner-outward face area vector
    and ``d`` the vector from the owner cell centre to the neighbour cell centre. 0 deg on
    a Cartesian mesh.
skewness
    ``t = ((Cf - C) . Sf) / (d . Sf)``; the cell-centre line ``C + t d`` meets the face
    plane at ``Pi = C + t d``; the reported number is ``|Pi - Cf| / |d|``. 0 when the line
    crosses the face at its centre, 1 when the crossing point is one cell-centre distance
    away from it. OpenFOAM's ``checkMesh`` warns above 4.

Only the standard library is used. See ``docs/mesh_quality_polymesh.md`` for the numbers
and ``examples/mesh_quality/run_polymesh_benchmark.py`` for the run that produces them.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from .report import Check, GateReport, summary_stats

__all__ = [
    "DEFAULT_POLYMESH_LIMITS",
    "Patch",
    "PolyMesh",
    "PolyMeshQuality",
    "cell_geometry",
    "check_polymesh",
    "face_geometry",
    "face_non_orthogonality",
    "face_skewness",
    "polymesh_from_hexes",
    "read_polymesh",
    "write_polymesh",
]

Vec3 = tuple[float, float, float]

#: ``checkMesh``'s own warning bands, as defaults: non-orthogonality above 70 deg is
#: "severely non-orthogonal", skewness above 4 is "highly skew". Both are policy, printed
#: with the measured margin by the benchmark.
DEFAULT_POLYMESH_LIMITS: dict[str, float] = {
    "max_non_orthogonality": 70.0,
    "max_skewness": 4.0,
}


# --------------------------------------------------------------------------------------
# container
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Patch:
    """One boundary patch: a contiguous run of boundary faces."""

    name: str
    kind: str
    n_faces: int
    start: int

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "type": self.kind, "nFaces": self.n_faces, "startFace": self.start}


@dataclass(frozen=True)
class PolyMesh:
    """An OpenFOAM polyMesh: shared points, polygonal faces, owner/neighbour addressing."""

    points: tuple[Vec3, ...]
    faces: tuple[tuple[int, ...], ...]
    owner: tuple[int, ...]
    neighbour: tuple[int, ...]
    patches: tuple[Patch, ...] = ()

    def __post_init__(self) -> None:
        if len(self.owner) != len(self.faces):
            raise ValueError(f"owner has {len(self.owner)} entries for {len(self.faces)} faces")
        if len(self.neighbour) > len(self.owner):
            raise ValueError("more neighbours than faces")
        for index, face in enumerate(self.faces):
            if len(face) < 3:
                raise ValueError(f"face {index} has {len(face)} points; a face needs at least 3")
            for point in face:
                if not 0 <= point < len(self.points):
                    raise ValueError(f"face {index} references point {point} of {len(self.points)}")
        for index, cell in enumerate(self.owner):
            if cell < 0:
                raise ValueError(f"face {index} has a negative owner {cell}")
        for index, cell in enumerate(self.neighbour):
            if cell < 0:
                raise ValueError(f"internal face {index} has a negative neighbour {cell}")
        covered = 0
        for patch in self.patches:
            if patch.start != len(self.neighbour) + covered:
                raise ValueError(
                    f"patch {patch.name!r} starts at face {patch.start}, expected "
                    f"{len(self.neighbour) + covered}; boundary faces must be contiguous and "
                    "cover the owner-only tail of the face list"
                )
            covered += patch.n_faces
        if self.patches and covered != len(self.owner) - len(self.neighbour):
            raise ValueError(
                f"patches cover {covered} boundary faces, the mesh has "
                f"{len(self.owner) - len(self.neighbour)}"
            )

    @property
    def n_cells(self) -> int:
        cells = set(self.owner)
        cells.update(self.neighbour)
        return max(cells) + 1 if cells else 0

    @property
    def n_internal_faces(self) -> int:
        return len(self.neighbour)

    @property
    def n_boundary_faces(self) -> int:
        return len(self.owner) - len(self.neighbour)

    def to_dict(self) -> dict[str, object]:
        return {
            "points": len(self.points),
            "faces": len(self.faces),
            "internal_faces": self.n_internal_faces,
            "boundary_faces": self.n_boundary_faces,
            "cells": self.n_cells,
            "patches": [p.to_dict() for p in self.patches],
        }


# --------------------------------------------------------------------------------------
# text plumbing
# --------------------------------------------------------------------------------------
def _strip_comments(text: str) -> str:
    """Remove ``//`` line comments and ``/* ... */`` blocks, keeping offsets readable."""
    out: list[str] = []
    index, length = 0, len(text)
    while index < length:
        if text.startswith("/*", index):
            end = text.find("*/", index + 2)
            index = length if end < 0 else end + 2
            out.append(" ")
        elif text.startswith("//", index):
            end = text.find("\n", index)
            index = length if end < 0 else end + 1
            out.append(" ")
        else:
            out.append(text[index])
            index += 1
    return "".join(out)


_TOKEN = re.compile(r"[()]|[^\s()]+")


class _Cursor:
    """Minimal cursor over the token stream of one OpenFOAM list."""

    def __init__(self, tokens: Sequence[str]) -> None:
        self._tokens = list(tokens)
        self._index = 0

    def peek(self) -> str | None:
        return self._tokens[self._index] if self._index < len(self._tokens) else None

    def take(self) -> str:
        token = self.peek()
        if token is None:
            raise ValueError("unexpected end of data")
        self._index += 1
        return token

    def expect(self, token: str) -> None:
        found = self.take()
        if found != token:
            raise ValueError(f"expected {token!r}, found {found!r}")


def _number(token: str) -> float:
    try:
        return float(token)
    except ValueError as error:  # pragma: no cover - defensive
        raise ValueError(f"expected a number, found {token!r}") from error


def _parse_value(cursor: _Cursor) -> object:
    token = cursor.peek()
    if token == "(":
        cursor.take()
        items: list[object] = []
        while cursor.peek() not in (")", None):
            items.append(_parse_value(cursor))
        cursor.expect(")")
        return items
    return _number(cursor.take())


def _parse_list(text: str) -> list[object]:
    """Parse ``count ( entry entry ... )`` and return the entries.

    OpenFOAM writes a list of faces in either ``( (0 1 2 3) (1 2 3 4) )`` form or the
    compact ``( 4(0 1 2 3) 4(1 2 3 4) )`` form; both collapse to one entry per item here.
    """
    start = text.find("(")
    if start < 0:
        raise ValueError("no list found")
    count_token = text[:start].split()
    if not count_token:
        raise ValueError("no item count before the list")
    declared = int(_number(count_token[-1]))
    body = _parse_value(_Cursor(_TOKEN.findall(text[start:])))
    if not isinstance(body, list):
        raise ValueError("expected a parenthesised list")
    if all(not isinstance(item, list) for item in body):
        # a flat list of scalars (owner, neighbour, a labelList)
        if len(body) != declared:
            raise ValueError(f"header declares {declared} entries, found {len(body)}")
        return body
    entries: list[object] = []
    pending: float | None = None
    for item in body:
        if isinstance(item, list):
            if pending is not None:
                if len(item) != int(pending):
                    raise ValueError(
                        f"compact face declares {int(pending)} points, lists {len(item)}"
                    )
                pending = None
            entries.append(item)
        else:
            if pending is not None:
                raise ValueError("two counts in a row")
            pending = item
    if pending is not None:
        raise ValueError("trailing count with no list")
    if len(entries) != declared:
        raise ValueError(f"header declares {declared} entries, found {len(entries)}")
    return entries


def _header(text: str) -> tuple[str, int]:
    """Return ``(format, offset)`` where ``offset`` is just past the ``FoamFile`` block."""
    match = re.search(r"FoamFile\s*\{", text)
    if not match:
        raise ValueError("no FoamFile header: this is not an OpenFOAM mesh file")
    end = text.find("}", match.end())
    if end < 0:
        raise ValueError("unterminated FoamFile header")
    block = text[match.end():end]
    fmt = re.search(r"format\s+(\w+)\s*;", block)
    if not fmt:
        raise ValueError("FoamFile header has no format entry")
    if fmt.group(1) != "ascii":
        raise ValueError(
            f"polyMesh format is {fmt.group(1)!r}, only 'ascii' is read here; convert the case "
            "with foamFormatConvert -ascii before handing it over"
        )
    return fmt.group(1), end + 1


def _read_file(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"missing polyMesh file {path}")
    return _strip_comments(path.read_text(encoding="utf-8", errors="replace"))


_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  (written by tensorcad.meshquality)    |
|   \\\\  /    A nd           |                                                 |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "constant/polyMesh";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def _write_list(path: Path, cls: str, obj: str, lines: Iterable[str], count: int) -> None:
    body = "\n".join(lines)
    path.write_text(_HEADER.format(cls=cls, obj=obj) + f"{count}\n(\n{body}\n)\n", encoding="utf-8")


def read_polymesh(path: str | Path) -> PolyMesh:
    """Read an ASCII ``polyMesh`` from a case directory or a ``polyMesh`` directory."""
    directory = Path(path)
    if (directory / "constant" / "polyMesh").is_dir():
        directory = directory / "constant" / "polyMesh"
    if not (directory / "points").is_file():
        raise ValueError(f"{directory} does not look like a polyMesh directory (no 'points')")

    points_text = _read_file(directory / "points")
    _header(points_text)
    points = tuple(
        (float(p[0]), float(p[1]), float(p[2]))  # type: ignore[index]
        for p in _parse_list(points_text)
    )

    faces_text = _read_file(directory / "faces")
    _header(faces_text)
    faces = tuple(tuple(int(v) for v in face) for face in _parse_list(faces_text))  # type: ignore[union-attr]

    owner_text = _read_file(directory / "owner")
    _header(owner_text)
    owner = tuple(int(v) for v in _parse_list(owner_text))

    neighbour_path = directory / "neighbour"
    neighbour: tuple[int, ...] = ()
    if neighbour_path.is_file():
        neighbour_text = _read_file(neighbour_path)
        _header(neighbour_text)
        neighbour = tuple(int(v) for v in _parse_list(neighbour_text))

    patches: tuple[Patch, ...] = ()
    boundary_path = directory / "boundary"
    if boundary_path.is_file():
        boundary_text = _read_file(boundary_path)
        _header(boundary_text)
        tail = boundary_text[boundary_text.find("(", boundary_text.find("}")) :]
        entries = re.finditer(r"(\w+)\s*\{([^}]*)\}", tail, re.S)
        found: list[Patch] = []
        for entry in entries:
            fields = dict(
                (key, value.strip()) for key, value in re.findall(r"(\w+)\s+([^;]+);", entry.group(2))
            )
            found.append(
                Patch(
                    name=entry.group(1),
                    kind=fields.get("type", "patch"),
                    n_faces=int(_number(fields["nFaces"])),
                    start=int(_number(fields["startFace"])),
                )
            )
        patches = tuple(found)

    return PolyMesh(points=points, faces=faces, owner=owner, neighbour=neighbour, patches=patches)


def write_polymesh(mesh: PolyMesh, directory: str | Path) -> Path:
    """Write ``mesh`` as an ASCII ``polyMesh`` (``points``/``faces``/``owner``/...)."""
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    _write_list(
        target / "points",
        "vectorField",
        "points",
        (f"({p[0]:.17g} {p[1]:.17g} {p[2]:.17g})" for p in mesh.points),
        len(mesh.points),
    )
    _write_list(
        target / "faces",
        "faceList",
        "faces",
        (f"{len(f)}({' '.join(str(i) for i in f)})" for f in mesh.faces),
        len(mesh.faces),
    )
    _write_list(
        target / "owner",
        "labelList",
        "owner",
        (str(c) for c in mesh.owner),
        len(mesh.owner),
    )
    _write_list(
        target / "neighbour",
        "labelList",
        "neighbour",
        (str(c) for c in mesh.neighbour),
        len(mesh.neighbour),
    )
    lines: list[str] = []
    for patch in mesh.patches:
        lines.append(f"{patch.name}")
        lines.append("{")
        lines.append(f"    type            {patch.kind};")
        lines.append(f"    nFaces          {patch.n_faces};")
        lines.append(f"    startFace       {patch.start};")
        lines.append("}")
    _write_list(target / "boundary", "polyBoundaryMesh", "boundary", lines, len(mesh.patches))
    return target


# --------------------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------------------
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def face_geometry(
    points: Sequence[Vec3], faces: Sequence[Sequence[int]]
) -> tuple[tuple[Vec3, ...], tuple[Vec3, ...], tuple[float, ...]]:
    """Face centres, area vectors and areas by the triangle fan about vertex 0.

    The fan is what ``Foam::face::centre``/``Foam::face::area`` use, so a non-planar face
    is scored the same way here as by the solver. A degenerate (zero area) face yields a
    zero area vector and a centre at vertex 0's fan average, and is counted downstream as
    degenerate rather than dropped.
    """
    centres: list[Vec3] = []
    area_vectors: list[Vec3] = []
    areas: list[float] = []
    for face in faces:
        rings = [points[i] for i in face]
        origin = rings[0]
        sum_centre = [0.0, 0.0, 0.0]
        sum_vector = [0.0, 0.0, 0.0]
        total = 0.0
        for index in range(1, len(rings) - 1):
            a = _sub(rings[index], origin)
            b = _sub(rings[index + 1], origin)
            vector = _cross(a, b)
            area = 0.5 * _norm(vector)
            if area <= 0.0:
                continue
            tri_centre = (
                (origin[0] + rings[index][0] + rings[index + 1][0]) / 3.0,
                (origin[1] + rings[index][1] + rings[index + 1][1]) / 3.0,
                (origin[2] + rings[index][2] + rings[index + 1][2]) / 3.0,
            )
            total += area
            for axis in range(3):
                sum_centre[axis] += area * tri_centre[axis]
                sum_vector[axis] += 0.5 * vector[axis]
        if total > 0.0:
            centre = (sum_centre[0] / total, sum_centre[1] / total, sum_centre[2] / total)
        else:
            n = float(len(rings))
            centre = (
                sum(p[0] for p in rings) / n,
                sum(p[1] for p in rings) / n,
                sum(p[2] for p in rings) / n,
            )
        area_vector = (sum_vector[0], sum_vector[1], sum_vector[2])
        centres.append(centre)
        area_vectors.append(area_vector)
        areas.append(_norm(area_vector))
    return tuple(centres), tuple(area_vectors), tuple(areas)


def _cell_faces(mesh: PolyMesh) -> list[list[tuple[int, int]]]:
    """Per cell: ``(face index, +1 if the face area vector points out of this cell)``."""
    fans: list[list[tuple[int, int]]] = [[] for _ in range(mesh.n_cells)]
    for index, cell in enumerate(mesh.owner):
        fans[cell].append((index, +1))
    for index, cell in enumerate(mesh.neighbour):
        fans[cell].append((index, -1))
    return fans


def cell_geometry(mesh: PolyMesh) -> tuple[tuple[Vec3, ...], tuple[float, ...]]:
    """Cell centres and volumes by OpenFOAM's pyramid decomposition.

    ``cEst`` is the area-weighted mean of the cell's face centres; the cell is then split
    into pyramids from ``cEst`` to each face, and volume and centroid are the signed sums.
    For a closed polyhedron with planar faces the volume is exact and the centre is the
    true volume centroid -- there is no truncation error to excuse, only a sign convention
    (outward face area vectors), which the orientation test below enforces.
    """
    centres, area_vectors, areas = face_geometry(mesh.points, mesh.faces)
    fans = _cell_faces(mesh)
    cell_centres: list[Vec3] = []
    cell_volumes: list[float] = []
    for fan in fans:
        if not fan:
            raise ValueError("a cell with no faces: owner/neighbour addressing is inconsistent")
        weight = sum(areas[face] for face, _ in fan)
        if weight <= 0.0:
            raise ValueError("a cell whose faces all have zero area")
        estimate = tuple(
            sum(areas[face] * centres[face][axis] for face, _ in fan) / weight
            for axis in range(3)
        )
        volume = 0.0
        accumulator = [0.0, 0.0, 0.0]
        for face, sign in fan:
            outward = tuple(sign * component for component in area_vectors[face])
            contribution = _dot(_sub(centres[face], estimate), outward)  # type: ignore[arg-type]
            pyramid_volume = contribution / 3.0
            pyramid_centre = (
                (estimate[0] + 3.0 * centres[face][0]) / 4.0,
                (estimate[1] + 3.0 * centres[face][1]) / 4.0,
                (estimate[2] + 3.0 * centres[face][2]) / 4.0,
            )
            volume += pyramid_volume
            for axis in range(3):
                accumulator[axis] += pyramid_volume * pyramid_centre[axis]
        if abs(volume) <= 0.0:
            cell_centres.append(estimate)  # type: ignore[arg-type]
            cell_volumes.append(0.0)
            continue
        cell_centres.append(
            (accumulator[0] / volume, accumulator[1] / volume, accumulator[2] / volume)
        )
        cell_volumes.append(abs(volume))
    return tuple(cell_centres), tuple(cell_volumes)


def face_non_orthogonality(mesh: PolyMesh) -> tuple[float, ...]:
    """``acos(Sf . d / (|Sf| |d|))`` in degrees for every internal face."""
    centres, _ = cell_geometry(mesh)
    _, area_vectors, areas = face_geometry(mesh.points, mesh.faces)
    angles: list[float] = []
    for index, neighbour in enumerate(mesh.neighbour):
        owner = mesh.owner[index]
        d = _sub(centres[neighbour], centres[owner])
        mag_d = _norm(d)
        if areas[index] <= 0.0 or mag_d <= 0.0:
            angles.append(math.nan)
            continue
        cosine = _dot(area_vectors[index], d) / (areas[index] * mag_d)
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return tuple(angles)


def face_skewness(mesh: PolyMesh) -> tuple[float, ...]:
    """``|Pi - Cf| / |d|`` for every internal face, ``Pi`` the cell-centre line crossing."""
    centres, _ = cell_geometry(mesh)
    face_centres, area_vectors, areas = face_geometry(mesh.points, mesh.faces)
    values: list[float] = []
    for index, neighbour in enumerate(mesh.neighbour):
        owner = mesh.owner[index]
        d = _sub(centres[neighbour], centres[owner])
        mag_d = _norm(d)
        denominator = _dot(d, area_vectors[index])
        if mag_d <= 0.0 or denominator == 0.0:
            values.append(math.nan)
            continue
        t = _dot(_sub(face_centres[index], centres[owner]), area_vectors[index]) / denominator
        crossing = (
            centres[owner][0] + t * d[0],
            centres[owner][1] + t * d[1],
            centres[owner][2] + t * d[2],
        )
        values.append(_norm(_sub(crossing, face_centres[index])) / mag_d)
    return tuple(values)


# --------------------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class PolyMeshQuality:
    """The ``checkMesh`` numbers of one polyMesh, as data."""

    mesh: PolyMesh
    non_orthogonality: tuple[float, ...]
    skewness: tuple[float, ...]
    volumes: tuple[float, ...]

    def statistics(self) -> dict[str, dict[str, object]]:
        return {
            "non_orthogonality": summary_stats("non_orthogonality", self.non_orthogonality).to_dict(),
            "skewness": summary_stats("skewness", self.skewness).to_dict(),
            "cell_volume": summary_stats("cell_volume", self.volumes).to_dict(),
        }

    def to_dict(self) -> dict[str, object]:
        return {
            "mesh": self.mesh.to_dict(),
            "total_volume": math.fsum(self.volumes),
            "statistics": self.statistics(),
        }

    def gate(self, limits: dict[str, float] | None = None) -> GateReport:
        limits = dict(DEFAULT_POLYMESH_LIMITS if limits is None else limits)
        ortho = summary_stats("max_non_orthogonality", self.non_orthogonality)
        skew = summary_stats("max_skewness", self.skewness)
        volume = summary_stats("cell_volume", self.volumes)
        non_positive = sum(1 for v in self.volumes if v <= 0.0)
        checks = [
            Check(
                name="polymesh.max_non_orthogonality",
                kind="polymesh",
                metric="max_non_orthogonality",
                passed=bool(ortho.count) and ortho.maximum <= limits["max_non_orthogonality"],
                value=ortho.maximum,
                limit=f"<= {limits['max_non_orthogonality']:g} deg",
                detail=(
                    f"max {ortho.maximum:.6g} deg, mean {ortho.mean:.6g} over "
                    f"{ortho.count} internal face(s)"
                ),
            ),
            Check(
                name="polymesh.max_skewness",
                kind="polymesh",
                metric="max_skewness",
                passed=bool(skew.count) and skew.maximum <= limits["max_skewness"],
                value=skew.maximum,
                limit=f"<= {limits['max_skewness']:g}",
                detail=f"max {skew.maximum:.6g}, mean {skew.mean:.6g} over {skew.count} internal face(s)",
            ),
            Check(
                name="polymesh.cell_volume",
                kind="polymesh",
                metric="cell_volume",
                passed=non_positive == 0,
                value=volume.minimum,
                limit="> 0 for every cell",
                detail=(
                    f"min {volume.minimum:.6g} m^3, max {volume.maximum:.6g} m^3, "
                    f"{non_positive} non-positive cell(s) of {volume.count}"
                ),
            ),
        ]
        return GateReport(tuple(checks))


def check_polymesh(mesh: PolyMesh) -> PolyMeshQuality:
    """Measure a polyMesh with the OpenFOAM definitions (non-orthogonality, skewness, volume)."""
    _, volumes = cell_geometry(mesh)
    return PolyMeshQuality(
        mesh=mesh,
        non_orthogonality=face_non_orthogonality(mesh),
        skewness=face_skewness(mesh),
        volumes=volumes,
    )


# --------------------------------------------------------------------------------------
# from this repository's own connectivity
# --------------------------------------------------------------------------------------
#: The six hexahedral faces with outward **cyclic** order for a cell whose nodes follow the
#: VTK hexahedron convention (0-3 bottom loop, 4-7 top loop). The orientation is checked
#: against the cell centre in :func:`polymesh_from_hexes`, not trusted.
HEX_FACES: tuple[tuple[int, ...], ...] = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def polymesh_from_hexes(
    points: Sequence[Sequence[float]],
    hexes: Sequence[Sequence[int]],
    *,
    patch_name: str = "walls",
) -> PolyMesh:
    """Build a ``PolyMesh`` from hexahedral connectivity, orienting faces outward.

    Faces are matched by their sorted node key, so conforming cells share one face and the
    first occurrence becomes the owner. A key shared by more than two cells means the input
    is not a conforming hex mesh and raises, rather than quietly picking an owner. The
    result follows OpenFOAM's layout: every internal face first (``neighbour`` aligned with
    it), then the boundary faces, which is what the ``boundary`` patches address.
    """
    vertices: tuple[Vec3, ...] = tuple(
        (float(p[0]), float(p[1]), float(p[2])) for p in points  # type: ignore[index]
    )
    records: list[tuple[tuple[int, ...], int, int | None]] = []
    seen: dict[tuple[int, ...], int] = {}
    for cell_index, cell in enumerate(hexes):
        if len(cell) != 8:
            raise ValueError(f"hex {cell_index} has {len(cell)} nodes, expected 8")
        centre = tuple(sum(vertices[i][axis] for i in cell) / 8.0 for axis in range(3))
        for template in HEX_FACES:
            ring = tuple(int(cell[i]) for i in template)
            centres_face, vectors_face, areas_face = face_geometry(vertices, (ring,))
            if areas_face[0] > 0.0:
                outward = _dot(vectors_face[0], _sub(centres_face[0], centre))
            else:
                # A collapsed column (zero breadth at keel/stem) leaves a zero-area face.
                # The mesh is not rejected -- measuring that defect is the reason to read
                # it -- so the ring is kept as given and the face measures as degenerate
                # (NaN non-orthogonality, zero cell volume) further down.
                outward = 1.0
            if outward < 0.0:
                ring = tuple(reversed(ring))
            key = tuple(sorted(ring))
            index = seen.get(key)
            if index is None:
                seen[key] = len(records)
                records.append((ring, cell_index, None))
                continue
            ring0, owner0, neighbour0 = records[index]
            if neighbour0 is not None:
                raise ValueError(f"face {key} is shared by more than two cells")
            if owner0 == cell_index:
                raise ValueError(f"cell {cell_index} uses face {key} twice")
            records[index] = (ring0, owner0, cell_index)

    internal = [record for record in records if record[2] is not None]
    boundary = [record for record in records if record[2] is None]
    ordered = internal + boundary
    mesh = PolyMesh(
        points=vertices,
        faces=tuple(record[0] for record in ordered),
        owner=tuple(record[1] for record in ordered),
        neighbour=tuple(int(record[2]) for record in internal),
    )
    patches = (
        (Patch(name=patch_name, kind="patch", n_faces=len(boundary), start=len(internal)),)
        if boundary
        else ()
    )
    return PolyMesh(
        points=mesh.points,
        faces=mesh.faces,
        owner=mesh.owner,
        neighbour=mesh.neighbour,
        patches=patches,
    )
