"""Backend-neutral geometry contracts.

Everything in this module is pure standard library on purpose: the core package
must install and run in constrained agent environments, while mature CAD kernels
(OCP/CadQuery, gmsh, ...) stay optional backends behind :mod:`tensorcad.adapters`.

The contract is intentionally small: a triangulated surface (positions plus
1-based-style triangle indices stored 0-based internally) carrying the units and
coordinate frame it was generated in, so that downstream consumers never have to
guess either.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

Vec3 = tuple[float, float, float]
Face = tuple[int, int, int]


def _canonical(point: Iterable[float]) -> Vec3:
    """Normalise a point so that ``-0.0`` and ``0.0`` are the same position."""
    values = tuple(float(value) for value in point)
    if len(values) != 3:
        raise ValueError(f"a point needs exactly 3 coordinates, got {len(values)}")
    result = (values[0] + 0.0, values[1] + 0.0, values[2] + 0.0)
    for axis, value in enumerate(result):
        if not math.isfinite(value):
            raise ValueError(f"vertex coordinate {axis} is not finite: {value!r}")
    return result


@dataclass(frozen=True)
class BoundingBox:
    """Axis-aligned bounds in the model's own units."""

    minimum: Vec3
    maximum: Vec3

    @property
    def sizes(self) -> Vec3:
        return (
            self.maximum[0] - self.minimum[0],
            self.maximum[1] - self.minimum[1],
            self.maximum[2] - self.minimum[2],
        )

    def to_dict(self) -> dict[str, tuple[float, float, float]]:
        return {"min": self.minimum, "max": self.maximum, "size": self.sizes}


@dataclass(frozen=True)
class Mesh:
    """Triangulated surface: vertex positions plus triangle index triples.

    ``faces`` holds 0-based indices into ``vertices``. Use :meth:`to_obj_faces`
    when writing formats that are 1-based (such as OBJ).
    """

    vertices: tuple[Vec3, ...]
    faces: tuple[Face, ...]
    units: str = "m"
    coordinate_frame: str = "unspecified right-handed frame"

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", tuple(self.vertices))
        object.__setattr__(self, "faces", tuple(tuple(f) for f in self.faces))
        for index, face in enumerate(self.faces):
            if len(face) != 3:
                raise ValueError(f"face {index} is not a triangle")
            for vertex in face:
                if not 0 <= vertex < len(self.vertices):
                    raise ValueError(f"face {index} references out-of-range vertex {vertex}")

    def __len__(self) -> int:
        return len(self.faces)

    def iter_vertices(self) -> Iterator[Vec3]:
        yield from self.vertices

    def to_obj_faces(self) -> Iterator[Face]:
        """Yield triangle indices offset by one, as OBJ expects."""
        for a, b, c in self.faces:
            yield (a + 1, b + 1, c + 1)

    def bounding_box(self) -> BoundingBox:
        if not self.vertices:
            raise ValueError("cannot bound an empty mesh")
        minimum = [math.inf, math.inf, math.inf]
        maximum = [-math.inf, -math.inf, -math.inf]
        for x, y, z in self.vertices:
            for axis, value in enumerate((x, y, z)):
                minimum[axis] = min(minimum[axis], value)
                maximum[axis] = max(maximum[axis], value)
        return BoundingBox(tuple(minimum), tuple(maximum))  # type: ignore[arg-type]

    def metadata(self) -> dict[str, object]:
        data: dict[str, object] = {
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "vertex_count": len(self.vertices),
            "triangle_count": len(self.faces),
        }
        if self.vertices:
            data["bounding_box"] = self.bounding_box().to_dict()
        return data


class MeshBuilder:
    """Accumulate a triangle mesh while deduplicating vertices by exact position.

    Deduplication is keyed on the exact coordinate tuple (with ``-0.0``
    normalised), never on a floating-point tolerance: two faces meeting at a
    shared corner get the same index by construction, which keeps the resulting
    topology reproducible across runs and platforms.
    """

    def __init__(self, units: str = "m", coordinate_frame: str = Mesh.coordinate_frame) -> None:
        self.units = units
        self.coordinate_frame = coordinate_frame
        self._vertices: list[Vec3] = []
        self._index: dict[Vec3, int] = {}
        self._faces: list[Face] = []
        self._degenerate_skipped = 0

    @property
    def vertex_count(self) -> int:
        return len(self._vertices)

    @property
    def degenerate_skipped(self) -> int:
        """Number of collapsed quads/triangles that were emitted as one triangle or dropped."""
        return self._degenerate_skipped

    def add_vertex(self, point: Iterable[float]) -> int:
        key = _canonical(point)
        existing = self._index.get(key)
        if existing is not None:
            return existing
        index = len(self._vertices)
        self._vertices.append(key)
        self._index[key] = index
        return index

    def index_of(self, point: Iterable[float]) -> int | None:
        return self._index.get(_canonical(point))

    def add_triangle(self, a: int, b: int, c: int) -> bool:
        """Add a triangle; return ``False`` when it collapsed to a line or point."""
        unique = {a, b, c}
        if len(unique) < 3:
            self._degenerate_skipped += 1
            return False
        self._faces.append((a, b, c))
        return True

    def add_quad(self, a: int, b: int, c: int, d: int) -> int:
        """Add ``(a, b, c)`` and ``(a, c, d)``, collapsing shared corners.

        A lofted quad whose adjacent corners are the same vertex (a section that
        degenerates to a line at a stem or transom) reduces to a single
        triangle; one whose *diagonal* corners coincide has zero area and is
        dropped. Either way no zero-area face ever reaches the output.
        """
        if a == c or b == d:
            self._degenerate_skipped += 1
            return 0
        if a == d:
            return int(self.add_triangle(a, b, c))
        if b == c:
            return int(self.add_triangle(a, b, d))
        self.add_triangle(a, b, c)
        self.add_triangle(a, c, d)
        return 2

    def build(self) -> Mesh:
        return Mesh(
            vertices=tuple(self._vertices),
            faces=tuple(self._faces),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )
