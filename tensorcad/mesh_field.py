"""Dense geometry fields over an *explicit triangle mesh*.

Why this exists
---------------
Every field TensorCAD shipped before this module is analytic: a plane, a sphere,
a parametric hull, a SUBOFF body. Those fields are cheap to sample (closed form
per point) and that is exactly what the voxelization ladder in
``tensorcad.perf_benchmark`` measured -- ``grid_resolution`` against wall clock,
with the triangle count column structurally zero.

Real pipelines do not work that way. A CAD interchange mesh, a repaired STL, a
scanned hull, a B-rep tessellation at 10^6--10^7 triangles has *no* closed form;
the solver-side SDF (LBM/FVM solid fraction) has to be derived from the triangles
themselves. That is the regime where voxelization becomes expensive, and until
there is a mesh-backed field it cannot even be measured.

:class:`MeshSDFField` closes that hole: it adapts a ``Mesh`` to the same
:class:`~tensorcad.geometry_field.GeometryField` protocol the analytic primitives
implement, so ``to_voxel_grid(field, resolution)``, ``to_mesh`` and everything
else that consumes a field works unchanged on triangle soup.

Definition
----------
``value(p) = -d(p)`` inside the closed surface, ``+d(p)`` outside, with ``d`` the
**exact Euclidean distance to the nearest point on the surface** (not a
nearest-vertex or k-ring approximation): :meth:`TriangleGrid.nearest_surface`
returns the closest point and prunes with a provable lower bound, so the value
agrees with a brute-force scan over every triangle to rounding error (pinned by
``tests/test_mesh_field.py``). The sign comes from ray-crossing parity
(:meth:`TriangleGrid.contains`), which requires a closed surface -- for an open
mesh the magnitude is still correct, the sign is meaningless, and the field
raises unless ``closed=False`` is passed explicitly.

Cost, measured
--------------
Per query the field costs one shell search plus one parity ray, so a grid of
``n`` cells costs ``n`` queries. That is O(cells) *queries* (each cheaper than a
full scan, but far from the closed-form cost of an analytic field), and it is
why voxelizing an explicit 10^7-triangle mesh is a different engineering problem
from sampling ``SuboffHullField``. ``docs/perf_scale_10m.md`` §12 records the
measured numbers.

Standard library only, like the rest of the core; NumPy is not needed even though
``to_voxel_grid`` uses it to hold the result.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any, ClassVar

from .geometry import BoundingBox, Mesh
from .geometry_field import (
    DEFAULT_FRAME,
    DEFAULT_UNITS,
    NEGATIVE_INSIDE,
    SIGNED_DISTANCE,
    FieldError,
    GeometryField,
    Point,
)
from .spatial_index import TriangleGrid

#: Field kind recorded in metadata / provenance.
MESH_SDF_KIND = "mesh_sdf"


class MeshSDFField(GeometryField):
    """Exact signed-distance field of a closed triangle mesh.

    Parameters
    ----------
    mesh:
        Source triangles. Vertices may be shared; nothing is mutated.
    grid:
        Pre-built :class:`TriangleGrid` for the mesh. Supplying one lets a caller
        pay the index build once and reuse it across grids/resolutions.
    resolution:
        Index resolution ``(nx, ny, nz)`` when ``grid`` is not given; ``None``
        lets the grid pick an occupancy-driven resolution.
    closed:
        Assert the surface is closed (the default). Pass ``False`` for an open
        mesh only when the *magnitude* is all that is wanted -- the sign of an
        open surface is not defined and the parity test will return arbitrary
        values there.
    units / coordinate_frame:
        Copied into the field metadata so a solver handoff is self-describing.
    band:
        Narrow-band truncation for the distance magnitude.  When set, the
        magnitude is capped at ``band`` (a deep-interior cell reports ``-band``
        instead of its true distance).  The sign is unaffected.  This is the
        narrow-band rule real voxelisers use: with ``band >= h/2`` (``h`` the
        largest cell edge) the partial-fill solid fraction is unchanged, so the
        voxelised volume is unchanged, but deep-interior cells stop paying for a
        full nearest-surface search.  ``None`` (the default) keeps the exact
        distance everywhere.
    """

    kind: ClassVar[str] = MESH_SDF_KIND
    quantity: ClassVar[str] = SIGNED_DISTANCE
    inside_sign: ClassVar[str] = NEGATIVE_INSIDE

    def __init__(
        self,
        mesh: Mesh,
        *,
        grid: TriangleGrid | None = None,
        resolution: Sequence[int] | None = None,
        closed: bool = True,
        units: str | None = None,
        coordinate_frame: str | None = None,
        band: float | None = None,
    ) -> None:
        if len(mesh.faces) == 0:
            raise FieldError("a mesh field needs at least one triangle")
        if band is not None and band <= 0:
            raise FieldError(f"band must be positive or None, got {band!r}")
        self.mesh = mesh
        if grid is not None:
            self.grid = grid
        elif resolution is not None:
            res = tuple(int(v) for v in resolution)
            if len(res) != 3:
                raise FieldError(f"resolution must have 3 entries, got {resolution!r}")
            self.grid = TriangleGrid(mesh, resolution=(res[0], res[1], res[2]))
        else:
            self.grid = TriangleGrid(mesh)
        self.closed = bool(closed)
        self.band = band
        self.units = units if units is not None else getattr(mesh, "units", DEFAULT_UNITS)
        self.coordinate_frame = (
            coordinate_frame
            if coordinate_frame is not None
            else getattr(mesh, "coordinate_frame", DEFAULT_FRAME)
        )

    # -- GeometryField contract -------------------------------------------- #

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        """Exact signed distance at each point (negative inside a closed mesh)."""
        out: list[float] = []
        nearest = self.grid.nearest_surface
        contains = self.grid.contains
        for point in points:
            p = (float(point[0]), float(point[1]), float(point[2]))
            hit = nearest(p, band=self.band)
            distance = 0.0 if hit is None else hit.distance
            out.append(-distance if contains(p) else distance)
        return out

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Central-difference gradient (the surface has no closed-form derivative).

        A distance field's gradient is the unit normal, but taking it from
        differences keeps the sign changes and the creases (where the distance
        field is not differentiable) honest instead of smoothing them over.
        """
        batch = [(float(p[0]), float(p[1]), float(p[2])) for p in points]
        return self.numerical_gradient(batch)

    def bounds(self) -> BoundingBox | None:
        vertices = self.mesh.vertices
        if not vertices:
            return None
        xs = [v[0] for v in vertices]
        ys = [v[1] for v in vertices]
        zs = [v[2] for v in vertices]
        return BoundingBox(
            minimum=(min(xs), min(ys), min(zs)),
            maximum=(max(xs), max(ys), max(zs)),
        )

    def state(self) -> dict[str, float]:
        """A mesh field has no editable scalars; its shape is the triangles."""
        return {}

    def with_state(self, **values: float) -> GeometryField:
        raise FieldError(
            "a mesh field has no editable state; edit the triangles instead"
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "triangles": len(self.mesh.faces),
            "vertices": len(self.mesh.vertices),
            "index_resolution": list(self.grid.resolution),
            "index_cells": self.grid.n_cells,
            "index_entries": self.grid.n_entries,
            "closed": self.closed,
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "triangle_mesh",
            "kind": self.kind,
            "triangles": len(self.mesh.faces),
            "index_build_seconds": self.grid.build_seconds,
        }

    # -- conveniences ------------------------------------------------------- #

    def nearest_distance(self, point: Sequence[float]) -> float:
        """Unsigned distance to the surface (one query, no sign test)."""
        hit = self.grid.nearest_surface((float(point[0]), float(point[1]), float(point[2])))
        return math.inf if hit is None else hit.distance

    def resolution_hint(self, spacing: float) -> tuple[int, int, int]:
        """Cells per axis that give approximately ``spacing`` cells around the body."""
        box = self.bounds()
        if box is None:
            raise FieldError("mesh has no vertices")
        if spacing <= 0:
            raise FieldError(f"spacing must be positive, got {spacing!r}")
        return tuple(max(1, int(math.ceil(size / spacing))) for size in box.sizes)


def voxelize_mesh(
    mesh: Mesh,
    resolution: int | Sequence[int] = 32,
    **kwargs: Any,
) -> Any:
    """Voxelise an explicit mesh: ``to_voxel_grid(MeshSDFField(mesh), resolution)``.

    Thin wrapper so callers do not have to import two modules to turn triangles
    into a solver-side SDF grid. NumPy is required (the grid is a NumPy array).
    """
    from .voxel import to_voxel_grid

    field_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in {"grid", "resolution_index", "closed", "units", "coordinate_frame", "band"}}
    field = MeshSDFField(
        mesh,
        grid=field_kwargs.get("grid"),
        resolution=field_kwargs.get("resolution_index"),
        closed=field_kwargs.get("closed", True),
        units=field_kwargs.get("units"),
        coordinate_frame=field_kwargs.get("coordinate_frame"),
        band=field_kwargs.get("band"),
    )
    return to_voxel_grid(field, resolution, **kwargs)