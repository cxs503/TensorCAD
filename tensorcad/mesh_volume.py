"""Closed-body volume and watertightness verification against analytic solids.

The question this module answers is the first one a mesh has to pass before any
solver or hydrostatics table can use it: *does this closed surface bound the
volume of the body it claims to be, and is it actually closed?*  A tessellated
ship hull has no analytic volume, so the question is asked where an exact answer
exists -- the sphere and the rectangular box of classical solid geometry -- and
the answer is decomposed into the two things that can be wrong:

* **the summation** (divergence theorem over triangles), which is *exact* for any
  closed polyhedron and must therefore reproduce an analytically representable
  body (a box) to floating-point round-off; and
* **the tessellation**, which for a curved body is an inscribed polyhedron and is
  therefore *below* the analytic volume by a discretisation deficit that must go
  to zero at the mesh's own order -- ``O(h**2)`` for the icosphere used here, and
  reported as an observed order rather than assumed.

Keeping the two apart is the point.  A single "sphere volume is 0.05 % low" row
cannot tell a reader whether the volume integrator has a bug or the mesh is
simply coarse; a box row (exact) plus a sphere ladder (deficit -> 0 at order 2)
can.

References are closed forms, evaluated in code, never transcribed:

========================  ==============================================
quantity                  reference (evaluated, not copied)
========================  ==============================================
sphere volume             ``4/3 * pi * r**3``
sphere area               ``4 * pi * r**2``
box volume                ``a * b * c``
box area                  ``2 * (a*b + b*c + c*a)``
regular icosahedron       ``V = (5/12)(3+sqrt(5)) * a**3`` with the
(circumradius ``R``)      edge ``a = 4R/sqrt(10 + 2*sqrt(5))``
========================  ==============================================

The icosahedron row matters: it is the ``subdivisions=0`` icosphere, i.e. the
*actual* 20-triangle mesh the ladder starts from, and its closed-form volume is
derived from the regular polyhedron rather than from any mesh summation.  If the
divergence-theorem integrator and that formula agree to round-off on 20 skewed
triangles, the integrator is verified on more than the six axis-aligned planes of
a box.

Tolerances are module constants so a caller cannot tune them after seeing a
result:

* :data:`RELATIVE_TOLERANCE` -- 0.1 %, the acceptance line of the benchmark case
  ``CASE-CAD-VOLUME-CLOSED``;
* :data:`EXACT_TOLERANCE` -- 1e-12 relative, for bodies a triangle mesh can
  represent exactly (the box): only floating-point round-off is allowed;
* :data:`ROUTE_TOLERANCE` -- 1e-10 relative, for the agreement of two summation
  routes over the same triangles.

What this module deliberately does **not** claim:

* It does not claim a coarse triangulation of a curved body has the analytic
  volume.  It reports the deficit and the observed order; only the refined level
  is gated against the 0.1 % line, and which level that is happens to be is read
  off the ladder, not asserted.
* It does not re-implement watertightness.  The closed/composite verdict and the
  defect vocabulary are :func:`tensorcad.geometry_audit.audit_watertightness`'s;
  this module only composes it with the volume number so that "closed" and "the
  right volume" are never reported apart.
* It does not repair anything.  Repair is
  :func:`tensorcad.geomquality.mesh_defects.repair_facet_defects`'s job, and a
  repair is only interesting *after* this check can see the damage.

The negative controls (:func:`defect_detection`) exist because a check that
passes everything proves nothing: a hole, a flipped facet, a duplicated facet and
an unwelded triangle soup (what an STL carries before welding) are injected into
a clean sphere and the report must name each one.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from .geometry import Mesh
from .geometry_audit import audit_watertightness, inject_duplicate_faces, inject_flip_faces, inject_hole
from .mesh_xval import TriMesh, cube, icosphere, ellipsoid as ellipsoid_mesh, torus as torus_mesh_gen

__all__ = [
    "EXACT_TOLERANCE",
    "FLOAT32_TOLERANCE",
    "RELATIVE_TOLERANCE",
    "ROUTE_TOLERANCE",
    "ClosedBodyCheck",
    "analytic_box_area",
    "analytic_box_volume",
    "analytic_cone_area",
    "analytic_cone_volume",
    "analytic_cylinder_area",
    "analytic_cylinder_volume",
    "analytic_icosahedron_volume",
    "analytic_ellipsoid_volume",
    "analytic_sphere_area",
    "analytic_sphere_volume",
    "analytic_sphere_volume_via_cylinder",
    "analytic_torus_area",
    "analytic_torus_volume",
    "as_geometry_mesh",
    "box_mesh",
    "check_closed_body",
    "cone_mesh",
    "convergence_study",
    "cylinder_convergence_study",
    "cylinder_mesh",
    "defect_detection",
    "ellipsoid_convergence_study",
    "mean_edge_length",
    "observed_order",
    "torus_convergence_study",
    "torus_mesh",
    "triangle_area",
    "volume_by_tet_centroid",
]

#: Acceptance line for the benchmark case: 0.1 % relative error.
RELATIVE_TOLERANCE = 1.0e-3
#: A triangle mesh represents a box exactly; only round-off is allowed.
EXACT_TOLERANCE = 1.0e-12
#: Two summation routes over the same closed surface must agree this closely.
ROUTE_TOLERANCE = 1.0e-10
#: A binary STL stores float32 coordinates (~1.2e-7 relative); three facets per
#: corner plus the volume sum allow a few ulp of that, so a 1e-6 line is the
#: honest statement "the file, not the integrator, is the error source".
FLOAT32_TOLERANCE = 1.0e-6


# ---------------------------------------------------------------------------
# analytic references (evaluated here, never typed in as decimal values)
# ---------------------------------------------------------------------------
def analytic_sphere_volume(radius: float) -> float:
    """``V = 4/3 pi r^3``."""
    return 4.0 / 3.0 * math.pi * float(radius) ** 3


def analytic_sphere_volume_via_cylinder(radius: float) -> float:
    """The same sphere volume from Euclid's circumscribed-cylinder relation.

    *Elements* XII.18 (exhaustion) proves the sphere is two thirds of the
    circumscribed cylinder, so ``V = (2/3)(pi r^2)(2r)``.  Evaluating that
    expression separately from ``4/3 pi r^3`` gives the reference its own
    provenance check: two derivations, no transcribed decimal.
    """
    r = float(radius)
    cylinder = math.pi * r**2 * (2.0 * r)
    return 2.0 / 3.0 * cylinder


def analytic_sphere_area(radius: float) -> float:
    """``A = 4 pi r^2``."""
    return 4.0 * math.pi * float(radius) ** 2


def analytic_box_volume(a: float, b: float, c: float) -> float:
    """``V = a b c``."""
    return float(a) * float(b) * float(c)


def analytic_box_area(a: float, b: float, c: float) -> float:
    """``A = 2 (ab + bc + ca)``."""
    a, b, c = float(a), float(b), float(c)
    return 2.0 * (a * b + b * c + c * a)


def analytic_icosahedron_volume(circumradius: float) -> float:
    """Volume of a regular icosahedron whose vertices lie on ``circumradius``.

    Edge ``a`` and circumradius ``R`` satisfy ``a = 4R / sqrt(10 + 2 sqrt(5))``;
    the volume of a regular icosahedron of edge ``a`` is ``(5/12)(3+sqrt(5)) a^3``.
    Both are evaluated from ``math.sqrt(5)`` so no table is transcribed, and this
    is the closed form for the ``subdivisions=0`` icosphere -- an independent
    check of the divergence-theorem sum on 20 non-axis-aligned triangles.
    """
    r = float(circumradius)
    edge = 4.0 * r / math.sqrt(10.0 + 2.0 * math.sqrt(5.0))
    return (5.0 / 12.0) * (3.0 + math.sqrt(5.0)) * edge**3


def analytic_cylinder_volume(radius: float, height: float) -> float:
    """``V = pi r^2 h``."""
    return math.pi * float(radius) ** 2 * float(height)


def analytic_cylinder_area(radius: float, height: float) -> float:
    """Total surface area: ``A = 2 pi r (r + h)`` (side + 2 caps)."""
    r, h = float(radius), float(height)
    return 2.0 * math.pi * r * (r + h)


def analytic_cone_volume(radius: float, height: float) -> float:
    """``V = 1/3 pi r^2 h``."""
    return math.pi * float(radius) ** 2 * float(height) / 3.0


def analytic_cone_area(radius: float, height: float) -> float:
    """Total surface area: ``A = pi r (r + sqrt(r^2 + h^2))`` (base + lateral)."""
    r, h = float(radius), float(height)
    slant = math.sqrt(r * r + h * h)
    return math.pi * r * (r + slant)


def analytic_torus_volume(major: float, minor: float) -> float:
    """Volume of a torus: ``V = 2 pi^2 R r^2`` (Pappus's centroid theorem).

    The torus is generated by revolving a circle of radius ``r`` (minor) whose
    centre is at distance ``R`` (major) from the axis of revolution.  By Pappus's
    second centroid theorem, the volume is the product of the generating circle's
    area (``pi r^2``) and the distance its centroid travels (``2 pi R``).
    """
    R, r = float(major), float(minor)
    return 2.0 * math.pi**2 * R * r**2


def analytic_torus_area(major: float, minor: float) -> float:
    """Surface area of a torus: ``A = 4 pi^2 R r`` (Pappus's first centroid theorem).

    The surface area is the product of the generating circle's circumference
    (``2 pi r``) and the distance its centroid travels (``2 pi R``).
    """
    R, r = float(major), float(minor)
    return 4.0 * math.pi**2 * R * r


# ---------------------------------------------------------------------------
# meshes
# ---------------------------------------------------------------------------

def analytic_ellipsoid_volume(semi_axes) -> float:
    """Volume of a triaxial ellipsoid: ``V = 4/3 pi a b c``.

    Delegates to :func:`tensorcad.analytic_solids.ellipsoid_volume` so the
    reference is evaluated by the same converged quadrature used elsewhere
    in the codebase.
    """
    from .analytic_solids import ellipsoid_volume
    return ellipsoid_volume(semi_axes)


def box_mesh(a: float, b: float, c: float) -> TriMesh:
    """Axis-aligned box ``a x b x c`` centred on the origin, as 12 triangles.

    Built by scaling the existing :func:`tensorcad.mesh_xval.cube` rather than by
    assembling eight new corners, so the box and the sphere share one generator
    and one orientation convention.
    """
    base = cube(2.0)
    scale = np.array([float(a), float(b), float(c)], dtype=np.float64) / 2.0
    return TriMesh(base.points * scale, base.faces, f"box{a:g}x{b:g}x{c:g}")


def cylinder_mesh(radius: float, height: float, n_radial: int = 32, n_axial: int = 8) -> TriMesh:
    """Closed cylinder mesh for volume verification.

    The mesh is an inscribed N-gonal prism with triangulated caps.  Because the
    side is flat-faced (not curved), the volume deficit comes only from the
    radial discretisation and converges at ``O(1/n_radial^2)``.
    """
    from .mesh_xval import cylinder
    return cylinder(float(radius), float(height), int(n_radial), int(n_axial))


def cone_mesh(radius: float, height: float, n_radial: int = 64) -> TriMesh:
    """Closed cone mesh for volume verification.

    The mesh is an inscribed N-gonal pyramid with a triangulated base.  The
    apex singularity is handled by the mesh generator; the volume deficit
    converges at ``O(1/n_radial^2)``.
    """
    from .mesh_xval import cone
    return cone(float(radius), float(height), int(n_radial))


def torus_mesh(major: float = 1.0, minor: float = 0.3, nu: int = 64, nv: int = 32) -> TriMesh:
    """Closed torus mesh for volume verification (genus 1).

    The torus is parameterised by major radius ``R`` (distance from centre of
    tube to centre of torus) and minor radius ``r`` (tube radius).  The mesh is
    a regular ``(nu x nv)`` grid on the torus surface, triangulated into
    ``2 * nu * nv`` triangles.

    Unlike the sphere/icosphere, the torus is **not** an inscribed polyhedron
    in the same sense -- the vertices lie *on* the torus surface, so the volume
    error can be positive or negative depending on how the triangles approximate
    the curved surface.  The error still converges at ``O(h^2)`` for a
    piecewise-linear approximation.
    """
    return torus_mesh_gen(float(major), float(minor), int(nu), int(nv))


def as_geometry_mesh(mesh: TriMesh) -> Mesh:
    """Convert the exchange container into the core :class:`tensorcad.geometry.Mesh`."""
    return Mesh(
        vertices=tuple(map(tuple, mesh.points.tolist())),
        faces=tuple(map(tuple, mesh.faces.tolist())),
    )


# ---------------------------------------------------------------------------
# measures
# ---------------------------------------------------------------------------
def triangle_area(mesh: TriMesh) -> float:
    """Total triangle area (not signed; area is orientation-blind)."""
    a, b, c = (mesh.points[mesh.faces[:, i]] for i in range(3))
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


def mean_edge_length(mesh: TriMesh) -> float:
    """Mean length of the distinct mesh edges -- the ``h`` of the convergence study."""
    edges = np.vstack([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    p, q = mesh.points[edges[:, 0]], mesh.points[edges[:, 1]]
    return float(np.linalg.norm(q - p, axis=1).mean())


def volume_by_tet_centroid(mesh: TriMesh) -> float:
    """Volume as the signed sum of tetrahedra from the mesh centroid.

    A second summation *route*, not a second theorem: ``V = sum det(a-C, b-C, c-C)/6``
    telescopes to the same identity as the divergence-theorem form, but the terms
    are different (``O(h**3)`` differences from an interior point instead of
    ``O(R**3)`` corner products), so a sign or indexing error in one does not
    cancel silently in the other.  Agreement is therefore evidence about the
    *implementation*, and this module does not present it as independent physics.
    """
    centroid = mesh.points[np.unique(mesh.faces)].mean(axis=0)
    a = mesh.points[mesh.faces[:, 0]] - centroid
    b = mesh.points[mesh.faces[:, 1]] - centroid
    c = mesh.points[mesh.faces[:, 2]] - centroid
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# ---------------------------------------------------------------------------
# one closed-body check
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ClosedBodyCheck:
    """Volume, area and topology of one closed body against its closed form."""

    name: str
    family: str
    faces: int
    vertices: int
    mean_edge_length: float
    reference_volume: float
    measured_volume: float
    tet_volume: float
    reference_area: float
    measured_area: float
    open_edges: int
    boundary_edges: int
    euler_characteristic: int
    audit_open_edges: int
    audit_nonmanifold_edges: int
    audit_duplicate_faces: int
    audit_winding_conflicts: int
    audit_signed_volume: float
    watertight: bool
    audit_defects: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    # -- derived ----------------------------------------------------------
    @property
    def volume_error(self) -> float:
        """Signed relative error ``(measured - reference) / reference``."""
        return (self.measured_volume - self.reference_volume) / self.reference_volume

    @property
    def area_error(self) -> float:
        return (self.measured_area - self.reference_area) / self.reference_area

    @property
    def route_difference(self) -> float:
        """Relative difference between the divergence and centroid-tet routes."""
        return abs(self.measured_volume - self.tet_volume) / abs(self.reference_volume)

    @property
    def topologically_closed(self) -> bool:
        """Closed, manifold, single-piece, consistently wound, outward, genus 0 or 1.

        Genus 0 (sphere-like): Euler characteristic χ = 2.
        Genus 1 (torus-like): Euler characteristic χ = 0.
        Both are valid closed surfaces for volume integration.
        """
        return (
            self.open_edges == 0
            and self.audit_open_edges == 0
            and self.audit_nonmanifold_edges == 0
            and self.audit_duplicate_faces == 0
            and self.audit_winding_conflicts == 0
            and self.audit_signed_volume > 0.0
            and self.euler_characteristic in (2, 0)
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "family": self.family,
            "faces": self.faces,
            "vertices": self.vertices,
            "mean_edge_length": self.mean_edge_length,
            "reference_volume": self.reference_volume,
            "measured_volume": self.measured_volume,
            "tet_volume": self.tet_volume,
            "volume_error": self.volume_error,
            "reference_area": self.reference_area,
            "measured_area": self.measured_area,
            "area_error": self.area_error,
            "route_difference": self.route_difference,
            "open_edges": self.open_edges,
            "boundary_edges": self.boundary_edges,
            "euler_characteristic": self.euler_characteristic,
            "audit": {
                "open_edges": self.audit_open_edges,
                "nonmanifold_edges": self.audit_nonmanifold_edges,
                "duplicate_faces": self.audit_duplicate_faces,
                "winding_conflicts": self.audit_winding_conflicts,
                "signed_volume": self.audit_signed_volume,
                "watertight": self.watertight,
                "defects": list(self.audit_defects),
            },
            "topologically_closed": self.topologically_closed,
            "notes": list(self.notes),
        }


def check_closed_body(
    mesh: TriMesh,
    *,
    name: str,
    family: str,
    reference_volume: float,
    reference_area: float,
    length: float | None = None,
    notes: Sequence[str] = (),
) -> ClosedBodyCheck:
    """Measure one closed body and compare it with its closed form.

    ``length`` is the audit's reference length for its relative gap tolerance; the
    default (bounding-box length) is what a caller who has no declared ``L`` gets.
    """
    surface = as_geometry_mesh(mesh)
    audit = audit_watertightness(surface, length=length)
    return ClosedBodyCheck(
        name=name,
        family=family,
        faces=len(mesh.faces),
        vertices=len(mesh.points),
        mean_edge_length=mean_edge_length(mesh),
        reference_volume=float(reference_volume),
        measured_volume=mesh.volume(),
        tet_volume=volume_by_tet_centroid(mesh),
        reference_area=float(reference_area),
        measured_area=triangle_area(mesh),
        open_edges=int(mesh.boundary_edges()),
        boundary_edges=int(audit.open_edges),
        euler_characteristic=int(mesh.euler()),
        audit_open_edges=int(audit.open_edges),
        audit_nonmanifold_edges=int(audit.nonmanifold_edges),
        audit_duplicate_faces=int(audit.duplicate_faces),
        audit_winding_conflicts=int(audit.winding_conflicts),
        audit_signed_volume=float(audit.signed_volume),
        watertight=bool(audit.passed),
        audit_defects=tuple(defect.kind for defect in audit.defects),
        notes=tuple(notes),
    )


# ---------------------------------------------------------------------------
# convergence: the tessellation error, separated from the integrator error
# ---------------------------------------------------------------------------
def convergence_study(
    radius: float = 1.0,
    levels: Sequence[int] = range(0, 7),
) -> list[ClosedBodyCheck]:
    """Volume/area deficit of the inscribed icosphere at each subdivision level.

    ``levels=0`` is the raw icosahedron, whose volume has a closed form
    (:func:`analytic_icosahedron_volume`); the remaining levels are checked
    against the sphere.  Topology is audited at *every* level even though it is
    expected not to change: "the chain stayed watertight through refinement" is a
    measurement, not an assumption.
    """
    volume = analytic_sphere_volume(radius)
    area = analytic_sphere_area(radius)
    rows: list[ClosedBodyCheck] = []
    for level in levels:
        mesh = icosphere(radius, int(level))
        row = check_closed_body(
            mesh,
            name=f"sphere_r{radius:g}_sub{level}",
            family="sphere",
            reference_volume=volume,
            reference_area=area,
            length=2.0 * float(radius),
        )
        if int(level) == 0:
            icosa = analytic_icosahedron_volume(radius)
            row = ClosedBodyCheck(
                **{
                    **row.__dict__,
                    "notes": (
                        f"icosahedron closed form {icosa:.12f}; divergence route "
                        f"{(row.measured_volume - icosa) / icosa:+.3e} relative",
                    ),
                }
            )
        rows.append(row)
    return rows


def cylinder_convergence_study(
    radius: float = 1.0,
    height: float = 2.0,
    radial_levels: Sequence[int] = (16, 32, 64, 128, 256),
) -> list[ClosedBodyCheck]:
    """Volume/area deficit of an inscribed cylinder at each radial resolution.

    The cylinder side is a flat-faced prism, so the volume deficit comes only
    from the radial discretisation (the caps are exact triangles).  The error
    converges at ``O(1/n^2)`` where n is the number of radial divisions.

    Note: ``mean_edge_length`` is not a good ``h`` here because the axial
    resolution is held fixed; the observed order should be computed against
    ``2 pi r / n_radial`` (the arc length between radial vertices).
    """
    volume = analytic_cylinder_volume(radius, height)
    area = analytic_cylinder_area(radius, height)
    rows: list[ClosedBodyCheck] = []
    for nr in radial_levels:
        mesh = cylinder_mesh(radius, height, n_radial=nr, n_axial=max(8, nr // 4))
        row = check_closed_body(
            mesh,
            name=f"cylinder_r{radius:g}_h{height:g}_nr{nr}",
            family="cylinder",
            reference_volume=volume,
            reference_area=area,
            length=2.0 * max(radius, height),
        )
        rows.append(row)
    return rows


def cone_convergence_study(
    radius: float = 1.0,
    height: float = 2.0,
    radial_levels: Sequence[int] = (16, 32, 64, 128, 256),
) -> list[ClosedBodyCheck]:
    """Volume/area deficit of an inscribed cone at each radial resolution.

    The cone side is a flat-faced pyramid, so the volume deficit comes only
    from the radial discretisation (the base is exact triangles).  The error
    converges at ``O(1/n^2)`` where n is the number of radial divisions.

    Note: ``mean_edge_length`` is not a good ``h`` here because the cone mesh
    has only one axial division; the observed order should be computed against
    ``2 pi r / n_radial`` (the arc length between radial vertices).
    """
    volume = analytic_cone_volume(radius, height)
    area = analytic_cone_area(radius, height)
    rows: list[ClosedBodyCheck] = []
    for nr in radial_levels:
        mesh = cone_mesh(radius, height, n_radial=nr)
        row = check_closed_body(
            mesh,
            name=f"cone_r{radius:g}_h{height:g}_nr{nr}",
            family="cone",
            reference_volume=volume,
            reference_area=area,
            length=2.0 * max(radius, height),
        )
        rows.append(row)
    return rows


def torus_convergence_study(
    major: float = 1.0,
    minor: float = 0.3,
    resolution_levels: Sequence[int] = (16, 32, 64, 128, 256),
) -> list[ClosedBodyCheck]:
    """Volume/area deficit of a torus mesh at each resolution.

    The torus mesh is a regular ``(nu x nv)`` grid on the surface, triangulated
    into ``2 * nu * nv`` triangles.  Both ``nu`` (poloidal) and ``nv`` (toroidal)
    are set to the same value from ``resolution_levels``.

    Unlike the icosphere (which is strictly inscribed), the torus mesh vertices
    lie *on* the surface, so the volume error can be positive or negative.
    The error converges at ``O(h^2)`` for a piecewise-linear approximation.

    This is the first **genus-1** body in the case pool, testing the volume
    integrator on a non-simply-connected surface.
    """
    volume = analytic_torus_volume(major, minor)
    area = analytic_torus_area(major, minor)
    rows: list[ClosedBodyCheck] = []
    for n in resolution_levels:
        mesh = torus_mesh(major, minor, nu=n, nv=max(8, n // 2))
        row = check_closed_body(
            mesh,
            name=f"torus_R{major:g}_r{minor:g}_n{n}",
            family="torus",
            reference_volume=volume,
            reference_area=area,
            length=2.0 * (major + minor),
        )
        rows.append(row)
    return rows


def observed_order(rows: Sequence[ClosedBodyCheck], *, metric: str = "volume_error") -> float:
    """Least-squares slope of ``log|error|`` against ``log h`` over the ladder.

    Only rows with a finite, non-zero error are used (a level whose error is at
    round-off carries no slope information).  The value is *observed*: nothing in
    the module assumes it is 2.
    """
    hs: list[float] = []
    es: list[float] = []
    for row in rows:
        err = abs(getattr(row, metric))
        if math.isfinite(err) and err > 0.0:
            hs.append(math.log(row.mean_edge_length))
            es.append(math.log(err))
    if len(hs) < 2:
        return float("nan")
    h = np.asarray(hs)
    e = np.asarray(es)
    return float(((h - h.mean()) * (e - e.mean())).sum() / ((h - h.mean()) ** 2).sum())


def first_level_within(
    rows: Sequence[ClosedBodyCheck], tolerance: float = RELATIVE_TOLERANCE
) -> str | None:
    """Name of the coarsest level whose ``|volume error|`` is inside ``tolerance``."""
    for row in rows:
        if abs(row.volume_error) <= tolerance:
            return row.name
    return None


# ---------------------------------------------------------------------------
# negative controls
# ---------------------------------------------------------------------------
def _soup(mesh: TriMesh) -> TriMesh:
    """The same facets with connectivity destroyed -- what an STL carries."""
    corners = mesh.points[mesh.faces].reshape(-1, 3)
    faces = np.arange(len(corners), dtype=np.int64).reshape(-1, 3)
    return TriMesh(corners, faces, mesh.name + "/soup")



def ellipsoid_convergence_study(
    semi_axes,
    levels: Sequence[int] = range(0, 7),
) -> list[ClosedBodyCheck]:
    """Volume/area deficit of an inscribed ellipsoid mesh at each subdivision level.

    The mesh is built by scaling a unit icosphere, so the topology is identical
    to the sphere ladder but the curvature is non-uniform.  This tests the
    volume integrator on skewed triangles that the sphere's symmetry hides.

    The reference volume is ``4/3 pi a b c`` and the reference area is computed
    by converged Gauss-Legendre quadrature (:func:`analytic_solids.ellipsoid_surface_area`).
    """
    from .analytic_solids import ellipsoid_surface_area, ellipsoid_volume

    a, b, c = (float(x) for x in semi_axes)
    volume = ellipsoid_volume(semi_axes)
    area = ellipsoid_surface_area(semi_axes)
    rows: list[ClosedBodyCheck] = []
    for level in levels:
        mesh = ellipsoid_mesh(semi_axes, int(level))
        row = check_closed_body(
            mesh,
            name=f"ellipsoid{a:g}x{b:g}x{c:g}_sub{level}",
            family="ellipsoid",
            reference_volume=volume,
            reference_area=area,
            length=2.0 * max(a, b, c),
        )
        rows.append(row)
    return rows


def defect_detection(
    mesh: TriMesh | None = None,
    *,
    radius: float = 1.0,
    level: int = 2,
) -> list[dict[str, Any]]:
    """Inject four defect families into a clean sphere and report what is named.

    Each row carries ``expected`` (the audit fields that must move) and
    ``detected`` (whether they did).  A check that returns ``clean`` for every one
    of these is worthless, which is why they are part of the module and not only
    of its tests.
    """
    clean = mesh if mesh is not None else icosphere(radius, level)
    surface = as_geometry_mesh(clean)
    reference = clean.volume()

    def row(name: str, defect: str, candidate: Mesh, expected: dict[str, int]) -> dict[str, Any]:
        audit = audit_watertightness(candidate, length=2.0 * float(radius))
        measured = {
            "open_edges": int(audit.open_edges),
            "nonmanifold_edges": int(audit.nonmanifold_edges),
            "duplicate_faces": int(audit.duplicate_faces),
            "winding_conflicts": int(audit.winding_conflicts),
            "watertight": 0 if audit.passed else 1,
        }
        detected = all(measured[key] >= want for key, want in expected.items())
        try:
            volume = TriMesh(
                np.array(candidate.vertices, dtype=np.float64),
                np.array(candidate.faces, dtype=np.int64),
            ).volume()
            volume_error = (volume - reference) / reference
        except Exception as exc:  # pragma: no cover - defensive
            volume_error = float("nan")
            del exc
        return {
            "name": name,
            "defect": defect,
            "expected": expected,
            "measured": measured,
            "defects_named": [defect.kind for defect in audit.defects],
            "detected": detected,
            "volume_error": volume_error,
        }

    return [
        row(
            "hole",
            "facets removed inside an x/z window",
            inject_hole(surface, x_centre=0.0, z_centre=0.0, x_size=0.3, z_size=0.3),
            {"open_edges": 1, "watertight": 1},
        ),
        row(
            "flipped_facet",
            "facet winding reversed on 4 facets",
            inject_flip_faces(surface, count=4),
            {"winding_conflicts": 1, "watertight": 1},
        ),
        row(
            "duplicate_facet",
            "3 facets repeated",
            inject_duplicate_faces(surface, count=3),
            {"duplicate_faces": 3, "watertight": 1},
        ),
        row(
            "unwelded_soup",
            "connectivity destroyed (STL before welding)",
            as_geometry_mesh(_soup(clean)),
            {"open_edges": 1, "watertight": 1},
        ),
    ]
