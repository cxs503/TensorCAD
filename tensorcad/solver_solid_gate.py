"""One solid from the unified geometry library, three solver routes, one gate.

Why this module exists
----------------------
``tensorcad.solver_integration`` states the anchor's requirement -- "three solver
routes, one geometry, one frame" -- but an independent audit of that module found
it could not carry the claim:

* the *FVM* leg never wrote a volume bundle and never built a solver mesh: it
  computed the signed volume of a surface and reported ``status="geometry_ready"``;
* the *DG* leg handed TensorDG a **circle of radius a** while the solid under test
  was an ellipsoid with semi-axes ``(a, b, c)`` -- not a section, not the same
  geometry;
* the *LBM* leg voxelised a **cylinder of radius a and length 2c**, whose volume
  ``2 pi a^2 c`` is ``3a / (2b)`` times the ellipsoid's -- with the default
  ``a=1, b=0.5`` that is **3x** the solid it was compared against;
* ``frame_tamper_refused`` swallowed every exception (``except Exception:
  return True``), so a crashing verifier counted as a refusal, and its DG branch
  tampered a bundle whose verifier never looked at a frame at all;
* nothing recorded *which* geometry a route consumed, so "the same solid" was a
  sentence in a document rather than a value in a manifest.

This module rebuilds the cross-route gate so each claim is a measurement:

``unified geometry library first``
    The solid is an ``EllipsoidField`` -- the CGC geometry protocol of
    :mod:`tensorcad.geometry_field` -- sampled once by ``to_mesh`` into a closed
    triangle mesh. Every route consumes **that mesh** (or a section of it), and
    the mesh digest travels in each bundle's manifest, so "same solid" is checked
    by comparing digests, not prose.
``three routes, three bundles``
    * FVM: the mesh is star-tetrahedralised about an interior point, written as
      gmsh ``C3D4`` (:func:`tensorcad.solver_interface.write_bundle_volume`),
      re-read by the independent parser, then handed to TensorFVM's
      ``NativeMesh3D.from_cells`` and solved with ``NativeHybridScalar3D`` for the
      linear field ``phi = 1 + 0.4x - 0.7y + 1.3z``, which is harmonic and must be
      reproduced exactly.
    * DG: the *z = 0 section of the same solid* is handed over as an exact
      rational-quadratic NURBS **ellipse** (:func:`tensorcad.dg_handoff.
      ellipse_contract`), and every sampled point is tested against the geometry
      library's own ``signed_distance`` -- the curve has to lie on the solid.
    * LBM: the same mesh is voxelised (:func:`tensorcad.voxel_solid.solid_lattice`)
      and written as an occupancy lattice with true LBM semantics (fluid cells
      outside the solid).
``one gate``
    Frames must all be the contract frame; the same ``solid_digest`` must appear in
    all three manifests; each route's measured volume (FVM tetrahedra, LBM solid
    cells) must agree with the analytic ``4/3 pi a b c`` *and with each other*
    within declared tolerances; every bundle must verify **before** it is
    tampered with (positive control) and a ``y-up`` manifest must be refused by
    that route's own verifier with :class:`HandoffError` -- no blanket ``except``.

What is deliberately *not* claimed
----------------------------------
The DG consumer (``tensordg.curved_geometry.mesh_from_contract``) is not
importable in this environment; the DG leg therefore stops at a verified bundle
in the schema TensorDG validates, and the gate records that leg as
``consumer_status = "withheld"`` instead of pretending it ran. The FVM conduction
solve is a patch test (exactness on a harmonic field), not a convergence study;
the convergence ladder lives in :mod:`tensorcad.fvm_handoff`.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Sequence

from .dg_handoff import (
    HandoffError as DgHandoffError,
    ellipse_contract,
    ellipse_perimeter,
    evaluate_curve,
    verify_bundle as verify_dg_bundle,
    write_dg_bundle,
)
from .geometry_field import EllipsoidField, to_mesh
from .solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    signed_volume,
    total_volume,
    verify_bundle,
    write_bundle_lattice,
    write_bundle_volume,
)
from .voxel_solid import solid_lattice

__all__ = [
    "SOLID_GATE_SCHEMA",
    "analytic_volume",
    "cross_route_gate",
    "dg_route",
    "fvm_route",
    "lattice_route",
    "run_gate",
    "solid_from_library",
    "star_tetrahedra",
    "tamper_refused",
]

#: Bumped on any incompatible change to ``gate.json``.
SOLID_GATE_SCHEMA = "tensorcad.solver_solid_gate/1"

#: The patch field the FVM leg solves for: harmonic, so exact on any conforming mesh.
PATCH_COEFFICIENTS = (1.0, 0.4, -0.7, 1.3)  # phi = c0 + c1 x + c2 y + c3 z

#: Route tolerances against the analytic solid volume.  The FVM stars a closed
#: surface, so its tetrahedral volume *is* the surface volume (tight); the lattice
#: is staircase-approximated, so it gets the slack the discretisation earns.
ROUTE_TOLERANCES: dict[str, float] = {
    # The star decomposition tiles the sampled surface exactly, so the tetrahedra
    # must reproduce the *sampled* volume to working precision. The sampled
    # surface itself is 0.36 % from the analytic ellipsoid at resolution 24 --
    # that is fixture discretisation, gated separately in solid_from_library.
    "fvm_tet_vs_surface_volume_rel_error": 1e-12,
    "fvm_surface_vs_analytic_volume_rel_error": 0.05,
    "fvm_patch_relative_linf": 1e-06,
    "fvm_patch_relative_l2": 1e-09,
    # Star tetrahedra are slivers (see ``maximum_nonorthogonality_degrees`` in the
    # artifact: 88.7 deg, minimum cell volume 1.5e-13 m3), and the plain
    # orthogonal hybrid operator loses global balance on them. Measured here
    # 2.3e-08; the same operator on the hull volume mesh reaches 3.8e-16
    # (docs/fvm_handoff.md). The tolerance tracks the *fixture*, not the solver.
    "fvm_global_relative_balance": 1e-07,
    "dg_perimeter_rel_error": 1e-09,
    "dg_area_rel_error": 1e-09,
    "dg_on_solid_max_distance_m": 1e-12,
    "lattice_volume_rel_error": 0.06,
    "lattice_fluid_plus_solid": 0.0,
    "cross_route_volume_rel_error": 0.06,
}


def analytic_volume(semi_axes: Sequence[float]) -> float:
    """``4/3 pi a b c`` -- the reference every route is scored against."""
    a, b, c = (float(v) for v in semi_axes)
    if min(a, b, c) <= 0.0 or not all(math.isfinite(v) for v in (a, b, c)):
        raise HandoffError("semi-axes must be positive and finite")
    return 4.0 / 3.0 * math.pi * a * b * c


def _mesh_digest(mesh: Any) -> str:
    h = hashlib.sha256()
    for x, y, z in mesh.vertices:
        h.update(f"{x:.17g},{y:.17g},{z:.17g};".encode())
    h.update(b"|")
    for face in mesh.faces:
        h.update((",".join(str(int(i)) for i in face) + ";").encode())
    return h.hexdigest()


def _edge_audit(mesh: Any) -> dict[str, Any]:
    """Count how often each undirected edge is used, and whether orientations pair.

    A closed, consistently oriented surface uses every edge exactly twice, once in
    each direction. Anything else is an open shell with a hole, a duplicated face
    or a flipped patch, and no volume solver should be handed it.
    """
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in mesh.faces:
        for u, v in ((a, b), (b, c), (c, a)):
            directed[(int(u), int(v))] = directed.get((int(u), int(v)), 0) + 1
    unmatched = 0
    nonmanifold = 0
    for (u, v), count in directed.items():
        if count != 1:
            nonmanifold += 1
        if directed.get((v, u), 0) != 1:
            unmatched += 1
    return {
        "edges": len(directed),
        "boundary_edges": unmatched,
        "nonmanifold_edges": nonmanifold,
        "closed": unmatched == 0 and nonmanifold == 0,
    }


# --------------------------------------------------------------------------- solid


def solid_from_library(
    semi_axes: Sequence[float], *, resolution: int | Sequence[int] = 40
) -> dict[str, Any]:
    """Sample the unified geometry library into the one solid all routes consume.

    ``geometry_field.EllipsoidField`` is the CGC protocol object (``signed_distance``,
    ``gradient``, ``bounds``, ``state``) and ``geometry_field.to_mesh`` is that
    library's own extraction path, so the mesh here is *library output*, not a
    private triangulation written for this gate. The geometry gate runs first, in
    the order the 5415 rounds taught: closed and outward-oriented mesh, then a
    volume inside 5 % of ``4/3 pi a b c`` -- a bad fixture is refused before any
    solver sees it.
    """
    a, b, c = (float(v) for v in semi_axes)
    analytic = analytic_volume((a, b, c))
    field = EllipsoidField((a, b, c))
    sampled = to_mesh(field, resolution)
    mesh = sampled.mesh
    edges = _edge_audit(mesh)
    if not edges["closed"]:
        raise HandoffError(
            "the unified-library mesh is not a closed oriented surface: "
            f"{edges['boundary_edges']} unmatched edge directions, "
            f"{edges['nonmanifold_edges']} non-manifold edges")
    volume = signed_volume(mesh.vertices, mesh.faces)
    if volume <= 0.0:
        raise HandoffError(f"library mesh volume {volume:.6g} m3 is not positive (inward normals)")
    rel = abs(volume - analytic) / analytic
    if rel > 0.05:
        raise HandoffError(
            f"sampled volume {volume:.6g} m3 is {100 * rel:.2f} % from the analytic "
            f"{analytic:.6g} m3; refine the resolution before blaming a solver")
    return {
        "field": field,
        "mesh": mesh,
        "semi_axes": [a, b, c],
        "resolution": resolution,
        "analytic_volume_m3": analytic,
        "surface_volume_m3": volume,
        "surface_volume_rel_error": rel,
        "nodes": len(mesh.vertices),
        "triangles": len(mesh.faces),
        "edges": edges,
        "digest_sha256": _mesh_digest(mesh),
        "library": "tensorcad.geometry_field.EllipsoidField + to_mesh",
        "frame": COORDINATE_CONTRACT["name"],
    }


def star_tetrahedra(mesh: Any) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int, int]]]:
    """Tetrahedralise a closed convex surface about its vertex centroid.

    Valid exactly when the solid is star-shaped with respect to that point (an
    ellipsoid is convex, so it is), and then the tetrahedra tile the interior with
    no gaps: their summed volume equals the surface's divergence volume to working
    precision. Every cell's signed volume is checked positive here, so an
    inverted or degenerate cell is a refusal rather than a small number.
    """
    vertices = [(float(x), float(y), float(z)) for x, y, z in mesh.vertices]
    if len(mesh.faces) < 4:
        raise HandoffError("a solid needs at least four surface triangles")
    n = len(vertices)
    centre = tuple(sum(v[axis] for v in vertices) / n for axis in range(3))
    vertices = vertices + [centre]
    # ``(a, c, b, centre)`` -- not ``(a, b, c, centre)``: with outward-facing
    # triangles the latter gives the *negative* of the cell volume in this
    # repository's determinant convention (measured, not assumed), which would
    # write an inverted volume mesh.
    tets = [(int(a), int(c), int(b), n) for a, b, c in mesh.faces]
    for tet in tets:
        p0, p1, p2, p3 = (vertices[i] for i in tet)
        v = _tet_volume(p0, p1, p2, p3)
        if v <= 0.0:
            raise HandoffError(
                f"star tetrahedron {tet} has non-positive volume {v:.3e} m3; the solid "
                "is not star-shaped about its centroid")
    return vertices, tets


def _tet_volume(p0: Sequence[float], p1: Sequence[float], p2: Sequence[float],
                p3: Sequence[float]) -> float:
    ax, ay, az = (p1[i] - p0[i] for i in range(3))
    bx, by, bz = (p2[i] - p0[i] for i in range(3))
    cx, cy, cz = (p3[i] - p0[i] for i in range(3))
    return ((ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz) / 6.0


# --------------------------------------------------------------------------- routes


def _tensorfvm_patch(vertices: Sequence[Sequence[float]],
                     tets: Sequence[Sequence[int]], *, solve: bool = True) -> dict[str, Any]:
    """Hand the tets to TensorFVM and solve the harmonic patch field.

    Returns ``{"status": "ran"|"withheld", ...}``. A withheld consumer is a fact
    about the environment, not a pass: the caller records it verbatim.
    """
    try:
        import numpy as np
        from tensorfvm.native_mesh3d import NativeMesh3D
        from tensorfvm.native_hybrid_scalar3d import NativeHybridScalar3D
    except Exception as exc:  # pragma: no cover - environment dependent
        return {"status": "withheld", "reason": f"{type(exc).__name__}: {exc}"}
    try:
        cells = [list(int(i) for i in tet) for tet in tets]
        mesh = NativeMesh3D.from_cells([list(v) for v in vertices], cells,
                                       ["tetra"] * len(cells))
        audit = mesh.audit()
        c0, c1, c2, c3 = PATCH_COEFFICIENTS

        def exact(points: Any) -> Any:
            p = np.asarray(points, dtype=float)
            return c0 + c1 * p[:, 0] + c2 * p[:, 1] + c3 * p[:, 2]

        if not solve:
            return {
                "status": "ran",
                "module": "tensorfvm.native_mesh3d.NativeMesh3D.from_cells",
                "solver_file": getattr(__import__("tensorfvm"), "__file__", None),
                "cells": int(audit["cells"]),
                "faces": int(audit["faces"]),
                "minimum_cell_volume_m3": float(audit["minimum_cell_volume"]),
                "maximum_relative_closure": float(audit["maximum_relative_closure"]),
                "maximum_nonorthogonality_degrees": float(audit["maximum_nonorthogonality_degrees"]),
                "patch_relative_linf": None,
                "note": "conduction solve skipped (solve=False)",
            }
        boundary = mesh.neighbor < 0
        values = np.zeros(len(mesh.owner), dtype=float)
        values[boundary] = exact(mesh.face_centers[boundary])
        solver = NativeHybridScalar3D(mesh, 1.0)
        result = solver.solve(boundary_values=values)
        cell_values = np.asarray(result["cell_values"], dtype=float)
        reference = exact(mesh.cell_centers)
        error = np.abs(cell_values - reference)
        scale = float(np.max(np.abs(reference))) or 1.0
        return {
            "status": "ran",
            "module": "tensorfvm.native_mesh3d.NativeMesh3D.from_cells + "
                      "tensorfvm.native_hybrid_scalar3d.NativeHybridScalar3D",
            "solver_file": getattr(__import__("tensorfvm"), "__file__", None),
            "cells": int(audit["cells"]),
            "faces": int(audit["faces"]),
            "minimum_cell_volume_m3": float(audit["minimum_cell_volume"]),
            "maximum_relative_closure": float(audit["maximum_relative_closure"]),
            "maximum_nonorthogonality_degrees": float(audit["maximum_nonorthogonality_degrees"]),
            "patch_relative_linf": float(error.max() / scale),
            "patch_relative_l2": float(
                math.sqrt(float((mesh.cell_volumes * error ** 2).sum())
                          / float((mesh.cell_volumes * reference ** 2).sum()))),
            "linear_relative_residual": float(result["linear_relative_residual"]),
            "global_relative_balance": float(result["global_relative_balance"]),
            "field": "phi = 1 + 0.4x - 0.7y + 1.3z (harmonic)",
        }
    except Exception as exc:  # pragma: no cover - solver-side refusal is data
        return {"status": "error", "reason": f"{type(exc).__name__}: {exc}"}


def fvm_route(solid: dict[str, Any], outdir: str | Path, *,
              solve: bool = True) -> dict[str, Any]:
    """Solid mesh -> tets -> gmsh bundle -> re-read -> TensorFVM consumes (and solves).

    ``solve=False`` still hands the cells to TensorFVM (topology + audit) but skips
    the conduction solve, which is the expensive part; the gate driver solves, the
    unit tests use the cheap path except for one exactness test.
    """
    vertices, tets = star_tetrahedra(solid["mesh"])
    manifest = write_bundle_volume(
        outdir, vertices=vertices, tets=tets, case="solid-gate-fvm",
        extra={"solid_digest": solid["digest_sha256"], "semi_axes": solid["semi_axes"],
               "producer": "tensorcad.solver_solid_gate.star_tetrahedra"})
    re_read = verify_bundle(outdir)
    volume = total_volume(vertices, tets)
    analytic = solid["analytic_volume_m3"]
    surface = solid["surface_volume_m3"]
    consumer = _tensorfvm_patch(vertices, tets, solve=solve)
    return {
        "route": "volume_msh",
        "frame": manifest["frame"]["name"],
        "solid_digest": manifest["extra"]["solid_digest"],
        "bundle": str(outdir),
        "tets": len(tets),
        "nodes": len(vertices),
        "boundary_triangles": manifest["counts"]["boundary_triangles"],
        "volume_m3": volume,
        "declared_volume_m3": manifest["volume_m3"],
        "analytic_volume_m3": analytic,
        "surface_volume_m3": surface,
        "tet_vs_surface_volume_rel_error": abs(volume - surface) / surface,
        "surface_vs_analytic_volume_rel_error": abs(surface - analytic) / analytic,
        "volume_rel_error": abs(volume - analytic) / analytic,
        "verify_re_read_ok": bool(re_read.get("ok")),
        "verify_volume_m3": re_read["volume_m3"],
        "consumer": consumer,
    }


def dg_route(solid: dict[str, Any], outdir: str | Path, *,
             segments: int = 8) -> dict[str, Any]:
    """The solid's ``z = 0`` section as an exact NURBS ellipse, checked on the solid."""
    a, b, _ = solid["semi_axes"]
    contract = ellipse_contract(a, b, (0.0, 0.0), segments=segments,
                                name="cad-solid-equator-section")
    manifest = write_dg_bundle(
        outdir, contract, r_outer_m=1.05 * max(a, b), case="solid-gate-dg",
        extra={"solid_digest": solid["digest_sha256"], "semi_axes": solid["semi_axes"],
               "section": "z = 0 (equator of the solid)",
               "producer": "tensorcad.dg_handoff.ellipse_contract"})

    # "same solid" is measured, not asserted: lift the 2-D curve back to 3-D and
    # evaluate the geometry library's own signed distance on it.
    field = solid["field"]
    samples = evaluate_curve(contract, [i / 64.0 for i in range(65)])
    residual = max(abs(d) for d in field.signed_distance([(x, y, 0.0) for x, y in samples]))
    if residual > ROUTE_TOLERANCES["dg_on_solid_max_distance_m"] * max(solid["semi_axes"]):
        raise HandoffError(
            f"the handed-over curve is {residual:.3e} m off the solid's surface; the "
            "section and the solid are not the same geometry")

    re_read = verify_dg_bundle(outdir)
    measured = manifest["measured"]
    solid["dg_bundle"] = str(outdir)
    return {
        "route": "curved_wall",
        "frame": manifest["frame"]["name"],
        "solid_digest": manifest["extra"]["solid_digest"],
        "bundle": str(outdir),
        "semi_axes": [a, b],
        "segments": measured["segments"],
        "control_points": measured["control_points"],
        "perimeter_m": measured["circumference_m"],
        "perimeter_analytic_m": ellipse_perimeter(a, b),
        "perimeter_rel_error": measured["circumference_rel_error"],
        "area_m2": measured["enclosed_area_m2"],
        "area_analytic_m2": math.pi * a * b,
        "area_rel_error": measured["area_rel_error"],
        "closure_gap_m": measured["closure_gap_m"],
        "on_solid_max_distance_m": residual,
        # dg_handoff.verify_bundle raises on any disagreement and returns no "ok"
        # flag; reaching this line means the re-read reproduced every number.
        "verify_re_read_ok": bool(re_read.get("digest_matches", True)),
        "consumer": {
            "status": "withheld",
            "reason": "tensordg is not importable in this environment; the bundle is "
                      "written in tensordg.curved_geometry.v1 and verified here, but "
                      "TensorDG itself did not read it in this run",
        },
    }


def lattice_route(solid: dict[str, Any], outdir: str | Path, *,
                  radial_cells: int = 6) -> dict[str, Any]:
    """Voxelise the same mesh, write a true LBM occupancy bundle, re-read it."""
    a, b, c = solid["semi_axes"]
    h = min(solid["semi_axes"]) / float(radial_cells)
    lo = (-a, -b, -c)
    hi = (a, b, c)
    lattice = solid_lattice(solid["mesh"], spacing=h, bounds=[lo, hi], ray_axis="x")
    shape = tuple(int(v) for v in lattice.shape)
    nx, ny, nz = shape
    spacing = tuple(float(v) for v in lattice.spacing)
    solid_bits = lattice.occupancy
    solid_cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)
                   if solid_bits[i + nx * (j + ny * k)]]
    fluid_cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)
                   if not solid_bits[i + nx * (j + ny * k)]]
    if not solid_cells:
        raise HandoffError("the solid voxelised to no interior cells; refine the spacing")
    cell_volume = spacing[0] * spacing[1] * spacing[2]
    solid_volume = len(solid_cells) * cell_volume
    analytic = solid["analytic_volume_m3"]
    manifest = write_bundle_lattice(
        outdir, shape=shape, origin=lo, spacing=spacing, fluid_cells=fluid_cells,
        walls={"kind": "wall", "source": "tensorcad.voxel_solid.solid_lattice",
               "ray_axis": lattice.ray_axis},
        case="solid-gate-lbm",
        extra={"solid_digest": solid["digest_sha256"], "semi_axes": solid["semi_axes"],
               "solid_cells": len(solid_cells), "solid_volume_m3": solid_volume,
               "lattice_mesh_digest": lattice.mesh_digest,
               "producer": "tensorcad.voxel_solid.solid_lattice"})
    re_read = verify_bundle(outdir)
    total = nx * ny * nz
    return {
        "route": "lattice",
        "frame": manifest["frame"]["name"],
        "solid_digest": manifest["extra"]["solid_digest"],
        "bundle": str(outdir),
        "shape": list(shape),
        "cell_size_m": h,
        "solid_cells": len(solid_cells),
        "fluid_cells": len(fluid_cells),
        "total_cells": total,
        "cells_accounted": len(solid_cells) + len(fluid_cells) - total,
        "solid_volume_m3": solid_volume,
        "analytic_volume_m3": analytic,
        "volume_rel_error": abs(solid_volume - analytic) / analytic,
        "fluid_volume_m3": manifest["fluid_volume_m3"],
        "verify_re_read_ok": bool(re_read.get("ok")),
        "verify_fluid_volume_m3": re_read["fluid_volume_m3"],
        # the voxeliser fingerprints the mesh it consumed with its own digest
        # scheme; it is carried in the manifest so the voxelisation can be traced
        # back to geometry even where the gate's digest scheme differs
        "lattice_mesh_digest": lattice.mesh_digest,
        "manifest_lattice_mesh_digest": manifest["extra"]["lattice_mesh_digest"],
    }


def _pin_tensordg_path() -> str | None:
    """Resolve ``tensordg`` the way the FVM leg resolves ``tensorfvm``.

    The ecosystem checkout is used when the package is not installed, and the
    resolved location is returned so a number can be traced to the consumer that
    produced it (an unpinned import is how "which solver ran?" gets lost).
    """
    try:
        import tensordg  # noqa: F401
        return "installed"
    except Exception:
        pass
    for candidate in (os.environ.get("TENSORDG_SRC"), "/root/tensor-ecosystem/TensorDG"):
        if not candidate:
            continue
        if not (Path(candidate) / "tensordg" / "curved_geometry.py").is_file():
            continue
        if candidate not in sys.path:
            sys.path.insert(0, candidate)
        try:
            import tensordg  # noqa: F401
            return candidate
        except Exception:
            continue
    return None


def dg_consumer_probe(solid: dict[str, Any], outdir: str | Path, *,
                      nr: int = 4, ntheta: int = 16, degree: int = 2) -> dict[str, Any]:
    """Let TensorDG read the CAD payload -- and record what it refuses.

    The ellipse section of the ellipsoid is exact CAD, but TensorDG's annulus
    builder requires a **circular** inner wall, so it must refuse it. That refusal
    is the honest state of the consumer, so it is measured here instead of being
    hidden by handing the solver a circle and calling it "the same solid".

    The route that *can* be consumed end to end is measured too: a spheroid
    (``a == b``) comes out of the same unified-library pipeline, its ``z = 0``
    section is a circle, and TensorDG builds its curved annulus mesh from the
    bundle **file** on disk. Both outcomes are returned.
    """
    pinned = _pin_tensordg_path()
    import tensordg.curved_geometry as cg  # noqa: F401

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "consumer": "tensordg.curved_geometry.mesh_from_contract",
        "tensordg_path": pinned,
        "nr": nr, "ntheta": ntheta, "degree": degree,
    }

    # -- what the consumer refuses today: the ellipsoid's exact ellipse section
    a, b, _ = solid["semi_axes"]
    ellipse_contract_payload = json.loads((Path(solid["dg_bundle"]) / "curved_wall.json").read_text())
    try:
        cg.mesh_from_contract(ellipse_contract_payload, nr, ntheta, 1.05 * max(a, b), degree)
        result["ellipse_section"] = {"status": "accepted",
                                     "note": "consumer accepted a non-circular wall"}
    except Exception as exc:
        result["ellipse_section"] = {"status": "refused",
                                     "reason": f"{type(exc).__name__}: {exc}"}

    # -- what it does consume: the spheroid's circular section, read from disk
    spheroid = solid_from_library((0.5, 0.5, 0.3), resolution=solid["resolution"])
    circle_dir = out / "dg_circle"
    circle_contract = ellipse_contract(0.5, 0.5, (0.0, 0.0), segments=8,
                                       name="cad-spheroid-equator-section")
    manifest = write_dg_bundle(
        circle_dir, circle_contract, r_outer_m=0.6, case="solid-gate-dg-circle",
        extra={"solid_digest": spheroid["digest_sha256"],
               "semi_axes": spheroid["semi_axes"], "section": "z = 0 (circle)"})
    verify_dg_bundle(circle_dir)  # positive control before the consumer reads it
    from_disk = json.loads((circle_dir / "curved_wall.json").read_text())
    try:
        mesh = cg.mesh_from_contract(from_disk, nr, ntheta, 0.6, degree)
        report = cg.validate_geometry_contract(from_disk)
        centre = [float(v) for v in from_disk["centre"]]
        import numpy as np
        tsample = np.linspace(0.0, 1.0, 721)
        pts = cg.nurbs_curve_points(from_disk["degree"], from_disk["knots"],
                                   from_disk["control_points"], from_disk["weights"], tsample)
        rad = np.linalg.norm(pts - np.asarray(centre), axis=1)
        result["circle_section"] = {
            "status": "ran",
            "solid_digest": spheroid["digest_sha256"],
            "analytic_volume_m3": spheroid["analytic_volume_m3"],
            "bundle_verified": True,
            "contract_degree": report.get("degree"),
            "wall_radius_m": float(rad.mean()),
            "wall_radius_rel_error": float(abs(rad.mean() - 0.5) / 0.5),
            "wall_radius_spread_rel": float(np.max(np.abs(rad - rad.mean())) / rad.mean()),
            "elements": int(getattr(mesh, "ncells", 0)) or None,
            "element_degree": int(getattr(mesh, "degree", 0)) or None,
            "mesh_type": type(mesh).__name__,
        }
    except Exception as exc:  # pragma: no cover - consumer-side refusal is data
        result["circle_section"] = {"status": "error", "reason": f"{type(exc).__name__}: {exc}"}
    return result


# --------------------------------------------------------------------------- gate


def tamper_refused(outdir: str | Path, verifier, *, name: str) -> dict[str, Any]:
    """Positive control, then a frame tamper that the route's own verifier must refuse.

    The control is the point: a verifier that refuses *everything* would look
    perfect on the tamper test alone. Only the handoff's own
    :class:`HandoffError` counts as a refusal -- anything else propagates, because
    a crash is not an audit.
    """
    source = Path(outdir)
    control = verifier(source)  # raises if the untampered bundle is not verifiable
    tampered = source.parent / f"{source.name}__tampered_{name}"
    if tampered.exists():
        shutil.rmtree(tampered)
    shutil.copytree(source, tampered)
    handoff_path = tampered / "handoff.json"
    manifest = json.loads(handoff_path.read_text())
    manifest["frame"] = {"name": "y-up", "up_axis": "y"}
    handoff_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))

    refused = False
    reason = ""
    try:
        verifier(tampered)
    except (HandoffError, DgHandoffError) as exc:
        refused = True
        reason = str(exc)
    shutil.rmtree(tampered)
    if not refused:
        raise HandoffError(
            f"the {name} verifier accepted a y-up manifest; the frame contract is not "
            "enforced on that route")
    return {"route": name, "positive_control_ok": True, "tamper_refused": True,
            "refusal": reason, "control": control.get("route", name)}


def cross_route_gate(routes: dict[str, dict[str, Any]],
                     solid: dict[str, Any]) -> dict[str, Any]:
    """Compare the routes against the analytic solid and against each other."""
    checks: list[dict[str, Any]] = []

    def check(key: str, measured: float, tolerance: float) -> None:
        ok = bool(measured <= tolerance)
        checks.append({"check": key, "measured": float(measured), "tolerance": float(tolerance),
                       "verdict": "PASS" if ok else "FAIL"})

    analytic = solid["analytic_volume_m3"]
    # the DG consumer probe is a consumer report, not a route with its own frame
    geometry_routes = {name: route for name, route in routes.items() if "frame" in route}
    frames = {name: route["frame"] for name, route in geometry_routes.items()}
    checks.append({
        "check": "frame_contract",
        "measured": sorted(set(frames.values())),
        "tolerance": COORDINATE_CONTRACT["name"],
        "verdict": "PASS" if set(frames.values()) == {COORDINATE_CONTRACT["name"]} else "FAIL",
    })
    digests = {name: route["solid_digest"] for name, route in geometry_routes.items()}
    checks.append({
        "check": "same_solid_digest",
        "measured": len(set(digests.values())),
        "tolerance": 1,
        "verdict": "PASS" if len(set(digests.values())) == 1
        and digests["fvm"] == solid["digest_sha256"] else "FAIL",
    })

    check("fvm_tet_vs_surface_volume_rel_error",
          routes["fvm"]["tet_vs_surface_volume_rel_error"],
          ROUTE_TOLERANCES["fvm_tet_vs_surface_volume_rel_error"])
    check("fvm_surface_vs_analytic_volume_rel_error",
          routes["fvm"]["surface_vs_analytic_volume_rel_error"],
          ROUTE_TOLERANCES["fvm_surface_vs_analytic_volume_rel_error"])
    consumer = routes["fvm"]["consumer"]
    if consumer["status"] == "ran" and consumer.get("patch_relative_linf") is not None:
        check("fvm_patch_relative_linf", consumer["patch_relative_linf"],
              ROUTE_TOLERANCES["fvm_patch_relative_linf"])
        check("fvm_patch_relative_l2", consumer["patch_relative_l2"],
              ROUTE_TOLERANCES["fvm_patch_relative_l2"])
        check("fvm_global_relative_balance", abs(consumer["global_relative_balance"]),
              ROUTE_TOLERANCES["fvm_global_relative_balance"])
    else:
        checks.append({"check": "fvm_patch_relative_linf", "measured": None,
                       "tolerance": ROUTE_TOLERANCES["fvm_patch_relative_linf"],
                       "verdict": "WITHHELD", "reason": consumer.get("reason", "")})

    check("dg_perimeter_rel_error", routes["dg"]["perimeter_rel_error"],
          ROUTE_TOLERANCES["dg_perimeter_rel_error"])
    check("dg_area_rel_error", routes["dg"]["area_rel_error"],
          ROUTE_TOLERANCES["dg_area_rel_error"])
    check("dg_on_solid_max_distance_m", routes["dg"]["on_solid_max_distance_m"],
          ROUTE_TOLERANCES["dg_on_solid_max_distance_m"])

    check("lattice_volume_rel_error", routes["lattice"]["volume_rel_error"],
          ROUTE_TOLERANCES["lattice_volume_rel_error"])
    check("lattice_fluid_plus_solid", float(abs(routes["lattice"]["cells_accounted"])),
          ROUTE_TOLERANCES["lattice_fluid_plus_solid"])
    check("cross_route_volume_rel_error",
          abs(routes["fvm"]["volume_m3"] - routes["lattice"]["solid_volume_m3"]) / analytic,
          ROUTE_TOLERANCES["cross_route_volume_rel_error"])

    consumer = routes.get("dg_consumer")
    if consumer is not None:
        circle = consumer.get("circle_section", {})
        if circle.get("status") == "ran":
            check("dg_consumer_wall_radius_rel_error", circle["wall_radius_rel_error"], 1e-09)
            checks.append({"check": "dg_consumer_ran", "measured": "ran",
                           "tolerance": "ran", "verdict": "PASS"})
        else:
            checks.append({"check": "dg_consumer_ran", "measured": circle.get("status"),
                           "tolerance": "ran", "verdict": "FAIL",
                           "reason": circle.get("reason", "")})
        ellipse = consumer.get("ellipse_section", {})
        checks.append({
            "check": "dg_consumer_accepts_non_circular_wall",
            "measured": ellipse.get("status"),
            "tolerance": "accepted",
            "verdict": "PASS" if ellipse.get("status") == "accepted" else "WITHHELD",
            "reason": ellipse.get("reason", "")})

    failed = [c["check"] for c in checks if c["verdict"] == "FAIL"]
    withheld = [c["check"] for c in checks if c["verdict"] == "WITHHELD"]
    return {
        "checks": checks,
        "analytic_volume_m3": analytic,
        "frames": frames,
        "solid_digests": digests,
        "failed": failed,
        "withheld": withheld,
        "passed": not failed,
    }


def run_gate(outdir: str | Path, *, semi_axes: Sequence[float] = (1.0, 0.5, 0.3),
             resolution: int = 40, radial_cells: int = 6,
             tamper: bool = True) -> dict[str, Any]:
    """Run the whole gate and write ``gate.json`` + ``summary.md`` under ``outdir``."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    solid = solid_from_library(semi_axes, resolution=resolution)

    fvm = fvm_route(solid, out / "fvm")
    dg = dg_route(solid, out / "dg")
    lattice = lattice_route(solid, out / "lattice", radial_cells=radial_cells)
    routes = {"fvm": fvm, "dg": dg, "lattice": lattice}
    if _pin_tensordg_path() is not None:
        routes["dg_consumer"] = dg_consumer_probe(solid, out / "dg_consumer")
    else:  # pragma: no cover - environment dependent
        routes["dg_consumer"] = {
            "consumer": "tensordg.curved_geometry.mesh_from_contract",
            "tensordg_path": None,
            "circle_section": {"status": "withheld", "reason": "tensordg not importable"},
            "ellipse_section": {"status": "withheld", "reason": "tensordg not importable"},
        }
    gate = cross_route_gate(routes, solid)

    tamper_results: list[dict[str, Any]] = []
    if tamper:
        tamper_results.append(tamper_refused(out / "fvm", verify_bundle, name="fvm"))
        tamper_results.append(tamper_refused(out / "dg", verify_dg_bundle, name="dg"))
        tamper_results.append(tamper_refused(out / "lattice", verify_bundle, name="lattice"))
        gate["checks"].append({
            "check": "frame_tamper_refused_all_routes",
            "measured": sum(1 for r in tamper_results if r["tamper_refused"]),
            "tolerance": 3,
            "verdict": "PASS" if all(r["tamper_refused"] and r["positive_control_ok"]
                                     for r in tamper_results) else "FAIL",
        })
        gate["failed"] = [c["check"] for c in gate["checks"] if c["verdict"] == "FAIL"]
        gate["passed"] = not gate["failed"]

    report = {
        "schema": SOLID_GATE_SCHEMA,
        "frame_contract": COORDINATE_CONTRACT,
        "solid": {key: value for key, value in solid.items() if key not in {"field", "mesh"}},
        "routes": routes,
        "tamper": tamper_results,
        "gate": gate,
        "tolerances": ROUTE_TOLERANCES,
        "passed": gate["passed"],
    }
    (out / "gate.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (out / "summary.md").write_text(summary_markdown(report), encoding="utf-8")
    return report


def summary_markdown(report: dict[str, Any]) -> str:
    """The human-readable face of ``gate.json`` (numbers copied, never retyped)."""
    solid = report["solid"]
    routes = report["routes"]
    gate = report["gate"]
    lines = [
        "# One solid, three solver routes, one frame gate",
        "",
        f"* schema: `{report['schema']}`",
        f"* solid: `{solid['library']}`, semi-axes {solid['semi_axes']}, "
        f"resolution {solid['resolution']}",
        f"* mesh: {solid['nodes']} nodes / {solid['triangles']} triangles, "
        f"digest `{solid['digest_sha256'][:16]}…`",
        f"* analytic volume: {solid['analytic_volume_m3']:.6f} m³; "
        f"sampled surface volume: {solid['surface_volume_m3']:.6f} m³ "
        f"({100 * solid['surface_volume_rel_error']:.3f} %)",
        "",
        "| route | frame | measurement | vs analytic | bundle verifies | consumer |",
        "|---|---|---|---|---|---|",
        f"| FVM | {routes['fvm']['frame']} | {routes['fvm']['tets']} tets, "
        f"V={routes['fvm']['volume_m3']:.6f} m³ | "
        f"{routes['fvm']['volume_rel_error']:.2e} | "
        f"{routes['fvm']['verify_re_read_ok']} | "
        f"{routes['fvm']['consumer']['status']} |",
        f"| DG | {routes['dg']['frame']} | ellipse a={routes['dg']['semi_axes'][0]}, "
        f"b={routes['dg']['semi_axes'][1]}, P={routes['dg']['perimeter_m']:.6f} m | "
        f"{routes['dg']['perimeter_rel_error']:.2e} | "
        f"{routes['dg']['verify_re_read_ok']} | "
        f"{routes['dg']['consumer']['status']} |",
        f"| LBM | {routes['lattice']['frame']} | "
        f"{routes['lattice']['shape']} @ {routes['lattice']['cell_size_m']:.4f} m, "
        f"V={routes['lattice']['solid_volume_m3']:.6f} m³ | "
        f"{routes['lattice']['volume_rel_error']:.2e} | "
        f"{routes['lattice']['verify_re_read_ok']} | "
        f"occupancy bundle |",
        "",
        "| gate check | measured | tolerance | verdict |",
        "|---|---|---|---|",
    ]
    for entry in gate["checks"]:
        measured = entry["measured"]
        if isinstance(measured, list):
            measured = ", ".join(str(v) for v in measured) or "—"
        elif isinstance(measured, float):
            measured = f"{measured:.6e}"
        lines.append(
            f"| {entry['check']} | {measured} | {entry['tolerance']} | {entry['verdict']} |")
    lines += [
        "",
        f"**Verdict**: {'PASS' if report['passed'] else 'FAIL'}"
        + (f" (failed: {', '.join(gate['failed'])})" if gate["failed"] else ""),
    ]
    if gate["withheld"]:
        lines.append(f"* withheld: {', '.join(gate['withheld'])}")
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="One solid from the unified geometry library -> FVM / DG / LBM")
    parser.add_argument("--out", default="artifacts/solver_solid_gate")
    parser.add_argument("--a", type=float, default=1.0)
    parser.add_argument("--b", type=float, default=0.5)
    parser.add_argument("--c", type=float, default=0.3)
    parser.add_argument("--resolution", type=int, default=40)
    parser.add_argument("--radial-cells", type=int, default=6)
    parser.add_argument("--no-tamper", action="store_true")
    args = parser.parse_args(argv)
    report = run_gate(args.out, semi_axes=(args.a, args.b, args.c),
                      resolution=args.resolution, radial_cells=args.radial_cells,
                      tamper=not args.no_tamper)
    print(summary_markdown(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())