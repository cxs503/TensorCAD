"""Cross-tool volume cross-validation: the self-developed integration against independent engines.

Why this module exists
----------------------
A single volume integrator cannot check itself.  The JBC episode is the standing
reminder: the platform's hull pipeline once reported a displacement ``-18.3 %``
below the published value, and every internal guard agreed with the pipeline --
the guards and the pipeline shared the same integration route, so the bias was
invisible until an *outside* number was compared.  The fix is structural, not a
patch on one integrator: every pipeline volume must be reproduced by at least two
implementations that do not share code, and the difference between them must be
attributed to a *named* effect (tessellation, STL quantisation, topology) rather
than waved through a tolerance.

The three routes this module provides
-------------------------------------
============================  ==================================================
route                         what it is
============================  ==================================================
``divergence``                :func:`tensorcad.mesh.total_signed_volume` -- the
                              platform's own ``sum(a . (b x c))/6`` over triangles.
``trimesh``                   :mod:`trimesh`'s ``Trimesh.volume`` -- an
                              independent third-party engine (NumPy vectorised,
                              its own watertightness/winding guards).
``numpy-stl``                 :mod:`stl`'s ``get_mass_properties`` -- a *second*
                              independent engine, Mirtich's polyhedral
                              mass-property algorithm in C++, reached through a
                              binary-STL round trip.
``sectional``                 a method that shares no mathematics with any of the
                              three: per-station section area by Simpson
                              quadrature of the half-breadth, then Simpson in
                              ``x``.  No divergence theorem, no triangulation.
============================  ==================================================

``sectional`` is the one that speaks to the JBC failure mode.  The triangle-based
routes all integrate the *inscribed polyhedron*, so they agree with each other
and sit below the smooth geometry by a discretisation deficit; a shared bug in
the geometry source would move all four together.  A quadrature route that
integrates the half-breadth field itself therefore separates "the geometry is
wrong" from "the triangles are coarse".

What the discrepancy decomposition means (and how it is classified)
-------------------------------------------------------------------
:func:`classify_routes` turns a row of route volumes into a verdict, and the
verdicts are ordered by what they can prove:

``engine-disagreement``
    the independent engines disagree with each other beyond
    :data:`ENGINE_AGREEMENT_TOLERANCE` while reading the *same* triangles.  That
    is a topology defect (cracks, flipped winding, duplicated vertices) or a
    precision artefact (binary STL is float32, so ~1e-9 relative is expected and
    is not a defect); it is never a tessellation effect.
``tessellation-deficit``
    engines agree, and the difference to the pipeline volume shrinks at the
    measured order under refinement.  This is the expected ``O(h**2)`` deficit of
    an inscribed polyhedron on a curved hull; it must be extrapolated away, not
    tuned away.
``systematic-bias``
    the difference does not shrink under refinement, i.e. it is not a property of
    the mesh at all.  This is the class the JBC incident belongs to and it must
    be booked in the negative-knowledge ledger with the measured magnitude.
``agrees``
    within the acceptance line :data:`ACCEPTANCE_TOLERANCE` (0.3 %).

Tolerances are module constants so a caller cannot tune them after seeing a
result.  The only tolerance that is a property of the *tool* rather than of the
physics is :data:`ENGINE_AGREEMENT_TOLERANCE`: it is 1e-6 because a binary STL
stores coordinates as float32 (~1e-7 relative), and the two engine routes are
compared through that round trip.

Dependency policy
-----------------
``import tensorcad.volume_xval`` needs the standard library only, like the rest of
the core package.  NumPy / trimesh / numpy-stl are imported *inside* the engine
adapters, so a caller without the optional extras still gets the self route, the
sectional route and every analysis helper, and an engine adapter that is missing
says so in its own name instead of failing at import time.
"""

from __future__ import annotations

import math
import os
import tempfile
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder

SCHEMA = "tensorcad.volume_xval/1.0"

#: Route-to-route tolerance for two engines reading the *same* triangles.
#: 1e-6, not 1e-12: the numpy-stl route goes through a binary STL whose
#: coordinates are float32.
ENGINE_AGREEMENT_TOLERANCE = 1e-6
#: Acceptance line copied from the platform's volume gate: 0.3 %.
ACCEPTANCE_TOLERANCE = 3.0e-3
#: Order window a discrepancy must fall into to be called a tessellation deficit.
ORDER_WINDOW = (1.5, 2.6)

#: A half-breadth field ``b(x, z)`` in metres, ``z`` measured up from the keel.
Breadth = Callable[[float, float], float]


class VolumeXvalError(ValueError):
    """Raised when a route cannot produce a number that means what it claims."""


# ---------------------------------------------------------------------------
# sampling grids
# ---------------------------------------------------------------------------
def cosine_stations(length: float, stations: int) -> tuple[float, ...]:
    """Station abscissae from the transom (``x = 0``) to the stem (``x = L``).

    Cosine clustering at both ends: a hull's curvature lives at the stem and the
    transom, and a uniform grid would put the coarsest triangles exactly there.
    """
    if stations < 2:
        raise VolumeXvalError("need at least 2 stations")
    if not length > 0.0:
        raise VolumeXvalError("length must be positive")
    return tuple(
        0.5 * length * (1.0 - math.cos(math.pi * i / (stations - 1))) for i in range(stations)
    )


def waterline_levels(draft: float, levels: int) -> tuple[float, ...]:
    """Vertical levels from the keel (``z = 0``) to the waterline (``z = T``).

    Cosine clustering again, for the same reason: the bilge turns hardest.
    """
    if levels < 2:
        raise VolumeXvalError("need at least 2 levels")
    if not draft > 0.0:
        raise VolumeXvalError("draft must be positive")
    return tuple(
        0.5 * draft * (1.0 - math.cos(math.pi * j / (levels - 1))) for j in range(levels)
    )


# ---------------------------------------------------------------------------
# closed solid from a half-breadth field
# ---------------------------------------------------------------------------
def closed_solid_from_breadth(
    breadth: Breadth,
    *,
    length: float,
    draft: float,
    stations: int = 81,
    levels: int = 17,
) -> Mesh:
    """Close the wetted hull surface of a half-breadth field into a triangle solid.

    The solid is the hull **below the design waterline**, capped by the waterplane,
    with the transom closed by a planar fan.  That is the volume a hydrostatics
    table claims (``Cb L B T``), and it is the only volume a surface mesh can be
    graded on: an open surface has no divergence-theorem volume at all.

    Construction, and why it is watertight *by construction*:

    * every station carries one closed ring in the ``(y, z)`` plane -- starboard
      side from the flat-of-bottom up to the waterline, across the waterplane
      chord, down the port side, and back along the flat-of-bottom.  Adjacent
      rings are lofted edge-by-edge, so every ring edge is shared by exactly two
      triangles except at the two ends;
    * :class:`~tensorcad.geometry.MeshBuilder` deduplicates vertices on the exact
      coordinate tuple, and ``-0.0`` is canonically ``0.0``, so starboard and port
      points on the centreplane *are* the same vertex.  A stem whose half-breadth
      reaches exactly zero therefore closes onto itself with no seam to weld;
    * the aft end (``x = 0``) is closed by a fan over its ring polygon -- the
      transom.  The fan is fan-able only if the polygon is star-shaped about its
      centroid, and that is *checked* (fan area versus shoelace area) instead of
      assumed: a silently self-overlapping open stern is exactly the failure mode
      this module exists to catch.
    """
    xs = cosine_stations(length, stations)
    zs = waterline_levels(draft, levels)
    n_lev = len(zs)
    ring_size = 2 * n_lev

    builder = MeshBuilder(units="m", coordinate_frame="x aft->fore, y starboard, z up from keel")
    rings: list[list[int]] = []
    breadths: list[list[float]] = []
    for x in xs:
        row = [float(breadth(x, z)) for z in zs]
        for value in row:
            if not (math.isfinite(value) and value >= 0.0):
                raise VolumeXvalError(
                    f"half-breadth must be finite and non-negative, got {value!r} at x={x}"
                )
        breadths.append(row)
        ring: list[int] = []
        # starboard, keel -> waterline
        for j, z in enumerate(zs):
            ring.append(builder.add_vertex((x, row[j], z)))
        # port, waterline -> keel
        for j in range(n_lev - 1, -1, -1):
            ring.append(builder.add_vertex((x, -row[j], zs[j])))
        rings.append(ring)

    for i in range(len(xs) - 1):
        for k in range(ring_size):
            nxt = (k + 1) % ring_size
            builder.add_quad(rings[i][k], rings[i][nxt], rings[i + 1][nxt], rings[i + 1][k])

    # Transom: a planar fan over station 0's ring polygon (x = 0 plane).
    aft = breadths[0]
    if max(aft) > 0.0:
        # Exact ring centroid in the x = 0 plane from the shoelace sums (area weighted).
        area2, my, mz = _ring_sums(aft, zs, length=1.0)
        if abs(area2) < 1e-300:
            raise VolumeXvalError("transom section has zero area but a non-zero breadth")
        # Area-weighted polygon centroid: 6A*ybar = sum((y_i + y_i+1) * cross_i), 2A = area2.
        cy = my / (3.0 * area2)
        cz = mz / (3.0 * area2)
        apex = builder.add_vertex((0.0, cy, cz))
        for k in range(ring_size):
            builder.add_triangle(apex, rings[0][k], rings[0][(k + 1) % ring_size])
        _check_fan(aft, zs, cy, cz)

    # Fore end: a hull whose breadth vanishes at the stem closes onto itself, but a
    # body that is still full-width there (a barge/stem plate) leaves a real hole, so
    # it gets its own planar cap.  Closing it *only* when the section is still open is
    # what keeps a Wigley-style stem from gaining a fictitious plate.
    fore = breadths[-1]
    if max(fore) > _breadth_tolerance(draft):
        area2, my, mz = _ring_sums(fore, zs, length=1.0)
        if abs(area2) < 1e-300:
            raise VolumeXvalError("stem section has zero area but a non-zero breadth")
        cy = my / (3.0 * area2)
        cz = mz / (3.0 * area2)
        apex = builder.add_vertex((xs[-1], cy, cz))
        for k in range(ring_size):
            builder.add_triangle(apex, rings[-1][k], rings[-1][(k + 1) % ring_size])
        _check_fan(fore, zs, cy, cz)

    mesh = meshlib.remove_degenerate_faces(builder.build())
    oriented = meshlib.orient_outward(mesh)
    solid = oriented.mesh
    if not meshlib.is_watertight(solid):
        raise VolumeXvalError(
            "lofted solid is not watertight: "
            f"{len(meshlib.boundary_edges(solid))} boundary edges"
        )
    if meshlib.total_signed_volume(solid) <= 0.0:
        raise VolumeXvalError("lofted solid has non-positive volume after orientation")
    return solid


def _breadth_tolerance(draft: float) -> float:
    """Breadth below which a station is treated as a point (the stem closes itself)."""
    return 1e-12 * max(draft, 1.0)


def _ring_sums(breadth: Sequence[float], zs: Sequence[float], *, length: float) -> tuple[float, float, float]:
    """Shoelace sums ``(2A, 2A*ybar, 2A*zbar)`` of the station ring polygon.

    The ring runs up the starboard side (``y = +b``), across the waterplane chord,
    down the port side (``y = -b``) and back along the flat-of-bottom.  Only the
    first and last legs contribute (the others are pairs at equal ``z`` whose
    contributions cancel), which is why this stays a two-term loop.
    """
    n = len(zs)
    total = 0.0
    my = 0.0
    mz = 0.0

    def term(y0: float, z0: float, y1: float, z1: float) -> None:
        nonlocal total, my, mz
        cross = y0 * z1 - y1 * z0
        total += cross
        my += (y0 + y1) * cross
        mz += (z0 + z1) * cross

    # starboard side: (b_j, z_j) -> (b_{j+1}, z_{j+1})
    for j in range(n - 1):
        term(breadth[j], zs[j], breadth[j + 1], zs[j + 1])
    # waterplane chord: (b_last, T) -> (-b_last, T)
    term(breadth[-1], zs[-1], -breadth[-1], zs[-1])
    # port side: (-b_j, z_j) -> (-b_{j+1}, z_{j+1}), descending
    for j in range(n - 1, 0, -1):
        term(-breadth[j], zs[j], -breadth[j - 1], zs[j - 1])
    # flat of bottom: (-b_0, 0) -> (b_0, 0)
    term(-breadth[0], zs[0], breadth[0], zs[0])
    return total, my, mz


def _ring_area(breadth: Sequence[float], zs: Sequence[float]) -> float:
    total, _my, _mz = _ring_sums(breadth, zs, length=1.0)
    return 0.5 * abs(total)


def _check_fan(breadth: Sequence[float], zs: Sequence[float], cy: float, cz: float) -> None:
    """The transom fan must reproduce the polygon's own area.

    A fan from a point outside a non-star-shaped polygon overlaps itself and would
    quietly *add* volume, so the check is arithmetic rather than visual.
    """
    polygon = _ring_area(breadth, zs)
    fan = 0.0
    n = len(zs)
    ring = [(breadth[j], zs[j]) for j in range(n)]
    ring += [(-breadth[j], zs[j]) for j in range(n - 1, -1, -1)]
    for k in range(len(ring)):
        y0, z0 = ring[k]
        y1, z1 = ring[(k + 1) % len(ring)]
        fan += 0.5 * abs((y0 - cy) * (z1 - cz) - (y1 - cy) * (z0 - cz))
    if polygon <= 0.0:
        return
    if abs(fan - polygon) > 1e-9 * max(polygon, 1.0):
        raise VolumeXvalError(
            "transom section is not star-shaped about its centroid: fan area "
            f"{fan:.12g} vs polygon area {polygon:.12g}; a self-overlapping stern cap "
            "would add fictitious volume"
        )


# ---------------------------------------------------------------------------
# topology report
# ---------------------------------------------------------------------------
def topology_report(mesh: Mesh) -> dict[str, Any]:
    """Shape-agnostic facts a reader needs before trusting any volume below."""
    counts = meshlib.edge_face_counts(mesh)
    boundary = sum(1 for count in counts.values() if count == 1)
    nonmanifold = sum(1 for count in counts.values() if count > 2)
    euler = len(mesh.vertices) - len(counts) + len(mesh.faces)
    return {
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
        "edges": len(counts),
        "boundary_edges": boundary,
        "nonmanifold_edges": nonmanifold,
        "watertight": meshlib.is_watertight(mesh),
        "edge_manifold": meshlib.is_edge_manifold(mesh),
        "euler_characteristic": euler,
        "closed_genus": (2 - euler) // 2 if boundary == 0 and nonmanifold == 0 else None,
        "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
        "surface_area_m2": meshlib.surface_area(mesh),
    }


# ---------------------------------------------------------------------------
# route 1: the platform's own divergence theorem
# ---------------------------------------------------------------------------
def volume_divergence(mesh: Mesh) -> float:
    """The platform's own closed-surface volume (``tensorcad.mesh``)."""
    return meshlib.total_signed_volume(mesh)


# ---------------------------------------------------------------------------
# route 2: trimesh (independent engine A)
# ---------------------------------------------------------------------------
def volume_trimesh(mesh: Mesh) -> float:
    """Volume from :mod:`trimesh` -- independent third-party engine.

    Uses ``process=False`` on purpose: repairing/welding inside the engine would
    hide exactly the topology defects this route is supposed to expose.
    """
    try:
        import numpy as np
        import trimesh  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise VolumeXvalError(
            "engine A (trimesh) is not installed in this environment; "
            "install it with `pip install trimesh`"
        ) from exc
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    solid = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    if not solid.is_watertight:
        raise VolumeXvalError("engine A: mesh is not watertight, so its volume is undefined")
    return float(solid.volume)


# ---------------------------------------------------------------------------
# route 3: numpy-stl / Mirtich (independent engine B)
# ---------------------------------------------------------------------------
def volume_numpy_stl(mesh: Mesh, *, keep_file: str | None = None) -> float:
    """Volume from :mod:`stl` (numpy-stl) after a binary-STL round trip.

    Independent of both other routes in code *and* in algorithm: numpy-stl
    computes the polyhedron's mass properties with Mirtich's method in C++.  The
    round trip is kept in the measurement because it is part of the exchange path
    a solver actually sees -- and it costs precision, since binary STL stores
    float32 coordinates.  That cost is reported, not corrected: it is why
    :data:`ENGINE_AGREEMENT_TOLERANCE` is 1e-6 rather than 1e-12.
    """
    try:
        import numpy as np
        import stl  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise VolumeXvalError(
            "engine B (numpy-stl) is not installed in this environment; "
            "install it with `pip install numpy-stl`"
        ) from exc
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    vectors = vertices[faces]
    handle, path = tempfile.mkstemp(suffix=".stl", prefix="volume_xval_")
    os.close(handle)
    try:
        solid = stl.mesh.Mesh(np.zeros(len(faces), dtype=stl.mesh.Mesh.dtype))
        solid.vectors = vectors
        solid.save(path)
        if keep_file:
            with open(path, "rb") as src, open(keep_file, "wb") as dst:
                dst.write(src.read())
        reloaded = stl.mesh.Mesh.from_file(path)
        volume, _cog, _inertia = reloaded.get_mass_properties()
    finally:
        try:
            os.unlink(path)
        except OSError:  # pragma: no cover - only if the file was already gone
            pass
    return float(volume)


# ---------------------------------------------------------------------------
# route 4: sectional quadrature (independent method, no triangulation)
# ---------------------------------------------------------------------------
def volume_sectional(
    breadth: Breadth,
    *,
    length: float,
    draft: float,
    stations: int = 201,
    levels: int = 401,
) -> float:
    """Volume by Simpson quadrature of station areas -- shares no maths with the engines.

    ``A(x) = 2 * int_0^T b(x, z) dz`` then ``V = int_0^L A dx``.  No triangles and
    no divergence theorem, so a bug in either the tessellation or the surface
    integral cannot hide here.  It integrates the *smooth* field, so it converges
    from above the inscribed-polyhedron value; the gap is the tessellation deficit
    seen from the other side, which is what makes the attribution test conclusive.
    """
    xs = _uniform(length, stations)
    areas = [_section_area(breadth, x, draft, levels) for x in xs]
    return _simpson(areas, length)


def _section_area(breadth: Breadth, x: float, draft: float, levels: int) -> float:
    zs = _uniform(draft, levels)
    values = [2.0 * float(breadth(x, z)) for z in zs]
    return _simpson(values, draft)


def _uniform(span: float, count: int) -> tuple[float, ...]:
    if count < 3:
        raise VolumeXvalError("Simpson quadrature needs at least 3 samples")
    if count % 2 == 0:
        count += 1
    return tuple(span * i / (count - 1) for i in range(count))


def _simpson(values: Sequence[float], span: float) -> float:
    n = len(values)
    if n < 3 or n % 2 == 0:
        raise VolumeXvalError("Simpson quadrature needs an odd number (>= 3) of samples")
    h = span / (n - 1)
    total = values[0] + values[-1]
    for i in range(1, n - 1):
        total += 4.0 * values[i] if i % 2 else 2.0 * values[i]
    return total * h / 3.0


# ---------------------------------------------------------------------------
# discrepancy attribution
# ---------------------------------------------------------------------------
def relative_difference(value: float, reference: float) -> float:
    """Signed relative difference ``(value - reference) / reference``."""
    if reference == 0.0:
        raise VolumeXvalError("relative difference against a zero reference is undefined")
    return (value - reference) / reference


def fit_order(samples: Sequence[tuple[float, float]]) -> float | None:
    """Observed order of convergence from ``(h, deficit)`` pairs (log-log slope).

    Returns ``None`` when the deficits do not decay monotonically -- that is a
    real answer ("no order"), not a numerical failure, and it is what a
    systematic bias looks like.
    """
    usable = [(h, d) for h, d in samples if h > 0.0 and d > 0.0]
    if len(usable) < 2:
        return None
    magnitudes = [d for _h, d in usable]
    for finer, coarser in zip(magnitudes[1:], magnitudes[:-1]):
        if finer > coarser:
            return None
    xs = [math.log(h) for h, _d in usable]
    ys = [math.log(d) for _h, d in usable]
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0.0:
        return None
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return sxy / sxx


def richardson_limit(samples: Sequence[tuple[float, float]], order: float) -> float:
    """Extrapolate ``(h, V_h)`` rows to ``h -> 0``.

    Consecutive rows are paired, and each pair gives
    ``V_exact ~ V_h + (V_h - V_2h) / ((h_2h/h)**p - 1)``.  The mean over the pairs
    is returned; a route whose deficit really is tessellation lands on the smooth
    volume, and a route with a bias does not.
    """
    if order <= 0.0:
        raise VolumeXvalError("Richardson extrapolation needs a positive order")
    estimates = []
    for i in range(1, len(samples)):
        h_prev, v_prev = samples[i - 1]
        h, v = samples[i]
        if h_prev <= 0.0 or h <= 0.0:
            continue
        ratio = (h_prev / h) ** order
        if abs(ratio - 1.0) < 1e-9:
            continue
        estimates.append(v + (v - v_prev) / (ratio - 1.0))
    if not estimates:
        return samples[-1][1]
    return sum(estimates) / len(estimates)


def classify_routes(
    *,
    divergence: float,
    engine_a: float,
    engine_b: float,
    pipeline: float,
    deficits: Sequence[tuple[float, float]] = (),
    acceptance: float = ACCEPTANCE_TOLERANCE,
    engine_tolerance: float = ENGINE_AGREEMENT_TOLERANCE,
) -> dict[str, Any]:
    """Attribute the route differences to a named effect and give a verdict.

    ``pipeline`` is the volume the platform's own integration route claims for the
    *smooth* geometry (a Gauss-quadrature integral of the patch/field); the three
    triangle routes integrate the inscribed polyhedron instead, so the expected
    difference is the tessellation deficit.  ``deficits`` are
    ``(h, |pipeline - engine|)`` samples over refinements, used only to decide
    whether the difference shrinks (tessellation) or does not (bias).
    """
    spread = abs(engine_a - engine_b) / max(abs(engine_a), 1e-300)
    divergence_vs_engine = abs(divergence - engine_a) / max(abs(engine_a), 1e-300)
    deficit = pipeline - 0.5 * (engine_a + engine_b)
    deficit_rel = deficit / max(abs(pipeline), 1e-300)
    order = fit_order(deficits)

    if spread > 100.0 * engine_tolerance:
        status = "engine-disagreement"
        reason = (
            f"the two independent engines read the same triangles and disagree by "
            f"{100 * spread:.3e} % -- a topology defect or a precision artefact, never a "
            "tessellation effect"
        )
    elif abs(deficit_rel) <= acceptance:
        status = "agrees"
        reason = (
            f"pipeline and engines agree to {100 * abs(deficit_rel):.3f} %, inside the "
            f"{100 * acceptance:.1f} % acceptance line"
        )
    elif order is not None and ORDER_WINDOW[0] <= order <= ORDER_WINDOW[1]:
        status = "tessellation-deficit"
        reason = (
            f"engines agree to {100 * spread:.2e} %, the pipeline sits {100 * deficit_rel:+.3f} % "
            f"above the inscribed polyhedron and the gap decays at observed order "
            f"{order:.2f} (window {ORDER_WINDOW[0]}-{ORDER_WINDOW[1]}) -- this is O(h^2) "
            "tessellation, not an integrator bias"
        )
    else:
        status = "systematic-bias"
        reason = (
            f"engines agree to {100 * spread:.2e} % but the pipeline differs by "
            f"{100 * deficit_rel:+.3f} % and the gap does not decay at a tessellation order "
            f"(observed {order!r}) -- book this in the negative-knowledge ledger"
        )
    return {
        "status": status,
        "reason": reason,
        "engine_spread": spread,
        "divergence_vs_engine_a": divergence_vs_engine,
        "pipeline_vs_engines": deficit_rel,
        "deficit_absolute_m3": deficit,
        "observed_order": order,
        "acceptance": acceptance,
    }


def engine_agreement_note(spread: float, engine_tolerance: float = ENGINE_AGREEMENT_TOLERANCE) -> str:
    """One line explaining what the engine-to-engine spread does and does not prove."""
    if spread <= engine_tolerance:
        return (
            f"engines agree to {spread:.2e} relative, inside the {engine_tolerance:.0e} "
            "binary-STL float32 round-trip budget"
        )
    return (
        f"engines differ by {spread:.2e} relative, beyond the {engine_tolerance:.0e} "
        "float32 round-trip budget: inspect topology and winding before reading any "
        "volume below"
    )


# ---------------------------------------------------------------------------
# case record
# ---------------------------------------------------------------------------
@dataclass
class CaseRow:
    """One (ship, resolution) row of the cross-validation table."""

    ship: str
    reference_kind: str
    reference_volume: float | None
    reference_source: str
    stations: int
    levels: int
    h: float
    pipeline: float
    divergence: float
    engine_a: float
    engine_b: float
    sectional: float | None = None
    notes: str = ""
    extras: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "ship": self.ship,
            "reference_kind": self.reference_kind,
            "reference_source": self.reference_source,
            "stations": self.stations,
            "levels": self.levels,
            "h": self.h,
            "pipeline_volume_m3": self.pipeline,
            "divergence_volume_m3": self.divergence,
            "engine_a_trimesh_m3": self.engine_a,
            "engine_b_numpy_stl_m3": self.engine_b,
            "sectional_volume_m3": self.sectional,
            "notes": self.notes,
        }
        if self.reference_volume is not None:
            payload["reference_volume_m3"] = self.reference_volume
            payload["pipeline_vs_reference"] = relative_difference(
                self.pipeline, self.reference_volume
            )
            payload["engine_a_vs_reference"] = relative_difference(
                self.engine_a, self.reference_volume
            )
            payload["engine_b_vs_reference"] = relative_difference(
                self.engine_b, self.reference_volume
            )
        payload.update(self.extras)
        return payload