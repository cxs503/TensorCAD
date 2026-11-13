"""Kovasznay (1948) steady two-dimensional Navier-Stokes exact solution.

Why this module exists
----------------------
The Kovasznay flow is the canonical *manufactured* solution for verifying
incompressible Navier-Stokes codes: it satisfies the steady momentum and
continuity equations **exactly**, so a solver's discretisation error is known in
closed form and the observed order of accuracy can be measured on a mesh family
(32/64/128 cells).  TensorCAD's role in the ecosystem is geometry and reference
data -- exactly the role it plays for the lid-driven cavity in
:mod:`tensorcad.cavity_mesh`; the flow solver itself lives in TensorFVM /
TensorLBM / TensorDG and consumes the artefacts written here.

Provided here
-------------
1. The exact solution with analytic derivatives (:class:`KovasznayFlow`), so
   residuals can be evaluated without finite-difference stencils.
2. A structured quadrilateral mesh family over the standard computational
   domain ``[-0.5, 1.5] x [-0.5, 1.5]`` with boundary-condition tags.
3. Exact reference fields (cell-centre and nodal values) exported as JSON,
   plus a probe table for spot checks.
4. Consistency (truncation-error) measurements: the second-order central
   difference Navier-Stokes operator applied to the exact solution must have a
   residual that decays as ``h^2`` -- the order the discrete scheme can deliver.
5. A solver-agnostic order-of-accuracy harness
   (:func:`convergence_study`): hand it a callback returning the cell-centre
   velocity a solver reports and it returns the observed log-log order ``p`` in
   ``||u_h - u||_2 = C h^p``.

Exact solution
--------------
With ``Re = U L / nu`` and ``U = L = 1``, the decay rate ``lam`` is the negative
root of the quadratic (the positive root is unbounded downstream)::

    lam^2 - Re*lam - 4*pi^2 = 0   =>   lam = (Re - sqrt(Re^2 + 16*pi^2)) / 2

    u(x, y) = 1 - exp(lam*x) * cos(2*pi*y)
    v(x, y) = lam/(2*pi) * exp(lam*x) * sin(2*pi*y)
    p(x, y) = p0 - exp(2*lam*x) / 2                    (p0 = pressure datum)

Stream function ``psi = y - exp(lam*x) sin(2*pi*y) / (2*pi)`` gives
``u = dpsi/dy`` and ``v = -dpsi/dx``; the vorticity obtained from ``psi`` and
from the velocity gradient are two algebraically independent routes to the same
quantity and are cross-checked in
:meth:`KovasznayFlow.vorticity_mismatch`.

Reference
---------
Kovasznay, L. I. G. (1948). "Laminar flow behind a two-dimensional grid."
*Proceedings of the Cambridge Philosophical Society*, 44(1), 58-62.
doi:10.1017/S0305004100023999
"""
from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "KOVASZNAY_1948",
    "KOVASZNAY_ORDER_ACCEPTANCE",
    "KovasznayFlow",
    "StructuredMesh",
    "case_spec",
    "cell_averages",
    "cell_centre_fields",
    "convergence_study",
    "fit_order",
    "generate_structured_mesh",
    "l2_norm",
    "lam_for_reynolds",
    "mesh_summary",
    "mesh_to_json",
    "midpoint_consistency_study",
    "nodal_fields",
    "probe_table",
    "truncation_error_study",
    "write_case_bundle",
    "write_mesh_json",
]

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
KOVASZNAY_1948: dict[str, Any] = {
    "citation": (
        'Kovasznay, L. I. G. (1948). "Laminar flow behind a two-dimensional '
        'grid." Proceedings of the Cambridge Philosophical Society, 44(1), 58-62.'
    ),
    "doi": "10.1017/S0305004100023999",
    "best_for": (
        "manufactured solution for incompressible Navier-Stokes: exact momentum "
        "and continuity balance, smooth fields, closed-form error, order study"
    ),
    "solution_form": (
        "u = 1 - exp(lam*x)cos(2*pi*y); v = lam/(2*pi)exp(lam*x)sin(2*pi*y); "
        "p = p0 - exp(2*lam*x)/2; lam^2 - Re*lam - 4*pi^2 = 0"
    ),
}

# Benchmark acceptance used by the case-pool entry FVM-KOVASZNAY-ORDER:
# observed L2 order >= 1.8 (nominal second order) from 32/64/128 meshes.
KOVASZNAY_ORDER_ACCEPTANCE: dict[str, Any] = {
    "metric": "log-log slope of velocity L2 error vs h on 32/64/128 meshes",
    "nominal_order": 2.0,
    "min_acceptable_order": 1.8,
    "levels": [32, 64, 128],
}

_TWO_PI = 2.0 * math.pi

# 4-point Gauss-Legendre rule on [-1, 1] (used for cell averages).
_GL4_X = (-0.8611363115940526, -0.3399810435848563,
          0.3399810435848563, 0.8611363115940526)
_GL4_W = (0.3478548451374538, 0.6521451548625461,
          0.6521451548625461, 0.3478548451374538)


def lam_for_reynolds(reynolds_number: float) -> float:
    """Decay rate of the Kovasznay solution (negative root of the quadratic)."""
    if reynolds_number <= 0.0:
        raise ValueError("reynolds_number must be positive")
    disc = reynolds_number * reynolds_number + 16.0 * math.pi * math.pi
    return 0.5 * (reynolds_number - math.sqrt(disc))


# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class KovasznayFlow:
    """Steady Kovasznay solution for a given Reynolds number and domain.

    Defaults follow the standard benchmark set-up: ``Re = 40`` on
    ``[-0.5, 1.5]^2`` with unit reference velocity and length.
    """

    reynolds_number: float = 40.0
    x_min: float = -0.5
    x_max: float = 1.5
    y_min: float = -0.5
    y_max: float = 1.5
    pressure_offset: float = 0.0

    def __post_init__(self) -> None:
        if self.reynolds_number <= 0.0:
            raise ValueError("reynolds_number must be positive")
        if not self.x_max > self.x_min:
            raise ValueError("x_max must exceed x_min")
        if not self.y_max > self.y_min:
            raise ValueError("y_max must exceed y_min")

    # -- parameters ---------------------------------------------------------
    @property
    def lam(self) -> float:
        """Decay rate ``lam`` (negative)."""
        return lam_for_reynolds(self.reynolds_number)

    @property
    def nu(self) -> float:
        """Kinematic viscosity ``U L / Re`` with ``U = L = 1``."""
        return 1.0 / self.reynolds_number

    @property
    def width(self) -> float:
        return self.x_max - self.x_min

    @property
    def height(self) -> float:
        return self.y_max - self.y_min

    # -- primary fields -----------------------------------------------------
    def velocity(self, x: float, y: float) -> tuple[float, float]:
        """Exact velocity ``(u, v)``."""
        e = math.exp(self.lam * x)
        return (1.0 - e * math.cos(_TWO_PI * y),
                self.lam / _TWO_PI * e * math.sin(_TWO_PI * y))

    def pressure(self, x: float, y: float) -> float:
        """Exact pressure (defined up to :attr:`pressure_offset`)."""
        return self.pressure_offset - 0.5 * math.exp(2.0 * self.lam * x)

    def streamfunction(self, x: float, y: float) -> float:
        """Stream function ``psi`` with ``u = dpsi/dy``, ``v = -dpsi/dx``."""
        e = math.exp(self.lam * x)
        return y - e * math.sin(_TWO_PI * y) / _TWO_PI

    def speed(self, x: float, y: float) -> float:
        u, v = self.velocity(x, y)
        return math.hypot(u, v)

    # -- analytic derivatives ----------------------------------------------
    def velocity_gradient(self, x: float, y: float) -> tuple[float, float, float, float]:
        """``(u_x, u_y, v_x, v_y)`` analytically."""
        e = math.exp(self.lam * x)
        c = math.cos(_TWO_PI * y)
        s = math.sin(_TWO_PI * y)
        lam = self.lam
        return (-lam * e * c,                    # u_x
                _TWO_PI * e * s,                 # u_y
                lam * lam / _TWO_PI * e * s,     # v_x
                lam * e * c)                     # v_y

    def vorticity_from_velocity(self, x: float, y: float) -> float:
        """``omega = v_x - u_y`` from the velocity gradient."""
        _, uy, vx, _ = self.velocity_gradient(x, y)
        return vx - uy

    def vorticity_from_streamfunction(self, x: float, y: float) -> float:
        """``omega = -(psi_xx + psi_yy)`` from the stream function."""
        e = math.exp(self.lam * x)
        s = math.sin(_TWO_PI * y)
        lam = self.lam
        psi_xx = -lam * lam / _TWO_PI * e * s
        psi_yy = _TWO_PI * e * s
        return -(psi_xx + psi_yy)

    def vorticity_mismatch(self, x: float, y: float) -> float:
        """Absolute difference between the two independent vorticity routes."""
        return abs(self.vorticity_from_velocity(x, y)
                   - self.vorticity_from_streamfunction(x, y))

    # -- verification residuals --------------------------------------------
    def ns_residual(self, x: float, y: float) -> tuple[float, float, float]:
        """Steady incompressible NS residual using analytic derivatives.

        Returns ``(R_u, R_v, div)`` where::

            R_u = u u_x + v u_y + p_x - (1/Re)(u_xx + u_yy)
            R_v = u v_x + v v_y + p_y - (1/Re)(v_xx + v_yy)
            div = u_x + v_y

        All three vanish for an exact solution; because the derivatives are
        analytic, only round-off remains.  That makes this a check of the
        implementation rather than of a finite-difference approximation.
        """
        lam = self.lam
        e = math.exp(lam * x)
        c = math.cos(_TWO_PI * y)
        s = math.sin(_TWO_PI * y)
        inv_re = 1.0 / self.reynolds_number
        u = 1.0 - e * c
        v = lam / _TWO_PI * e * s
        u_x = -lam * e * c
        u_y = _TWO_PI * e * s
        v_x = lam * lam / _TWO_PI * e * s
        v_y = lam * e * c
        u_xx = -lam * lam * e * c
        u_yy = 4.0 * math.pi * math.pi * e * c
        v_xx = lam ** 3 / _TWO_PI * e * s
        v_yy = -_TWO_PI * lam * e * s
        p_x = -lam * math.exp(2.0 * lam * x)   # p = p0 - exp(2 lam x)/2
        p_y = 0.0
        r_u = u * u_x + v * u_y + p_x - inv_re * (u_xx + u_yy)
        r_v = u * v_x + v * v_y + p_y - inv_re * (v_xx + v_yy)
        div = u_x + v_y
        return r_u, r_v, div

    # -- serialisation ------------------------------------------------------
    def as_dict(self) -> dict[str, Any]:
        return {
            "model": "kovasznay_1948",
            "reynolds_number": self.reynolds_number,
            "nu": self.nu,
            "lam": self.lam,
            "domain": {"x": [self.x_min, self.x_max],
                       "y": [self.y_min, self.y_max]},
            "pressure_offset": self.pressure_offset,
            "reference": KOVASZNAY_1948,
        }


# ---------------------------------------------------------------------------
# Structured mesh family
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class StructuredMesh:
    """Uniform structured quadrilateral mesh on the flow domain.

    Node ``(i, j)`` has index ``j * (nx + 1) + i``; quad ``(ci, cj)`` is
    ``(p(i,j), p(i+1,j), p(i+1,j+1), p(i,j+1))`` -- counter-clockwise, positive
    Jacobian.
    """

    flow: KovasznayFlow
    nx: int
    ny: int
    nodes: list[tuple[float, float]] = field(default_factory=list)
    quads: list[tuple[int, int, int, int]] = field(default_factory=list)

    @property
    def hx(self) -> float:
        return self.flow.width / self.nx

    @property
    def hy(self) -> float:
        return self.flow.height / self.ny

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def node_index(self, i: int, j: int) -> int:
        return j * (self.nx + 1) + i

    def node(self, i: int, j: int) -> tuple[float, float]:
        return self.nodes[self.node_index(i, j)]

    def centroids(self) -> list[tuple[float, float]]:
        """Cell-centre coordinates in quad order."""
        hx, hy = self.hx, self.hy
        x0, y0 = self.flow.x_min, self.flow.y_min
        return [(x0 + (ci + 0.5) * hx, y0 + (cj + 0.5) * hy)
                for cj in range(self.ny) for ci in range(self.nx)]

    def areas(self) -> list[float]:
        """Cell areas via the shoelace formula (sign exposes bad connectivity)."""
        out = []
        for quad in self.quads:
            s = 0.0
            for k in range(4):
                xa, ya = self.nodes[quad[k]]
                xb, yb = self.nodes[quad[(k + 1) % 4]]
                s += xa * yb - xb * ya
            out.append(0.5 * s)
        return out

    def boundary_sets(self) -> dict[str, list[int]]:
        """Node indices per side: ``left``/``right``/``bottom``/``top``."""
        nx, ny = self.nx, self.ny
        return {
            "left": [self.node_index(0, j) for j in range(ny + 1)],
            "right": [self.node_index(nx, j) for j in range(ny + 1)],
            "bottom": [self.node_index(i, 0) for i in range(nx + 1)],
            "top": [self.node_index(i, ny) for i in range(nx + 1)],
        }

    def summary(self) -> dict[str, Any]:
        areas = self.areas()
        hx, hy = self.hx, self.hy
        domain_area = self.flow.width * self.flow.height
        return {
            "nx": self.nx,
            "ny": self.ny,
            "n_nodes": self.n_nodes,
            "n_cells": self.n_cells,
            "hx": hx,
            "hy": hy,
            "min_cell_area": min(areas),
            "max_cell_area": max(areas),
            "total_area": sum(areas),
            "domain_area": domain_area,
            "area_error": sum(areas) - domain_area,
            "aspect_ratio": max(hx, hy) / min(hx, hy),
            "negative_jacobians": sum(1 for a in areas if a <= 0.0),
        }


def generate_structured_mesh(flow: KovasznayFlow, nx: int,
                             ny: int | None = None) -> StructuredMesh:
    """Build the uniform quad mesh with ``nx`` x ``ny`` cells (square default)."""
    if ny is None:
        ny = nx
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be >= 1")
    x0, y0 = flow.x_min, flow.y_min
    hx = flow.width / nx
    hy = flow.height / ny
    nodes = [(x0 + i * hx, y0 + j * hy)
             for j in range(ny + 1) for i in range(nx + 1)]
    quads = []
    for cj in range(ny):
        for ci in range(nx):
            p00 = cj * (nx + 1) + ci
            quads.append((p00, p00 + 1, p00 + nx + 2, p00 + nx + 1))
    return StructuredMesh(flow=flow, nx=nx, ny=ny, nodes=nodes, quads=quads)


def mesh_summary(mesh: StructuredMesh) -> dict[str, Any]:
    """Convenience wrapper around :meth:`StructuredMesh.summary`."""
    return mesh.summary()


# ---------------------------------------------------------------------------
# Exact reference fields
# ---------------------------------------------------------------------------
def cell_centre_fields(mesh: StructuredMesh) -> dict[str, list[float]]:
    """Exact ``u``, ``v``, ``p`` at every cell centroid."""
    flow = mesh.flow
    us, vs, ps = [], [], []
    for x, y in mesh.centroids():
        u, v = flow.velocity(x, y)
        us.append(u)
        vs.append(v)
        ps.append(flow.pressure(x, y))
    return {"u": us, "v": vs, "p": ps}


def nodal_fields(mesh: StructuredMesh) -> dict[str, list[float]]:
    """Exact ``u``, ``v``, ``p`` at every mesh node."""
    flow = mesh.flow
    us, vs, ps = [], [], []
    for x, y in mesh.nodes:
        u, v = flow.velocity(x, y)
        us.append(u)
        vs.append(v)
        ps.append(flow.pressure(x, y))
    return {"u": us, "v": vs, "p": ps}


def cell_averages(mesh: StructuredMesh) -> dict[str, list[float]]:
    """Cell averages of the exact solution via a 4x4 Gauss-Legendre rule.

    A finite-volume solver carries *cell averages*; comparing those against the
    point value at the centroid isolates the ``O(h^2)`` midpoint-rule error and
    is what :func:`midpoint_consistency_study` measures.
    """
    flow = mesh.flow
    hx, hy = mesh.hx, mesh.hy
    us, vs, ps = [], [], []
    for xc, yc in mesh.centroids():
        su = sv = sp = 0.0
        for xi, wi in zip(_GL4_X, _GL4_W, strict=False):
            x = xc + 0.5 * hx * xi
            for eta, wj in zip(_GL4_X, _GL4_W, strict=False):
                y = yc + 0.5 * hy * eta
                u, v = flow.velocity(x, y)
                su += wi * wj * u
                sv += wi * wj * v
                sp += wi * wj * flow.pressure(x, y)
        us.append(su / 4.0)
        vs.append(sv / 4.0)
        ps.append(sp / 4.0)
    return {"u": us, "v": vs, "p": ps}


# ---------------------------------------------------------------------------
# Norms and order fitting
# ---------------------------------------------------------------------------
def l2_norm(values: Sequence[float],
            weights: Sequence[float] | None = None) -> float:
    """Weighted discrete L2 norm ``sqrt(sum w v^2 / sum w)``."""
    if not values:
        raise ValueError("values must be non-empty")
    if weights is None:
        return math.sqrt(sum(v * v for v in values) / len(values))
    if len(weights) != len(values):
        raise ValueError("weights and values must have the same length")
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("weights must sum to a positive value")
    return math.sqrt(sum(w * v * v for w, v in zip(weights, values, strict=False)) / total)


def fit_order(hs: Sequence[float], errors: Sequence[float]) -> dict[str, Any]:
    """Least-squares fit of ``log(err) = p log(h) + c``; returns ``p`` and ``R^2``.

    ``err ~ h**p`` means the fitted slope *is* the order of accuracy.  Fewer than
    two samples, or non-positive values, raise ``ValueError`` so a caller cannot
    silently accept a meaningless order.
    """
    if len(hs) != len(errors):
        raise ValueError("hs and errors must have the same length")
    if len(hs) < 2:
        raise ValueError("at least two grid levels are required")
    if any(h <= 0.0 for h in hs):
        raise ValueError("grid spacings must be positive")
    if any(e <= 0.0 for e in errors):
        raise ValueError("errors must be positive to fit a log-log slope")
    xs = [math.log(h) for h in hs]
    ys = [math.log(e) for e in errors]
    n = float(len(xs))
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys, strict=False))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        raise ValueError("grid spacings are degenerate; cannot fit a slope")
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    mean_y = sy / n
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys, strict=False))
    r2 = 1.0 if ss_tot < 1e-30 else 1.0 - ss_res / ss_tot
    # err ~ h**p  =>  log(err) = p*log(h) + c, so the fitted slope *is* the order.
    return {"order": slope, "intercept": intercept, "r2": r2}


_ROUNDOFF_FLOOR = 1e-13


# ---------------------------------------------------------------------------
# Solver-agnostic order-of-accuracy harness
# ---------------------------------------------------------------------------
def convergence_study(
    flow: KovasznayFlow,
    levels: Sequence[int] = (32, 64, 128),
    discrete_velocity: Callable[[StructuredMesh], Sequence[tuple[float, float]]] | None = None,
) -> dict[str, Any]:
    """Measure the observed order of a solver on the 32/64/128 mesh family.

    ``discrete_velocity(mesh)`` must return the solver's velocity at every cell
    centroid, in :meth:`StructuredMesh.centroids` order, as ``(u, v)`` pairs.
    The returned dict holds the per-level L2 error and the fitted log-log order
    -- i.e. exactly the measurement the FVM case asks for.

    Round-off is flagged, not fitted: if the finest-level error drops below
    ``1e-13`` the order is reported as ``None`` with a ``note``, because a slope
    fitted through round-off carries no information.
    """
    if discrete_velocity is None:
        raise ValueError("discrete_velocity callback is required")
    rows: list[dict[str, Any]] = []
    hs: list[float] = []
    errors: list[float] = []
    for n in levels:
        mesh = generate_structured_mesh(flow, n)
        centroids = mesh.centroids()
        reported = list(discrete_velocity(mesh))
        if len(reported) != len(centroids):
            raise ValueError(
                f"discrete_velocity returned {len(reported)} values for "
                f"{len(centroids)} cells at level {n}")
        defects = []
        for (x, y), (u, v) in zip(centroids, reported, strict=False):
            ue, ve = flow.velocity(x, y)
            # l2_norm squares its input, so pass the defect magnitude itself.
            defects.append(math.hypot(u - ue, v - ve))
        err = l2_norm(defects, mesh.areas())
        h = math.sqrt(flow.width * flow.height / mesh.n_cells)
        rows.append({"n": n, "h": h, "l2_error": err,
                     "n_cells": mesh.n_cells})
        hs.append(h)
        errors.append(err)
    result: dict[str, Any] = {
        "levels": rows,
        "h": hs,
        "l2_error": errors,
        "acceptance": KOVASZNAY_ORDER_ACCEPTANCE,
    }
    if min(errors) <= _ROUNDOFF_FLOOR:
        result["order"] = None
        result["note"] = (
            "finest-level error at round-off "
            f"({min(errors):.3e}); order not meaningful")
    else:
        fit = fit_order(hs, errors)
        result.update(fit)
        result["passes_acceptance"] = fit["order"] >= KOVASZNAY_ORDER_ACCEPTANCE["min_acceptable_order"]
    return result


# ---------------------------------------------------------------------------
# Consistency (truncation error) measurements
# ---------------------------------------------------------------------------
def midpoint_consistency_study(
    flow: KovasznayFlow,
    levels: Sequence[int] = (32, 64, 128),
) -> dict[str, Any]:
    """Order of the cell-average <-> centroid-value (midpoint rule) difference.

    A cell-centred FVM stores cell averages; interpreting them as point values at
    the centroid costs ``h^2 (u_xx + u_yy)/24``.  Measuring that difference on
    the mesh family proves the ``O(h^2)`` representation error independently of
    any solver, and pins the asymptotic regime in which a second-order method can
    actually be observed.
    """
    rows: list[dict[str, Any]] = []
    hs: list[float] = []
    err_u: list[float] = []
    err_v: list[float] = []
    err_uv: list[float] = []
    for n in levels:
        mesh = generate_structured_mesh(flow, n)
        point = cell_centre_fields(mesh)
        avg = cell_averages(mesh)
        du = [a - b for a, b in zip(avg["u"], point["u"], strict=False)]
        dv = [a - b for a, b in zip(avg["v"], point["v"], strict=False)]
        areas = mesh.areas()
        eu = l2_norm(du, areas)
        ev = l2_norm(dv, areas)
        euv = math.sqrt(sum(w * (a * a + b * b)
                            for w, a, b in zip(areas, du, dv, strict=False)) / sum(areas))
        h = math.sqrt(flow.width * flow.height / mesh.n_cells)
        rows.append({"n": n, "h": h, "l2_u": eu, "l2_v": ev,
                     "l2_velocity": euv})
        hs.append(h)
        err_u.append(eu)
        err_v.append(ev)
        err_uv.append(euv)
    result: dict[str, Any] = {"levels": rows, "h": hs}
    for name, errs in (("u", err_u), ("v", err_v), ("velocity", err_uv)):
        if min(errs) > _ROUNDOFF_FLOOR:
            fit = fit_order(hs, errs)
            result[f"order_{name}"] = fit["order"]
            result[f"r2_{name}"] = fit["r2"]
        else:
            result[f"order_{name}"] = None
    return result


def _nodal(mesh: StructuredMesh, field: Sequence[float],
           i: int, j: int) -> float:
    """Nodal value of ``field`` at grid index ``(i, j)``."""
    return field[mesh.node_index(i, j)]


def truncation_error_study(
    flow: KovasznayFlow,
    levels: Sequence[int] = (32, 64, 128),
) -> dict[str, Any]:
    """Second-order central-difference NS operator applied to the exact solution.

    The nodal exact values are frozen and the discrete steady momentum and
    continuity operators are evaluated with standard 3-point central
    differences.  If the implementation of the operators is consistent, the
    residual decays as ``h^2``; the fitted log-log slopes are the consistency
    order of the discretisation (not the error order of a converged solution).
    """
    inv_re = 1.0 / flow.reynolds_number
    rows: list[dict[str, Any]] = []
    hs: list[float] = []
    res_x: list[float] = []
    res_y: list[float] = []
    res_c: list[float] = []
    for n in levels:
        mesh = generate_structured_mesh(flow, n)
        f = nodal_fields(mesh)
        h = mesh.hx
        nx, ny = mesh.nx, mesh.ny
        rx: list[float] = []
        ry: list[float] = []
        rc: list[float] = []
        for j in range(1, ny):
            for i in range(1, nx):
                u = _nodal(mesh, f["u"], i, j)
                v = _nodal(mesh, f["v"], i, j)
                _p = _nodal(mesh, f["p"], i, j)
                u_x = (_nodal(mesh, f["u"], i + 1, j) - _nodal(mesh, f["u"], i - 1, j)) / (2 * h)
                u_y = (_nodal(mesh, f["u"], i, j + 1) - _nodal(mesh, f["u"], i, j - 1)) / (2 * h)
                v_x = (_nodal(mesh, f["v"], i + 1, j) - _nodal(mesh, f["v"], i - 1, j)) / (2 * h)
                v_y = (_nodal(mesh, f["v"], i, j + 1) - _nodal(mesh, f["v"], i, j - 1)) / (2 * h)
                p_x = (_nodal(mesh, f["p"], i + 1, j) - _nodal(mesh, f["p"], i - 1, j)) / (2 * h)
                p_y = (_nodal(mesh, f["p"], i, j + 1) - _nodal(mesh, f["p"], i, j - 1)) / (2 * h)
                u_xx = (_nodal(mesh, f["u"], i + 1, j) - 2 * u + _nodal(mesh, f["u"], i - 1, j)) / (h * h)
                u_yy = (_nodal(mesh, f["u"], i, j + 1) - 2 * u + _nodal(mesh, f["u"], i, j - 1)) / (h * h)
                v_xx = (_nodal(mesh, f["v"], i + 1, j) - 2 * v + _nodal(mesh, f["v"], i - 1, j)) / (h * h)
                v_yy = (_nodal(mesh, f["v"], i, j + 1) - 2 * v + _nodal(mesh, f["v"], i, j - 1)) / (h * h)
                rx.append(u * u_x + v * u_y + p_x - inv_re * (u_xx + u_yy))
                ry.append(u * v_x + v * v_y + p_y - inv_re * (v_xx + v_yy))
                rc.append(u_x + v_y)
        ex, ey, ec = l2_norm(rx), l2_norm(ry), l2_norm(rc)
        rows.append({"n": n, "h": h, "l2_momentum_x": ex,
                     "l2_momentum_y": ey, "l2_continuity": ec})
        hs.append(h)
        res_x.append(ex)
        res_y.append(ey)
        res_c.append(ec)
    result: dict[str, Any] = {"levels": rows, "h": hs}
    for name, errs in (("momentum_x", res_x), ("momentum_y", res_y),
                       ("continuity", res_c)):
        if min(errs) > _ROUNDOFF_FLOOR:
            fit = fit_order(hs, errs)
            result[f"order_{name}"] = fit["order"]
            result[f"r2_{name}"] = fit["r2"]
        else:
            result[f"order_{name}"] = None
    return result


# ---------------------------------------------------------------------------
# Reference probe table and case specification
# ---------------------------------------------------------------------------
DEFAULT_PROBE_X = (-0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5)
DEFAULT_PROBE_Y = (-0.5, -0.25, 0.0, 0.25, 0.5)


def probe_table(flow: KovasznayFlow,
                xs: Sequence[float] = DEFAULT_PROBE_X,
                ys: Sequence[float] = DEFAULT_PROBE_Y) -> list[dict[str, float]]:
    """Exact solution on a probe lattice (numbers generated, never transcribed)."""
    rows: list[dict[str, float]] = []
    for y in ys:
        for x in xs:
            u, v = flow.velocity(x, y)
            ru, rv, div = flow.ns_residual(x, y)
            rows.append({
                "x": x, "y": y, "u": u, "v": v,
                "p": flow.pressure(x, y),
                "vorticity": flow.vorticity_from_velocity(x, y),
                "ns_residual_u": ru, "ns_residual_v": rv, "continuity": div,
                "vorticity_mismatch": flow.vorticity_mismatch(x, y),
            })
    return rows


def case_spec(flow: KovasznayFlow,
              levels: Sequence[int] = (32, 64, 128)) -> dict[str, Any]:
    """Machine-readable specification a solver repo can consume."""
    spec = flow.as_dict()
    spec.update({
        "case_code": "FVM-KOVASZNAY-ORDER",
        "acceptance": KOVASZNAY_ORDER_ACCEPTANCE,
        "levels": list(levels),
        "meshes": [generate_structured_mesh(flow, n).summary() for n in levels],
        "boundary_conditions": {
            "left": "Dirichlet inflow: u, v from the exact solution",
            "right": "Dirichlet outflow: u, v from the exact solution",
            "bottom": "Dirichlet wall: u, v from the exact solution",
            "top": "Dirichlet wall: u, v from the exact solution",
        },
        "exact_solution": KOVASZNAY_1948["solution_form"],
        "error_measure": "||u_h - u_exact||_2 / ||u_exact||_2 over cell centroids",
    })
    return spec


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------
def _sig(value: float, digits: int = 12) -> float:
    """Round to ``digits`` significant digits (keeps bundles small and stable).

    12 significant digits is ~1e-12 relative, eight orders of magnitude below
    any solver error a mesh family at these resolutions can produce, so the
    reference fields lose nothing measurable.
    """
    return float(f"{value:.{digits}g}")


def mesh_to_json(mesh: StructuredMesh, include_fields: bool = True) -> dict[str, Any]:
    """Serialise the mesh (nodes, quads, boundary tags, optional exact fields)."""
    payload: dict[str, Any] = {
        "summary": mesh.summary(),
        "nodes": [[_sig(x), _sig(y)] for x, y in mesh.nodes],
        "quads": [list(q) for q in mesh.quads],
        "boundary_sets": mesh.boundary_sets(),
        "cell_centroids": [[_sig(x), _sig(y)] for x, y in mesh.centroids()],
    }
    if include_fields:
        centre = cell_centre_fields(mesh)
        nodal = nodal_fields(mesh)
        payload["cell_centre_exact"] = {k: [_sig(v) for v in vals]
                                        for k, vals in centre.items()}
        payload["nodal_exact"] = {k: [_sig(v) for v in vals]
                                  for k, vals in nodal.items()}
    return payload


def _write_json(path: Path, payload: Any, compact: bool = False) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compact:
        text = json.dumps(payload, separators=(",", ":"), sort_keys=True)
    else:
        text = json.dumps(payload, indent=2, sort_keys=True)
    path.write_text(text + "\n", encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_mesh_json(mesh: StructuredMesh, path: str | Path,
                    include_fields: bool = True) -> Path:
    """Write one mesh as JSON and return the path."""
    out = Path(path)
    _write_json(out, mesh_to_json(mesh, include_fields=include_fields))
    return out


def write_case_bundle(flow: KovasznayFlow,
                      outdir: str | Path,
                      levels: Sequence[int] = (32, 64, 128)) -> dict[str, Any]:
    """Write the full solver-consumable bundle; returns a manifest with hashes."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {"outdir": str(out), "files": {}, "case": "FVM-KOVASZNAY-ORDER"}
    manifest["files"]["case_spec.json"] = _write_json(
        out / "case_spec.json", case_spec(flow, levels))
    for n in levels:
        mesh = generate_structured_mesh(flow, n)
        manifest["files"][f"mesh_{n}.json"] = _write_json(
            out / f"mesh_{n}.json", mesh_to_json(mesh), compact=True)
    manifest["files"]["probe_table.json"] = _write_json(
        out / "probe_table.json", probe_table(flow))
    return manifest
