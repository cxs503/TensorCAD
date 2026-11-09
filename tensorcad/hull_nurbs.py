"""NURBS hull surface -> structured body mesh -> geometry consistency loop.

This module is the geometry *source* for the ship-performance workflow: it turns a
handful of principal dimensions (``L``, ``B``, ``T``, target block coefficient ``Cb``)
into a tensor-product NURBS hull patch, integrates the hydrostatic quantities of that
patch with high-order Gauss quadrature, generates a body-fitted hexahedral volume mesh
from the same patch, and finally re-integrates the quantities **from the mesh** so the
two independent routes can be compared.

The point of the module is not the hull form itself -- it is the loop. A geometry source
is only usable by a solver if the numbers a solver measures on the mesh agree with the
numbers the design claims on the surface. Everything here is therefore stated twice:
once on the NURBS patch and once on the mesh.

Design of the surface
---------------------
The patch is a graph over ``(x, z)``::

    S(u, v) = ( L u,  y(u, v),  D v )

``x`` and ``z`` are carried by the control net at the Greville abscissae, so B-spline
linear precision makes ``x(u) = L u`` and ``z(v) = D v`` **exactly**, while ``y`` is a
genuine free-form B-spline coefficient field. The geometry therefore stays a real NURBS
patch (:class:`tensorcad.nurbs.NURBSSurface`, rational weights supported) while the
hydrostatics become an ordinary 2-D integral::

    volume = 2 * L * D * int_0^1 du int_0^{v_wl} dv  y(u, v)

The factor 2 is the mirror plane: the hull is described on the starboard side only.

Hull form
---------
The half-breadth factorises into a plan-form ratio ``f(xi)`` (waterline breadth over
``B/2``) and a section shape ``psi(zeta)`` (local half-breadth over the local waterline
half-breadth)::

    y(x, z) = (B/2) * f(xi) * psi(zeta; xi),    xi = x/L,  zeta = z/T

* ``f`` is 1 over a parallel middle body, a smoothstep run aft down to the transom
  ratio, and a superelliptical entrance forward (``(1 - s^a)^(1/b)``), which gives the
  stem a vertical tangent as a real bow has.
* ``psi = r_b + (1 - r_b) (1 - (1-zeta)^p)^(1/p)`` is a superelliptical section family:
  ``p -> 1`` is fine and V-shaped, large ``p`` is wall-sided with a wide flat of bottom.
  The local fullness ``p(xi)`` relaxes from the midship value to ``p_ends`` at the ends,
  so the entrance sections are finer than the midship one -- which is why the waterplane
  coefficient comes out above the prismatic one, as on real hulls.

Only the parallel-body length ``parallel_body_ratio`` is *not* a design input: it is
calibrated by bisection so that the integrated block coefficient of the NURBS patch
equals the requested ``Cb``. Nothing else is tuned to an acceptance number.

Deliberate limitations (stated, not hidden)
-------------------------------------------
* Single tensor-product patch, no trimming: no bulbous bow, no sonar dome, no
  appendages, no flare above the design waterline (``psi = 1`` for ``zeta > 1``).
* No bilge-radius *knuckle* and no separate flat-of-bottom plane: the flat of bottom is
  the ``psi(0) = r_b`` breadth carried by the side patch, and the bottom is the plane
  ``z = 0``.
* The body mesh is a body-fitted curvilinear grid, not an unstructured tet mesh. The bow
  column and the shell edges are reaches where a structured grid must collapse; those
  cells are counted as ``collapsed`` and reported, not smoothed away. A six-tet
  subdivision of every hex is provided as well.
* Hydrostatics are for the bare hull at even keel, in infinite-depth water, with no
  trim: no free-surface correction, no sinkage, no heel, no appendages.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Any

from .nurbs import NURBSError, NURBSSurface

__all__ = [
    "ASPECT_RATIO_MAX",
    "HULL_NURBS_SCHEMA",
    "SCALED_JACOBIAN_GOOD",
    "SCALED_JACOBIAN_MIN",
    "HexQuality",
    "HullBodyMesh",
    "HullForm",
    "HullHydrostatics",
    "HullNurbsError",
    "HullNurbsParameters",
    "HullNurbsSurface",
    "MeshIntegrals",
    "NURBSError",
    "acceptance_report",
    "build_body_mesh",
    "build_nurbs_hull",
    "calibrate_form",
    "clamped_knots",
    "gauss_legendre",
    "greville_abscissae",
    "hex_scaled_jacobian",
    "hex_to_tets",
    "hex_volume",
    "holtrop_wetted_area",
    "kcs_lightweight_parameters",
    "mesh_integrals",
    "nurbs_hydrostatics",
    "realised_dimensions",
    "schneekluth_lcb_ratio",
]

#: Schema tag of the acceptance artifact produced by :func:`acceptance_report`.
HULL_NURBS_SCHEMA = "tensorcad.hull_nurbs/1.0"

#: Mesh-quality admission thresholds. The scaled-Jacobian threshold follows common CFD
#: practice (ANSYS Meshing: below 0.2 is a "bad" element, above 0.5 a good one); the
#: aspect-ratio cap is the bulk (boundary-layer-free) limit used in the same tradition.
#: Both are declared here, before any measurement, so the pass fraction is not a number
#: tuned to the answer.
SCALED_JACOBIAN_MIN = 0.2
SCALED_JACOBIAN_GOOD = 0.5
ASPECT_RATIO_MAX = 100.0


class HullNurbsError(ValueError):
    """Raised when a hull form or a mesh request is geometrically inadmissible."""


# ---------------------------------------------------------------------------
# quadrature (stdlib only: the core package must install with no dependencies)
# ---------------------------------------------------------------------------
def gauss_legendre(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Nodes and weights of the ``order``-point Gauss-Legendre rule on ``[-1, 1]``.

    Newton iteration on the Legendre polynomial from the ``cos`` initial guess, which
    converges to machine precision in a handful of steps at these orders.
    """
    if order < 1:
        raise HullNurbsError(f"quadrature order must be >= 1, got {order}")
    nodes: list[float] = []
    weights: list[float] = []
    for i in range(1, order + 1):
        x = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        for _ in range(64):
            p0, p1 = 1.0, x
            for k in range(2, order + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = order * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / dp
            x -= step
            if abs(step) < 1e-16:
                break
        p0, p1 = 1.0, x
        for k in range(2, order + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = order * (x * p1 - p0) / (x * x - 1.0)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return tuple(nodes), tuple(weights)


def _mapped_rule(lo: float, hi: float, order: int) -> tuple[list[float], list[float]]:
    nodes, weights = gauss_legendre(order)
    half = 0.5 * (hi - lo)
    mid = 0.5 * (hi + lo)
    return [mid + half * n for n in nodes], [half * w for w in weights]


# ---------------------------------------------------------------------------
# B-spline scaffolding
# ---------------------------------------------------------------------------
def clamped_knots(degree: int, count: int) -> tuple[float, ...]:
    """Uniform clamped knot vector on ``[0, 1]`` for ``count`` control points."""
    if count < degree + 1:
        raise HullNurbsError(f"{count} control points cannot carry degree {degree}")
    inner = count - degree - 1
    knots = [0.0] * (degree + 1)
    knots.extend(k / (inner + 1) for k in range(1, inner + 1))
    knots.extend([1.0] * (degree + 1))
    return tuple(knots)


def greville_abscissae(knots: Sequence[float], degree: int) -> tuple[float, ...]:
    """Greville abscissae: the parameters whose B-spline reproduces linear functions."""
    count = len(knots) - degree - 1
    return tuple(sum(knots[i + 1 : i + degree + 1]) / degree for i in range(count))


# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HullNurbsParameters:
    """Principal dimensions plus the shape knobs of the parametric hull.

    The defaults are the reference case of this verification surface: a KCS-derived
    main-dimension triple at model scale -- ``L = 3.0 m``, ``B = 0.4 m``, ``T = 0.15 m``
    -- with the (lighter) block coefficient ``Cb = 0.555``. The full-scale KCS has
    ``Cb = 0.651``; this form is calibrated to the requested ``Cb`` and is therefore a
    different hull, which the report says outright instead of implying KCS identity.
    """

    length: float = 3.0
    beam: float = 0.4
    draft: float = 0.15
    block_coefficient: float = 0.555

    freeboard_ratio: float = 0.25
    transom_breadth_ratio: float = 0.15
    aft_body_share: float = 0.47
    bottom_breadth_ratio: float = 0.22
    section_exponent_midship: float = 3.6
    section_exponent_ends: float = 0.8
    entrance_exponent: float = 1.4
    entrance_root: float = 1.4
    bottom_taper_power: float = 0.5
    section_relaxation_power: float = 1.0

    parallel_body_ratio: float | None = None
    degree: tuple[int, int] = (3, 3)
    control_points: tuple[int, int] = (21, 13)
    quadrature_order: tuple[int, int] = (24, 16)
    design_froude_number: float = 0.20

    # -- derived ------------------------------------------------------------
    @property
    def deck_height(self) -> float:
        """Depth to the deck: the draft plus the freeboard."""
        return self.draft * (1.0 + self.freeboard_ratio)

    @property
    def design_displacement(self) -> float:
        """The displacement the block coefficient asks for, ``Cb L B T`` (m^3)."""
        return self.block_coefficient * self.length * self.beam * self.draft

    @property
    def length_to_beam(self) -> float:
        return self.length / self.beam if self.beam else float("inf")

    @property
    def beam_to_draft(self) -> float:
        return self.beam / self.draft if self.draft else float("inf")

    def as_dict(self) -> dict[str, Any]:
        return {
            "length_m": self.length,
            "beam_m": self.beam,
            "draft_m": self.draft,
            "block_coefficient": self.block_coefficient,
            "deck_height_m": self.deck_height,
            "design_displacement_m3": self.design_displacement,
            "length_to_beam": self.length_to_beam,
            "beam_to_draft": self.beam_to_draft,
            "freeboard_ratio": self.freeboard_ratio,
            "transom_breadth_ratio": self.transom_breadth_ratio,
            "aft_body_share": self.aft_body_share,
            "bottom_breadth_ratio": self.bottom_breadth_ratio,
            "section_exponent_midship": self.section_exponent_midship,
            "section_exponent_ends": self.section_exponent_ends,
            "entrance_exponent": self.entrance_exponent,
            "entrance_root": self.entrance_root,
            "bottom_taper_power": self.bottom_taper_power,
            "section_relaxation_power": self.section_relaxation_power,
            "parallel_body_ratio": self.parallel_body_ratio,
            "degree": list(self.degree),
            "control_points": list(self.control_points),
            "quadrature_order": list(self.quadrature_order),
            "design_froude_number": self.design_froude_number,
        }

    # -- validation ---------------------------------------------------------
    def problems(self) -> tuple[list[str], list[str]]:
        """Structural errors (refuse to build) and physical warnings (report)."""
        errors: list[str] = []
        warnings: list[str] = []
        for name, value in (
            ("length", self.length),
            ("beam", self.beam),
            ("draft", self.draft),
        ):
            if not math.isfinite(value) or value <= 0.0:
                errors.append(f"{name} must be finite and positive, got {value!r}")
        if errors:
            return errors, warnings
        if not 0.30 < self.block_coefficient < 0.90:
            errors.append(
                f"block_coefficient must lie in (0.30, 0.90), got {self.block_coefficient!r}"
            )
        if not 0.0 <= self.freeboard_ratio <= 1.0:
            errors.append(f"freeboard_ratio must lie in [0, 1], got {self.freeboard_ratio!r}")
        if not 0.0 < self.transom_breadth_ratio < 1.0:
            errors.append(
                f"transom_breadth_ratio must lie in (0, 1), got {self.transom_breadth_ratio!r}"
            )
        if not 0.0 < self.bottom_breadth_ratio < 1.0:
            errors.append(
                f"bottom_breadth_ratio must lie in (0, 1), got {self.bottom_breadth_ratio!r}"
            )
        if self.section_exponent_midship <= 0.0 or self.section_exponent_ends <= 0.0:
            errors.append("section exponents must be positive")
        if self.entrance_exponent <= 0.0 or self.entrance_root <= 0.0:
            errors.append("entrance exponents must be positive")
        if self.bottom_taper_power <= 0.0:
            errors.append("bottom_taper_power must be positive")
        if self.section_relaxation_power <= 0.0:
            errors.append("section_relaxation_power must be positive")
        if not 0.0 < self.aft_body_share < 1.0:
            errors.append(f"aft_body_share must lie in (0, 1), got {self.aft_body_share!r}")
        if self.parallel_body_ratio is not None and not 0.0 <= self.parallel_body_ratio < 0.9:
            errors.append(
                f"parallel_body_ratio must lie in [0, 0.9), got {self.parallel_body_ratio!r}"
            )
        if self.degree[0] < 1 or self.degree[1] < 1:
            errors.append(f"degrees must be >= 1, got {self.degree!r}")
        if self.control_points[0] < self.degree[0] + 2 or self.control_points[1] < self.degree[1] + 2:
            errors.append(f"control_points {self.control_points!r} too small for {self.degree!r}")
        if errors:
            return errors, warnings

        # Physical sanity band of the principal dimensions: a parametric patch can be
        # built outside it, but it would not be a displacement hull.
        if not 4.0 <= self.length_to_beam <= 12.0:
            warnings.append(
                f"L/B = {self.length_to_beam:.3f} is outside the displacement-hull band [4, 12]"
            )
        if not 1.5 <= self.beam_to_draft <= 5.5:
            warnings.append(
                f"B/T = {self.beam_to_draft:.3f} is outside the displacement-hull band [1.5, 5.5]"
            )
        if not 0.45 <= self.block_coefficient <= 0.75:
            warnings.append(
                f"Cb = {self.block_coefficient:.3f} is outside the calibration band [0.45, 0.75]"
            )
        return errors, warnings

    def validate(self) -> None:
        errors, _ = self.problems()
        if errors:
            raise HullNurbsError("; ".join(errors))

    def replaced(self, **changes: Any) -> HullNurbsParameters:
        return replace(self, **changes)


def kcs_lightweight_parameters(**changes: Any) -> HullNurbsParameters:
    """The reference case: ``L = 3.0 m``, ``B = 0.4 m``, ``T = 0.15 m``, ``Cb = 0.555``."""
    return HullNurbsParameters(**changes)


# ---------------------------------------------------------------------------
# hull form
# ---------------------------------------------------------------------------
def _smoothstep(s: float) -> float:
    return s * s * (3.0 - 2.0 * s)


@dataclass(frozen=True)
class HullForm:
    """Resolved plan form: where the run, the parallel body and the entrance meet."""

    parameters: HullNurbsParameters
    parallel_body_ratio: float

    @property
    def aft_taper_ratio(self) -> float:
        return (1.0 - self.parallel_body_ratio) * self.parameters.aft_body_share

    @property
    def entrance_ratio(self) -> float:
        return (1.0 - self.parallel_body_ratio) * (1.0 - self.parameters.aft_body_share)

    @property
    def fore_start(self) -> float:
        return self.aft_taper_ratio + self.parallel_body_ratio

    def waterline_ratio(self, xi: float) -> float:
        """Half-breadth at the waterline over ``B/2`` as a function of ``xi = x/L``."""
        p = self.parameters
        aft = self.aft_taper_ratio
        fore = self.fore_start
        if xi <= 0.0:
            return p.transom_breadth_ratio
        if xi >= 1.0:
            return 0.0
        if aft > 0.0 and xi < aft:
            s = xi / aft
            return p.transom_breadth_ratio + (1.0 - p.transom_breadth_ratio) * _smoothstep(s)
        if xi <= fore:
            return 1.0
        s = (xi - fore) / (1.0 - fore) if fore < 1.0 else 1.0
        return float(max(0.0, (1.0 - s**p.entrance_exponent) ** (1.0 / p.entrance_root)))

    def section_exponent(self, xi: float) -> float:
        """Local section fullness: midship value weighted by the local breadth ratio.

        The weight is ``f^2``, not ``(1 - f)^2``: a section at half the maximum breadth is
        already a fine, V-shaped section on a real hull, so the fullness must fall away
        as soon as the plan form narrows. Relaxing only at the very ends would leave the
        entrance sections as full as the midship one and push ``Cb`` far above the
        requested value.
        """
        p = self.parameters
        weight = self.waterline_ratio(xi) ** p.section_relaxation_power
        return float(
            p.section_exponent_ends
            + (p.section_exponent_midship - p.section_exponent_ends) * weight
        )

    def bottom_ratio(self, xi: float) -> float:
        """Flat-of-bottom half-breadth, as a fraction of the local waterline breadth.

        The flat of bottom exists in the body of the hull and vanishes towards the ends:
        a stem or a post has no flat of bottom. Without this the ``zeta = 0`` row would
        keep a finite breadth all the way to the stem, which no hull does.
        """
        p = self.parameters
        return float(p.bottom_breadth_ratio * self.waterline_ratio(xi) ** p.bottom_taper_power)

    def section_ratio(self, zeta: float, xi: float) -> float:
        """Half-breadth at ``zeta = z/T`` over the local waterline half-breadth."""
        r = self.bottom_ratio(xi)
        if zeta <= 0.0:
            return r
        if zeta >= 1.0:
            # Wall-sided above the design waterline: this single patch has no flare.
            return 1.0
        n = self.section_exponent(xi)
        return float(r + (1.0 - r) * (1.0 - (1.0 - zeta) ** n) ** (1.0 / n))

    def half_breadth(self, xi: float, zeta: float) -> float:
        """The form function the NURBS control net is sampled from, in metres."""
        p = self.parameters
        return 0.5 * p.beam * self.waterline_ratio(xi) * self.section_ratio(zeta, xi)

    def as_dict(self) -> dict[str, Any]:
        return {
            "parallel_body_ratio": self.parallel_body_ratio,
            "aft_taper_ratio": self.aft_taper_ratio,
            "entrance_ratio": self.entrance_ratio,
            "fore_start_ratio": self.fore_start,
        }


# ---------------------------------------------------------------------------
# NURBS patch
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HullNurbsSurface:
    """A NURBS hull patch plus the form data it was built from."""

    parameters: HullNurbsParameters
    form: HullForm
    surface: NURBSSurface

    @property
    def deck_height(self) -> float:
        return self.parameters.deck_height

    @property
    def waterline_parameter(self) -> float:
        """``v`` at the design waterline (``z = T``)."""
        return self.parameters.draft / self.deck_height

    def patch_half_breadth(self, xi: float, zeta: float) -> float:
        """Half-breadth of the *patch* at ``(xi, zeta)`` (not of the form function)."""
        p = self.parameters
        return self.surface.evaluate(xi, zeta * p.draft / self.deck_height)[1]

    def form_half_breadth(self, xi: float, zeta: float) -> float:
        """Half-breadth of the underlying form function at ``(xi, zeta)``."""
        return self.form.half_breadth(xi, zeta)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": HULL_NURBS_SCHEMA,
            "parameters": self.parameters.as_dict(),
            "form": self.form.as_dict(),
            "degree": list(self.surface.degree),
            "control_net_shape": [
                len(self.surface.control_net),
                len(self.surface.control_net[0]),
            ],
            "rational": self.surface.rational,
            "realised_dimensions": realised_dimensions(self),
        }


def realised_dimensions(hull: HullNurbsSurface, *, samples: int = 2001) -> dict[str, float]:
    """Measure L, B and the waterline breadth of the *patch* against the request.

    The geometry gate comes before any integral: if the surface does not realise the
    principal dimensions it was asked for, every number integrated from it is describing
    a different hull. ``x`` and ``z`` are exact by construction (linear precision), but
    the beam is a spline maximum and is therefore measured, not assumed.
    """
    p = hull.parameters
    v_wl = hull.waterline_parameter
    peak = 0.0
    for k in range(samples):
        u = k / (samples - 1)
        peak = max(peak, hull.surface.evaluate(u, v_wl)[1])
    deck = 0.0
    for k in range(samples):
        u = k / (samples - 1)
        for m in range(11):
            deck = max(deck, hull.surface.evaluate(u, m / 10)[1])
    return {
        "length_m": p.length,
        "beam_at_waterline_m": 2.0 * peak,
        "beam_requested_m": p.beam,
        "beam_relative_deficit": (p.beam - 2.0 * peak) / p.beam,
        "beam_max_over_depth_m": 2.0 * deck,
        "draft_m": p.draft,
        "deck_height_m": p.deck_height,
        "waterline_area_m2": p.length * p.beam,
    }


def _build_patch(form: HullForm) -> NURBSSurface:
    p = form.parameters
    nu, nv = p.control_points
    pu, pv = p.degree
    knots_u = clamped_knots(pu, nu)
    knots_v = clamped_knots(pv, nv)
    xs = greville_abscissae(knots_u, pu)
    vs = greville_abscissae(knots_v, pv)
    deck = p.deck_height
    net: list[tuple[tuple[float, float, float], ...]] = []
    for xi in xs:
        row: list[tuple[float, float, float]] = []
        for v in vs:
            zeta = v * deck / p.draft
            row.append((p.length * xi, form.half_breadth(xi, zeta), deck * v))
        net.append(tuple(row))
    return NURBSSurface(degree=(pu, pv), knots=(knots_u, knots_v), control_net=tuple(net))


def _quadrature_volume(form: HullForm, order: tuple[int, int]) -> float:
    """Volume of the patch by Gauss quadrature on the graph ``y(u, v)``."""
    p = form.parameters
    patch = _build_patch(form)
    v_wl = p.draft / p.deck_height
    us, wu = _mapped_rule(0.0, 1.0, order[0])
    vs, wv = _mapped_rule(0.0, v_wl, order[1])
    total = 0.0
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            total += wu[i] * wv[j] * patch.evaluate(u, v)[1]
    # 2 (mirror plane) * L*D (constant Jacobian of x = L u, z = D v)
    return 2.0 * p.length * p.deck_height * total


def calibrate_form(
    parameters: HullNurbsParameters | None = None,
    *,
    tolerance: float = 1e-6,
    max_iterations: int = 80,
) -> HullForm:
    """Find the parallel-body length whose integrated ``Cb`` matches the request.

    ``Cb`` increases monotonically with the parallel-body length, so bisection is both
    safe and cheap. The calibration integral is the same Gauss quadrature used for the
    acceptance numbers, at reduced order; the tests verify that the order dependence is
    far below the calibration tolerance.
    """
    p = parameters or HullNurbsParameters()
    p.validate()
    if p.parallel_body_ratio is not None:
        return HullForm(parameters=p, parallel_body_ratio=p.parallel_body_ratio)

    target = p.design_displacement
    order = (10, 6)
    lo, hi = 0.0, 0.85
    lo_volume = _quadrature_volume(HullForm(parameters=p, parallel_body_ratio=lo), order)
    hi_volume = _quadrature_volume(HullForm(parameters=p, parallel_body_ratio=hi), order)
    if lo_volume > target:
        raise HullNurbsError(
            f"target displacement {target:.6g} m3 is below the finest form this family can "
            f"build ({lo_volume:.6g} m3): Cb is too small for these shape parameters"
        )
    if hi_volume < target:
        raise HullNurbsError(
            f"target displacement {target:.6g} m3 is above the fullest form this family can "
            f"build ({hi_volume:.6g} m3): raise the section fullness exponents"
        )
    form = HullForm(parameters=p, parallel_body_ratio=0.5 * (lo + hi))
    for _ in range(max_iterations):
        mid = 0.5 * (lo + hi)
        form = HullForm(parameters=p, parallel_body_ratio=mid)
        volume = _quadrature_volume(form, order)
        if abs(volume - target) <= tolerance * target:
            return form
        if volume < target:
            lo = mid
        else:
            hi = mid
    return form


def build_nurbs_hull(
    parameters: HullNurbsParameters | None = None,
    *,
    tolerance: float = 1e-6,
) -> HullNurbsSurface:
    """Build the calibrated NURBS hull patch."""
    form = calibrate_form(parameters, tolerance=tolerance)
    return HullNurbsSurface(parameters=form.parameters, form=form, surface=_build_patch(form))


# ---------------------------------------------------------------------------
# hydrostatics on the patch
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HullHydrostatics:
    """Hydrostatics of the NURBS patch, integrated with Gauss quadrature."""

    displaced_volume: float
    lcb_from_ap: float
    lcb_ratio_from_midship: float
    vertical_centre_of_buoyancy: float
    waterplane_area: float
    waterplane_coefficient: float
    maximum_section_area: float
    midship_section_coefficient: float
    prismatic_coefficient: float
    wetted_surface_area: float
    side_area: float
    bottom_area: float
    transom_area_below_waterline: float
    quadrature_order: tuple[int, int]

    @property
    def block_coefficient(self) -> float:
        return self.prismatic_coefficient * self.midship_section_coefficient

    def as_dict(self) -> dict[str, Any]:
        return {
            "displaced_volume_m3": self.displaced_volume,
            "lcb_from_ap_m": self.lcb_from_ap,
            "lcb_ratio_from_midship": self.lcb_ratio_from_midship,
            "vertical_centre_of_buoyancy_m": self.vertical_centre_of_buoyancy,
            "waterplane_area_m2": self.waterplane_area,
            "waterplane_coefficient": self.waterplane_coefficient,
            "maximum_section_area_m2": self.maximum_section_area,
            "midship_section_coefficient": self.midship_section_coefficient,
            "prismatic_coefficient": self.prismatic_coefficient,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "side_area_m2": self.side_area,
            "bottom_area_m2": self.bottom_area,
            "transom_area_below_waterline_m2": self.transom_area_below_waterline,
            "quadrature_order": list(self.quadrature_order),
        }


def _side_area(patch: NURBSSurface, v_lo: float, v_hi: float, order: tuple[int, int]) -> float:
    us, wu = _mapped_rule(0.0, 1.0, order[0])
    vs, wv = _mapped_rule(v_lo, v_hi, order[1])
    total = 0.0
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            d = patch.derivatives(u, v, order=1)
            du = d[(1, 0)]
            dv = d[(0, 1)]
            cx = du[1] * dv[2] - du[2] * dv[1]
            cy = du[2] * dv[0] - du[0] * dv[2]
            cz = du[0] * dv[1] - du[1] * dv[0]
            total += wu[i] * wv[j] * math.sqrt(cx * cx + cy * cy + cz * cz)
    return total


def nurbs_hydrostatics(
    hull: HullNurbsSurface,
    *,
    quadrature_order: tuple[int, int] | None = None,
) -> HullHydrostatics:
    """Hydrostatics of the NURBS patch: volume, LCB, KB, coefficients, wetted area.

    Every quantity is a Gauss-quadrature integral of the *patch*, not of the form
    function it was sampled from, so what is verified is the geometry a consumer will
    actually evaluate.
    """
    p = hull.parameters
    patch = hull.surface
    order = quadrature_order or p.quadrature_order
    v_wl = hull.waterline_parameter
    deck = hull.deck_height
    length, beam, draft = p.length, p.beam, p.draft
    jacobian_xz = length * deck  # dx dz = L D du dv

    us, wu = _mapped_rule(0.0, 1.0, order[0])
    vs, wv = _mapped_rule(0.0, v_wl, order[1])
    volume_local = 0.0
    moment_x = 0.0
    moment_z = 0.0
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            weight = wu[i] * wv[j]
            y = patch.evaluate(u, v)[1]
            volume_local += weight * y
            moment_x += weight * y * (length * u)
            moment_z += weight * y * (deck * v)
    volume = 2.0 * jacobian_xz * volume_local
    lcb = 2.0 * jacobian_xz * moment_x / volume
    vcb = 2.0 * jacobian_xz * moment_z / volume

    waterplane = 2.0 * length * sum(
        w * patch.evaluate(u, v_wl)[1] for u, w in zip(us, wu, strict=True)
    )

    # Maximum section area over a dense station grid: the maximum of a spline section
    # area is not at a quadrature node, so it is searched rather than assumed.
    vs_s, wv_s = _mapped_rule(0.0, v_wl, order[1])
    best_area = 0.0
    station_samples = 401
    for k in range(station_samples + 1):
        u = k / station_samples
        area = 2.0 * deck * sum(
            w * patch.evaluate(u, v)[1] for v, w in zip(vs_s, wv_s, strict=True)
        )
        best_area = max(best_area, area)
    midship_coefficient = best_area / (beam * draft)
    prismatic = volume / (best_area * length)

    side = _side_area(patch, 0.0, v_wl, (max(order[0], 16), max(order[1], 12)))
    bottom = 2.0 * length * sum(
        w * patch.evaluate(u, 0.0)[1] for u, w in zip(us, wu, strict=True)
    )
    transom = 2.0 * deck * sum(
        w * patch.evaluate(0.0, v)[1] for v, w in zip(vs, wv, strict=True)
    )
    return HullHydrostatics(
        displaced_volume=volume,
        lcb_from_ap=lcb,
        lcb_ratio_from_midship=(lcb - 0.5 * length) / length,
        vertical_centre_of_buoyancy=vcb,
        waterplane_area=waterplane,
        waterplane_coefficient=waterplane / (length * beam),
        maximum_section_area=best_area,
        midship_section_coefficient=midship_coefficient,
        prismatic_coefficient=prismatic,
        wetted_surface_area=2.0 * side + bottom + transom,
        side_area=side,
        bottom_area=bottom,
        transom_area_below_waterline=transom,
        quadrature_order=order,
    )


# ---------------------------------------------------------------------------
# empirical baselines
# ---------------------------------------------------------------------------
def holtrop_wetted_area(
    parameters: HullNurbsParameters,
    *,
    midship_section_coefficient: float,
    waterplane_coefficient: float,
) -> float:
    """Holtrop-Mennen (1984) statistical estimate of the bare-hull wetted area.

    ``S = L (2T + B) sqrt(CM) (0.453 + 0.4425 Cb - 0.2862 CM - 0.003467 B/T
    + 0.3696 CWP) + 2.38 A_BT / Cb``, with no bulbous-bow transverse area
    (``A_BT = 0``), which is what this single patch represents. The formula is a
    regression over a ship series that starts around ``Cb = 0.55``, so a 10 % band is
    the honest way to read it, not a 2 % one.
    """
    p = parameters
    bracket = (
        0.453
        + 0.4425 * p.block_coefficient
        - 0.2862 * midship_section_coefficient
        - 0.003467 * p.beam_to_draft
        + 0.3696 * waterplane_coefficient
    )
    return p.length * (2.0 * p.draft + p.beam) * math.sqrt(midship_section_coefficient) * bracket


def schneekluth_lcb_ratio(parameters: HullNurbsParameters) -> float:
    """Schneekluth-Bertram regression for the LCB position, as a fraction of ``L``.

    ``LCB[%] = 9.7 - 22.5 Fn - 13.5 Cb`` measured from midship, positive forward
    (Schneekluth & Bertram, *Ship Design for Efficiency and Economy*, 2nd ed., 1998).
    It is negative (aft of midship) across the block coefficients of interest.
    """
    p = parameters
    return (9.7 - 22.5 * p.design_froude_number - 13.5 * p.block_coefficient) / 100.0


# ---------------------------------------------------------------------------
# structured body mesh
# ---------------------------------------------------------------------------
def _cosine_stations(count: int) -> tuple[float, ...]:
    """Station distribution clustering at both ends, where the plan form turns."""
    return tuple(0.5 * (1.0 - math.cos(math.pi * i / count)) for i in range(count + 1))


def _graded_levels(count: int, total: float, exponent: float = 1.0) -> tuple[float, ...]:
    """``count`` intervals on ``[0, total]``, clustered towards ``0`` when ``exponent > 1``.

    ``exponent = 1`` is a uniform distribution. The grading is used to put cells where the
    section curvature is: at the bilge, i.e. near ``zeta = 0``.
    """
    if count <= 0:
        return (0.0,)
    if exponent < 1.0:
        raise HullNurbsError(f"level grading exponent must be >= 1, got {exponent!r}")
    return tuple(total * (i / count) ** exponent for i in range(count + 1))


@dataclass(frozen=True)
class HexQuality:
    """Per-cell hexahedral quality, reduced to the numbers a report quotes."""

    cells: int
    collapsed: int
    min_scaled_jacobian: float
    mean_scaled_jacobian: float
    min_aspect_ratio: float
    max_aspect_ratio: float
    pass_count: int
    good_count: int

    @property
    def valid_cells(self) -> int:
        return self.cells - self.collapsed

    @property
    def pass_fraction(self) -> float:
        return self.pass_count / self.cells if self.cells else 0.0

    @property
    def good_fraction(self) -> float:
        return self.good_count / self.cells if self.cells else 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "cells": self.cells,
            "collapsed_cells": self.collapsed,
            "valid_cells": self.valid_cells,
            "min_scaled_jacobian": self.min_scaled_jacobian,
            "mean_scaled_jacobian": self.mean_scaled_jacobian,
            "min_aspect_ratio": self.min_aspect_ratio,
            "max_aspect_ratio": self.max_aspect_ratio,
            "pass_count": self.pass_count,
            "pass_fraction": self.pass_fraction,
            "good_count": self.good_count,
            "good_fraction": self.good_fraction,
            "thresholds": {
                "scaled_jacobian_min": SCALED_JACOBIAN_MIN,
                "scaled_jacobian_good": SCALED_JACOBIAN_GOOD,
                "aspect_ratio_max": ASPECT_RATIO_MAX,
            },
        }


@dataclass(frozen=True)
class HullBodyMesh:
    """A body-fitted hexahedral mesh of the starboard half hull."""

    parameters: HullNurbsParameters
    nodes: tuple[tuple[float, float, float], ...]
    hexes: tuple[tuple[int, ...], ...]
    shape: tuple[int, int, int]
    waterline_level: int
    quality: HexQuality

    def node(self, i: int, j: int, k: int) -> tuple[float, float, float]:
        _, nz, nw = self.shape
        return self.nodes[(i * (nz + 1) + j) * (nw + 1) + k]

    def as_dict(self) -> dict[str, Any]:
        return {
            "shape_cells": list(self.shape),
            "nodes": len(self.nodes),
            "hex_cells": len(self.hexes),
            "waterline_level": self.waterline_level,
            "quality": self.quality.as_dict(),
        }


#: Hexahedron corners in the ``(di, dj, dk)`` order used throughout this module.
_HEX_CORNERS = (
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
)

#: Six-tetrahedron subdivision of a hexahedron along its 0-6 diagonal. The face
#: diagonals are consistent across shared faces, so the subdivision is conforming and
#: the tet volumes sum to the hexahedron volume.
_HEX_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)

#: For each corner, the three neighbours whose edge vectors span the corner trihedron.
_CORNER_NEIGHBOURS: tuple[tuple[int, int, int], ...] = (
    (1, 4, 3),
    (2, 5, 0),
    (3, 6, 1),
    (0, 7, 2),
    (7, 0, 5),
    (4, 1, 6),
    (5, 2, 7),
    (6, 3, 4),
)


def _sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def hex_scaled_jacobian(points: Sequence[Sequence[float]]) -> float:
    """Minimum corner scaled Jacobian of a hexahedron (1 = cube, 0 = degenerate).

    A corner whose trihedron is left-handed or flat makes the whole cell score 0, so an
    inverted corner cannot hide behind the absolute value in :func:`hex_volume`. The
    global orientation is taken from corner 0, so the metric is independent of whether a
    caller numbers the corners in the ``(di, dj, dk)`` or the ``(di, dk, dj)`` convention.
    """
    reference = _CORNER_NEIGHBOURS[0]
    origin = points[0]
    orientation = _dot(
        _sub(points[reference[0]], origin),
        _cross(_sub(points[reference[1]], origin), _sub(points[reference[2]], origin)),
    )
    sign = 1.0 if orientation >= 0.0 else -1.0
    worst = 1.0
    for corner, neighbours in enumerate(_CORNER_NEIGHBOURS):
        origin = points[corner]
        e1 = _sub(points[neighbours[0]], origin)
        e2 = _sub(points[neighbours[1]], origin)
        e3 = _sub(points[neighbours[2]], origin)
        lengths = (_norm(e1), _norm(e2), _norm(e3))
        if min(lengths) <= 0.0:
            return 0.0
        triple = sign * _dot(e1, _cross(e2, e3))
        if triple <= 0.0:
            return 0.0
        worst = min(worst, triple / (lengths[0] * lengths[1] * lengths[2]))
    return max(worst, 0.0)


def hex_characteristic_lengths(points: Sequence[Sequence[float]]) -> tuple[float, float]:
    """Shortest and longest of the twelve hexahedron edges."""
    edges = (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    )
    lengths = [_norm(_sub(points[b], points[a])) for a, b in edges]
    return min(lengths), max(lengths)


def hex_volume(points: Sequence[Sequence[float]]) -> float:
    """Volume of a hexahedron via its six-tet subdivision (absolute value)."""
    total = 0.0
    for a, b, c, d in _HEX_TETS:
        pa, pb, pc, pd = points[a], points[b], points[c], points[d]
        total += abs(_dot(_sub(pb, pa), _cross(_sub(pc, pa), _sub(pd, pa)))) / 6.0
    return total


def build_body_mesh(
    hull: HullNurbsSurface,
    *,
    stations: int = 32,
    levels_below_waterline: int = 12,
    levels_above_waterline: int = 4,
    layers: int = 6,
    vertical_grading: float = 1.2,
) -> HullBodyMesh:
    """Body-fitted hexahedral mesh of the starboard half hull.

    ``i`` runs fore-and-aft on a cosine distribution (clustered at the stem and the
    transom), ``j`` runs vertically from the bottom to the deck with the design
    waterline an *exact* grid level (so the wetted area is clipped exactly rather than by
    centroid tests), and ``k`` runs from the centreplane to the shell. Every node on the
    ``k = layers`` surface is the NURBS patch itself, so the mesh boundary *is* the
    geometry and not an approximation of it.
    """
    p = hull.parameters
    if stations < 2 or levels_below_waterline < 1 or layers < 1:
        raise HullNurbsError(
            "mesh resolution needs at least 2 stations, 1 waterline level and 1 layer"
        )
    deck = hull.deck_height
    draft = p.draft

    xs = _cosine_stations(stations)
    zetas = _graded_levels(levels_below_waterline, draft, vertical_grading)
    if levels_above_waterline > 0:
        above = _graded_levels(levels_above_waterline, deck - draft, 1.0)
        zetas = zetas + tuple(draft + z for z in above[1:])
    nz = len(zetas) - 1

    nodes: list[tuple[float, float, float]] = []
    for xi in xs:
        for z in zetas:
            y_shell = hull.surface.evaluate(xi, z / deck)[1] if deck > 0 else 0.0
            for layer in range(layers + 1):
                nodes.append((p.length * xi, y_shell * layer / layers, z))

    def node_index(i: int, j: int, k: int) -> int:
        return (i * (nz + 1) + j) * (layers + 1) + k

    hexes: list[tuple[int, ...]] = []
    for i in range(stations):
        for j in range(nz):
            for k in range(layers):
                hexes.append(
                    tuple(node_index(i + di, j + dj, k + dk) for di, dj, dk in _HEX_CORNERS)
                )

    return HullBodyMesh(
        parameters=p,
        nodes=tuple(nodes),
        hexes=tuple(hexes),
        shape=(stations, nz, layers),
        waterline_level=levels_below_waterline,
        quality=_mesh_quality(nodes, hexes),
    )


def _mesh_quality(
    nodes: Sequence[Sequence[float]], hexes: Sequence[Sequence[int]]
) -> HexQuality:
    sj_values: list[float] = []
    ar_values: list[float] = []
    collapsed = 0
    pass_count = 0
    good_count = 0
    for cell in hexes:
        points = [nodes[index] for index in cell]
        volume = hex_volume(points)
        sj = hex_scaled_jacobian(points)
        short, long = hex_characteristic_lengths(points)
        if volume <= 1e-18 or short <= 1e-15 or sj <= 0.0:
            collapsed += 1
            sj_values.append(0.0)
            continue
        aspect = long / short
        sj_values.append(sj)
        ar_values.append(aspect)
        if sj >= SCALED_JACOBIAN_MIN and aspect <= ASPECT_RATIO_MAX:
            pass_count += 1
        if sj >= SCALED_JACOBIAN_GOOD and aspect <= ASPECT_RATIO_MAX:
            good_count += 1
    return HexQuality(
        cells=len(hexes),
        collapsed=collapsed,
        min_scaled_jacobian=min(sj_values) if sj_values else 0.0,
        mean_scaled_jacobian=(sum(sj_values) / len(sj_values)) if sj_values else 0.0,
        min_aspect_ratio=min(ar_values) if ar_values else float("inf"),
        max_aspect_ratio=max(ar_values) if ar_values else float("inf"),
        pass_count=pass_count,
        good_count=good_count,
    )


def hex_to_tets(mesh: HullBodyMesh) -> tuple[tuple[tuple[int, int, int, int], ...], dict[str, float]]:
    """Six-tetrahedron subdivision of every hexahedron, with its volume agreement.

    The report dictionary carries the tet count and the relative volume mismatch between
    the tet mesh and the hex mesh, which must be zero for a conforming subdivision.
    """
    tets: list[tuple[int, int, int, int]] = []
    hex_total = 0.0
    tet_total = 0.0
    worst = 0.0
    for cell in mesh.hexes:
        points = [mesh.nodes[index] for index in cell]
        hex_total += hex_volume(points)
        local_hex = 0.0
        for a, b, c, d in _HEX_TETS:
            tets.append((cell[a], cell[b], cell[c], cell[d]))
            pa, pb, pc, pd = points[a], points[b], points[c], points[d]
            volume = abs(_dot(_sub(pb, pa), _cross(_sub(pc, pa), _sub(pd, pa)))) / 6.0
            local_hex += volume
            tet_total += volume
        local = hex_volume(points)
        if local > 0.0:
            worst = max(worst, abs(local_hex - local) / local)
    return tuple(tets), {
        "tetrahedra": float(len(tets)),
        "hexahedra": float(len(mesh.hexes)),
        "hex_volume_m3": hex_total,
        "tet_volume_m3": tet_total,
        "relative_volume_mismatch": abs(tet_total - hex_total) / hex_total if hex_total else 0.0,
        "worst_cell_relative_mismatch": worst,
    }


# ---------------------------------------------------------------------------
# integration from the mesh
# ---------------------------------------------------------------------------
def _triangle_area(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def _quad_area(points: Sequence[Sequence[float]], order: Sequence[int]) -> float:
    p0, p1, p2, p3 = (points[index] for index in order)
    return _triangle_area(p0, p1, p2) + _triangle_area(p0, p2, p3)


@dataclass(frozen=True)
class MeshIntegrals:
    """Hydrostatics re-integrated **from the body mesh** rather than from the patch.

    All quantities are for the full hull (the starboard half times the mirror factor 2),
    so they are directly comparable with :class:`HullHydrostatics`.
    """

    cells: int
    buoyant_cells: int
    displaced_volume: float
    lcb_from_ap: float
    lcb_ratio_from_midship: float
    vertical_centre_of_buoyancy: float
    wetted_surface_area: float
    side_area: float
    bottom_area: float
    transom_area: float
    deck_area: float
    symmetry_area: float
    wetted_faces: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "cells": self.cells,
            "buoyant_cells": self.buoyant_cells,
            "displaced_volume_m3": self.displaced_volume,
            "lcb_from_ap_m": self.lcb_from_ap,
            "lcb_ratio_from_midship": self.lcb_ratio_from_midship,
            "vertical_centre_of_buoyancy_m": self.vertical_centre_of_buoyancy,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "side_area_m2": self.side_area,
            "bottom_area_m2": self.bottom_area,
            "transom_area_m2": self.transom_area,
            "deck_area_m2": self.deck_area,
            "symmetry_plane_area_m2": self.symmetry_area,
            "wetted_faces": self.wetted_faces,
        }


def mesh_integrals(mesh: HullBodyMesh) -> MeshIntegrals:
    """Volume, LCB, KB and wetted area measured on the hexahedral mesh itself.

    Wetted area is the sum of the boundary faces that see water: the outer shell below
    the design waterline, the bottom (``z = 0``) and the transom below the waterline.
    The deck and the symmetry plane are not wetted and are reported separately, so the
    classification can be audited.
    """
    stations, nz, layers = mesh.shape
    wetted_levels = mesh.waterline_level
    volume = 0.0
    moment_x = 0.0
    moment_z = 0.0
    buoyant_cells = 0
    for i in range(stations):
        for j in range(wetted_levels):
            for k in range(layers):
                cell = [
                    mesh.node(i + di, j + dj, k + dk) for di, dj, dk in _HEX_CORNERS
                ]
                points = cell
                cell_volume = hex_volume(points)
                if cell_volume <= 0.0:
                    continue
                buoyant_cells += 1
                cx = sum(point[0] for point in points) / 8.0
                cz = sum(point[2] for point in points) / 8.0
                volume += cell_volume
                moment_x += cell_volume * cx
                moment_z += cell_volume * cz

    draft = mesh.parameters.draft
    tolerance = draft * (1.0 + 1e-12)
    side = 0.0
    bottom = 0.0
    transom = 0.0
    deck = 0.0
    symmetry = 0.0
    wetted_faces = 0
    for i in range(stations):
        for j in range(nz):
            for k in range(layers):
                cell = [mesh.node(i + di, j + dj, k + dk) for di, dj, dk in _HEX_CORNERS]
                if k == layers - 1 and max(cell[4][2], cell[5][2], cell[6][2], cell[7][2]) <= tolerance:
                    side += _quad_area(cell, (4, 5, 6, 7))
                    wetted_faces += 1
                if j == 0:
                    bottom += _quad_area(cell, (0, 1, 5, 4))
                    wetted_faces += 1
                if i == 0 and max(cell[0][2], cell[3][2], cell[7][2], cell[4][2]) <= tolerance:
                    transom += _quad_area(cell, (0, 3, 7, 4))
                    wetted_faces += 1
                if j == nz - 1:
                    deck += _quad_area(cell, (4, 5, 6, 7))
                if k == 0:
                    symmetry += _quad_area(cell, (0, 1, 2, 3))

    mirror = 2.0
    full_volume = mirror * volume
    lcb = mirror * moment_x / full_volume
    return MeshIntegrals(
        cells=len(mesh.hexes),
        buoyant_cells=buoyant_cells,
        displaced_volume=full_volume,
        lcb_from_ap=lcb,
        lcb_ratio_from_midship=(lcb - 0.5 * mesh.parameters.length) / mesh.parameters.length,
        vertical_centre_of_buoyancy=mirror * moment_z / full_volume,
        wetted_surface_area=mirror * (side + bottom + transom),
        side_area=mirror * side,
        bottom_area=mirror * bottom,
        transom_area=mirror * transom,
        deck_area=mirror * deck,
        symmetry_area=2.0 * symmetry,
        wetted_faces=2 * wetted_faces,
    )


# ---------------------------------------------------------------------------
# acceptance harness
# ---------------------------------------------------------------------------
def _criterion(name: str, value: float, limit: float, unit: str, detail: str) -> dict[str, Any]:
    return {
        "criterion": name,
        "value": value,
        "limit": limit,
        "unit": unit,
        "passed": bool(value <= limit),
        "detail": detail,
    }


def acceptance_report(
    parameters: HullNurbsParameters | None = None,
    *,
    refinement_factor: int = 2,
) -> dict[str, Any]:
    """Evaluate the geometry-consistency acceptance criteria C1-C7 on one hull.

    ======  =========================================================================
    C1      ``|V - Cb L B T| / (Cb L B T) <= 2 %``       displaced volume vs design
    C2      ``|S - S_holtrop| / S_holtrop <= 10 %``       wetted area vs empirical fit
    C3      ``|LCB - LCB_schneekluth| / L <= 5 %``        LCB vs empirical fit
    C4      hex cells passing the quality gate ``> 90 %``
    C5      NURBS route vs mesh route: volume, wetted area, KB ``<= 1 %`` and
            ``|dLCB| / L <= 1 %``
    C6      ``Cb`` sweep: displacement strictly increasing, calibration residual small
    C7      mesh refinement by ``refinement_factor`` in every direction changes the
            measured volume and wetted area by ``<= 2 %``
    ======  =========================================================================

    C8 (the pytest count) is not a geometric measurement; the report tool runs the suite
    and writes it into the artifact next to this block.
    """
    p = parameters or HullNurbsParameters()
    p.validate()
    hull = build_nurbs_hull(p)
    hydro = nurbs_hydrostatics(hull)

    criteria: dict[str, Any] = {}
    design = p.design_displacement
    c1 = abs(hydro.displaced_volume - design) / design
    criteria["C1"] = _criterion(
        "displaced_volume_vs_design",
        c1,
        0.02,
        "relative",
        f"V = {hydro.displaced_volume:.6f} m3 vs Cb*L*B*T = {design:.6f} m3",
    )

    s_empirical = holtrop_wetted_area(
        p,
        midship_section_coefficient=hydro.midship_section_coefficient,
        waterplane_coefficient=hydro.waterplane_coefficient,
    )
    c2 = abs(hydro.wetted_surface_area - s_empirical) / s_empirical
    criteria["C2"] = _criterion(
        "wetted_area_vs_holtrop",
        c2,
        0.10,
        "relative",
        f"S = {hydro.wetted_surface_area:.6f} m2 vs Holtrop-Mennen {s_empirical:.6f} m2 "
        f"(CM = {hydro.midship_section_coefficient:.4f}, CWP = {hydro.waterplane_coefficient:.4f})",
    )

    lcb_empirical = schneekluth_lcb_ratio(p)
    lcb_deviation = abs(hydro.lcb_ratio_from_midship - lcb_empirical)
    criteria["C3"] = _criterion(
        "lcb_vs_schneekluth",
        lcb_deviation,
        0.05,
        "fraction_of_L",
        f"LCB = {hydro.lcb_ratio_from_midship * 100:+.3f} %L vs Schneekluth "
        f"{lcb_empirical * 100:+.3f} %L (Fn = {p.design_froude_number:.3f})",
    )

    mesh = build_body_mesh(hull)
    quality = mesh.quality
    criteria["C4"] = {
        "criterion": "hex_mesh_quality",
        "value": 1.0 - quality.pass_fraction,
        "limit": 0.10,
        "unit": "failing_fraction",
        "passed": bool(quality.pass_fraction > 0.90),
        "detail": (
            f"{quality.pass_count}/{quality.cells} cells pass (sj >= {SCALED_JACOBIAN_MIN}, "
            f"AR <= {ASPECT_RATIO_MAX}); {quality.collapsed} collapsed at the stem column; "
            f"mean sj = {quality.mean_scaled_jacobian:.4f}"
        ),
        "quality": quality.as_dict(),
    }

    integrals = mesh_integrals(mesh)
    volume_error = abs(integrals.displaced_volume - hydro.displaced_volume) / hydro.displaced_volume
    area_error = abs(integrals.wetted_surface_area - hydro.wetted_surface_area) / (
        hydro.wetted_surface_area
    )
    kb_error = abs(integrals.vertical_centre_of_buoyancy - hydro.vertical_centre_of_buoyancy) / (
        hydro.vertical_centre_of_buoyancy
    )
    lcb_error = abs(integrals.lcb_from_ap - hydro.lcb_from_ap) / p.length
    worst = max(volume_error, area_error, kb_error, lcb_error)
    criteria["C5"] = _criterion(
        "nurbs_vs_mesh_consistency",
        worst,
        0.01,
        "relative",
        (
            f"volume {volume_error * 100:.4f} %, wetted area {area_error * 100:.4f} %, "
            f"KB {kb_error * 100:.4f} %, |dLCB|/L {lcb_error * 100:.4f} %"
        ),
    )
    criteria["C5"]["components"] = {
        "volume_relative": volume_error,
        "wetted_area_relative": area_error,
        "kb_relative": kb_error,
        "lcb_over_length": lcb_error,
    }

    # The calibration can only reach block coefficients the form family can build, so the
    # sweep is placed inside that band and the band itself is reported.
    band_low = _quadrature_volume(
        HullForm(parameters=p.replaced(parallel_body_ratio=0.0), parallel_body_ratio=0.0), (10, 6)
    ) / p.design_displacement * p.block_coefficient
    band_high = _quadrature_volume(
        HullForm(parameters=p.replaced(parallel_body_ratio=0.85), parallel_body_ratio=0.85), (10, 6)
    ) / p.design_displacement * p.block_coefficient
    sweep_cb = (p.block_coefficient, 0.575, 0.60, 0.625, 0.65)
    sweep = []
    monotone = True
    previous = -1.0
    for cb in sweep_cb:
        variant = p.replaced(block_coefficient=cb)
        variant_hull = build_nurbs_hull(variant)
        variant_hydro = nurbs_hydrostatics(variant_hull)
        residual = abs(variant_hydro.displaced_volume - variant.design_displacement) / (
            variant.design_displacement
        )
        monotone = monotone and variant_hydro.displaced_volume > previous
        previous = variant_hydro.displaced_volume
        sweep.append(
            {
                "block_coefficient": cb,
                "displaced_volume_m3": variant_hydro.displaced_volume,
                "design_displacement_m3": variant.design_displacement,
                "calibration_residual": residual,
                "parallel_body_ratio": variant_hull.form.parallel_body_ratio,
                "prismatic_coefficient": variant_hydro.prismatic_coefficient,
            }
        )
    # dV/dCb is analytic (V = Cb L B T), so the measured slope has a reference value.
    analytic_slope = p.length * p.beam * p.draft
    measured_slope = (sweep[-1]["displaced_volume_m3"] - sweep[0]["displaced_volume_m3"]) / (
        sweep[-1]["block_coefficient"] - sweep[0]["block_coefficient"]
    )
    slope_error = abs(measured_slope - analytic_slope) / analytic_slope
    criteria["C6"] = {
        "criterion": "block_coefficient_monotonicity",
        "value": 0.0 if monotone else 1.0,
        "limit": 0.0,
        "unit": "boolean",
        "passed": bool(monotone),
        "detail": "displacement strictly increases over the Cb sweep: "
        + ", ".join(
            f"{item['block_coefficient']:.3f}->{item['displaced_volume_m3']:.6f}" for item in sweep
        ),
        "sweep": sweep,
        "reachable_block_coefficient_band": [band_low, band_high],
        "analytic_dv_dcb_m3": analytic_slope,
        "measured_dv_dcb_m3": measured_slope,
        "slope_relative_error": slope_error,
    }

    levels_above = mesh.shape[1] - mesh.waterline_level
    refined = build_body_mesh(
        hull,
        stations=mesh.shape[0] * refinement_factor,
        levels_below_waterline=mesh.waterline_level * refinement_factor,
        levels_above_waterline=max(levels_above * refinement_factor, 1),
        layers=mesh.shape[2] * refinement_factor,
    )
    refined_integrals = mesh_integrals(refined)
    refined_volume_change = abs(
        refined_integrals.displaced_volume - integrals.displaced_volume
    ) / integrals.displaced_volume
    refined_area_change = abs(
        refined_integrals.wetted_surface_area - integrals.wetted_surface_area
    ) / integrals.wetted_surface_area
    refined_vs_nurbs = abs(
        refined_integrals.displaced_volume - hydro.displaced_volume
    ) / hydro.displaced_volume
    criteria["C7"] = _criterion(
        "mesh_refinement_stability",
        max(refined_volume_change, refined_area_change),
        0.02,
        "relative",
        (
            f"{mesh.shape} -> {refined.shape}: volume {refined_volume_change * 100:.4f} %, "
            f"wetted area {refined_area_change * 100:.4f} %; refined-vs-NURBS volume "
            f"{refined_vs_nurbs * 100:.4f} %"
        ),
    )
    criteria["C7"]["components"] = {
        "volume_relative_change": refined_volume_change,
        "wetted_area_relative_change": refined_area_change,
        "refined_volume_m3": refined_integrals.displaced_volume,
        "refined_wetted_area_m2": refined_integrals.wetted_surface_area,
        "coarse_volume_m3": integrals.displaced_volume,
        "coarse_wetted_area_m2": integrals.wetted_surface_area,
    }

    _, tet_report = hex_to_tets(mesh)
    passed = all(criteria[key]["passed"] for key in ("C1", "C2", "C3", "C4", "C5", "C6", "C7"))
    return {
        "schema": HULL_NURBS_SCHEMA,
        "parameters": p.as_dict(),
        "hull": hull.as_dict(),
        "hydrostatics": hydro.as_dict(),
        "mesh": mesh.as_dict(),
        "mesh_integrals": integrals.as_dict(),
        "tetrahedra": tet_report,
        "criteria": criteria,
        "all_passed": passed,
    }
