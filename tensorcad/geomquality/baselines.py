"""Reference hull forms as CAD patch sets: Wigley III, Series 60, a KVLCC2 class.

Each baseline exists twice, and that duplication is the verification: an **analytic**
definition (closed-form half-breadth plus :func:`reference_integrals`, an independent
composite-Gauss quadrature) and a **CAD** definition (a tensor-product B-spline patch
set with sewn topology -- what a design office hands to a solver).  The pipeline
reports numbers from the CAD definition; the analytic one is the reference it is
judged against.

Vertical parameterisation
-------------------------
The Series 60 / KVLCC2 section shape is a power law ``y ~ zeta**(1/m)`` whose
derivative at the keel is unbounded, which no spline in ``z`` can represent
cheaply.  The patch therefore uses the parameter ``v`` with ``z = -T + T*v**m`` and
``y = a(x)*v``: ``z(v)`` is a polynomial of degree ``m`` (exactly one Bezier span)
and ``y`` is *linear* in ``v``, so the vertical direction is exact and all
representation error sits in the one-dimensional plan form ``a(x)``, where it can be
measured directly.  Wigley III is polynomial in both directions and is exact.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

from ..nurbs import NURBSSurface, gauss_legendre
from .brep import BRep, Face, mirror, sew

Vec3 = tuple[float, float, float]


# --------------------------------------------------------------------------- #
# small dense linear algebra
# --------------------------------------------------------------------------- #
def _solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting."""
    n = len(rhs)
    a = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(a[r][column]))
        if abs(a[pivot][column]) < 1e-300:
            raise ValueError("singular interpolation matrix")
        a[column], a[pivot] = a[pivot], a[column]
        for row in range(column + 1, n):
            factor = a[row][column] / a[column][column]
            if factor == 0.0:
                continue
            for k in range(column, n + 1):
                a[row][k] -= factor * a[column][k]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        x[row] = (a[row][n] - sum(a[row][k] * x[k] for k in range(row + 1, n))) / a[row][row]
    return x


def greville_abscissae(knots: Sequence[float], degree: int) -> list[float]:
    count = len(knots) - degree - 1
    return [sum(knots[i + 1 : i + degree + 1]) / degree for i in range(count)]


def basis_row(knots: Sequence[float], degree: int, u: float) -> list[float]:
    """Nonzero basis values at ``u``, padded to the full control-point count.

    Piegl & Tiller A2.2: the triangular scheme needs no recursion table and cannot
    index outside the ``degree + 1``-wide window it is given.
    """
    count = len(knots) - degree - 1
    values = [0.0] * count
    if u >= knots[-1]:
        span = count - 1
    else:
        span = degree
        for i in range(degree, count):
            if knots[i] <= u < knots[i + 1]:
                span = i
                break
    n = [0.0] * (degree + 1)
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    n[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = n[r] / denominator if denominator != 0.0 else 0.0
            n[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        n[j] = saved
    for i in range(degree + 1):
        values[span - degree + i] = n[i]
    return values


def _knot_vector(points: int, degree: int, params: Sequence[float] | None = None) -> list[float]:
    if params is None:
        interior = [i / (points - degree) for i in range(1, points - degree)]
    else:
        interior = list(params)
    return [0.0] * (degree + 1) + interior + [1.0] * (degree + 1)


def fit_plan_form(
    plan_form: Callable[[float], float], *, degree: int = 3, spans: int = 48
) -> tuple[list[float], list[float]]:
    """Interpolate ``plan_form(u)`` on ``[0, 1]``; returns ``(knots, controls)``.

    Knots are graded toward midship (``u = 0.5``), where the Series 60 plan form has
    a ``|x|**n`` kink in its second derivative.
    """
    points = spans + degree
    half = spans // 2
    graded = [0.5 * (i / half) ** 1.6 for i in range(1, half)]
    params = [0.5 - value for value in reversed(graded)] + [0.5] + [
        0.5 + value for value in graded
    ]
    knots = _knot_vector(points, degree, params)
    greville = greville_abscissae(knots, degree)
    matrix = [basis_row(knots, degree, u) for u in greville]
    return knots, _solve(matrix, [plan_form(u) for u in greville])


# --------------------------------------------------------------------------- #
# patch construction
# --------------------------------------------------------------------------- #
def _bezier_coefficients(
    function: Callable[[float], float], degree: int
) -> list[float]:
    """Bernstein coefficients of a polynomial ``function(v)`` of degree <= ``degree``.

    Solved by collocation at the Greville abscissae of the single Bezier span, which
    is exact for polynomials of that degree (and the reason every vertical shape below
    is represented exactly, not fitted).
    """
    knots = [0.0] * (degree + 1) + [1.0] * (degree + 1)
    greville = greville_abscissae(knots, degree)
    matrix = [basis_row(knots, degree, v) for v in greville]
    return _solve(matrix, [function(v) for v in greville])


def hull_patch(
    plan_form: Callable[[float], float],
    *,
    half_length: float,
    depth: float,
    depth_power: int,
    profile_of_v: Callable[[float], float],
    z_of_v: Callable[[float], float],
    profile_degree: int,
    degree_u: int = 3,
    spans_u: int = 48,
) -> tuple[NURBSSurface, list[float], list[float], list[float]]:
    """Starboard hull patch: ``x = L(u-1/2)``, ``y = a(x) h(v)``, ``z = z(v)``.

    ``h`` is the section ordinate normalised to 1 at the waterline and ``z(v)`` the
    depth mapping; both are represented exactly as one Bezier span whose degree covers
    the larger of the two (``depth_power`` for a power-law section, 2 for Wigley).
    """
    def plan_form_u(u: float) -> float:
        return plan_form(half_length * (2.0 * u - 1.0))

    knots_u, controls_u = fit_plan_form(plan_form_u, degree=degree_u, spans=spans_u)
    degree_v = max(depth_power, profile_degree)
    knots_v = [0.0] * (degree_v + 1) + [1.0] * (degree_v + 1)
    w = _bezier_coefficients(profile_of_v, degree_v)
    z_controls = _bezier_coefficients(z_of_v, degree_v)
    greville_u = greville_abscissae(knots_u, degree_u)
    net = []
    for index, gu in enumerate(greville_u):
        x = half_length * (2.0 * gu - 1.0)
        net.append(
            [(x, controls_u[index] * w[j], z_controls[j]) for j in range(degree_v + 1)]
        )
    surface = NURBSSurface(
        degree=(degree_u, degree_v),
        control_net=net,
        weights=[[1.0] * (degree_v + 1) for _ in net],
        knots=[knots_u, knots_v],
    )
    return surface, knots_u, controls_u, greville_u


def deck_patch(
    controls_u: Sequence[float],
    knots_u: Sequence[float],
    *,
    degree_u: int,
    half_length: float,
) -> NURBSSurface:
    """Waterplane half patch ``y = t*a(x)`` at ``z = 0``; closes the model."""
    greville_u = greville_abscissae(knots_u, degree_u)
    net = []
    for index, gu in enumerate(greville_u):
        x = half_length * (2.0 * gu - 1.0)
        net.append([(x, 0.0, 0.0), (x, controls_u[index], 0.0)])
    return NURBSSurface(
        degree=(degree_u, 1),
        control_net=net,
        weights=[[1.0, 1.0] for _ in net],
        knots=[list(knots_u), [0.0, 0.0, 1.0, 1.0]],
    )


def build_hull(
    plan_form: Callable[[float], float],
    *,
    length: float,
    depth: float,
    depth_power: int,
    profile_of_v: Callable[[float], float],
    profile_degree: int = 1,
    degree_u: int = 3,
    spans_u: int = 48,
) -> tuple[BRep, dict[str, Any]]:
    """Half model (hull + deck) mirrored into a closed body, plus its build settings."""
    half_length = 0.5 * length
    def z_of_v(v: float) -> float:
        return -depth + depth * v**depth_power

    hull, knots_u, controls_u, _ = hull_patch(
        plan_form,
        half_length=half_length,
        depth=depth,
        depth_power=depth_power,
        profile_of_v=profile_of_v,
        z_of_v=z_of_v,
        profile_degree=profile_degree,
        degree_u=degree_u,
        spans_u=spans_u,
    )
    deck = deck_patch(controls_u, knots_u, degree_u=degree_u, half_length=half_length)
    # The starboard hull patch has S_u x S_v pointing inboard (-y), so its face
    # carries same_sense = False; the deck patch's natural normal already points up.
    half = sew(
        [hull, deck],
        names=["hull_half", "deck_half"],
        same_sense=[False, True],
        source="baseline-half",
    )
    faces = tuple(
        Face(
            index=face.index,
            surface=face.surface,
            boundary=face.boundary,
            same_sense=face.same_sense,
            name=face.name,
            wetted=(face.name != "deck_half"),
        )
        for face in half.faces
    )
    half = BRep(faces=faces, edges=half.edges, source="baseline-half")
    full = mirror(half, axis=1, plane=0.0)
    settings = {
        "length": length,
        "depth": depth,
        "depth_power": depth_power,
        "profile_degree": max(depth_power, profile_degree),
        "degree_u": degree_u,
        "spans_u": spans_u,
        "faces": len(full.faces),
        "edges": len(full.edges),
    }
    return full, settings


# --------------------------------------------------------------------------- #
# analytic reference
# --------------------------------------------------------------------------- #
def _composite_gauss(
    function: Callable[[float], float], low: float, high: float, panels: int = 12, order: int = 12
) -> float:
    nodes, weights = gauss_legendre(order)
    total = 0.0
    for panel in range(panels):
        a = low + (high - low) * panel / panels
        b = low + (high - low) * (panel + 1) / panels
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        total += half * sum(
            weight * function(mid + half * node) for node, weight in zip(nodes, weights)
        )
    return total


def reference_integrals(
    y_of_x_zeta: Callable[[float, float], float],
    *,
    length: float,
    depth: float,
    depth_power: int,
    panels_x: int = 12,
    order: int = 12,
) -> dict[str, Any]:
    """Half-model displacement, centres and waterplane area from the analytic form.

    ``zeta`` runs from 0 at the keel to 1 at the waterline; the substitution
    ``zeta = u**depth_power`` turns the power-law section into a polynomial in ``u``,
    so the quadrature converges fast and is independent of any mesh or patch set.
    """
    half_length = 0.5 * length
    power = depth_power
    nodes, weights = gauss_legendre(order)

    def integrate(weight: Callable[[float, float], float]) -> float:
        total = 0.0
        for panel in range(panels_x):
            a = -half_length + length * panel / panels_x
            b = -half_length + length * (panel + 1) / panels_x
            half = 0.5 * (b - a)
            mid = 0.5 * (a + b)
            for node, weight_x in zip(nodes, weights):
                x = mid + half * node
                for node_u, weight_u in zip(nodes, weights):
                    u = 0.5 * (node_u + 1.0)
                    total += weight_x * weight_u * half * 0.5 * weight(x, u)
        return total

    def volume_integrand(x: float, u: float) -> float:
        return y_of_x_zeta(x, u**power) * (power * u ** (power - 1) * depth)

    half_volume = integrate(volume_integrand)
    half_moment_x = integrate(lambda x, u: volume_integrand(x, u) * x)
    half_moment_z = integrate(
        lambda x, u: volume_integrand(x, u) * (-depth + depth * u**power)
    )
    half_waterplane = _composite_gauss(
        lambda x: y_of_x_zeta(x, 1.0), -half_length, half_length, panels=12, order=order
    )
    volume = 2.0 * half_volume
    return {
        "volume": volume,
        "lcb_x": 2.0 * half_moment_x / volume,
        "vcb_z": 2.0 * half_moment_z / volume,
        "waterplane_area": 2.0 * half_waterplane,
        "method": "analytic half-breadth + composite Gauss (independent of the patch set)",
    }


# --------------------------------------------------------------------------- #
# hull families
# --------------------------------------------------------------------------- #
@dataclass
class HullSpecimen:
    """A reference hull: its CAD model, its analytic twin and the deviation between."""

    name: str
    brep: BRep
    particulars: dict[str, Any]
    reference: dict[str, Any]
    plan_form: Callable[[float], float]
    y_of_x_zeta: Callable[[float, float], float]
    depth_power: int
    settings: dict[str, Any] = field(default_factory=dict)
    notes: str = ""

    @property
    def length(self) -> float:
        return self.particulars["length"]

    @property
    def beam(self) -> float:
        return self.particulars["beam"]

    @property
    def draft(self) -> float:
        return self.particulars["draft"]

    def deviation(self, samples: int = 61, zeta_samples: int = 41) -> dict[str, Any]:
        """Max/RMS of ``|y_patch - y_analytic|``, in metres and as a fraction of ``L``."""
        hull = [face for face in self.brep.faces if face.name.endswith("hull_half")][0]
        debug = 0.0
        worst = None
        half_length = 0.5 * self.length
        power = self.depth_power
        squares = 0.0
        for i in range(samples):
            u = i / (samples - 1)
            x = half_length * (2.0 * u - 1.0)
            for j in range(zeta_samples):
                zeta = j / (zeta_samples - 1)
                v = zeta ** (1.0 / power)
                y_patch = hull.surface.evaluate(u, v)[1]
                y_exact = self.y_of_x_zeta(x, zeta)
                error = abs(y_patch - y_exact)
                squares += error * error
                if error > debug:
                    debug = error
                    worst = {"x": x, "zeta": zeta, "y_patch": y_patch, "y_exact": y_exact}
        count = samples * zeta_samples
        return {
            "max_abs": debug,
            "max_fraction_of_L": debug / self.length,
            "rms": math.sqrt(squares / count),
            "worst_at": worst,
            "samples": count,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "particulars": dict(self.particulars),
            "reference": dict(self.reference),
            "depth_power": self.depth_power,
            "settings": dict(self.settings),
            "notes": self.notes,
            "brep": {
                "faces": len(self.brep.faces),
                "edges": len(self.brep.edges),
                "control_points": sum(
                    len(face.surface.control_net) * len(face.surface.control_net[0])
                    for face in self.brep.faces
                ),
            },
        }


def wigley_iii(
    *, length: float = 10.0, beam: float = 1.0, draft: float = 0.625, spans_u: int = 8
) -> HullSpecimen:
    """Wigley III: ``y = (B/2)(1-(2x/L)^2)(1-(z/T)^2)``, Cb = 4/9, Cwp = 2/3."""
    half_length, half_beam = 0.5 * length, 0.5 * beam

    def plan_form(x: float) -> float:
        return half_beam * (1.0 - (x / half_length) ** 2)

    def y_of_x_zeta(x: float, zeta: float) -> float:
        return half_beam * (1.0 - (x / half_length) ** 2) * (1.0 - (zeta - 1.0) ** 2)

    brep, settings = build_hull(
        plan_form,
        length=length,
        depth=draft,
        depth_power=1,
        profile_of_v=lambda v: 2.0 * v - v * v,
        profile_degree=2,
        degree_u=2,
        spans_u=spans_u,
    )
    reference = reference_integrals(y_of_x_zeta, length=length, depth=draft, depth_power=1)
    reference["closed_form"] = {
        "volume": 4.0 / 9.0 * length * beam * draft,
        "waterplane_area": 2.0 / 3.0 * length * beam,
        # 3/8 of the draft below the waterline (z measured up from the waterplane).
        "vcb_z": -3.0 / 8.0 * draft,
        "lcb_x": 0.0,
    }
    reference["block_coefficient"] = 4.0 / 9.0
    particulars = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "block_coefficient": 4.0 / 9.0,
        "waterplane_coefficient": 2.0 / 3.0,
        "midship_coefficient": 2.0 / 3.0,
        "lcb_fraction": 0.0,
        "source": "Wigley III parabolic hull (L/B = 10, B/T = 1.6), Cb = 4/9",
    }
    return HullSpecimen(
        name="wigley_iii",
        brep=brep,
        particulars=particulars,
        reference=reference,
        plan_form=plan_form,
        y_of_x_zeta=y_of_x_zeta,
        depth_power=1,
        settings=settings,
        notes="polynomial in x and z: the B-spline patch set is an exact representation",
    )


_SERIES60_PARENTS: dict[float, float] = {0.60: 5.0 / 3.0, 0.65: 0.670 / (1.0 - 0.670)}


def series60(
    cb: float = 0.60,
    *,
    length: float = 10.0,
    beam: float = 1.0,
    draft: float = 0.625,
    spans_u: int = 40,
) -> HullSpecimen:
    """Series 60 parent form ``y = (B/2)(1-|2xi-1|^n) zeta**(1/m)`` with closed forms.

    ``Cp = n/(n+1)``, ``Cm = m/(m+1)``, ``Cb = Cp*Cm``, ``Cwp = Cp`` and
    ``VCB/T = (m+1)/(2m+1)`` hold exactly for this parameterisation.
    """
    if cb not in _SERIES60_PARENTS:
        raise ValueError(f"Series 60 parent forms available: {sorted(_SERIES60_PARENTS)}")
    exponent = _SERIES60_PARENTS[cb]
    cp = exponent / (exponent + 1.0)
    # Cb = Cp*Cm with Cm = m/(m+1), so m = Cb/(Cp - Cb); the section exponent must be an
    # integer for z = -T + T*v**m to be an exact Bezier span, so it is rounded and the
    # realised Cb (Cp*m/(m+1)) is what the reference below reports.
    power = int(round(cb / (cp - cb)))
    cm = power / (power + 1.0)
    realised_cb = cp * cm
    half_beam = 0.5 * beam

    def plan_form(x: float) -> float:
        xi = 0.5 + x / length
        return half_beam * (1.0 - abs(2.0 * xi - 1.0) ** exponent)

    def y_of_x_zeta(x: float, zeta: float) -> float:
        return plan_form(x) * zeta ** (1.0 / power)

    brep, settings = build_hull(
        plan_form,
        length=length,
        depth=draft,
        depth_power=int(power),
        profile_of_v=lambda v: v,
        profile_degree=1,
        degree_u=3,
        spans_u=spans_u,
    )
    reference = reference_integrals(
        y_of_x_zeta, length=length, depth=draft, depth_power=int(power)
    )
    reference["closed_form"] = {
        "volume": realised_cb * length * beam * draft,
        "waterplane_area": cp * length * beam,
        "vcb_z": -draft + (power + 1.0) / (2.0 * power + 1.0) * draft,
        "lcb_x": 0.0,
        "cp": cp,
        "cm": cm,
    }
    reference["block_coefficient"] = realised_cb
    particulars = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "block_coefficient": realised_cb,
        "waterplane_coefficient": cp,
        "midship_coefficient": cm,
        "lcb_fraction": 0.0,
        "series60_exponents": {
            "n": exponent,
            "m": power,
            "requested_cb": cb,
            "realised_cb": realised_cb,
        },
        "source": f"Series 60 parent form Cb = {cb:.2f} (Todd 1951 closed forms)",
    }
    return HullSpecimen(
        name=f"series60_cb{cb:.2f}",
        brep=brep,
        particulars=particulars,
        reference=reference,
        plan_form=plan_form,
        y_of_x_zeta=y_of_x_zeta,
        depth_power=int(power),
        settings=settings,
        notes="vertical direction exact (z = -T + T*v**m, y linear in v)",
    )


def kvlcc2_class(
    *,
    length: float = 320.0,
    beam: float = 58.0,
    draft: float = 20.8,
    block_coefficient: float = 0.8098,
    lcb_fraction: float = -0.0348,
    depth_power: int = 10,
    spans_u: int = 96,
) -> HullSpecimen:
    """A KVLCC2-class tanker form calibrated to the published particulars.

    The official KVLCC2 offset table is not in this repository, so the form is a math
    hull: plan form ``1-(1-2xi)**a_aft`` aft of midship and ``1-(2xi-1)**a_fore``
    forward, section ``zeta**(1/m)``.  The two exponents are solved so the analytic Cb
    and LCB equal the published KVLCC2 values (0.8098, -3.48 % L).  The Cb check is
    therefore a consistency check of the *representation*, not of KVLCC2 itself.
    """
    power = float(depth_power)
    section_mean = power / (power + 1.0)
    area_target = block_coefficient / section_mean
    moment_target = (lcb_fraction + 0.5) * area_target

    def area(a_aft: float, a_fore: float) -> float:
        return 1.0 - 0.5 / (a_aft + 1.0) - 0.5 / (a_fore + 1.0)

    def moment(a_aft: float, a_fore: float) -> float:
        return 0.25 * (0.5 - 1.0 / (a_aft + 1.0) + 1.0 / (a_aft + 2.0)) + 0.25 * (
            1.5 - 1.0 / (a_fore + 1.0) - 1.0 / (a_fore + 2.0)
        )

    def fore_from_aft(a_aft: float) -> float:
        denominator = 2.0 * (1.0 - area_target - 0.5 / (a_aft + 1.0))
        return 1.0 / denominator - 1.0 if denominator > 0 else float("inf")

    low, high = 1.001, 200.0
    for _ in range(200):
        mid = 0.5 * (low + high)
        a_fore = fore_from_aft(mid)
        if not math.isfinite(a_fore) or a_fore <= 1.0:
            low = mid
            continue
        if moment(mid, a_fore) > moment_target:
            low = mid
        else:
            high = mid
    a_aft = 0.5 * (low + high)
    a_fore = fore_from_aft(a_aft)
    half_beam = 0.5 * beam

    def plan_form(x: float) -> float:
        xi = 0.5 + x / length
        if xi <= 0.5:
            return half_beam * (1.0 - (1.0 - 2.0 * xi) ** a_aft)
        return half_beam * (1.0 - (2.0 * xi - 1.0) ** a_fore)

    def y_of_x_zeta(x: float, zeta: float) -> float:
        return plan_form(x) * zeta ** (1.0 / power)

    brep, settings = build_hull(
        plan_form,
        length=length,
        depth=draft,
        depth_power=depth_power,
        profile_of_v=lambda v: v,
        profile_degree=1,
        degree_u=3,
        spans_u=spans_u,
    )
    reference = reference_integrals(
        y_of_x_zeta, length=length, depth=draft, depth_power=depth_power
    )
    realised_cb = reference["volume"] / (length * beam * draft)
    reference["closed_form"] = {
        "volume": block_coefficient * length * beam * draft,
        "block_coefficient": block_coefficient,
        "lcb_fraction": lcb_fraction,
        "note": "published KVLCC2 particulars used to calibrate the plan-form exponents",
    }
    particulars = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "block_coefficient": block_coefficient,
        "lcb_fraction": lcb_fraction,
        "calibrated_exponents": {"a_aft": a_aft, "a_fore": a_fore, "m": depth_power},
        "realised_block_coefficient": realised_cb,
        "source": (
            "KVLCC2-class math hull calibrated to published L/B/T/Cb/LCB "
            "(no official offset table in repo)"
        ),
    }
    return HullSpecimen(
        name="kvlcc2_class",
        brep=brep,
        particulars=particulars,
        reference=reference,
        plan_form=plan_form,
        y_of_x_zeta=y_of_x_zeta,
        depth_power=depth_power,
        settings=settings,
        notes=(
            "|Cb and LCB are reproduced by construction: the block-coefficient check is "
            "a representation consistency test, not an independent KVLCC2 validation"
        ),
    )


def specimens() -> dict[str, Callable[[], HullSpecimen]]:
    """The three reference hulls the pipeline must handle without tuning."""
    return {"wigley_iii": wigley_iii, "series60_cb0.60": series60, "kvlcc2_class": kvlcc2_class}
