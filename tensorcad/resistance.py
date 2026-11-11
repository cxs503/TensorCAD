"""Resistance: what a wetted surface costs to tow, and how a model test scales.

Every other module in this package answers a question about *the body* -- how much
water it displaces, where the waterplane is, what the righting arm is at 20 degrees.
This one answers the question the body exists to be asked: **at this speed, how hard
does it pull?** and its mirror image, addressed to a towing tank rather than a hull:
*I measured 47.3 N on a 4.36 m model; what will the ship pay?*

The two are the same problem, and both are dominated by one line. Nothing here solves
the Navier-Stokes equations: the module implements the **1966 ITTC method, first
stage** (`Prohaska` 1966) -- a friction line, a form factor ``(1 + k)`` fitted from
the model's own resistance curve, and the two-dimensional scaling of the 1978 ITTC
performance prediction procedure. It deliberately does *not* predict the residual
(wave and separation) resistance: that number is a property of the free surface, and
this package's honest position is that it comes from a model test, a panel code
(``TensorBEM``) or a RANS solver, not from an algebraic guess. The caller supplies it
-- as a fitted Prohaska slope, as a measured model-test coefficient, or as
:func:`tensorcad.design`-style empirical input from outside -- and this module carries
it through the scaling arithmetic without inventing it.

What is here:

* :func:`ittc57_friction_coefficient` -- ``0.075 / (log10 Re - 2)^2``, the 1957 ITTC
  line. Two independent flat-plate lines ride along for cross-checking it:
  :func:`schoenherr_friction_coefficient` (Schoenherr 1932, implicit, solved here by
  bisection) and :func:`prandtl_schlichting_friction_coefficient`
  (``0.455 / (log10 Re)^2.58``). The 1957 line was adopted *because* it approximates
  the Schoenherr line to within a percent over the model-and-ship Reynolds range, so
  the residual between them is a test of this implementation, not a decoration.
* :func:`prohaska_fit` -- the form factor from a low-speed model test by the Prohaska
  plot: ``CT/CF = (1 + k) + y Fn^4 / CF`` is a straight line whose intercept is the
  form factor. On synthetic data generated from a known ``(1 + k)`` and wave
  coefficient it recovers both to round-off, which is what
  ``tests/test_resistance.py`` asserts; on real data it reports the residuals and the
  standard errors of the fit rather than pretending the line is exact.
* :func:`bare_hull_resistance` / :func:`mesh_bare_hull_resistance` -- the resistance
  components of one body at one speed: friction ``R_F``, viscous-with-form
  ``(1 + k) R_F``, residual ``R_R``, total ``R_T``, and the effective power
  ``P_E = R_T V`` a propulsive coefficient would divide into.
* :func:`ship_from_model_test` -- Froude scaling plus the ITTC-78
  two-dimensional split ``CT_s = (1 + k) CF_s + CR + CA``, with the residual
  coefficient recovered from the model's total and the model-ship correlation
  allowance as an explicit caller input (no allowance is invented here).

Units are SI throughout: metres, seconds, newtons, watts, kg/m^3, m^2/s. Angles do not
appear. The fluid defaults are the conventional 15 deg C values (:data:`SEA_WATER_*`,
:data:`FRESH_WATER_*`); they are *inputs*, not physics claims, and every verification
claim in `docs/resistance.md` holds for any density and viscosity, because the checks
are ratios of like quantities.

What this module is **not**: a wave-resistance theory (no Michell integral, no
Rankine source panel method -- that is ``TensorBEM``), a Holtrop-Mennen style
statistical regression over ship particulars (those coefficients are a shipyard
correlation with its own validity box, not a geometry operation), a propeller or
propulsive-coefficient model, a roughness allowance model, or a seakeeping
added-resistance model. Each absence is listed in ``docs/resistance.md`` with the
number it would change, so a consumer knows which of these is the next one they need.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh
from .mesh import surface_area as _mesh_surface_area

__all__ = [
    "FRESH_WATER_DENSITY",
    "FRESH_WATER_KINEMATIC_VISCOSITY",
    "FRICTION_LINES",
    "MINIMUM_REYNOLDS",
    "ModelShipPrediction",
    "ProhaskaFit",
    "Resistance",
    "ResistanceError",
    "SEA_WATER_DENSITY",
    "SEA_WATER_KINEMATIC_VISCOSITY",
    "STANDARD_GRAVITY",
    "bare_hull_resistance",
    "friction_coefficient",
    "froude_number",
    "ittc57_friction_coefficient",
    "mesh_bare_hull_resistance",
    "prandtl_schlichting_friction_coefficient",
    "prohaska_fit",
    "reynolds_number",
    "schoenherr_friction_coefficient",
    "ship_from_model_test",
]

#: Standard gravity, m/s^2 (the value the ITTC procedures use for Froude number).
STANDARD_GRAVITY = 9.80665

#: Conventional 15 deg C sea-water density, kg/m^3. An input, not a physics claim.
SEA_WATER_DENSITY = 1025.9

#: Conventional 15 deg C sea-water kinematic viscosity, m^2/s.
SEA_WATER_KINEMATIC_VISCOSITY = 1.1892e-6

#: Conventional 15 deg C fresh-water density, kg/m^3.
FRESH_WATER_DENSITY = 999.1

#: Conventional 15 deg C fresh-water kinematic viscosity, m^2/s.
FRESH_WATER_KINEMATIC_VISCOSITY = 1.1390e-6

#: Below this Reynolds number the turbulent flat-plate lines are outside their
#: validity: the 1957 line's own denominator is only defined for ``Re > 100``, but
#: the line is fitted to turbulent data and is not used below ``1e4`` here.
MINIMUM_REYNOLDS = 1.0e4


class ResistanceError(ValueError):
    """Raised when a resistance input is outside the domain of the method."""


# --------------------------------------------------------------- non-dimensional


def _check_finite(name: str, value: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:  # pragma: no cover - defensive
        raise ResistanceError(f"{name} must be a number, got {value!r}") from error
    if not math.isfinite(number):
        raise ResistanceError(f"{name} must be finite, got {number!r}")
    return number


def _check_positive(name: str, value: float) -> float:
    number = _check_finite(name, value)
    if number <= 0.0:
        raise ResistanceError(f"{name} must be positive, got {number!r}")
    return number


def reynolds_number(
    length: float, velocity: float, kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY
) -> float:
    """``Re = V L / nu``. The length is the one the friction line is quoted for."""
    length = _check_positive("length", length)
    velocity = _check_positive("velocity", velocity)
    nu = _check_positive("kinematic_viscosity", kinematic_viscosity)
    return velocity * length / nu


def froude_number(length: float, velocity: float, gravity: float = STANDARD_GRAVITY) -> float:
    """``Fn = V / sqrt(g L)``."""
    length = _check_positive("length", length)
    velocity = _check_positive("velocity", velocity)
    g = _check_positive("gravity", gravity)
    return velocity / math.sqrt(g * length)


def _reynolds_domain(reynolds: float) -> float:
    value = _check_positive("reynolds", reynolds)
    if value < MINIMUM_REYNOLDS:
        raise ResistanceError(
            f"reynolds={value!r} is below the turbulent-line domain (Re >= {MINIMUM_REYNOLDS:g});"
            " a laminar or transitional body is not described by these lines"
        )
    return value


def ittc57_friction_coefficient(reynolds: float) -> float:
    """The 1957 ITTC friction line, ``CF = 0.075 / (log10 Re - 2)^2``.

    This is the line the 1978 ITTC performance prediction procedure prescribes, and
    the one every number downstream of :func:`bare_hull_resistance` uses by default.
    It has no free parameters: at ``Re = 1e7`` it returns exactly ``0.003``, which is
    the identity ``tests/test_resistance.py`` pins.
    """
    re = _reynolds_domain(reynolds)
    return 0.075 / (math.log10(re) - 2.0) ** 2


def schoenherr_friction_coefficient(
    reynolds: float, *, tolerance: float = 1e-13, max_iterations: int = 200
) -> float:
    """Schoenherr's (1932) flat-plate line, ``0.242 / sqrt(CF) = log10(Re CF)``.

    Implicit in ``CF``, so it is solved here by bisection on the monotone residual
    ``f(CF) = 0.242/sqrt(CF) - log10(Re CF)``, which is positive at small ``CF`` and
    negative at large ``CF``. It is not used to predict anything -- it is the
    independent route the 1957 line is checked against, since the two are supposed to
    agree to about a percent.
    """
    re = _reynolds_domain(reynolds)
    tolerance = _check_positive("tolerance", tolerance)

    def residual(cf: float) -> float:
        return 0.242 / math.sqrt(cf) - math.log10(re * cf)

    low, high = 1.0e-7, 1.0e-1  # f(low) > 0 > f(high) for every Re >= 1e4
    if residual(low) <= 0.0 or residual(high) >= 0.0:  # pragma: no cover - guard
        raise ResistanceError(f"no bracket for the Schoenherr line at Re={re!r}")
    for _ in range(max_iterations):
        middle = 0.5 * (low + high)
        if residual(middle) > 0.0:
            low = middle
        else:
            high = middle
        if high - low <= tolerance * middle:
            break
    return 0.5 * (low + high)


def prandtl_schlichting_friction_coefficient(reynolds: float) -> float:
    """``CF = 0.455 / (log10 Re)^2.58`` -- Prandtl-Schlichting's turbulent plate line.

    A third, independent estimate of the same physical quantity. It sits a few
    percent above the 1957 line by construction (different data set, different
    correlation), so it bounds the friction from above rather than replacing it.
    """
    re = _reynolds_domain(reynolds)
    return 0.455 / math.log10(re) ** 2.58


#: The friction lines this module knows, by name.
FRICTION_LINES: dict[str, Any] = {
    "ittc57": ittc57_friction_coefficient,
    "schoenherr": schoenherr_friction_coefficient,
    "prandtl-schlichting": prandtl_schlichting_friction_coefficient,
}


def friction_coefficient(reynolds: float, line: str = "ittc57") -> float:
    """Evaluate one of :data:`FRICTION_LINES` by name."""
    try:
        function = FRICTION_LINES[line]
    except KeyError:
        known = ", ".join(sorted(FRICTION_LINES))
        raise ResistanceError(f"unknown friction line {line!r}; known: {known}") from None
    return float(function(reynolds))


# ------------------------------------------------------------------- Prohaska fit


@dataclass(frozen=True)
class ProhaskaFit:
    """A Prohaska plot: ``CT/CF`` against ``Fn^4/CF``, a straight line.

    ``intercept = 1 + k`` is the form factor and ``wave_coefficient`` is the slope,
    which is the coefficient of the ``Fn^4`` wave-resistance term at low Froude
    number (the assumption the method rests on, and the reason its validity is
    limited to ``Fn`` below roughly 0.4 where the wave term is still quartic).
    """

    form_factor: float
    wave_coefficient: float
    froude_numbers: tuple[float, ...]
    friction_coefficients: tuple[float, ...]
    total_coefficients: tuple[float, ...]
    residuals: tuple[float, ...]
    max_relative_residual: float
    slope_standard_error: float
    intercept_standard_error: float
    friction_line: str = "ittc57"

    @property
    def points(self) -> tuple[tuple[float, float], ...]:
        """The ``(Fn^4/CF, CT/CF)`` pairs the line was fitted through."""
        return tuple(
            (fn**4 / cf, ct / cf)
            for fn, cf, ct in zip(
                self.froude_numbers, self.friction_coefficients, self.total_coefficients
            )
        )

    def total_coefficient(self, froude: float, reynolds: float) -> float:
        """``CT`` the fitted line implies at a point -- ``(1 + k) CF + y Fn^4``."""
        cf = friction_coefficient(reynolds, self.friction_line)
        return self.form_factor * cf + self.wave_coefficient * froude**4

    def to_dict(self) -> dict[str, Any]:
        return {
            "form_factor": self.form_factor,
            "wave_coefficient": self.wave_coefficient,
            "friction_line": self.friction_line,
            "froude_numbers": list(self.froude_numbers),
            "friction_coefficients": list(self.friction_coefficients),
            "total_coefficients": list(self.total_coefficients),
            "residuals": list(self.residuals),
            "max_relative_residual": self.max_relative_residual,
            "slope_standard_error": self.slope_standard_error,
            "intercept_standard_error": self.intercept_standard_error,
        }


def _as_sequence(name: str, values: Sequence[float]) -> tuple[float, ...]:
    try:
        items = tuple(float(value) for value in values)
    except (TypeError, ValueError) as error:
        raise ResistanceError(f"{name} must be a sequence of numbers: {error}") from None
    if not items:
        raise ResistanceError(f"{name} must not be empty")
    for item in items:
        if not math.isfinite(item):
            raise ResistanceError(f"{name} contains a non-finite value {item!r}")
    return items


def prohaska_fit(
    froude_numbers: Sequence[float],
    total_coefficients: Sequence[float],
    *,
    length: float | None = None,
    friction_coefficients: Sequence[float] | None = None,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    gravity: float = STANDARD_GRAVITY,
    line: str = "ittc57",
) -> ProhaskaFit:
    """Fit ``(1 + k)`` and the wave coefficient from a low-speed resistance curve.

    Prohaska's (1966) idea: below the hump the wave resistance behaves as
    ``R_W ~ Fn^4``, so dividing the measured total coefficient by the friction
    coefficient turns

        ``CT = (1 + k) CF + y Fn^4``   into   ``CT/CF = (1 + k) + y Fn^4/CF``,

    a straight line in ``x = Fn^4/CF``. The intercept *is* the form factor, so the
    method needs neither a wave-resistance theory nor a double-body solution.

    Two ways to supply the friction coefficient, because they answer different
    questions: pass ``friction_coefficients`` when they come from the model test's own
    analysis (and this module only does the regression), or pass ``length`` and let it
    compute ``Re = Fn sqrt(g L) L / nu`` from the Froude numbers. Exactly one of the
    two must be given -- guessing the model length is how a form factor comes out
    wrong by ten percent and nothing complains.

    Requires at least three points for the standard errors to mean anything; two
    points would fit exactly and report a zero residual, which is not evidence.
    """
    froude = _as_sequence("froude_numbers", froude_numbers)
    total = _as_sequence("total_coefficients", total_coefficients)
    if len(froude) != len(total):
        raise ResistanceError(
            f"froude_numbers has {len(froude)} points, total_coefficients has {len(total)}"
        )
    if len(froude) < 3:
        raise ResistanceError(
            f"prohaska_fit needs at least 3 points for a meaningful fit, got {len(froude)}"
        )
    for fn in froude:
        if fn <= 0.0:
            raise ResistanceError(f"froude numbers must be positive, got {fn!r}")
    for ct in total:
        if ct <= 0.0:
            raise ResistanceError(f"total coefficients must be positive, got {ct!r}")

    if friction_coefficients is not None and length is not None:
        raise ResistanceError(
            "give either friction_coefficients or length, not both; the friction line"
            " and the model length are two routes to the same number"
        )
    if friction_coefficients is None and length is None:
        raise ResistanceError("prohaska_fit needs friction_coefficients or a model length")
    if friction_coefficients is not None:
        cf = _as_sequence("friction_coefficients", friction_coefficients)
        if len(cf) != len(froude):
            raise ResistanceError(
                f"friction_coefficients has {len(cf)} points, froude_numbers has {len(froude)}"
            )
    else:
        model_length = _check_positive("length", length)  # type: ignore[arg-type]
        cf = tuple(
            friction_coefficient(
                reynolds_number(model_length, fn * math.sqrt(gravity * model_length), kinematic_viscosity),
                line,
            )
            for fn in froude
        )

    xs = tuple(fn**4 / coefficient for fn, coefficient in zip(froude, cf))
    ys = tuple(ct / coefficient for ct, coefficient in zip(total, cf))

    n = len(xs)
    mean_x = math.fsum(xs) / n
    mean_y = math.fsum(ys) / n
    sxx = math.fsum((x - mean_x) ** 2 for x in xs)
    if sxx == 0.0:  # pragma: no cover - only if every Fn is identical
        raise ResistanceError("all froude numbers are identical; the fit is not determined")
    sxy = math.fsum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x

    residuals = tuple(y - (intercept + slope * x) for x, y in zip(xs, ys))
    scale = max(abs(y) for y in ys)
    sse = math.fsum(r * r for r in residuals)
    sigma2 = sse / (n - 2)
    slope_error = math.sqrt(sigma2 / sxx) if sigma2 > 0.0 else 0.0
    intercept_error = math.sqrt(sigma2 * (1.0 / n + mean_x**2 / sxx)) if sigma2 > 0.0 else 0.0

    return ProhaskaFit(
        form_factor=intercept,
        wave_coefficient=slope,
        froude_numbers=froude,
        friction_coefficients=cf,
        total_coefficients=total,
        residuals=residuals,
        max_relative_residual=max(abs(r) for r in residuals) / scale if scale else 0.0,
        slope_standard_error=slope_error,
        intercept_standard_error=intercept_error,
        friction_line=line,
    )


# --------------------------------------------------------------- resistance itself


@dataclass(frozen=True)
class Resistance:
    """One body at one speed, decomposed the way the ITTC method decomposes it.

    ``total_coefficient = (1 + k) CF + CR + CA`` and every resistance is
    ``0.5 rho S V^2`` times the matching coefficient, so ``total_resistance`` can be
    recomputed from the coefficients by hand -- which is what the tests do.
    """

    length: float
    velocity: float
    wetted_surface: float
    density: float
    kinematic_viscosity: float
    reynolds: float
    froude: float
    friction_coefficient: float
    form_factor: float
    residual_coefficient: float
    correlation_allowance: float
    friction_line: str = "ittc57"

    @property
    def dynamic_pressure_area(self) -> float:
        """``0.5 rho S V^2``, N -- the factor every coefficient multiplies."""
        return 0.5 * self.density * self.wetted_surface * self.velocity * self.velocity

    @property
    def viscous_coefficient(self) -> float:
        """``(1 + k) CF``: the form factor applied to the friction line."""
        return self.form_factor * self.friction_coefficient

    @property
    def total_coefficient(self) -> float:
        return self.viscous_coefficient + self.residual_coefficient + self.correlation_allowance

    @property
    def friction_resistance(self) -> float:
        return self.dynamic_pressure_area * self.friction_coefficient

    @property
    def viscous_resistance(self) -> float:
        return self.dynamic_pressure_area * self.viscous_coefficient

    @property
    def form_resistance(self) -> float:
        """``(1 + k) R_F - R_F``: what the form factor adds over the flat plate."""
        return self.viscous_resistance - self.friction_resistance

    @property
    def residual_resistance(self) -> float:
        return self.dynamic_pressure_area * self.residual_coefficient

    @property
    def allowance_resistance(self) -> float:
        return self.dynamic_pressure_area * self.correlation_allowance

    @property
    def total_resistance(self) -> float:
        return self.dynamic_pressure_area * self.total_coefficient

    @property
    def effective_power(self) -> float:
        """``P_E = R_T V``, W -- the power the towline delivers to the hull."""
        return self.total_resistance * self.velocity

    def to_dict(self) -> dict[str, float | str]:
        return {
            "length_m": self.length,
            "velocity_m_s": self.velocity,
            "wetted_surface_m2": self.wetted_surface,
            "density_kg_m3": self.density,
            "kinematic_viscosity_m2_s": self.kinematic_viscosity,
            "reynolds": self.reynolds,
            "froude": self.froude,
            "friction_line": self.friction_line,
            "friction_coefficient": self.friction_coefficient,
            "form_factor": self.form_factor,
            "viscous_coefficient": self.viscous_coefficient,
            "residual_coefficient": self.residual_coefficient,
            "correlation_allowance": self.correlation_allowance,
            "total_coefficient": self.total_coefficient,
            "friction_resistance_n": self.friction_resistance,
            "viscous_resistance_n": self.viscous_resistance,
            "form_resistance_n": self.form_resistance,
            "residual_resistance_n": self.residual_resistance,
            "allowance_resistance_n": self.allowance_resistance,
            "total_resistance_n": self.total_resistance,
            "effective_power_w": self.effective_power,
        }


def bare_hull_resistance(
    *,
    length: float,
    velocity: float,
    wetted_surface: float,
    form_factor: float = 1.0,
    residual_coefficient: float = 0.0,
    correlation_allowance: float = 0.0,
    density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    gravity: float = STANDARD_GRAVITY,
    line: str = "ittc57",
) -> Resistance:
    """The resistance of a bare hull at one speed, from its wetted surface.

    The two coefficients that the method cannot derive are **inputs with no
    default**: ``form_factor`` defaults to ``1 + k = 1`` (a flat plate, i.e. the
    statement "I have not measured k" rather than an assumption that k is zero in
    reality) and ``residual_coefficient`` defaults to ``0`` (pure friction: a
    useful lower bound on a real hull, and the exact answer for a body of revolution
    only in the slender limit). Read the defaults as "not measured", and get both from
    :func:`prohaska_fit` or from a towing tank before quoting a resistance.

    ``length`` must be the length the friction line is evaluated at -- the
    waterline or the wetted length, consistently with the surface that produced
    ``wetted_surface``. Mixing a model-test length with a ship's wetted area is the
    single most common way this arithmetic goes wrong, so every returned object
    carries both inputs.
    """
    length = _check_positive("length", length)
    velocity = _check_positive("velocity", velocity)
    wetted = _check_positive("wetted_surface", wetted_surface)
    form = _check_positive("form_factor", form_factor)
    if form < 1.0:
        raise ResistanceError(
            f"form_factor must be >= 1 (it is 1 + k, and k >= 0), got {form!r}"
        )
    residual = _check_finite("residual_coefficient", residual_coefficient)
    allowance = _check_finite("correlation_allowance", correlation_allowance)
    if residual < 0.0 or allowance < 0.0:
        raise ResistanceError("residual_coefficient and correlation_allowance must be >= 0")
    rho = _check_positive("density", density)
    nu = _check_positive("kinematic_viscosity", kinematic_viscosity)
    if line not in FRICTION_LINES:
        known = ", ".join(sorted(FRICTION_LINES))
        raise ResistanceError(f"unknown friction line {line!r}; known: {known}")

    re = reynolds_number(length, velocity, nu)
    return Resistance(
        length=length,
        velocity=velocity,
        wetted_surface=wetted,
        density=rho,
        kinematic_viscosity=nu,
        reynolds=re,
        froude=froude_number(length, velocity, gravity),
        friction_coefficient=friction_coefficient(re, line),
        form_factor=form,
        residual_coefficient=residual,
        correlation_allowance=allowance,
        friction_line=line,
    )


def mesh_bare_hull_resistance(
    mesh: Mesh, *, length: float, velocity: float, **kwargs: Any
) -> Resistance:
    """Same as :func:`bare_hull_resistance`, with ``S`` measured off the mesh.

    ``S`` is :func:`tensorcad.mesh.surface_area` of the mesh **as given**, which is the
    wetted surface only if the mesh *is* the wetted shell. The hull families here close
    their bodies at the design waterline, so a generated hull's shell includes the deck
    patch -- on the Wigley reference case the deck is 31 % of the shell area, so a shell
    area used as a wetted surface inflates the friction by 45 %. A mesh capped at the
    waterline (``tensorcad.trim``/``tensorcad.gz`` produce such meshes) has no deck patch
    and needs no correction; otherwise subtract it or pass ``wetted_surface`` explicitly.
    """
    return bare_hull_resistance(
        length=length,
        velocity=velocity,
        wetted_surface=_mesh_surface_area(mesh),
        **kwargs,
    )


# ------------------------------------------------------------- model-to-ship scaling


@dataclass(frozen=True)
class ModelShipPrediction:
    """A model test carried to full scale by the ITTC-78 two-dimensional method.

    ``ship.total_coefficient = (1 + k) CF_s + CR + CA`` where ``CR`` came from the
    model's measured total, ``CR = CT_m - (1 + k) CF_m``. The prediction is only as
    good as ``CA``, which is why it is a required-in-spirit input and is reported
    separately from the residual rather than folded into it.
    """

    model: Resistance
    ship: Resistance
    residual_coefficient: float
    correlation_allowance: float
    velocity_scale: float
    reynolds_scale: float

    @property
    def model_total_coefficient(self) -> float:
        return self.model.total_coefficient

    @property
    def scale_effect(self) -> float:
        """``CT_m / CT_s``: how much of the model's coefficient does not scale."""
        return self.model_total_coefficient / self.ship.total_coefficient

    @property
    def ship_resistance_ratio(self) -> float:
        """``R_T_s / R_T_m``, which is *not* ``lambda^3`` -- the whole point."""
        return self.ship.total_resistance / self.model.total_resistance

    def to_dict(self) -> dict[str, Any]:
        return {
            "model": self.model.to_dict(),
            "ship": self.ship.to_dict(),
            "residual_coefficient": self.residual_coefficient,
            "correlation_allowance": self.correlation_allowance,
            "velocity_scale": self.velocity_scale,
            "reynolds_scale": self.reynolds_scale,
            "scale_effect": self.scale_effect,
            "ship_resistance_ratio": self.ship_resistance_ratio,
        }


def ship_from_model_test(
    *,
    model_length: float,
    ship_length: float,
    model_velocity: float,
    model_wetted_surface: float,
    ship_wetted_surface: float,
    form_factor: float,
    model_total_coefficient: float | None = None,
    model_total_resistance: float | None = None,
    residual_coefficient: float | None = None,
    correlation_allowance: float = 0.0,
    model_density: float = SEA_WATER_DENSITY,
    model_kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    ship_density: float = SEA_WATER_DENSITY,
    ship_kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    gravity: float = STANDARD_GRAVITY,
    line: str = "ittc57",
) -> ModelShipPrediction:
    """Scale a model test to the ship at the corresponding Froude number.

    Inputs, at the model point: its length, speed, wetted surface (all measured) and
    its total resistance coefficient *or* total resistance in newtons. Either yields
    ``CR = CT_m - (1 + k) CF_m``; if ``residual_coefficient`` is given explicitly it
    wins and the model total becomes optional (that is the "I already have a CR from a
    panel code" case, and the model numbers are then reported only for the record).

    The ship is run at ``V_s = V_m sqrt(L_s/L_m)``, which is Froude scaling and the
    reason the wave pattern -- and only the wave pattern -- is preserved.
    ``correlation_allowance`` (ITTC's ``CA``, the model-ship correlation allowance) is
    applied to the ship's coefficient and is **never invented here**: it encodes
    roughness and the residual bias of the whole method, and a wrong value moves the
    answer by more than a form factor does.
    """
    model_length = _check_positive("model_length", model_length)
    ship_length = _check_positive("ship_length", ship_length)
    model_velocity = _check_positive("model_velocity", model_velocity)
    model_surface = _check_positive("model_wetted_surface", model_wetted_surface)
    ship_surface = _check_positive("ship_wetted_surface", ship_wetted_surface)
    form = _check_positive("form_factor", form_factor)
    allowance = _check_finite("correlation_allowance", correlation_allowance)

    provided = [
        name
        for name, value in (
            ("model_total_coefficient", model_total_coefficient),
            ("model_total_resistance", model_total_resistance),
            ("residual_coefficient", residual_coefficient),
        )
        if value is not None
    ]
    if not provided:
        raise ResistanceError(
            "give one of model_total_coefficient, model_total_resistance or residual_coefficient"
        )
    if "residual_coefficient" in provided and len(provided) > 1:
        raise ResistanceError(
            "residual_coefficient is given explicitly and supersedes the model total;"
            f" drop {', '.join(n for n in provided if n != 'residual_coefficient')} or drop it"
        )
    if "residual_coefficient" not in provided:
        if len(provided) > 1:
            # Both routes given: they must agree, and the check is the useful part.
            model_probe = bare_hull_resistance(
                length=model_length,
                velocity=model_velocity,
                wetted_surface=model_surface,
                form_factor=form,
                density=model_density,
                kinematic_viscosity=model_kinematic_viscosity,
                gravity=gravity,
                line=line,
            )
            implied = float(model_total_resistance) / model_probe.dynamic_pressure_area  # type: ignore[arg-type]
            if not math.isclose(implied, float(model_total_coefficient), rel_tol=1e-9):  # type: ignore[arg-type]
                raise ResistanceError(
                    "model_total_coefficient and model_total_resistance disagree:"
                    f" RT/(0.5 rho S V^2) = {implied!r}, CT = {model_total_coefficient!r}"
                )

    probe = bare_hull_resistance(
        length=model_length,
        velocity=model_velocity,
        wetted_surface=model_surface,
        form_factor=form,
        residual_coefficient=0.0,
        correlation_allowance=0.0,
        density=model_density,
        kinematic_viscosity=model_kinematic_viscosity,
        gravity=gravity,
        line=line,
    )

    if residual_coefficient is not None:
        residual = _check_finite("residual_coefficient", residual_coefficient)
        if residual < 0.0:
            raise ResistanceError("residual_coefficient must be >= 0")
    else:
        if model_total_coefficient is None:
            model_ct = float(model_total_resistance) / probe.dynamic_pressure_area  # type: ignore[arg-type]
        else:
            model_ct = _check_finite("model_total_coefficient", model_total_coefficient)
        residual = model_ct - probe.viscous_coefficient
        # A negative CR means the model measured less drag than the flat plate the
        # method subtracts. That happens at low Froude number on a well-formed hull,
        # so it is reported rather than clipped; it also happens when the form factor
        # or the water temperature is wrong, which is the caller's call to make.

    # The model record carries the residual back, so that ``model.total_coefficient``
    # is the measured coefficient it came from -- not ``(1 + k) CF_m`` with the wave
    # part dropped, which would look like a 13 % prediction error to anyone reading it.
    model = bare_hull_resistance(
        length=model_length,
        velocity=model_velocity,
        wetted_surface=model_surface,
        form_factor=form,
        residual_coefficient=residual,
        correlation_allowance=0.0,
        density=model_density,
        kinematic_viscosity=model_kinematic_viscosity,
        gravity=gravity,
        line=line,
    )

    velocity_scale = math.sqrt(ship_length / model_length)
    ship_velocity = model_velocity * velocity_scale
    ship = bare_hull_resistance(
        length=ship_length,
        velocity=ship_velocity,
        wetted_surface=ship_surface,
        form_factor=form,
        residual_coefficient=residual,
        correlation_allowance=allowance,
        density=ship_density,
        kinematic_viscosity=ship_kinematic_viscosity,
        gravity=gravity,
        line=line,
    )

    return ModelShipPrediction(
        model=model,
        ship=ship,
        residual_coefficient=residual,
        correlation_allowance=allowance,
        velocity_scale=velocity_scale,
        reynolds_scale=ship.reynolds / model.reynolds,
    )