"""Savitsky 1964 planing prediction: lift, drag and running trim for planing craft.

This module implements the classical Savitsky empirical method for predicting the
hydrodynamic forces and running attitude of planing hulls with constant deadrise.
The method is written for the prismatic planing hull family that
:mod:`tensorcad.planing` generates, and consumes its geometric parameters directly.

Reference
---------
Savitsky, D. (1964). "Hydrodynamic design of planing hulls." *Marine Technology*,
1(1), 71-95. The equations here are transcribed from that paper and from the
subsequent summary in Blount & Codega (1992) and Clement & Blount (1963).

What is here
------------
* :func:`running_trim` -- the equilibrium trim angle τ at a given speed and weight,
  solved from the lift equation CL = Δ / (0.5 ρ V² B²) where CL is a function of τ,
  the wetted length-beam ratio λ, and the deadrise β.
* :func:`lift_coefficient` -- the Savitsky CL(τ, λ, β) correlation.
* :func:`resistance_components` -- friction drag (ITTC-57 line on the wetted area),
  pressure drag (Savitsky correlation), and total resistance.
* :func:`planing_prediction` -- the complete prediction at one speed: trim, resistance,
  wetted length, and the individual drag components.
* :func:`speed_sweep` -- a speed range scan that finds the planing onset and tracks
  the equilibrium through the speed range.

What this is not
----------------
* A seakeeping or dynamic stability analysis (no heave/pitch RAOs, no porpoising).
* A method for non-prismatic hulls (variable deadrise, warped hulls, or hulls with
  spray rails/steps need a different approach).
* A validation against model test data beyond the Savitsky original correlation band
  (±10% on resistance is the published accuracy; this module reproduces the correlation,
  not the physics).

Units are SI throughout: metres, seconds, newtons, watts, kg/m^3. Angles are degrees
in the API but converted internally to radians for the trig.

The method is empirical and has a validity box:
* Speed range: Fr = V / sqrt(g L) typically 1.0 to 4.0 (planing regime)
* Deadrise: 10° to 30° (the correlation was fitted in this range)
* Length-beam ratio: L/B = 3 to 6 (the test matrix)
* Weight displacement: the hull must be planing, not displacement mode

Outside this box the predictions are extrapolations and the errors grow.

Audit status (m16954)
---------------------
An adversarial audit of this module (``docs/savitsky_audit_m16954.md``) found that
the original single-equation trim solve could not produce a wetted length ratio at
all: ``lambda`` was pinned at its 3.0 initial guess in every case, so the "equilibrium
wetted length" reported by :func:`planing_prediction` was an assumption, not a result.
The audit also found that the module's :func:`lift_coefficient` correlation carries no
speed coefficient ``Cv`` while Savitsky's lift equation does, so it cannot be combined
with Savitsky's centre-of-pressure condition to close the two-equation system.

Both are now explicit:

* :func:`trim_for_wetted_length` -- exact closed-form trim for an *assumed* ``λ``
  (what the old ``running_trim`` actually computed; kept for backwards compatibility).
* :func:`equilibrium_trim` -- the classical two-equation Savitsky closure (lift plus
  centre-of-pressure moment about the LCG) using :func:`classic_lift_coefficient`,
  which does carry ``Cv``. This is the only path in the module that *solves* ``λ``.
* :func:`planing_prediction` reports which of the two paths produced its numbers via
  ``solution_method`` and ``wetted_length_ratio_solved``, and never silently clamps.

The ``Cv``-dependent coefficients in :func:`classic_lift_coefficient` and
:func:`center_of_pressure_ratio` are the classical statement of Savitsky (1964)
eqs. (1)-(3); they could not be checked against the primary paper in the m16954
environment (no offline copy), so they carry evidence level WITHHELD and the two
correlations are cross-tabulated rather than silently merged. See the audit document.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

__all__ = [
    "PlaningPrediction",
    "PlaningPredictionError",
    "TrimSolution",
    "center_of_pressure_ratio",
    "classic_lift_coefficient",
    "equilibrium_trim",
    "lift_coefficient",
    "planing_prediction",
    "resistance_components",
    "running_trim",
    "speed_coefficient",
    "speed_sweep",
    "trim_for_wetted_length",
    "validity_warnings",
]

#: Standard gravity, m/s^2.
STANDARD_GRAVITY = 9.80665

#: Sea water density at 15 deg C, kg/m^3 (the ITTC standard).
SEA_WATER_DENSITY = 1025.9

#: Fresh water density at 15 deg C, kg/m^3.
FRESH_WATER_DENSITY = 999.1

#: Sea water kinematic viscosity at 15 deg C, m^2/s.
SEA_WATER_KINEMATIC_VISCOSITY = 1.1892e-6


class PlaningPredictionError(ValueError):
    """Raised when inputs are outside the Savitsky validity box or the solver fails."""


@dataclass
class PlaningPrediction:
    """Result of a Savitsky planing prediction at one speed.

    Attributes
    ----------
    speed_ms : float
        Forward speed, m/s.
    froude_number : float
        Froude number V / sqrt(g L).
    trim_deg : float
        Running trim angle τ, degrees.
    wetted_length_m : float
        Mean wetted length along the keel, m.
    wetted_length_ratio : float
        Wetted length / beam, λ = L_w / B (dimensionless).
    lift_coefficient : float
        Savitsky CL (dimensionless).
    friction_drag_n : float
        Friction drag (ITTC-57 line on wetted area), N.
    pressure_drag_n : float
        Pressure drag (Savitsky correlation), N.
    total_drag_n : float
        Total resistance R_T = R_F + R_P, N.
    effective_power_w : float
        Effective power P_E = R_T × V, W.
    weight_n : float
        Hull weight (displacement), N.
    reynolds_number : float
        Reynolds number V L_w / ν (based on wetted length).
    wetted_length_ratio_solved : bool
        ``True`` only when ``λ`` came out of :func:`equilibrium_trim` (the two-equation
        closure). ``False`` means λ was an input assumption or a hull-length clamp.
    solution_method : str
        ``'classic_equilibrium'`` or ``'assumed_wetted_length'``.
    warnings : tuple of str
        Honesty channel: every clamp, extrapolation or assumption applied to this
        result. Empty only when λ was solved and the case is inside the validity box.
    """

    speed_ms: float
    froude_number: float
    trim_deg: float
    wetted_length_m: float
    wetted_length_ratio: float
    lift_coefficient: float
    friction_drag_n: float
    pressure_drag_n: float
    total_drag_n: float
    effective_power_w: float
    weight_n: float
    reynolds_number: float
    wetted_length_ratio_solved: bool = False
    solution_method: str = "assumed_wetted_length"
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable dictionary."""
        return {
            "speed_ms": self.speed_ms,
            "froude_number": self.froude_number,
            "trim_deg": self.trim_deg,
            "wetted_length_m": self.wetted_length_m,
            "wetted_length_ratio": self.wetted_length_ratio,
            "lift_coefficient": self.lift_coefficient,
            "friction_drag_n": self.friction_drag_n,
            "pressure_drag_n": self.pressure_drag_n,
            "total_drag_n": self.total_drag_n,
            "effective_power_w": self.effective_power_w,
            "weight_n": self.weight_n,
            "reynolds_number": self.reynolds_number,
            "wetted_length_ratio_solved": self.wetted_length_ratio_solved,
            "solution_method": self.solution_method,
            "warnings": list(self.warnings),
        }


def lift_coefficient(trim_deg: float, wetted_length_ratio: float, deadrise_deg: float) -> float:
    """Savitsky lift coefficient CL(τ, λ, β).

    The correlation from Savitsky (1964) equation (1):

        CL = 0.00222 τ^0.9 λ^0.67 / (β/1000)^0.33

    where τ is in degrees, λ = L_w / B, and β is in degrees.

    Parameters
    ----------
    trim_deg : float
        Running trim angle τ, degrees (typically 2° to 8°).
    wetted_length_ratio : float
        Wetted length / beam, λ (typically 2 to 6).
    deadrise_deg : float
        Deadrise angle β, degrees (typically 10° to 30°).

    Returns
    -------
    float
        Lift coefficient CL (dimensionless).

    Raises
    ------
    PlaningPredictionError
        If inputs are outside the Savitsky validity box.
    """
    if trim_deg <= 0.0 or trim_deg > 15.0:
        raise PlaningPredictionError(
            f"trim_deg must be in (0, 15], got {trim_deg}"
        )
    if wetted_length_ratio <= 0.0 or wetted_length_ratio > 10.0:
        raise PlaningPredictionError(
            f"wetted_length_ratio must be in (0, 10], got {wetted_length_ratio}"
        )
    if deadrise_deg < 0.0 or deadrise_deg > 45.0:
        raise PlaningPredictionError(
            f"deadrise_deg must be in [0, 45], got {deadrise_deg}"
        )

    # Savitsky (1964) equation (1): CL = 0.00222 τ^0.9 λ^0.67 / (β/1000)^0.33
    # Note: β is in degrees in the denominator, so β/1000 is a small number.
    # The original paper writes this as (β/1000)^0.33 to keep the coefficient 0.00222.
    tau = trim_deg
    lam = wetted_length_ratio
    beta = deadrise_deg

    cl = 0.00222 * (tau ** 0.9) * (lam ** 0.67) / ((beta / 1000.0) ** 0.33)
    return cl



# ---------------------------------------------------------------------------
# Constants of the lift correlation and the trim window
# ---------------------------------------------------------------------------

#: Coefficient of the Savitsky lift correlation CL = k τ^a λ^b (β/1000)^-c.
_LIFT_K = 0.00222
_LIFT_TRIM_EXPONENT = 0.9
_LIFT_LAMBDA_EXPONENT = 0.67
_LIFT_DEADRISE_EXPONENT = 0.33

#: Trim window the lift correlation is used in, degrees. Documented in the module
#: docstring as the validity box; enforced rather than clamped (m16954).
TRIM_MIN_DEG = 0.1
TRIM_MAX_DEG = 15.0

#: Wetted length ratio assumed when the caller does not supply one.
DEFAULT_WETTED_LENGTH_RATIO = 3.0


def validity_warnings(
    froude_number: float,
    deadrise_deg: float,
    length_m: float | None = None,
    beam_m: float | None = None,
    trim_deg: float | None = None,
) -> tuple[str, ...]:
    """Return one message per documented validity-box violation (empty if inside)."""
    messages: list[str] = []
    if trim_deg is not None and (trim_deg < 2.0 or trim_deg > 15.0):
        messages.append(
            f"trim={trim_deg:.2f} deg is outside the fitted band 2-15 deg"
        )
    if froude_number < 1.0 or froude_number > 4.0:
        messages.append(
            f"Fr={froude_number:.3f} is outside the planing box 1.0-4.0: the prediction "
            "is an extrapolation of the correlation"
        )
    if deadrise_deg < 10.0 or deadrise_deg > 30.0:
        messages.append(
            f"deadrise={deadrise_deg:.1f} deg is outside the fitted band 10-30 deg"
        )
    if length_m is not None and beam_m is not None:
        lb = length_m / beam_m
        if lb < 3.0 or lb > 6.0:
            messages.append(f"L/B={lb:.2f} is outside the test matrix 3-6")
    return tuple(messages)


def trim_for_wetted_length(
    weight_n: float,
    speed_ms: float,
    beam_m: float,
    deadrise_deg: float,
    wetted_length_ratio: float,
    *,
    water_density: float = SEA_WATER_DENSITY,
) -> float:
    """Equilibrium running trim for an *assumed* wetted length ratio λ.

    The lift equation CL(τ, λ, β) × 0.5 ρ V² B² = Δ is linear in no unknown when λ
    is given, so it inverts exactly:

        τ = ( CL_target (β/1000)^c / (k λ^b) )^(1/a),   CL_target = Δ / (0.5 ρ V² B²)

    Parameters
    ----------
    weight_n, speed_ms, beam_m, deadrise_deg : float
        Hull weight (N), speed (m/s), beam (m) and deadrise (degrees).
    wetted_length_ratio : float
        λ = L_w / B, an input here -- this function does not solve for it.
    water_density : float, optional
        Water density, kg/m^3.

    Returns
    -------
    float
        Trim angle τ in degrees.

    Raises
    ------
    PlaningPredictionError
        If an input is non-positive or the required trim lies outside the
        implemented window ``[TRIM_MIN_DEG, TRIM_MAX_DEG]``. The message reports the
        trim the lift equation actually asks for, instead of a spurious
        "failed to converge" (m16954).
    """
    if weight_n <= 0.0:
        raise PlaningPredictionError(f"weight_n must be positive, got {weight_n}")
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")
    if wetted_length_ratio <= 0.0:
        raise PlaningPredictionError(
            f"wetted_length_ratio must be positive, got {wetted_length_ratio}"
        )

    cl_target = weight_n / (0.5 * water_density * speed_ms ** 2 * beam_m ** 2)
    inner = (
        cl_target
        * (deadrise_deg / 1000.0) ** _LIFT_DEADRISE_EXPONENT
        / (_LIFT_K * wetted_length_ratio ** _LIFT_LAMBDA_EXPONENT)
    )
    trim_deg = inner ** (1.0 / _LIFT_TRIM_EXPONENT)

    if trim_deg < TRIM_MIN_DEG or trim_deg > TRIM_MAX_DEG:
        raise PlaningPredictionError(
            f"required trim {trim_deg:.2f} deg is outside the implemented window "
            f"{TRIM_MIN_DEG}-{TRIM_MAX_DEG} deg: the craft is not planing at this "
            "weight/speed/beam combination"
        )
    return trim_deg


def speed_coefficient(speed_ms: float, beam_m: float) -> float:
    """Savitsky speed coefficient Cv = V / sqrt(g B)."""
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")
    return speed_ms / math.sqrt(STANDARD_GRAVITY * beam_m)


def classic_lift_coefficient(
    trim_deg: float,
    wetted_length_ratio: float,
    cv: float,
    deadrise_deg: float = 0.0,
) -> float:
    """Savitsky lift coefficient in its classical, speed-dependent statement.

        CL_beta = τ^1.1 ( 0.0120 λ^0.5 + 0.0055 λ^2.5 / Cv² )
        CL_0    = CL_beta - 0.0065 β CL_beta^0.6        (β in degrees, 0 when β=0)

    Written with τ in degrees, as in the original correlation. Unlike
    :func:`lift_coefficient` this form depends on the speed coefficient ``Cv``, which
    is what makes the two-equation closure of :func:`equilibrium_trim` possible.

    Evidence level WITHHELD: transcribed from the standard statement of Savitsky
    (1964) eqs. (1)-(2); not checked against the primary paper in this environment.
    """
    if trim_deg <= 0.0:
        raise PlaningPredictionError(f"trim_deg must be positive, got {trim_deg}")
    if wetted_length_ratio <= 0.0:
        raise PlaningPredictionError(
            f"wetted_length_ratio must be positive, got {wetted_length_ratio}"
        )
    if cv <= 0.0:
        raise PlaningPredictionError(f"cv must be positive, got {cv}")
    if deadrise_deg < 0.0:
        raise PlaningPredictionError(f"deadrise_deg must be non-negative, got {deadrise_deg}")

    cl_beta = trim_deg ** 1.1 * (
        0.0120 * wetted_length_ratio ** 0.5
        + 0.0055 * wetted_length_ratio ** 2.5 / cv ** 2
    )
    if deadrise_deg == 0.0:
        return cl_beta
    return cl_beta - 0.0065 * deadrise_deg * cl_beta ** 0.6


def center_of_pressure_ratio(wetted_length_ratio: float, cv: float) -> float:
    """Centre of pressure as a fraction of the wetted length, from the transom.

        l_cp / (λ B) = 0.75 - 1 / (5.21 Cv² / λ² + 2.39)

    Evidence level WITHHELD: standard statement of the Savitsky (1964) centre-of-
    pressure correlation; not checked against the primary paper in this environment.
    """
    if wetted_length_ratio <= 0.0:
        raise PlaningPredictionError(
            f"wetted_length_ratio must be positive, got {wetted_length_ratio}"
        )
    if cv <= 0.0:
        raise PlaningPredictionError(f"cv must be positive, got {cv}")
    return 0.75 - 1.0 / (5.21 * cv ** 2 / wetted_length_ratio ** 2 + 2.39)


@dataclass
class TrimSolution:
    """Equilibrium solution of the classical two-equation Savitsky closure.

    Attributes
    ----------
    trim_deg : float
        Running trim angle τ, degrees.
    wetted_length_ratio : float
        λ = L_w / B, *solved*, not assumed.
    speed_coefficient : float
        Cv = V / sqrt(g B).
    lift_residual : float
        |CL_0 q B² - Δ| / Δ, dimensionless.
    moment_residual_m : float
        |l_cp - lcg|, m.
    iterations : int
        Bisection iterations spent on the outer (trim) loop.
    method : str
        Provenance tag.
    """

    trim_deg: float
    wetted_length_ratio: float
    speed_coefficient: float
    lift_residual: float
    moment_residual_m: float
    iterations: int
    method: str = "classic_savitsky_lift_and_moment"

    def as_dict(self) -> dict[str, Any]:
        return {
            "trim_deg": self.trim_deg,
            "wetted_length_ratio": self.wetted_length_ratio,
            "speed_coefficient": self.speed_coefficient,
            "lift_residual": self.lift_residual,
            "moment_residual_m": self.moment_residual_m,
            "iterations": self.iterations,
            "method": self.method,
        }


def equilibrium_trim(
    weight_n: float,
    speed_ms: float,
    beam_m: float,
    deadrise_deg: float,
    lcg_m: float,
    *,
    water_density: float = SEA_WATER_DENSITY,
    max_wetted_length_m: float | None = None,
    tolerance: float = 1e-10,
    max_iterations: int = 200,
) -> TrimSolution:
    """Solve the classical two-equation Savitsky system for (τ, λ).

    Equations
    ---------
    1. Lift balances weight:      CL_0(τ, λ, Cv, β) 0.5 ρ V² B² = Δ
    2. Moment about the LCG:      l_cp(λ, Cv) = lcg

    Both lengths are measured along the keel from the transom, which is the standard
    statement of the method: the centre of pressure of the planing force must sit
    under the centre of gravity for the craft to hold a steady trim. Equation 1 gives
    λ(τ) (CL_0 is monotone increasing in λ), equation 2 then fixes τ by bisection
    (l_cp decreases monotonically with τ because more trim needs less wetted length).

    Parameters
    ----------
    weight_n, speed_ms, beam_m, deadrise_deg : float
        Hull weight (N), speed (m/s), beam (m), deadrise (degrees).
    lcg_m : float
        Longitudinal centre of gravity, m forward of the transom, along the keel.
    max_wetted_length_m : float, optional
        Cap on λ B (typically the hull length). Raises if the solution needs more.
    water_density : float, optional
        Water density, kg/m^3.
    tolerance : float
        Residual tolerance asked of the bisection.
    max_iterations : int
        Maximum outer bisection iterations.

    Returns
    -------
    TrimSolution

    Raises
    ------
    PlaningPredictionError
        If an input is invalid, or if no solution exists inside the trim window
        (message reports the range of centre-of-pressure positions that *is*
        reachable, instead of silently clamping).
    """
    if weight_n <= 0.0:
        raise PlaningPredictionError(f"weight_n must be positive, got {weight_n}")
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")
    if lcg_m <= 0.0:
        raise PlaningPredictionError(f"lcg_m must be positive, got {lcg_m}")

    cv = speed_coefficient(speed_ms, beam_m)
    cl_target = weight_n / (0.5 * water_density * speed_ms ** 2 * beam_m ** 2)

    lam_max = 60.0 if max_wetted_length_m is None else max_wetted_length_m / beam_m
    if lam_max <= 0.0:
        raise PlaningPredictionError(
            f"max_wetted_length_m={max_wetted_length_m} leaves no room for a wetted length"
        )

    def lift(trim_deg: float, lam: float) -> float:
        return classic_lift_coefficient(trim_deg, lam, cv, deadrise_deg)

    def lambda_for_lift(trim_deg: float) -> float | None:
        """Bisect the monotone CL_0(λ) for the λ that balances weight, None if capped."""
        if lift(trim_deg, 1e-6) > cl_target:
            return None
        if lift(trim_deg, lam_max) < cl_target:
            return None
        lo, hi = 1e-6, lam_max
        for _ in range(400):
            mid = 0.5 * (lo + hi)
            if lift(trim_deg, mid) < cl_target:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-14 * max(1.0, hi):
                break
        return 0.5 * (lo + hi)

    def cp_position(trim_deg: float) -> float | None:
        lam = lambda_for_lift(trim_deg)
        if lam is None:
            return None
        return center_of_pressure_ratio(lam, cv) * lam * beam_m

    hi_trim = TRIM_MAX_DEG
    if lift(hi_trim, lam_max) < cl_target:
        raise PlaningPredictionError(
            f"no planing equilibrium: at the maximum trim {TRIM_MAX_DEG:.1f} deg and the "
            f"maximum wetted length ratio {lam_max:.2f} the lift coefficient reaches only "
            f"{lift(hi_trim, lam_max):.4f} against the {cl_target:.4f} needed "
            f"(Cv={cv:.3f}): overloaded or too slow to plane"
        )
    # Lowest trim that can still reach the target lift at the longest wetted length.
    # Below it the λ(τ) branch of the lift equation simply does not exist (the
    # deadrise correction can even drive CL_0 negative), so it is not a solution
    # window and must not be bisected over.
    lo_trim = TRIM_MIN_DEG
    if lift(lo_trim, lam_max) < cl_target:
        scan_lo, scan_hi = TRIM_MIN_DEG, TRIM_MAX_DEG
        for _ in range(200):
            mid = 0.5 * (scan_lo + scan_hi)
            if lift(mid, lam_max) < cl_target:
                scan_lo = mid
            else:
                scan_hi = mid
            if scan_hi - scan_lo <= 1e-13 * max(1.0, scan_hi):
                break
        lo_trim = scan_hi

    cp_lo, cp_hi = cp_position(lo_trim), cp_position(hi_trim)
    assert cp_lo is not None and cp_hi is not None  # both reachable by construction
    if cp_lo < lcg_m or cp_hi > lcg_m:
        raise PlaningPredictionError(
            f"no planing equilibrium for lcg={lcg_m:.3f} m: the centre of pressure can "
            f"only be placed between {cp_hi:.3f} m and {cp_lo:.3f} m from the transom "
            f"over the trim window {TRIM_MIN_DEG}-{TRIM_MAX_DEG} deg "
            f"(Cv={cv:.3f}); move the LCG or change the speed"
        )

    # l_cp decreases with trim: bisect for l_cp(τ) = lcg.
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        mid_trim = 0.5 * (lo_trim + hi_trim)
        cp_mid = cp_position(mid_trim)
        assert cp_mid is not None
        if cp_mid > lcg_m:
            lo_trim = mid_trim
        else:
            hi_trim = mid_trim
        if hi_trim - lo_trim <= tolerance * max(1.0, mid_trim):
            break

    trim_deg = 0.5 * (lo_trim + hi_trim)
    lam = lambda_for_lift(trim_deg)
    assert lam is not None
    cl0 = classic_lift_coefficient(trim_deg, lam, cv, deadrise_deg)
    lift_residual = abs(
        cl0 * 0.5 * water_density * speed_ms ** 2 * beam_m ** 2 - weight_n
    ) / weight_n
    cp = center_of_pressure_ratio(lam, cv) * lam * beam_m

    return TrimSolution(
        trim_deg=trim_deg,
        wetted_length_ratio=lam,
        speed_coefficient=cv,
        lift_residual=lift_residual,
        moment_residual_m=abs(cp - lcg_m),
        iterations=iterations,
    )


def running_trim(
    weight_n: float,
    speed_ms: float,
    beam_m: float,
    deadrise_deg: float,
    *,
    wetted_length_ratio: float = DEFAULT_WETTED_LENGTH_RATIO,
    water_density: float = SEA_WATER_DENSITY,
    max_iterations: int = 50,
    tolerance: float = 1e-6,
) -> tuple[float, float]:
    """Equilibrium running trim for an *assumed* wetted length ratio; returns (τ, λ).

    Signature note (m16954 audit, ``docs/savitsky_audit_m16954.md``)
    --------------------------------------------------------------
    This function does **not** solve λ, and never did: the lift equation
    ``CL(τ, λ, β) 0.5 ρ V² B² = Δ`` is a single equation in two unknowns, and the
    original implementation carried ``λ = 3.0`` from the initial guess to the return
    value in every case. Rather than pretend, λ is now an explicit input (default 3.0)
    and is echoed back unchanged, and the trim follows in closed form from
    :func:`trim_for_wetted_length`.

    Callers that need the solved equilibrium -- trim *and* wetted length from weight,
    speed and LCG -- must use :func:`equilibrium_trim`, which closes the system with
    the Savitsky centre-of-pressure condition.

    Parameters
    ----------
    weight_n, speed_ms, beam_m, deadrise_deg : float
        Hull weight (N), speed (m/s), beam (m), deadrise (degrees).
    wetted_length_ratio : float, optional
        λ = L_w / B assumed for the solve (default 3.0).
    water_density : float, optional
        Water density, kg/m^3.
    max_iterations, tolerance : int, float
        Accepted for backwards compatibility and ignored: the lift equation inverts
        in closed form, so there is no iteration to converge.

    Returns
    -------
    tuple of float
        ``(trim_deg, wetted_length_ratio)``; the second element is the input λ.

    Raises
    ------
    PlaningPredictionError
        If an input is non-positive, or if the trim the lift equation asks for lies
        outside ``[TRIM_MIN_DEG, TRIM_MAX_DEG]`` (the message reports the required
        trim, not a spurious convergence failure).
    """
    del max_iterations, tolerance  # closed-form solve; kept in the signature only
    trim_deg = trim_for_wetted_length(
        weight_n,
        speed_ms,
        beam_m,
        deadrise_deg,
        wetted_length_ratio,
        water_density=water_density,
    )
    return (trim_deg, wetted_length_ratio)


def resistance_components(
    speed_ms: float,
    wetted_length_m: float,
    beam_m: float,
    deadrise_deg: float,
    trim_deg: float,
    weight_n: float,
    *,
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
) -> dict[str, float]:
    """Compute the resistance components using the Savitsky method.

    The total resistance is the sum of friction drag and pressure drag:

        R_T = R_F + R_P

    where:
    * R_F = 0.5 ρ V² S_wet CF(Re) is the friction drag on the wetted area, with
      CF from the ITTC-57 friction line: CF = 0.075 / (log10 Re - 2)^2.
    * R_P is the pressure drag from the Savitsky correlation, which depends on
      the lift coefficient, trim angle, and deadrise.

    Parameters
    ----------
    speed_ms : float
        Forward speed, m/s.
    wetted_length_m : float
        Mean wetted length along the keel, m.
    beam_m : float
        Beam (width) of the planing surface, m.
    deadrise_deg : float
        Deadrise angle β, degrees.
    trim_deg : float
        Running trim angle τ, degrees.
    weight_n : float
        Hull weight (displacement), N.
    water_density : float, optional
        Water density, kg/m^3 (default: sea water at 15°C).
    kinematic_viscosity : float, optional
        Kinematic viscosity, m^2/s (default: sea water at 15°C).

    Returns
    -------
    dict[str, float]
        Dictionary with keys:
        * 'reynolds_number': Reynolds number V L_w / ν
        * 'friction_coefficient': ITTC-57 CF
        * 'wetted_area_m2': wetted surface area S_wet
        * 'friction_drag_n': R_F
        * 'pressure_drag_n': R_P
        * 'total_drag_n': R_T = R_F + R_P
        * 'lift_coefficient': CL
    """
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if wetted_length_m <= 0.0:
        raise PlaningPredictionError(f"wetted_length_m must be positive, got {wetted_length_m}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")

    # Reynolds number based on wetted length
    re = speed_ms * wetted_length_m / kinematic_viscosity

    # ITTC-57 friction line
    cf = 0.075 / (math.log10(re) - 2.0) ** 2

    # Wetted area: for a prismatic hull with constant deadrise, the wetted area is
    # approximately the wetted length × beam / cos(beta) (the bottom is inclined).
    # More precisely, the wetted area is the planing surface area.
    beta_rad = math.radians(deadrise_deg)
    wetted_area = wetted_length_m * beam_m / math.cos(beta_rad)

    # Friction drag
    q = 0.5 * water_density * speed_ms ** 2
    friction_drag = cf * q * wetted_area

    # Lift coefficient
    lam = wetted_length_m / beam_m
    cl = lift_coefficient(trim_deg, lam, deadrise_deg)

    # Pressure drag from Savitsky correlation
    # The pressure drag is related to the lift and the trim angle:
    # R_P = Δ tan(τ) + correction for the pressure distribution
    # Savitsky gives: R_P / Δ = tan(τ) + 0.008 / sqrt(τ)
    # where τ is in degrees.
    tau_rad = math.radians(trim_deg)
    pressure_drag = weight_n * (math.tan(tau_rad) + 0.008 / math.sqrt(trim_deg))

    total_drag = friction_drag + pressure_drag

    return {
        "reynolds_number": re,
        "friction_coefficient": cf,
        "wetted_area_m2": wetted_area,
        "friction_drag_n": friction_drag,
        "pressure_drag_n": pressure_drag,
        "total_drag_n": total_drag,
        "lift_coefficient": cl,
    }


def planing_prediction(
    weight_n: float,
    speed_ms: float,
    length_m: float,
    beam_m: float,
    deadrise_deg: float,
    *,
    lcg_m: float | None = None,
    wetted_length_ratio: float | None = None,
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
) -> PlaningPrediction:
    """Complete Savitsky planing prediction at one speed.

    Two modes, and the result says which one produced it
    ----------------------------------------------------
    ``lcg_m given`` (recommended)
        :func:`equilibrium_trim` solves τ *and* λ from the two-equation closure
        (lift = weight, centre of pressure under the LCG). ``solution_method`` is
        ``'classic_equilibrium'`` and ``wetted_length_ratio_solved`` is ``True``.
    ``lcg_m omitted`` (legacy behaviour)
        λ is an input assumption (``wetted_length_ratio``, default 3.0) and only the
        trim is solved. ``solution_method`` is ``'assumed_wetted_length'`` and
        ``wetted_length_ratio_solved`` is ``False``.

    In both modes the result carries a ``warnings`` tuple: every clamp, extrapolation
    or assumption is reported there instead of being applied silently (m16954 audit).

    Parameters
    ----------
    weight_n : float
        Hull weight (displacement), N.
    speed_ms : float
        Forward speed, m/s.
    length_m : float
        Hull length on the waterline, m.
    beam_m : float
        Beam (width) of the planing surface, m.
    deadrise_deg : float
        Deadrise angle β, degrees.
    lcg_m : float, optional
        Longitudinal centre of gravity, m forward of the transom along the keel.
        Supplying it switches on the two-equation equilibrium solve.
    wetted_length_ratio : float, optional
        Assumed λ, used only when ``lcg_m`` is omitted.
    water_density, kinematic_viscosity : float, optional
        Fluid properties.

    Returns
    -------
    PlaningPrediction

    Raises
    ------
    PlaningPredictionError
        If inputs are invalid, or (equilibrium mode) if no planing equilibrium exists
        for the requested LCG/speed/beam.
    """
    if weight_n <= 0.0:
        raise PlaningPredictionError(f"weight_n must be positive, got {weight_n}")
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if length_m <= 0.0:
        raise PlaningPredictionError(f"length_m must be positive, got {length_m}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")

    warnings: list[str] = []

    # Froude number
    fn = speed_ms / math.sqrt(STANDARD_GRAVITY * length_m)

    if lcg_m is not None:
        solution = equilibrium_trim(
            weight_n,
            speed_ms,
            beam_m,
            deadrise_deg,
            lcg_m,
            water_density=water_density,
            max_wetted_length_m=length_m,
        )
        trim_deg = solution.trim_deg
        lam = solution.wetted_length_ratio
        solved = True
        method = "classic_equilibrium"
        lift_coefficient_reported = classic_lift_coefficient(
            trim_deg, lam, solution.speed_coefficient, deadrise_deg
        )
    else:
        lam = (
            DEFAULT_WETTED_LENGTH_RATIO
            if wetted_length_ratio is None
            else wetted_length_ratio
        )
        trim_deg = trim_for_wetted_length(
            weight_n, speed_ms, beam_m, deadrise_deg, lam, water_density=water_density
        )
        solved = False
        method = "assumed_wetted_length"
        lift_coefficient_reported = lift_coefficient(trim_deg, lam, deadrise_deg)
        warnings.append(
            f"wetted length ratio assumed lambda={lam:.3f}; the trim solve is one "
            "equation in two unknowns, so lambda is NOT solved -- pass lcg_m to "
            "equilibrium_trim()/planing_prediction() for the solved equilibrium"
        )

    # Physical cap: the wetted length cannot exceed the hull length. Instead of
    # clamping silently (which used to break lift equilibrium and misreport CL), the
    # clamp is reported and the trim is re-solved at the capped lambda.
    wetted_length_m = lam * beam_m
    if wetted_length_m > length_m:
        lam_clamped = length_m / beam_m
        warnings.append(
            f"wetted length {wetted_length_m:.3f} m exceeds the hull length "
            f"{length_m:.3f} m: lambda clamped from {lam:.3f} to {lam_clamped:.3f}, "
            "the craft is fully wetted / in the hump, outside the planing box"
        )
        lam = lam_clamped
        solved = False
        if lcg_m is None:
            trim_deg = trim_for_wetted_length(
                weight_n,
                speed_ms,
                beam_m,
                deadrise_deg,
                lam,
                water_density=water_density,
            )
            lift_coefficient_reported = lift_coefficient(trim_deg, lam, deadrise_deg)
        wetted_length_m = lam * beam_m

    # Resistance components
    resistance = resistance_components(
        speed_ms, wetted_length_m, beam_m, deadrise_deg, trim_deg, weight_n,
        water_density=water_density,
        kinematic_viscosity=kinematic_viscosity,
    )

    # In equilibrium mode the CL that closes the lift equation is the classic one.
    if solved:
        resistance["lift_coefficient"] = lift_coefficient_reported

    # Effective power
    effective_power = resistance["total_drag_n"] * speed_ms

    # Reynolds number
    re = resistance["reynolds_number"]

    warnings.extend(
        validity_warnings(fn, deadrise_deg, length_m, beam_m, trim_deg)
    )

    return PlaningPrediction(
        speed_ms=speed_ms,
        froude_number=fn,
        trim_deg=trim_deg,
        wetted_length_m=wetted_length_m,
        wetted_length_ratio=lam,
        lift_coefficient=resistance["lift_coefficient"],
        friction_drag_n=resistance["friction_drag_n"],
        pressure_drag_n=resistance["pressure_drag_n"],
        total_drag_n=resistance["total_drag_n"],
        effective_power_w=effective_power,
        weight_n=weight_n,
        reynolds_number=re,
        wetted_length_ratio_solved=solved,
        solution_method=method,
        warnings=tuple(warnings),
    )


def speed_sweep(
    weight_n: float,
    length_m: float,
    beam_m: float,
    deadrise_deg: float,
    speed_range_ms: tuple[float, ...],
    *,
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    lcg_m: float | None = None,
    failures: list[tuple[float, str]] | None = None,
) -> list[PlaningPrediction]:
    """Run the Savitsky prediction over a speed range.

    Parameters
    ----------
    weight_n : float
        Hull weight (displacement), N.
    length_m : float
        Hull length on the waterline, m.
    beam_m : float
        Beam (width) of the planing surface, m.
    deadrise_deg : float
        Deadrise angle β, degrees.
    speed_range_ms : tuple[float, ...]
        Speeds to evaluate, m/s.
    water_density : float, optional
        Water density, kg/m^3 (default: sea water at 15°C).
    kinematic_viscosity : float, optional
        Kinematic viscosity, m^2/s (default: sea water at 15°C).

    Returns
    -------
    list[PlaningPrediction]
        List of predictions, one per speed. Speeds that fail are skipped -- pass
        ``failures`` to see them, because a silently shortened sweep is not evidence
        (m16954 audit).
    """
    results = []
    for speed in speed_range_ms:
        try:
            prediction = planing_prediction(
                weight_n, speed, length_m, beam_m, deadrise_deg,
                water_density=water_density,
                kinematic_viscosity=kinematic_viscosity,
                lcg_m=lcg_m,
            )
            results.append(prediction)
        except PlaningPredictionError as exc:
            # If the prediction fails at this speed (e.g., too slow for planing),
            # skip it -- but record why when the caller asked for an audit trail.
            if failures is not None:
                failures.append((speed, str(exc)))
            continue
    return results
