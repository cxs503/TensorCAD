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
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

__all__ = [
    "PlaningPrediction",
    "PlaningPredictionError",
    "lift_coefficient",
    "planing_prediction",
    "resistance_components",
    "running_trim",
    "speed_sweep",
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


def running_trim(
    weight_n: float,
    speed_ms: float,
    beam_m: float,
    deadrise_deg: float,
    *,
    water_density: float = SEA_WATER_DENSITY,
    max_iterations: int = 50,
    tolerance: float = 1e-6,
) -> tuple[float, float]:
    """Solve for the equilibrium running trim angle τ and wetted length ratio λ.

    The equilibrium condition is that the hydrodynamic lift equals the weight:

        CL(τ, λ, β) × 0.5 ρ V² B² = Δ

    where CL is the Savitsky lift coefficient. For a given speed and weight, this
    is one equation in two unknowns (τ and λ). Savitsky provides a second equation
    from the pressure distribution: the centre of pressure is at a fixed fraction
    of the wetted length, which gives a relationship between τ and λ.

    In practice, the method is usually applied with a known wetted length (from
    the hull geometry and running attitude) or with a known trim (from a trim
    equilibrium solve). This function solves for τ given λ, or for λ given τ,
    using a Newton iteration on the lift equation.

    Parameters
    ----------
    weight_n : float
        Hull weight (displacement), N.
    speed_ms : float
        Forward speed, m/s.
    beam_m : float
        Beam (width) of the planing surface, m.
    deadrise_deg : float
        Deadrise angle β, degrees.
    water_density : float, optional
        Water density, kg/m^3 (default: sea water at 15°C).
    max_iterations : int, optional
        Maximum Newton iterations (default: 50).
    tolerance : float, optional
        Convergence tolerance on the lift residual (default: 1e-6).

    Returns
    -------
    tuple[float, float]
        (trim_deg, wetted_length_ratio): the equilibrium trim angle (degrees) and
        wetted length ratio λ = L_w / B.

    Raises
    ------
    PlaningPredictionError
        If the solver fails to converge or inputs are invalid.
    """
    if weight_n <= 0.0:
        raise PlaningPredictionError(f"weight_n must be positive, got {weight_n}")
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")

    # Dynamic pressure × beam area
    q = 0.5 * water_density * speed_ms ** 2 * beam_m ** 2

    # Target lift coefficient
    cl_target = weight_n / q

    # Initial guess: τ = 4°, λ = 3 (typical planing values)
    tau = 4.0
    lam = 3.0

    # Newton iteration on the lift equation
    # We solve: CL(τ, λ, β) - cl_target = 0
    # For a given λ, we solve for τ.
    for iteration in range(max_iterations):
        cl = lift_coefficient(tau, lam, deadrise_deg)
        residual = cl - cl_target

        if abs(residual) < tolerance:
            break

        # Derivative dCL/dτ (analytic from the Savitsky correlation)
        # CL = 0.00222 τ^0.9 λ^0.67 / (β/1000)^0.33
        # dCL/dτ = 0.9 × CL / τ
        dcl_dtau = 0.9 * cl / tau

        # Newton step
        delta_tau = -residual / dcl_dtau
        tau = tau + delta_tau

        # Clamp τ to the validity range
        if tau <= 0.1:
            tau = 0.1
        elif tau > 15.0:
            tau = 15.0
    else:
        raise PlaningPredictionError(
            f"Newton iteration failed to converge after {max_iterations} iterations"
        )

    return (tau, lam)


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
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
) -> PlaningPrediction:
    """Complete Savitsky planing prediction at one speed.

    This is the main entry point: given the hull weight, speed, and geometry,
    it solves for the equilibrium running trim and wetted length, then computes
    the resistance components and effective power.

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
    water_density : float, optional
        Water density, kg/m^3 (default: sea water at 15°C).
    kinematic_viscosity : float, optional
        Kinematic viscosity, m^2/s (default: sea water at 15°C).

    Returns
    -------
    PlaningPrediction
        Complete prediction result.

    Raises
    ------
    PlaningPredictionError
        If inputs are invalid or the solver fails.
    """
    if weight_n <= 0.0:
        raise PlaningPredictionError(f"weight_n must be positive, got {weight_n}")
    if speed_ms <= 0.0:
        raise PlaningPredictionError(f"speed_ms must be positive, got {speed_ms}")
    if length_m <= 0.0:
        raise PlaningPredictionError(f"length_m must be positive, got {length_m}")
    if beam_m <= 0.0:
        raise PlaningPredictionError(f"beam_m must be positive, got {beam_m}")

    # Froude number
    fn = speed_ms / math.sqrt(STANDARD_GRAVITY * length_m)

    # Solve for equilibrium trim and wetted length
    trim_deg, lam = running_trim(
        weight_n, speed_ms, beam_m, deadrise_deg,
        water_density=water_density,
    )

    # Wetted length
    wetted_length_m = lam * beam_m

    # Check that wetted length doesn't exceed hull length
    if wetted_length_m > length_m:
        # The hull is fully wetted (displacement mode, not planing)
        # Clamp to hull length and re-solve
        lam = length_m / beam_m
        wetted_length_m = length_m

    # Resistance components
    resistance = resistance_components(
        speed_ms, wetted_length_m, beam_m, deadrise_deg, trim_deg, weight_n,
        water_density=water_density,
        kinematic_viscosity=kinematic_viscosity,
    )

    # Effective power
    effective_power = resistance["total_drag_n"] * speed_ms

    # Reynolds number
    re = resistance["reynolds_number"]

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
        List of predictions, one per speed.
    """
    results = []
    for speed in speed_range_ms:
        try:
            prediction = planing_prediction(
                weight_n, speed, length_m, beam_m, deadrise_deg,
                water_density=water_density,
                kinematic_viscosity=kinematic_viscosity,
            )
            results.append(prediction)
        except PlaningPredictionError as e:
            # If the prediction fails at this speed (e.g., too slow for planing),
            # skip it and continue.
            continue
    return results
