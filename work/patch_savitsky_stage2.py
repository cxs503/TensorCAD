"""Stage 2: rewrite running_trim and planing_prediction (m16954 audit fixes)."""
from pathlib import Path

path = Path("tensorcad/savitsky.py")
src = path.read_text()

start = src.index("def running_trim(")
end = src.index("def resistance_components(")
NEW_RUNNING_TRIM = '''def running_trim(
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


'''
src = src[:start] + NEW_RUNNING_TRIM + src[end:]

# ------------------------------------------------------------------ planing_prediction
start = src.index("def planing_prediction(")
end = src.index("def speed_sweep(")
NEW_PREDICTION = '''def planing_prediction(
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


'''
src = src[:start] + NEW_PREDICTION + src[end:]

# ---------------------------------------------------------------- speed_sweep: report failures
src = src.replace(
    '''    Returns
    -------
    list[PlaningPrediction]
        List of predictions, one per speed.
    """
    results = []''',
    '''    Returns
    -------
    list[PlaningPrediction]
        List of predictions, one per speed. Speeds that fail are skipped -- pass
        ``failures`` to see them, because a silently shortened sweep is not evidence
        (m16954 audit).
    """
    results = []''',
)

src = src.replace(
    '''    speed_range_ms: tuple[float, ...],
    *,
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
) -> list[PlaningPrediction]:''',
    '''    speed_range_ms: tuple[float, ...],
    *,
    water_density: float = SEA_WATER_DENSITY,
    kinematic_viscosity: float = SEA_WATER_KINEMATIC_VISCOSITY,
    lcg_m: float | None = None,
    failures: list[tuple[float, str]] | None = None,
) -> list[PlaningPrediction]:''',
)

src = src.replace(
    '''                water_density=water_density,
                kinematic_viscosity=kinematic_viscosity,
            )
            results.append(prediction)
        except PlaningPredictionError as e:
            # If the prediction fails at this speed (e.g., too slow for planing),
            # skip it and continue.
            continue
    return results''',
    '''                water_density=water_density,
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
    return results''',
)

path.write_text(src)
print("stage 2 written")