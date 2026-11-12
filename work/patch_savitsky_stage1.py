"""Apply the m16954 audit fixes to tensorcad/savitsky.py.

Anchored, idempotent-free (fails loudly on missing anchors), minimal-diff patch.
"""
from pathlib import Path

path = Path("tensorcad/savitsky.py")
src = path.read_text()


def replace_once(text: str, old: str, new: str, label: str) -> str:
    assert text.count(old) == 1, f"anchor {label!r} found {text.count(old)} times"
    return text.replace(old, new)


# ---------------------------------------------------------------- __all__
src = replace_once(
    src,
    '''__all__ = [
    "PlaningPrediction",
    "PlaningPredictionError",
    "lift_coefficient",
    "planing_prediction",
    "resistance_components",
    "running_trim",
    "speed_sweep",
]''',
    '''__all__ = [
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
]''',
    "__all__",
)

# ---------------------------------------------------------------- docstring
src = replace_once(
    src,
    '''Outside this box the predictions are extrapolations and the errors grow.
"""''',
    '''Outside this box the predictions are extrapolations and the errors grow.

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
"""''',
    "module docstring",
)

# ---------------------------------------------------------------- dataclass
src = replace_once(
    src,
    '''    reynolds_number : float
        Reynolds number V L_w / ν (based on wetted length).
    """''',
    '''    reynolds_number : float
        Reynolds number V L_w / ν (based on wetted length).
    wetted_length_ratio_solved : bool
        ``True`` only when ``λ`` came out of :func:`equilibrium_trim` (the two-equation
        closure). ``False`` means λ was an input assumption or a hull-length clamp.
    solution_method : str
        ``'classic_equilibrium'`` or ``'assumed_wetted_length'``.
    warnings : tuple of str
        Honesty channel: every clamp, extrapolation or assumption applied to this
        result. Empty only when λ was solved and the case is inside the validity box.
    """''',
    "dataclass docstring",
)

src = replace_once(
    src,
    '''    weight_n: float
    reynolds_number: float

    def as_dict(self) -> dict[str, Any]:''',
    '''    weight_n: float
    reynolds_number: float
    wetted_length_ratio_solved: bool = False
    solution_method: str = "assumed_wetted_length"
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:''',
    "dataclass fields",
)

src = replace_once(
    src,
    '''            "weight_n": self.weight_n,
            "reynolds_number": self.reynolds_number,
        }''',
    '''            "weight_n": self.weight_n,
            "reynolds_number": self.reynolds_number,
            "wetted_length_ratio_solved": self.wetted_length_ratio_solved,
            "solution_method": self.solution_method,
            "warnings": list(self.warnings),
        }''',
    "as_dict",
)

# ------------------------------------------------ new classic-closure code
NEW_BLOCK = '''
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

    lo_trim, hi_trim = TRIM_MIN_DEG, TRIM_MAX_DEG
    cp_lo, cp_hi = cp_position(lo_trim), cp_position(hi_trim)
    if cp_lo is None or cp_hi is None:
        raise PlaningPredictionError(
            "no planing equilibrium: the lift equation cannot be satisfied for any "
            f"wetted length ratio up to {lam_max:.2f} inside the trim window "
            f"{TRIM_MIN_DEG}-{TRIM_MAX_DEG} deg (CL needed {cl_target:.4f}, "
            f"Cv={cv:.3f})"
        )
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


'''

src = replace_once(
    src,
    "def running_trim(\n",
    NEW_BLOCK + "def running_trim(\n",
    "insert classic closure",
)

path.write_text(src)
print("stage 1 written")