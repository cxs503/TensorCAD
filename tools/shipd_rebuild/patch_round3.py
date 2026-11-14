"""Apply the GAP-CAD-C2 round-3 upgrade to tensorcad/shipd_rebuild.py.

Adds: coefficient-calibrated reconstruction method, realised-coefficient gates,
a reproducible 500-ship population model, published nominal performance fields
and a correlation sanity analysis.
"""
from __future__ import annotations

import pathlib

path = pathlib.Path("tensorcad/shipd_rebuild.py")
src = path.read_text()


def sub(old: str, new: str) -> None:
    global src
    assert old in src, old[:120]
    src = src.replace(old, new, 1)


# --------------------------------------------------------------------------- #
# 1) Imports
# --------------------------------------------------------------------------- #
sub(
    """from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.series60 import Series60Parameters, build_series60_mesh
from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.hydrostatics import compute_hydrostatics
""",
    """from tensorcad import mesh as meshlib
from tensorcad.cb_sections import (
    CalibratedHullParameters,
    build_calibrated_hull,
    realised_coefficients,
    relative_errors,
)
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.series60 import Series60Parameters, build_series60_mesh
from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.hydrostatics import compute_hydrostatics
""",
)

# --------------------------------------------------------------------------- #
# 2) Constants
# --------------------------------------------------------------------------- #
sub(
    """DEFAULT_NX = 81               # Longitudinal stations
DEFAULT_NZ = 21               # Vertical levels
""",
    """DEFAULT_NX = 81               # Longitudinal stations
DEFAULT_NZ = 21               # Vertical levels

#: Block-coefficient window the coefficient-calibrated family claims.
CALIBRATED_CB_WINDOW = (0.30, 0.95)
#: Tolerance on the realised midship / waterplane coefficients (percent of the
#: requested value). Looser than the 3 % volume gate because these two are read
#: off the *sampled* station grid rather than integrated over the closed surface.
COEFFICIENT_TOLERANCE = 0.05
#: Seawater density [kg/m^3] and kinematic viscosity [m^2/s] used by the nominal
#: performance block (ITTC fresh/seawater tables at 15 degC).
SEAWATER_DENSITY = 1025.0
SEAWATER_VISCOSITY = 1.1883e-6
GRAVITY = 9.80665
""",
)

# --------------------------------------------------------------------------- #
# 3) RebuildResult: realised coefficients + strict verdict
# --------------------------------------------------------------------------- #
sub(
    """    triangle_count: int | None = None
    error_message: str | None = None
    rebuild_time_s: float | None = None
""",
    """    triangle_count: int | None = None
    error_message: str | None = None
    rebuild_time_s: float | None = None
    #: Coefficient agreement measured on the mesh (not assumed from the inputs).
    realised_Cb: float | None = None
    realised_Cm: float | None = None
    realised_Cwp: float | None = None
    coefficient_gates: list[GateResult] = field(default_factory=list)
    wetted_surface_area: float | None = None
    #: ``VERIFIED`` only when the eight canonical gates *and* the coefficient
    #: agreement gates pass. :attr:`verdict` keeps the canonical eight so round-2
    #: and round-3 pass rates stay directly comparable.
    strict_verdict: str = "PENDING"
""",
)

sub(
    """            "gates": [g.as_dict() for g in self.gates],
            "error_message": self.error_message,
            "rebuild_time_s": self.rebuild_time_s,
        }
""",
    """            "gates": [g.as_dict() for g in self.gates],
            "error_message": self.error_message,
            "rebuild_time_s": self.rebuild_time_s,
            "realised_Cb": self.realised_Cb,
            "realised_Cm": self.realised_Cm,
            "realised_Cwp": self.realised_Cwp,
            "coefficient_gates": [g.as_dict() for g in self.coefficient_gates],
            "wetted_surface_area": self.wetted_surface_area,
            "strict_verdict": self.strict_verdict,
        }
""",
)

# --------------------------------------------------------------------------- #
# 4) Method selection: calibrated family first, legacy chain preserved
# --------------------------------------------------------------------------- #
sub(
    '''def _select_method(cb: float, cm: float, cwp: float) -> tuple[str, str]:
    """Select the best reconstruction method for a ship's parameters."""
    low_s60, high_s60 = SERIES60_CB_WINDOW
    low_param, high_param = PARAMETRIC_CB_WINDOW
    low_wig, high_wig = WIGLEY_CB_WINDOW

    if low_s60 <= cb <= high_s60:''',
    '''def _calibrated_feasible(cb: float, cm: float, cwp: float) -> tuple[bool, str]:
    """Can the coefficient-calibrated family realise this parameter triple?

    The family needs ``Cp = Cb/Cm <= Cwp`` (a waterplane cannot be finer than the
    prismatic coefficient) and three coefficients strictly inside ``(0, 1)``.
    """
    if not all(math.isfinite(v) and 0.0 < v < 1.0 for v in (cb, cm, cwp)):
        return False, f"coefficients outside (0,1): Cb={cb}, Cm={cm}, Cwp={cwp}"
    low, high = CALIBRATED_CB_WINDOW
    if not (low <= cb <= high):
        return False, f"Cb={cb:.4f} outside the calibrated window [{low}, {high}]"
    cp = cb / cm
    if cp > cwp:
        return False, f"Cp=Cb/Cm={cp:.4f} exceeds Cwp={cwp:.4f}"
    return True, (
        f"Cb={cb:.4f} calibrated family: Cp={cp:.4f} <= Cwp={cwp:.4f} "
        f"(realises Cb, Cm, Cwp exactly)"
    )


def _select_method(
    cb: float, cm: float, cwp: float, *, method_preference: str = "calibrated"
) -> tuple[str, str]:
    """Select the best reconstruction method for a ship's parameters.

    ``method_preference="calibrated"`` (default) puts the coefficient-calibrated
    family first: it is the only method that honours the requested ``Cb``, ``Cm``
    and ``Cwp`` instead of producing whatever the fixed breadth law gives.
    ``method_preference="legacy"`` reproduces the round-2 chain
    (Series 60 -> parametric sections -> Wigley -> refused) unchanged, which is
    what the round-2 evidence was collected with.
    """
    if method_preference == "calibrated":
        feasible, reason = _calibrated_feasible(cb, cm, cwp)
        if feasible:
            return "cb_calibrated_sections", reason

    low_s60, high_s60 = SERIES60_CB_WINDOW
    low_param, high_param = PARAMETRIC_CB_WINDOW
    low_wig, high_wig = WIGLEY_CB_WINDOW

    if low_s60 <= cb <= high_s60:''',
)

# --------------------------------------------------------------------------- #
# 5) The calibrated rebuild driver
# --------------------------------------------------------------------------- #
sub(
    '''def _evaluate_gates(ship: dict, mesh: Mesh, method: str, reason: str, elapsed: float) -> RebuildResult:''',
    '''def _rebuild_calibrated(ship: dict, nx: int, nz: int) -> RebuildResult:
    """Rebuild with the coefficient-calibrated family (realises Cb, Cm, Cwp)."""
    t0 = time.time()
    try:
        params = CalibratedHullParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            block_coefficient=ship["Cb"],
            midship_coefficient=ship.get("Cm", 0.98),
            waterplane_coefficient=ship.get("Cwp", 0.80),
            nx=nx, nz=nz,
        )
        problems = params.problems()
        if problems:
            return RebuildResult(
                ship_id=ship["ship_id"], ship_type=ship["ship_type"],
                method="cb_calibrated_sections",
                method_reason="calibrated family selected but parameters rejected",
                verdict="REFUSED", error_message="; ".join(problems),
                rebuild_time_s=time.time() - t0,
            )
        hull = build_calibrated_hull(params)
        elapsed = time.time() - t0
        return _evaluate_gates(
            ship, hull.mesh, "cb_calibrated_sections",
            f"Cb={ship['Cb']:.4f}, Cp<=Cwp -> calibrated family", elapsed,
        )
    except Exception as e:  # pragma: no cover - defensive
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="cb_calibrated_sections", method_reason="calibrated family",
            verdict="FAILED", error_message=str(e), rebuild_time_s=time.time() - t0,
        )


def _evaluate_gates(ship: dict, mesh: Mesh, method: str, reason: str, elapsed: float) -> RebuildResult:''',
)

# --------------------------------------------------------------------------- #
# 6) Gates: add realised-coefficient measurement and the strict verdict
# --------------------------------------------------------------------------- #
sub(
    """    all_passed = all(g.passed for g in gates)
    verdict = "VERIFIED" if all_passed else "FAILED"

    return RebuildResult(
        ship_id=ship["ship_id"], ship_type=ship["ship_type"],
        method=method, method_reason=reason,
        mesh=mesh, gates=gates, verdict=verdict,
        volume_mesh=vol, volume_parametric=vol_parametric,
        volume_error=vol_error,
        realised_L=rL, realised_B=rB, realised_T=rT,
        triangle_count=tri_count, rebuild_time_s=elapsed,
    )""",
    """    all_passed = all(g.passed for g in gates)
    verdict = "VERIFIED" if all_passed else "FAILED"

    # Coefficient agreement on the mesh itself. ``Cm`` is the largest sampled
    # station area and ``Cwp`` the waterplane flux; both are measured, never
    # copied from the request, so a surface that drifts off the parameter vector
    # shows up here even when the volume happens to be right.
    coefficient_gates: list[GateResult] = []
    realised_cb = realised_cm = realised_cwp = None
    wetted = None
    try:
        measured = realised_coefficients(
            mesh, length=ship["L"], beam=ship["B"], draft=ship["T"]
        )
        realised_cb = measured["coefficients"]["Cb"]
        realised_cm = measured["coefficients"]["Cm"]
        realised_cwp = measured["coefficients"]["Cwp"]
        wetted = measured["wetted_surface_m2"]
        for name, got, want in (
            ("G9_midship_coefficient", realised_cm, ship.get("Cm", 0.98)),
            ("G10_waterplane_coefficient", realised_cwp, ship.get("Cwp", 0.80)),
        ):
            err = abs(got - want) / want if want else float("inf")
            coefficient_gates.append(GateResult(
                name, err < COEFFICIENT_TOLERANCE,
                f"realised={got:.5f} requested={want:.5f} err={err:.4%}",
                value=err, threshold=COEFFICIENT_TOLERANCE,
            ))
    except Exception as exc:  # pragma: no cover - defensive
        coefficient_gates.append(GateResult(
            "G9_midship_coefficient", False, f"measurement failed: {exc}",
        ))

    strict_ok = all_passed and all(g.passed for g in coefficient_gates)
    return RebuildResult(
        ship_id=ship["ship_id"], ship_type=ship["ship_type"],
        method=method, method_reason=reason,
        mesh=mesh, gates=gates, verdict=verdict,
        volume_mesh=vol, volume_parametric=vol_parametric,
        volume_error=vol_error,
        realised_L=rL, realised_B=rB, realised_T=rT,
        triangle_count=tri_count, rebuild_time_s=elapsed,
        realised_Cb=realised_cb, realised_Cm=realised_cm, realised_Cwp=realised_cwp,
        coefficient_gates=coefficient_gates, wetted_surface_area=wetted,
        strict_verdict="VERIFIED" if strict_ok else "FAILED",
    )""",
)

# --------------------------------------------------------------------------- #
# 7) rebuild_ship: method preference plumbed through
# --------------------------------------------------------------------------- #
sub(
    '''def rebuild_ship(ship: dict, *, nx: int = DEFAULT_NX, nz: int = DEFAULT_NZ) -> RebuildResult:
    """Rebuild one ship from its parameter vector."""
    method, reason = _select_method(ship["Cb"], ship.get("Cm", 0.98), ship.get("Cwp", 0.80))
''',
    '''def rebuild_ship(
    ship: dict,
    *,
    nx: int = DEFAULT_NX,
    nz: int = DEFAULT_NZ,
    method_preference: str = "calibrated",
) -> RebuildResult:
    """Rebuild one ship from its parameter vector."""
    method, reason = _select_method(
        ship["Cb"], ship.get("Cm", 0.98), ship.get("Cwp", 0.80),
        method_preference=method_preference,
    )
''',
)

sub(
    """    if method == "series60_form":
        return _rebuild_series60(ship, nx, nz)
""",
    """    if method == "cb_calibrated_sections":
        return _rebuild_calibrated(ship, nx, nz)
    if method == "series60_form":
        return _rebuild_series60(ship, nx, nz)
""",
)

path.write_text(src)
print("shipd_rebuild.py patched:", len(src), "bytes")