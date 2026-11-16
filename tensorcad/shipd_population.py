"""Reproducible ship-parameter population model and nominal performance fields.

Two jobs, both needed by the GAP-CAD-C2 rebuild campaign:

1. :func:`build_shipd_database` builds the 500-ship parameter table the rebuild
   pipeline consumes, including the longitudinal centre of buoyancy ``lcb`` (a
   signed fraction of ``L``) -- a principal particular every real hull has and the
   five scalar coefficients cannot imply. Round 2 of the campaign cited
   ``artifacts/shipd_rebuild/shipd_500.json`` but shipped no generator -- the file
   was not in the tree, so its 59.4 % pass rate could not be re-derived. This
   module *is* the generator: seeded, deterministic, no external file.

2. :func:`nominal_performance_fields` turns a parameter vector into the
   performance quantities a concept design quotes -- design speed, Reynolds
   number, ITTC-57 friction coefficient, Holtrop & Mennen form factor and
   wetted-surface area, hence a nominal (viscous-only) resistance. Those are the
   "nominal fields" the rebuild is cross-checked against: if the rebuilt surface
   carries the same displacement, waterplane and wetted area as the parameter
   vector implies, the geometry is a usable input to a resistance or stability
   study.

Honest provenance
=================
The MIT ShipD distribution is *not* used: there is no network access and no local
copy, so the population below is a **synthetic model built from the design-space
ranges published in the naval-architecture literature** (Watson, *Practical Ship
Design*, ch. 2 principal-dimension tables; Schneekluth & Bertram, *Ship Design for
Efficiency and Economy*, ch. 2 coefficient relations). Ranges are typical service
values per ship type, not fitted to any single vessel, and the sampling code is the
only ground truth. Coefficients follow two textbook approximations:

``Cwp ~= (1 + 2*Cb) / 3``      Schneekluth & Bertram coefficient relation
``Cm  ~= 1 - 0.20*(1 - Cb)**1.8``  fineness trend, ``Cm`` 0.93 (Cb 0.45) to 0.995 (Cb 0.85)

Both carry a small uniform scatter so the population is not degenerate on the two
regression lines. A sampled ``Cwp`` is raised when it would violate
``Cp = Cb/Cm <= Cwp`` (a waterplane cannot be finer than the prismatic
coefficient); those ships are flagged ``cwp_adjusted`` rather than silently
repaired.
"""
from __future__ import annotations

import math
import random
from typing import Any, Sequence

GRAVITY = 9.80665
SEAWATER_DENSITY = 1025.0
SEAWATER_VISCOSITY = 1.1883e-6

#: type -> (L_min [m], L_max [m], L/B, B/T, Cb_min, Cb_max, Fn_min, Fn_max)
#: Ranges are typical service values; see module docstring for provenance.
SHIP_TYPE_POPULATION: dict[str, tuple[float, float, float, float, float, float, float, float]] = {
    "container":     (150.0, 400.0, 7.0, 2.6, 0.55, 0.68, 0.20, 0.26),
    "tanker":        (180.0, 340.0, 6.0, 3.2, 0.78, 0.86, 0.12, 0.16),
    "bulk_carrier":  (150.0, 300.0, 6.2, 2.9, 0.72, 0.84, 0.12, 0.16),
    "general_cargo": (90.0, 180.0, 6.4, 2.7, 0.65, 0.78, 0.14, 0.19),
    "roro":          (120.0, 240.0, 7.4, 3.0, 0.60, 0.72, 0.18, 0.24),
    "cruise":        (200.0, 340.0, 8.2, 3.4, 0.62, 0.74, 0.20, 0.26),
    "ferry":         (60.0, 160.0, 6.6, 3.0, 0.55, 0.68, 0.28, 0.38),
    "tug":           (25.0, 45.0, 4.2, 2.6, 0.50, 0.60, 0.28, 0.36),
    "naval_escort":  (100.0, 180.0, 8.6, 3.4, 0.45, 0.55, 0.30, 0.42),
    "fishing":       (16.0, 60.0, 4.0, 2.2, 0.45, 0.60, 0.22, 0.32),
}


#: type -> (LCB_min, LCB_max) as a signed fraction of ``L``, positive forward of
#: mid-length. Typical *service* positions of the longitudinal centre of buoyancy,
#: recollected offline from the principal-dimension / coefficient tables in Watson,
#: *Practical Ship Design* (ch. 2) and Schneekluth & Bertram, *Ship Design for
#: Efficiency and Economy* (ch. 2, LCB vs Fn plots): full hulls carry their centroid
#: forward of amidships, fine hulls slightly aft.
#:
#: Provenance is ``offline_recollection, verified=false`` -- the same grade as the
#: fleet manifest. These ranges are wide enough to bracket every value quoted in
#: the campaign (KCS -1.48 %L, KVLCC2 +3.48 %L) and are *targets for the rebuild
#: control loop*, not a fit to any particular vessel. They must not be read as a
#: claim about any individual ship.
LCB_FRACTION_RANGE: dict[str, tuple[float, float]] = {
    "container":     (-0.020, -0.005),
    "tanker":        (+0.015, +0.035),
    "bulk_carrier":  (+0.010, +0.030),
    "general_cargo": (+0.000, +0.020),
    "roro":          (-0.010, +0.015),
    "cruise":        (-0.015, +0.005),
    "ferry":         (-0.010, +0.010),
    "tug":           (-0.015, +0.015),
    "naval_escort":  (-0.015, +0.010),
    "fishing":       (-0.020, +0.010),
}


def sample_lcb_fraction(ship_type: str, rng: random.Random) -> float:
    """One LCB target for ``ship_type``, uniform over :data:`LCB_FRACTION_RANGE`."""
    try:
        low, high = LCB_FRACTION_RANGE[ship_type]
    except KeyError as exc:  # pragma: no cover - table and population are paired
        raise KeyError(f"no LCB range for ship type {ship_type!r}") from exc
    return rng.uniform(low, high)


def _midship_from_block(cb: float) -> float:
    """Textbook fineness trend ``Cm = 1 - 0.20*(1-Cb)**1.8``."""
    return 1.0 - 0.20 * (1.0 - cb) ** 1.8


def _waterplane_from_block(cb: float) -> float:
    """Schneekluth & Bertram approximation ``Cwp = (1 + 2*Cb)/3``."""
    return (1.0 + 2.0 * cb) / 3.0


def build_shipd_database(n_ships: int = 500, *, seed: int = 42) -> list[dict[str, Any]]:
    """Deterministic parameter table: ``n_ships`` vessels over ten ship types.

    Type assignment is round robin so every type is represented at any fleet size;
    within a type the principal dimensions are sampled log-uniformly over the
    type's ``L`` range and the coefficients uniformly over their ranges, then
    perturbed by a small scatter. Identical ``(n_ships, seed)`` always yields the
    same table.
    """
    if n_ships <= 0:
        raise ValueError("n_ships must be positive")
    rng = random.Random(seed)
    # The LCB column has its own stream: drawing it from ``rng`` would shift every
    # later ship's dimensions and coefficients, and the round-3/5 pass rates would
    # then describe a different fleet than the one reported here. Separate streams
    # keep the historical columns bit-for-bit identical while adding the new one.
    lcb_rng = random.Random(seed * 7919 + 12345)
    types = list(SHIP_TYPE_POPULATION)
    ships: list[dict[str, Any]] = []
    for index in range(n_ships):
        ship_type = types[index % len(types)]
        (l_min, l_max, lb, bt, cb_min, cb_max, fn_min, fn_max) = SHIP_TYPE_POPULATION[ship_type]
        length = math.exp(rng.uniform(math.log(l_min), math.log(l_max)))
        beam = length / (lb * rng.uniform(0.92, 1.08))
        draft = beam / (bt * rng.uniform(0.92, 1.08))
        cb = rng.uniform(cb_min, cb_max)
        cm = _midship_from_block(cb) + rng.uniform(-0.010, 0.010)
        cm = min(max(cm, 0.80), 0.995)
        cwp = _waterplane_from_block(cb) + rng.uniform(-0.010, 0.010)
        cp = cb / cm
        cwp_adjusted = False
        if cp > 0.97 * cwp:  # keep a feasibility margin for the calibrated family
            cwp = cp / 0.97
            cwp_adjusted = True
        cwp = min(cwp, 0.95)
        ships.append({
            "ship_id": index + 1,
            "ship_type": ship_type,
            "L": round(length, 2),
            "B": round(beam, 2),
            "T": round(draft, 2),
            "Cb": round(cb, 4),
            "Cm": round(cm, 4),
            "Cwp": round(cwp, 4),
            "design_froude": round(rng.uniform(fn_min, fn_max), 4),
            "cwp_adjusted": cwp_adjusted,
            # Separate stream (see ``lcb_rng`` above): the historical columns of this
            # table are unchanged by the addition of the LCB column.
            "lcb": round(sample_lcb_fraction(ship_type, lcb_rng), 5),
            "seed": seed,
        })
    return ships


# --------------------------------------------------------------------------- #
# Nominal performance fields
# --------------------------------------------------------------------------- #
def holtrop_mennen_wetted_surface(ship: dict[str, Any]) -> float:
    """Holtrop & Mennen (1982) wetted surface, ``A_BT = 0`` (no bulb).

    ``S = L*(2T + B)*sqrt(Cm)*(0.453 + 0.4425*Cb - 0.2862*Cm
         - 0.003467*B/T + 0.3696*Cwp) + 2.38*A_BT/Cb``
    """
    length, beam, draft = ship["L"], ship["B"], ship["T"]
    cb, cm, cwp = ship["Cb"], ship["Cm"], ship["Cwp"]
    return length * (2.0 * draft + beam) * math.sqrt(cm) * (
        0.453 + 0.4425 * cb - 0.2862 * cm - 0.003467 * beam / draft + 0.3696 * cwp
    )


def holtrop_mennen_form_factor(ship: dict[str, Any], *, c_stern: float = 0.0) -> float:
    """Holtrop & Mennen (1982) form factor ``1 + k1`` (no bulbous bow, no stern shape).

    ``1 + k1 = 0.93 + 0.487118*c14*(B/L)^1.06806*(T/L)^0.46106
                     *(L/nabla^(1/3))^0.121563*(L**3/nabla)^0.36486
                     *(1-Cp)^(-0.604247)``
    with ``c14 = 1 + 0.011*c_stern``.
    """
    length, beam, draft = ship["L"], ship["B"], ship["T"]
    cb, cm = ship["Cb"], ship["Cm"]
    cp = cb / cm
    nabla = cb * length * beam * draft
    c14 = 1.0 + 0.011 * c_stern
    return 0.93 + 0.487118 * c14 * (beam / length) ** 1.06806 * (draft / length) ** 0.46106 * (
        length / nabla ** (1.0 / 3.0)
    ) ** 0.121563 * (length ** 3.0 / nabla) ** 0.36486 * (1.0 - cp) ** (-0.604247)


def ittc57_friction_coefficient(reynolds: float) -> float:
    """ITTC-57 line ``Cf = 0.075 / (log10 Re - 2)**2``."""
    if reynolds <= 0.0:
        raise ValueError("Reynolds number must be positive")
    return 0.075 / (math.log10(reynolds) - 2.0) ** 2


def nominal_performance_fields(ship: dict[str, Any]) -> dict[str, Any]:
    """Nominal performance block implied by the parameter vector alone.

    Nothing here sees the rebuilt mesh: the block is the *target* the geometry must
    be able to reproduce. Resistance is viscous only (``Rt = (1+k1)*Rf``) -- no
    wave, appendage or air resistance -- so it is a shape-similarity proxy, not a
    power prediction.
    """
    length, beam, draft = ship["L"], ship["B"], ship["T"]
    cb, cm, cwp = ship["Cb"], ship["Cm"], ship["Cwp"]
    # Benchmark hulls carry ``design_froude: null`` when the manifest has no service
    # Froude number; fall back to the mid-range 0.20 rather than raising, so the
    # geometry-side fields (volume, wetted surface) stay reportable for them.
    froude = ship.get("design_froude") or 0.20
    volume = cb * length * beam * draft
    speed = froude * math.sqrt(GRAVITY * length)      # m/s
    reynolds = speed * length / SEAWATER_VISCOSITY
    cf = ittc57_friction_coefficient(reynolds)
    form_factor = holtrop_mennen_form_factor(ship)
    wetted = holtrop_mennen_wetted_surface(ship)
    resistance = 0.5 * SEAWATER_DENSITY * speed**2 * wetted * cf * form_factor
    return {
        "design_froude": froude,
        "design_speed_mps": speed,
        "design_speed_knots": speed * 1.9438444924406046,
        "displacement_volume_m3": volume,
        "displacement_t": SEAWATER_DENSITY * volume / 1000.0,
        "reynolds": reynolds,
        "friction_coefficient": cf,
        "form_factor_1pk": form_factor,
        "wetted_surface_m2": wetted,
        "viscous_resistance_kN": resistance / 1000.0,
        "prismatic_coefficient": cb / cm,
    }


# --------------------------------------------------------------------------- #
# Correlation sanity: rebuilt mesh vs the nominal fields
# --------------------------------------------------------------------------- #
def _pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def correlation_sanity(
    ships: Sequence[dict[str, Any]],
    results: Sequence[Any],
    *,
    min_wetted_relative_tolerance: float = 0.15,
) -> dict[str, Any]:
    """Compare rebuilt-surface quantities against the nominal performance fields.

    For every ship whose rebuild produced a surface, the mesh's own wetted area and
    displacement are compared with the values the parameter vector implies
    (Holtrop & Mennen wetted surface, ``Cb*L*B*T`` displacement). The point of the
    check is not that the two agree to the last percent -- the generated family has
    no bulb, transom or parallel body -- but that they agree *systematically*:
    a correlation slope near 1 and an intercept near 0 across 500 hulls, with a
    bias small enough that a resistance estimate taken off the mesh is not
    misleading.

    ``min_wetted_relative_tolerance`` splits the fleet into "wetted surface inside
    the tolerance" for accounting; it is a reporting threshold, not a gate.
    """
    nominal_s: list[float] = []
    mesh_s: list[float] = []
    nominal_v: list[float] = []
    mesh_v: list[float] = []
    wetted_errors: list[float] = []
    displacement_errors: list[float] = []
    per_type: dict[str, list[float]] = {}
    worst: list[dict[str, Any]] = []

    for ship, result in zip(ships, results):
        if result.mesh is None or result.wetted_surface_area is None:
            continue
        nominal = nominal_performance_fields(ship)
        s_nom = nominal["wetted_surface_m2"]
        v_nom = nominal["displacement_volume_m3"]
        s_mesh = result.wetted_surface_area
        v_mesh = result.volume_mesh
        if not (s_nom > 0.0 and v_nom > 0.0) or v_mesh is None:
            continue
        nominal_s.append(s_nom)
        mesh_s.append(s_mesh)
        nominal_v.append(v_nom)
        mesh_v.append(v_mesh)
        s_err = (s_mesh - s_nom) / s_nom
        v_err = (v_mesh - v_nom) / v_nom
        wetted_errors.append(s_err)
        displacement_errors.append(v_err)
        per_type.setdefault(ship["ship_type"], []).append(s_err)
        worst.append({
            "ship_id": ship["ship_id"], "ship_type": ship["ship_type"],
            "Cb": ship["Cb"], "wetted_nominal_m2": s_nom, "wetted_mesh_m2": s_mesh,
            "wetted_relative_error": s_err,
        })

    if not nominal_s:
        return {"n_pairs": 0, "note": "no rebuilt surfaces to correlate"}

    def _stats(values: list[float]) -> dict[str, float]:
        ordered = sorted(values)
        n = len(ordered)
        return {
            "mean": sum(values) / n,
            "median": ordered[n // 2],
            "min": ordered[0],
            "max": ordered[-1],
            "mean_abs": sum(abs(v) for v in values) / n,
            "within_tolerance_fraction": sum(
                1 for v in values if abs(v) < min_wetted_relative_tolerance
            ) / n,
        }

    r_wetted = _pearson(nominal_s, mesh_s)
    r_disp = _pearson(nominal_v, mesh_v)
    slope = sum(s * m for s, m in zip(nominal_s, mesh_s)) / sum(n * n for n in nominal_s)
    worst.sort(key=lambda entry: abs(entry["wetted_relative_error"]), reverse=True)
    return {
        "n_pairs": len(nominal_s),
        "wetted_surface": {
            "reference": "Holtrop & Mennen (1982) wetted-surface regression, A_BT=0",
            "pearson_r": r_wetted,
            "r_squared": r_wetted**2 if r_wetted == r_wetted else float("nan"),
            "slope_through_origin": slope,
            "relative_error": _stats(wetted_errors),
            "by_ship_type_mean_relative_error": {
                key: sum(vals) / len(vals) for key, vals in sorted(per_type.items())
            },
            "tolerance_used": min_wetted_relative_tolerance,
        },
        "displacement": {
            "reference": "Cb*L*B*T from the parameter vector",
            "pearson_r": r_disp,
            "relative_error": _stats(displacement_errors),
        },
        "worst_wetted_offsets": worst[:10],
    }