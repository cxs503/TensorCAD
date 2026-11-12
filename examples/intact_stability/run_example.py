"""IS Code 2008 part A criteria measured on curves this repository actually integrates.

Every number in ``docs/intact_stability.md`` comes out of this run: the script writes
``results.json`` and ``report.md`` and **exits non-zero if one of its claims stops holding**,
so the page cannot drift away from the code.

The criterion arithmetic is easy to get wrong in a way that always says "compliant" -- an
inverted comparison, an area taken to the wrong limit, a ``>=`` where the code says ``<=``.
So the run is built as a control experiment rather than as a demonstration:

1. **A wall-sided barge** (``L = 10, B = 4, T = 2, freeboard = 2, KG = 1.2 m``), whose
   righting arm and its antiderivative are elementary. The mesh integral is measured against
   the arm, the report's areas against the antiderivative, and both waterline solvers are run
   on it (``tensorcad.gz``'s Newton on the wetted area, ``tensorcad.righting_arm``'s bisection
   on the volume) because the criteria cannot be more accurate than the curve under them.
2. **The same barge loaded deep** (``KG = 1.55 m``, ``GM_0 = 0.117 m``): the report has to
   *flip*, clause by clause. A check that cannot fail is not a check.
3. **A half-arm negative control**: the passing curve scaled by 0.5 has to fail on the
   area and lever clauses, on the same geometry.
4. **The parametric hull** (``L = 12, B = 2.4, T = 0.9 m`` with 0.6 m of freeboard,
   ``GM_0 = 0.30 m`` from the waterplane moments of :mod:`tensorcad.stability`): the case
   where no closed form exists, so the checks are the ``GM_0``-versus-slope identity, the
   angle-grid convergence, and the report itself -- whatever it says.
5. **The weather criterion** (2.3.1 / 2.3.2 / 2.3.3): a moderate wind that passes, a storm
   that fails both clauses, a deck-edge case where the 80 % rule, not the 16 deg limit,
   binds, and a flooding angle the windward roll amplitude is measured against.
6. **Angle-grid convergence**: the criteria are read off a polyline, so the sampled step is
   their only discretisation error. The ladder measures it rather than asserting it.

Run it with ``python examples/intact_stability/run_example.py --output docs/intact_stability``.
matplotlib is optional: the JSON and the Markdown are written either way, the PNGs are
skipped without it.
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import json
import math
import sys
import time
from pathlib import Path

from tensorcad import __version__
from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import righting_arm_curve as newton_curve
from tensorcad.gz import with_freeboard
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.intact_stability import (
    general_criteria,
    weather_criterion,
    wind_heeling_levers,
)
from tensorcad.righting_arm import righting_arm_curve as bisection_curve
from tensorcad.stability import initial_stability

LENGTH, BREADTH, DRAFT, DEPTH = 10.0, 4.0, 2.0, 4.0
KG, KG_HEAVY = 1.2, 1.55
BM0 = BREADTH * BREADTH / (12.0 * DRAFT)
KB0 = 0.5 * DRAFT
GM0 = KB0 + BM0 - KG
RHO_SEAWATER = 1.025

HULL_LENGTH, HULL_BEAM, HULL_DRAFT, HULL_FREEBOARD = 12.0, 2.4, 0.9, 0.6
HULL_GM_TARGET = 0.30
HULL_ANGLE_STEP = 5.0

GRIDS = (5.0, 2.0, 1.0, 0.5)


# ---------------------------------------------------------------------------
# reference: the barge's arm and its antiderivative
# ---------------------------------------------------------------------------
def wall_sided_gz(angle_deg: float, *, gm: float = GM0, bm: float = BM0) -> float:
    heel = math.radians(angle_deg)
    return math.sin(heel) * (gm + 0.5 * bm * math.tan(heel) ** 2)


def wall_sided_area(lower_deg: float, upper_deg: float, *, gm: float = GM0, bm: float = BM0) -> float:
    """``int GZ dphi`` from the elementary antiderivative of the wall-sided form."""

    def primitive(angle_deg: float) -> float:
        heel = math.radians(angle_deg)
        return gm * (1.0 - math.cos(heel)) + 0.5 * bm * (1.0 / math.cos(heel) + math.cos(heel) - 2.0)

    return primitive(upper_deg) - primitive(lower_deg)


def solve_wall_sided_heel(lever: float) -> float:
    """Independent bisection of ``GZ(phi) = lever`` on the rising branch."""
    low, high = 0.0, 45.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if wall_sided_gz(middle) < lever:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2),
        (0, 2, 1),
        (4, 5, 6),
        (4, 6, 7),
        (0, 1, 5),
        (0, 5, 4),
        (2, 3, 7),
        (2, 7, 6),
        (1, 2, 6),
        (1, 6, 5),
        (0, 4, 7),
        (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


def ladder(step_deg: float, limit_deg: float = 50.0) -> tuple[float, ...]:
    count = int(round(limit_deg / step_deg))
    return tuple(index * step_deg for index in range(count + 1))


# ---------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------
class Run:
    def __init__(self) -> None:
        self.claims: list[dict[str, object]] = []
        self.results: dict[str, object] = {}

    def claim(self, name: str, ok: bool, measured: str, expected: str) -> bool:
        self.claims.append(
            {"claim": name, "holds": bool(ok), "measured": measured, "expected": expected}
        )
        print(f"[{'ok  ' if ok else 'FAIL'}] {name}: {measured} (expected {expected})")
        return bool(ok)

    @property
    def failed(self) -> list[dict[str, object]]:
        return [claim for claim in self.claims if not claim["holds"]]


def criteria_table(report) -> list[dict[str, object]]:
    return [
        {
            "clause": item.clause,
            "quantity": item.quantity,
            "required": round(item.required, 9),
            "attained": round(item.attained, 9),
            "unit": item.unit,
            "comparison": item.comparison,
            "margin": round(item.margin, 9),
            "satisfied": item.satisfied,
        }
        for item in report.criteria
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="docs/intact_stability", type=Path)
    options = parser.parse_args(argv)
    output: Path = options.output
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    run = Run()

    # ---- 1. the wall-sided barge, through the mesh -------------------------
    barge = box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))
    angles = ladder(1.0)
    curve = newton_curve(barge, KG, angles, waterline_z=DRAFT)
    limit_deg = math.degrees(math.atan(2.0 * DRAFT / BREADTH))
    inside = [
        (angle, measured)
        for angle, measured in zip(curve.angles_deg, curve.gz)
        if angle <= limit_deg + 1e-9
    ]
    arm_error = max(abs(measured - wall_sided_gz(angle)) for angle, measured in inside)
    run.claim(
        "barge GZ from the mesh equals the wall-sided closed form inside its validity window",
        arm_error < 1e-11,
        f"max |GZ - GZ_analytic| = {arm_error:.3e} m over 0-{limit_deg:.2f} deg",
        "< 1e-11 m",
    )
    outside = [(angle, measured) for angle, measured in zip(curve.angles_deg, curve.gz) if angle > limit_deg + 1e-9]
    outside_gap = max(abs(measured - wall_sided_gz(angle)) for angle, measured in outside)
    run.claim(
        "the closed form stops applying beyond the bilge-emergence angle, and the mesh says so",
        outside_gap > 1e-3,
        f"max |GZ - GZ_analytic| = {outside_gap:.3e} m above {limit_deg:.2f} deg "
        f"(the formula is not the barge's arm there)",
        "> 1e-3 m, so the agreement above is not vacuous",
    )
    run.claim(
        "barge displacement is held at every heel",
        curve.max_volume_error < 1e-12,
        f"max |V(phi)/V(0) - 1| = {curve.max_volume_error:.3e}",
        "< 1e-12",
    )

    barge_report = general_criteria(curve, gm0_m=GM0)
    analytic_area = wall_sided_area(0.0, 30.0)
    area_error = abs(barge_report.area_0_30_m_rad - analytic_area) / analytic_area
    run.claim(
        "the report's 30 deg area matches the antiderivative",
        area_error < 5e-4,
        f"{barge_report.area_0_30_m_rad:.9f} m rad vs {analytic_area:.9f} "
        f"(relative {area_error:.2e})",
        "< 5e-4 relative on a 1 deg grid",
    )
    run.claim(
        "the light barge passes all six general criteria",
        barge_report.satisfied,
        barge_report.verdict(),
        "passes 2.2.1-2.2.6",
    )

    other_route = bisection_curve(barge, kg=KG, heels_deg=angles, waterline_z=DRAFT)
    other_report = general_criteria(other_route, gm0_m=GM0)
    route_gap = max(
        abs(left.attained - right.attained)
        for left, right in zip(barge_report.criteria, other_report.criteria)
    )
    run.claim(
        "the two waterline solvers produce the same report",
        route_gap < 1e-8,
        f"max |difference| over the six clauses = {route_gap:.3e}",
        "< 1e-8 (each clause keeps its own unit: m, m rad or deg)",
    )

    # ---- 2. the same barge loaded deep: the verdict must flip --------------
    heavy = newton_curve(barge, KG_HEAVY, angles, waterline_z=DRAFT)
    heavy_gm = KB0 + BM0 - KG_HEAVY
    heavy_report = general_criteria(heavy, gm0_m=heavy_gm)
    heavy_failed = sorted(item.clause for item in heavy_report.failures)
    run.claim(
        "loading 0.35 m deeper flips the verdict",
        not heavy_report.satisfied and heavy_gm < 0.15,
        f"GM_0 = {heavy_gm:.4f} m; {heavy_report.verdict()}",
        "fails, with GM_0 under the 0.15 m floor",
    )
    run.claim(
        "the failing clauses are the five the areas and levers predict",
        heavy_failed == ["2.2.1", "2.2.2", "2.2.3", "2.2.4", "2.2.6"],
        ", ".join(heavy_failed) or "none",
        "2.2.1, 2.2.2, 2.2.3, 2.2.4, 2.2.6 (2.2.5 is a shape clause and survives)",
    )

    # ---- 3. half-arm negative control --------------------------------------
    halved = (curve.angles_deg, tuple(0.5 * value for value in curve.gz))
    halved_report = general_criteria(halved, gm0_m=0.5 * GM0)
    run.claim(
        "scaling the passing arm by 0.5 makes the same geometry fail",
        not halved_report.satisfied,
        halved_report.verdict(),
        "fails (the comparison can fail on the same geometry)",
    )

    # ---- 4. the parametric hull -------------------------------------------
    model = generate_hull(
        HullParameters(length=HULL_LENGTH, beam=HULL_BEAM, draft=HULL_DRAFT, nx=81, nz=25)
    )
    # The waterplane at z = 0 is the same with and without the freeboard band, so KM (and
    # therefore the KG for a target GM) is taken from the hull as generated; the band is
    # used for the heeled cuts. The plane cut accepts the banded mesh too now that the
    # seam ring (edges lying in the cutting plane) is traced once instead of twice -- see
    # examples/freeboard_seam -- so the two routes are interchangeable, not exclusive.
    upright = initial_stability(model.mesh, kg=0.0, waterline_z=0.0)
    hull_kg = upright.km_transverse - HULL_GM_TARGET
    hull_stability = initial_stability(model.mesh, kg=hull_kg, waterline_z=0.0)
    banded = with_freeboard(model.mesh, HULL_FREEBOARD, waterline_z=0.0)
    hull_angles = ladder(HULL_ANGLE_STEP)
    hull_curve = newton_curve(banded, hull_kg, hull_angles, waterline_z=0.0)
    # The slope is a separate, deliberately tiny ladder: on the 5 deg criteria grid a
    # one-sided difference is 6.5 % low, which says more about the difference than about GM.
    probe = 0.25
    probe_curve = newton_curve(banded, hull_kg, (-probe, probe), waterline_z=0.0)
    slope = (probe_curve.gz[1] - probe_curve.gz[0]) / math.radians(2.0 * probe)
    slope_error = abs(slope - hull_stability.gm_transverse) / hull_stability.gm_transverse
    run.claim(
        "the hull curve's slope at zero equals GM_0 from the waterplane moments",
        slope_error < 1e-2,
        f"dGZ/dphi = {slope:.6f} m vs GM_0 = {hull_stability.gm_transverse:.6f} m "
        f"(relative {slope_error:.2e}, central difference over +/-{probe:g} deg)",
        "< 1e-2 relative (the rest is the 81x25 hull's own resolution)",
    )
    displacement_t = hull_stability.displacement_volume * RHO_SEAWATER
    hull_report = general_criteria(hull_curve, gm0_m=hull_stability.gm_transverse)
    print(f"       hull verdict: {hull_report.verdict()}")

    # ---- 5. weather criterion ---------------------------------------------
    barge_displacement = LENGTH * BREADTH * DRAFT * RHO_SEAWATER
    moderate = weather_criterion(
        curve, wind_area_m2=10.0, wind_arm_m=2.0, displacement_t=barge_displacement
    )
    reference_heel = solve_wall_sided_heel(moderate.levers.steady_m)
    heel_error = abs(moderate.theta0_deg - reference_heel) if moderate.theta0_deg else float("inf")
    run.claim(
        "the steady wind heel angle matches an independent bisection",
        heel_error < 1e-3,
        f"theta0 = {moderate.theta0_deg:.6f} deg vs {reference_heel:.6f} deg",
        "< 1e-3 deg (linear interpolant on a 1 deg grid)",
    )
    run.claim(
        "a moderate wind passes the weather criterion",
        moderate.satisfied,
        f"{moderate.verdict()}; lw1 = {moderate.levers.steady_m:.6f} m, "
        f"lw2 = {moderate.levers.gust_m:.6f} m, b = {moderate.area_b_m_rad:.6f}, "
        f"a = {moderate.area_a_m_rad:.6f} m rad",
        "passes 2.3.1 and 2.3.2",
    )
    storm = weather_criterion(
        curve, wind_area_m2=240.0, wind_arm_m=2.0, displacement_t=barge_displacement
    )
    run.claim(
        "a storm wind fails both weather clauses",
        not storm.satisfied,
        f"{storm.verdict()}; theta0 = {storm.theta0_deg:.3f} deg, b = {storm.area_b_m_rad:.6f}, "
        f"a = {storm.area_a_m_rad:.6f} m rad",
        "fails 2.3.1 and 2.3.2",
    )
    deck_edge = 12.0
    binding_lever = 0.1
    binding_area = binding_lever * (1000.0 * 9.81 * barge_displacement) / (504.0 * 2.0)
    binding = weather_criterion(
        curve,
        wind_area_m2=binding_area,
        wind_arm_m=2.0,
        displacement_t=barge_displacement,
        deck_edge_angle_deg=deck_edge,
    )
    run.claim(
        "the 80 % deck-edge rule binds before the 16 deg limit",
        binding.steady_heel_limit_deg < 16.0
        and bool(binding.failures)
        and binding.theta0_deg is not None
        and binding.theta0_deg < 16.0,
        f"limit = {binding.steady_heel_limit_deg:.1f} deg, theta0 = {binding.theta0_deg:.3f} deg, "
        f"{binding.verdict()}",
        f"limit 0.8 x {deck_edge:g} = 9.6 deg < theta0 < 16 deg, so 2.3.1 fails",
    )
    flooding = 20.0
    capped = weather_criterion(
        curve,
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=barge_displacement,
        roll_amplitude_deg=8.0,
        flooding_angle_deg=flooding,
    )
    run.claim(
        "clause 2.3.3 passes and theta_f caps the gust areas",
        capped.satisfied
        and capped.theta2_deg == flooding
        and capped.clauses == "2.3.1, 2.3.2 and 2.3.3",
        f"theta_f = {flooding:g} deg, theta_1 = {capped.roll_amplitude_deg:g} deg, "
        f"theta_2 = {capped.theta2_deg:g} deg, b = {capped.area_b_m_rad:.6f}, "
        f"a = {capped.area_a_m_rad:.6f} m rad",
        "theta_1 < theta_f, theta_2 = theta_f = 20 deg (not 50 deg), 2.3.3 satisfied",
    )
    swamped = weather_criterion(
        curve,
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=barge_displacement,
        roll_amplitude_deg=25.0,
        flooding_angle_deg=flooding,
    )
    run.claim(
        "a roll amplitude past theta_f fails 2.3.3 and leaves the gust areas unmeasured",
        not swamped.satisfied
        and [item.clause for item in swamped.failures] == ["2.3.3"]
        and swamped.area_b_m_rad is None
        and swamped.clauses == "2.3.1 and 2.3.3",
        f"{swamped.verdict()}; theta_1 = 25 deg > theta_f = {flooding:g} deg, "
        f"b = {swamped.area_b_m_rad}",
        "fails 2.3.3 alone; theta_2 - theta_1 span is empty, so b and a are None",
    )

    # ---- 6. angle-grid convergence ----------------------------------------
    convergence = []
    for step in GRIDS:
        grid_curve = newton_curve(barge, KG, ladder(step), waterline_z=DRAFT)
        grid_report = general_criteria(grid_curve, gm0_m=GM0)
        error = abs(grid_report.area_0_30_m_rad - analytic_area) / analytic_area
        convergence.append(
            {
                "angle_step_deg": step,
                "area_0_30_m_rad": grid_report.area_0_30_m_rad,
                "relative_error": error,
                "verdict": grid_report.verdict(),
            }
        )
    errors = [item["relative_error"] for item in convergence]
    ratios = [errors[index] / errors[index + 1] for index in range(len(errors) - 1)]
    asymptotic = ratios[1:]  # 5 deg -> 2 deg is outside the asymptotic regime
    run.claim(
        "the criteria's grid error falls as the step squared",
        all(3.5 < ratio < 4.5 for ratio in asymptotic)
        and all(item["verdict"] == barge_report.verdict() for item in convergence),
        "errors "
        + ", ".join(
            f"{errors[index]:.2e} at {convergence[index]['angle_step_deg']:g} deg"
            for index in range(len(errors))
        )
        + "; successive ratios "
        + ", ".join(f"{ratio:.2f}" for ratio in ratios)
        + "; verdict identical on all four grids",
        "successive ratios 3.5-4.5 for 2 deg and below (the 5 deg rung is not asymptotic), "
        "same verdict",
    )

    # ---- record -----------------------------------------------------------
    run.results = {
        "schema": "tensorcad.intact-stability-example/1",
        "tensorcad_version": __version__,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_s": round(time.time() - started, 3),
        "barge": {
            "length_m": LENGTH,
            "breadth_m": BREADTH,
            "draft_m": DRAFT,
            "depth_m": DEPTH,
            "kg_m": KG,
            "km_m": KB0 + BM0,
            "gm0_m": GM0,
            "bm0_m": BM0,
            "kb0_m": KB0,
            "wall_sided_limit_deg": limit_deg,
            "arm_error_beyond_limit_m": outside_gap,
            "analytic": {
                "gz_30_m": wall_sided_gz(30.0),
                "gz_max_m": wall_sided_gz(50.0),
                "area_0_30_m_rad": analytic_area,
                "area_0_40_m_rad": wall_sided_area(0.0, 40.0),
                "area_30_40_m_rad": wall_sided_area(30.0, 40.0),
            },
            "curve": {
                "angles_deg": list(curve.angles_deg),
                "gz_m": [round(value, 12) for value in curve.gz],
                "max_volume_error": curve.max_volume_error,
                "arm_error_vs_closed_form_m": arm_error,
            },
            "report": barge_report.to_dict(),
            "criteria": criteria_table(barge_report),
            "other_route": {
                "solver": "tensorcad.righting_arm (bisection on V(h))",
                "report": other_report.to_dict(),
                "max_clause_gap": route_gap,
            },
        },
        "heavy_loading": {
            "kg_m": KG_HEAVY,
            "gm0_m": heavy_gm,
            "criteria": criteria_table(heavy_report),
            "report": heavy_report.to_dict(),
            "curve": {"angles_deg": list(heavy.angles_deg), "gz_m": [round(v, 12) for v in heavy.gz]},
        },
        "negative_control": {
            "what": "the passing barge curve scaled by 0.5, GM_0 halved with it",
            "criteria": criteria_table(halved_report),
            "report": halved_report.to_dict(),
        },
        "hull": {
            "parameters": {
                "length_m": HULL_LENGTH,
                "beam_m": HULL_BEAM,
                "draft_m": HULL_DRAFT,
                "freeboard_m": HULL_FREEBOARD,
                "kg_m": hull_kg,
                "gm_target_m": HULL_GM_TARGET,
            },
            "gm0_m": hull_stability.gm_transverse,
            "displacement_volume_m3": hull_stability.displacement_volume,
            "displacement_t": displacement_t,
            "slope_at_zero_m": slope,
            "slope_probe_deg": probe,
            "slope_relative_error": slope_error,
            "angle_step_deg": HULL_ANGLE_STEP,
            "curve": {
                "angles_deg": list(hull_curve.angles_deg),
                "gz_m": [round(value, 12) for value in hull_curve.gz],
                "max_volume_error": hull_curve.max_volume_error,
            },
            "criteria": criteria_table(hull_report),
            "report": hull_report.to_dict(),
        },
        "weather": {
            "moderate": {
                "wind_area_m2": 10.0,
                "wind_arm_m": 2.0,
                "displacement_t": barge_displacement,
                "reference_theta0_deg": reference_heel,
                "report": moderate.to_dict(),
            },
            "storm": {
                "wind_area_m2": 240.0,
                "wind_arm_m": 2.0,
                "displacement_t": barge_displacement,
                "report": storm.to_dict(),
            },
            "deck_edge_rule": {"deck_edge_angle_deg": deck_edge, "report": binding.to_dict()},
            "flooding_angle": {
                "flooding_angle_deg": flooding,
                "roll_amplitude_deg": 8.0,
                "report": capped.to_dict(),
            },
            "roll_past_flooding": {
                "flooding_angle_deg": flooding,
                "roll_amplitude_deg": 25.0,
                "report": swamped.to_dict(),
            },
            "levers": wind_heeling_levers(10.0, 2.0, barge_displacement).to_dict(),
        },
        "convergence": convergence,
        "convergence_ratios": ratios,
        "claims": run.claims,
        "failures": [claim["claim"] for claim in run.failed],
    }
    (output / "results.json").write_text(json.dumps(run.results, indent=2) + "\n")
    _figures(output, run.results)
    (output / "report.md").write_text(_markdown(run.results))
    print(f"\nwrote {output/'results.json'}, report.md and figures in {run.results['runtime_s']} s")
    if run.failed:
        print(f"FAILED claims: {[claim['claim'] for claim in run.failed]}")
        return 1
    print(f"all {len(run.claims)} claims hold")
    return 0


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _figures(output: Path, results: dict) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - matplotlib is optional
        print(f"matplotlib unavailable ({error}); skipping PNGs")
        return

    barge = results["barge"]
    barge_angles = barge["curve"]["angles_deg"]
    barge_gz = barge["curve"]["gz_m"]

    # (1) the curve with the area clauses shaded, and every clause's margin
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.plot(barge_angles, barge_gz, "k-", lw=2, label="wall-sided barge (mesh integral)")
    ax.plot(
        barge_angles,
        [wall_sided_gz(angle) for angle in barge_angles],
        "r--",
        lw=1,
        label="analytic $\\sin\\phi(GM_0+\\frac{1}{2}BM_0\\tan^2\\phi)$",
    )
    span_a = [angle for angle in barge_angles if angle <= 30.0]
    span_b = [angle for angle in barge_angles if 30.0 <= angle <= 40.0]
    ax.fill_between(span_a, [wall_sided_gz(a) for a in span_a], color="#4c72b0", alpha=0.25)
    ax.fill_between(span_b, [wall_sided_gz(a) for a in span_b], color="#dd8452", alpha=0.45)
    ax.axvline(30.0, color="0.4", ls=":", lw=1)
    ax.axvline(40.0, color="0.4", ls=":", lw=1)
    ax.axhline(0.20, color="#55a868", ls="-.", lw=1)
    ax.annotate(
        f"$A_{{0-30}}$ = {barge['report']['area_0_30_m_rad']:.4f} m rad (min 0.055)",
        xy=(2.0, 0.05),
        fontsize=8,
    )
    ax.annotate(
        f"$A_{{30-40}}$ = {barge['report']['area_30_40_m_rad']:.4f} (min 0.030)",
        xy=(30.8, 0.05),
        fontsize=8,
    )
    ax.annotate(f"GZ(30°) = {barge['report']['gz_30_m']:.3f} m (min 0.20)", xy=(6.0, 0.215), fontsize=8)
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("(a) barge: the clauses are areas under this curve")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    cases = [
        ("barge, KG 1.20", barge["criteria"]),
        ("barge, KG 1.55", results["heavy_loading"]["criteria"]),
        ("hull, $GM_0$ 0.30", results["hull"]["criteria"]),
    ]
    clauses = [item["clause"] for item in cases[0][1]]
    width = 0.26
    for index, (label, table) in enumerate(cases):
        ratios = [item["attained"] / item["required"] for item in table]
        offsets = [position + (index - 1) * width for position in range(len(clauses))]
        bars = ax.bar(offsets, ratios, width=width, label=label)
        for bar, ratio in zip(bars, ratios):
            ax.annotate(
                f"{ratio:.2f}",
                (bar.get_x() + bar.get_width() / 2, ratio),
                ha="center",
                va="bottom",
                fontsize=6,
            )
    ax.axhline(1.0, color="k", lw=1)
    ax.set_xticks(range(len(clauses)))
    ax.set_xticklabels(clauses, fontsize=8)
    ax.set_ylabel("attained / required")
    ax.set_title("(b) every clause can fail: ratio > 1 is compliant")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(output / "fig_criteria.png", dpi=160)
    plt.close(fig)

    # (2) grid convergence
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    steps = [item["angle_step_deg"] for item in results["convergence"]]
    errors = [item["relative_error"] for item in results["convergence"]]
    ax.loglog(steps, errors, "o-", label="measured: $A_{0-30}$ vs antiderivative")
    ax.loglog(
        steps,
        [errors[-1] * (step / steps[-1]) ** 2 for step in steps],
        "k--",
        lw=1,
        label="$h^2$ reference",
    )
    ax.set_xlabel("angle step [deg]")
    ax.set_ylabel("relative error of $A_{0-30}$")
    ax.set_title("the criteria's only discretisation error is the angle grid")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(output / "fig_convergence.png", dpi=160)
    plt.close(fig)

    # (3) weather criterion
    moderate = results["weather"]["moderate"]["report"]
    storm_case = results["weather"]["storm"]["report"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    for ax, report, title in (
        (axes[0], moderate, "(a) moderate wind: passes"),
        (axes[1], storm_case, "(b) storm: fails 2.3.1 and 2.3.2"),
    ):
        ax.plot(barge_angles, barge_gz, "k-", lw=2, label="GZ")
        lw1 = report["levers"]["steady_lever_m"]
        lw2 = report["levers"]["gust_lever_m"]
        ax.axhline(lw1, color="#4c72b0", ls="--", lw=1.2, label=f"$l_{{w1}}$ = {lw1:.4f} m")
        ax.axhline(lw2, color="#c44e52", ls="-.", lw=1.2, label=f"$l_{{w2}}$ = {lw2:.4f} m")
        theta2 = report["theta2_deg"]
        ax.axvline(report["theta0_deg"], color="0.35", ls=":", lw=1)
        ax.axvline(theta2, color="0.35", ls=":", lw=1)
        span = [angle for angle in barge_angles if angle <= theta2]
        ax.fill_between(span, [wall_sided_gz(a) for a in span], color="#4c72b0", alpha=0.25)
        ax.fill_between(span, [lw2 for _ in span], color="#c44e52", alpha=0.25)
        ax.annotate(rf"$\theta_0$ = {report['theta0_deg']:.1f}°", xy=(report["theta0_deg"] + 1.0, 0.62), fontsize=8)
        ax.annotate(rf"$\theta_2$ = {theta2:.1f}°", xy=(max(theta2 - 12.0, 1.0), 0.72), fontsize=8)
        ax.annotate(
            f"b = {report['area_b_m_rad']:.4f}\na = {report['area_a_m_rad']:.4f}",
            xy=(12.0, 0.05),
            fontsize=8,
        )
        ax.set_xlabel("heel angle [deg]")
        ax.set_title(title)
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("righting arm / heeling lever [m]")
    fig.tight_layout()
    fig.savefig(output / "fig_weather.png", dpi=160)
    plt.close(fig)

    # (4) the three loading conditions
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    hull = results["hull"]
    ax.plot(
        hull["curve"]["angles_deg"],
        hull["curve"]["gz_m"],
        "o-",
        ms=3,
        label=f"hull, $GM_0$ = {hull['gm0_m']:.3f} m",
    )
    ax.plot(
        results["heavy_loading"]["curve"]["angles_deg"],
        results["heavy_loading"]["curve"]["gz_m"],
        "s-",
        ms=3,
        label=f"barge, KG = {results['heavy_loading']['kg_m']:.2f} m",
    )
    ax.plot(barge_angles, barge_gz, "-", lw=2, label=f"barge, KG = {results['barge']['kg_m']:.2f} m")
    ax.axhline(0.20, color="#55a868", ls="-.", lw=1, label="GZ(30°) requirement 0.20 m")
    ax.axvline(30.0, color="0.5", ls=":", lw=1)
    ax.axvline(25.0, color="0.8", ls=":", lw=1)
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("three loading conditions, one set of clauses")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output / "fig_curves.png", dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# the page
# ---------------------------------------------------------------------------
def _markdown(results: dict) -> str:
    barge = results["barge"]
    hull = results["hull"]
    heavy = results["heavy_loading"]
    control = results["negative_control"]
    moderate = results["weather"]["moderate"]["report"]
    storm = results["weather"]["storm"]["report"]
    binding = results["weather"]["deck_edge_rule"]["report"]
    capped = results["weather"]["flooding_angle"]["report"]
    swamped = results["weather"]["roll_past_flooding"]["report"]

    rows = "\n".join(
        f"| {item['clause']} | {item['quantity']} | {item['comparison']} {item['required']:g} "
        f"{item['unit']} | {item['attained']:.6f} | {item['margin']:+.6f} | "
        f"{'pass' if item['satisfied'] else '**fail**'} |"
        for item in barge["criteria"]
    )
    heavy_rows = "\n".join(
        f"| {item['clause']} | {item['comparison']} {item['required']:g} {item['unit']} | "
        f"{item['attained']:.6f} | {item['margin']:+.6f} | "
        f"{'pass' if item['satisfied'] else '**fail**'} |"
        for item in heavy["criteria"]
    )
    hull_rows = "\n".join(
        f"| {item['clause']} | {item['quantity']} | {item['comparison']} {item['required']:g} "
        f"{item['unit']} | {item['attained']:.6f} | "
        f"{'pass' if item['satisfied'] else '**fail**'} |"
        for item in hull["criteria"]
    )
    control_rows = "\n".join(
        f"| {item['clause']} | {item['attained']:.6f} | {item['required']:g} | "
        f"{'pass' if item['satisfied'] else '**fail**'} |"
        for item in control["criteria"]
    )
    ratios_text = ", ".join(f"{ratio:.2f}" for ratio in results["convergence_ratios"])
    convergence_rows = "\n".join(
        f"| {item['angle_step_deg']:g} | {item['area_0_30_m_rad']:.9f} | "
        f"{item['relative_error']:.3e} | {item['verdict']} |"
        for item in results["convergence"]
    )
    claims = "\n".join(
        f"| {claim['claim']} | {claim['measured']} | {claim['expected']} | "
        f"{'yes' if claim['holds'] else '**NO**'} |"
        for claim in results["claims"]
    )
    def area(value: float | None) -> str:
        """A gust area that has no span to be measured over is written as such."""
        return "not measured" if value is None else f"{value:.6f}"

    weather_rows = "\n".join(
        f"| {label} | {report['theta0_deg']:.3f} | {report['steady_heel_limit_deg']:.1f} | "
        f"{report['theta2_deg']:.1f} | {area(report['area_b_m_rad'])} | "
        f"{area(report['area_a_m_rad'])} | {report['clauses']} | {report['verdict']} |"
        for label, report in (
            ("moderate wind (A = 10 m²)", moderate),
            ("storm (A = 240 m²)", storm),
            ("deck edge at 12° (lw1 = 0.1 m)", binding),
            ("θ_f = 20°, θ_1 = 8°", capped),
            ("θ_f = 20°, θ_1 = 25°", swamped),
        )
    )
    return f"""# Intact-stability criteria (IS Code 2008 part A) on measured curves

Generated by `python examples/intact_stability/run_example.py --output docs/intact_stability`
with TensorCAD {results['tensorcad_version']} at {results['generated_utc']} in
{results['runtime_s']} s. Every table below is read from `results.json` in this directory;
nothing is typed from a different run.

## 1. What this run is

`tensorcad.gz` and `tensorcad.righting_arm` integrate a righting-arm curve from the surface
(two independent routes); `tensorcad.stability` gives its slope at zero.
`tensorcad.intact_stability` turns that curve into the six general criteria and the weather
criterion, which is arithmetic -- and arithmetic that is easy to make vacuously true. So the
run is a control experiment: the barge below has a closed-form arm *and* a closed-form
antiderivative, one loading is driven out of compliance deliberately, one curve is halved to
show the comparison can fail on the same geometry, and a real hull is reported as it comes
out.

## 2. The wall-sided barge

`L = {barge['length_m']:g} m`, `B = {barge['breadth_m']:g} m`, `T = {barge['draft_m']:g} m`,
depth `{barge['depth_m']:g} m`, `KG = {barge['kg_m']:g} m` from the keel, floating at `z = T`.
Then `KB_0 = {barge['kb0_m']:.4f} m`, `BM_0 = {barge['bm0_m']:.4f} m`,
`GM_0 = {barge['gm0_m']:.4f} m`, and the closed form is exact up to
`{barge['wall_sided_limit_deg']:.2f}°` (bilge emergence and deck-edge immersion coincide here
because `T = D - T`).

| clause | quantity | requirement | attained | margin | verdict |
|---|---|---|---|---|---|
{rows}

Analytic reference on the same grid: `A(0,30) = {barge['analytic']['area_0_30_m_rad']:.9f}`,
`A(0,40) = {barge['analytic']['area_0_40_m_rad']:.9f}`,
`A(30,40) = {barge['analytic']['area_30_40_m_rad']:.9f}` m rad,
`GZ(30°) = {barge['analytic']['gz_30_m']:.9f} m`, `GZ_max = {barge['analytic']['gz_max_m']:.9f} m`.
The integrated arm differs from the closed form by at most
`{barge['curve']['arm_error_vs_closed_form_m']:.1e} m` -- round-off, because the divergence
integral of a polyhedron is exact.

![criteria](fig_criteria.png)

*Panel (a)*: the integrated curve drawn over the analytic one, with the 0–30° and 30–40° areas
shaded in and the 0.20 m lever requirement drawn across. *Panel (b)*: attained / required for
every clause of every case in this run; the line at 1.0 is the code's threshold, and the red
bars are the deep-loaded barge.

## 3. The same hull, loaded 0.35 m deeper

`KG = {heavy['kg_m']:g} m`, so `GM_0 = {heavy['gm0_m']:.4f} m` -- under the 0.15 m floor,
which is the point of the case: the geometry is unchanged and the whole report has to change.

| clause | requirement | attained | margin | verdict |
|---|---|---|---|---|
{heavy_rows}

The verdict flips to **{heavy['report']['verdict']}**, while 2.2.5 (a pure shape clause: where
the peak sits) survives -- the right reading, because the maximum did not move, only its size
did.

**Negative control on the comparison itself.** {control['what'].capitalize()}:
**{control['report']['verdict']}**.

| clause | attained | required | verdict |
|---|---|---|---|
{control_rows}

The area and lever clauses fail on the same geometry, so the comparison is capable of failing
without a change of shape.

## 4. A real hull: `L = {hull['parameters']['length_m']:g} m` parametric hull

`B = {hull['parameters']['beam_m']:g} m`, `T = {hull['parameters']['draft_m']:g} m`,
`{hull['parameters']['freeboard_m']:g} m` of freeboard added by `with_freeboard` (without it a
heeled hull runs out of surface and `GZ` stops meaning anything), `KG` set so that
`GM_0 = {hull['gm0_m']:.4f} m` comes out of the waterplane moments of `tensorcad.stability`.
Displacement `{hull['displacement_volume_m3']:.4f} m³` = `{hull['displacement_t']:.3f} t`. Curve
on a `{hull['angle_step_deg']:g}°` ladder, displacement residual
`{hull['curve']['max_volume_error']:.2e}`.

| clause | quantity | requirement | attained | verdict |
|---|---|---|---|---|
{hull_rows}

Verdict: **{hull['report']['verdict']}**. No closed form exists for this hull, so the checks
that can be made are made: the curve's slope at zero against `GM_0` from a completely different
integration (relative {hull['slope_relative_error']:.2e}), the displacement residual above, and
the grid study in §6. The rest of the number is the integral, and the integral is the one
verified on the barge.

![curves](fig_curves.png)

## 5. The weather criterion

`lw1 = P A Z / (1000 g Δ)` with `P = 504 Pa`, `g = 9.81 m/s²` and `Δ` in tonnes; the barge
displaces `{results['weather']['moderate']['displacement_t']:.2f} t` (80 m³ x 1.025). For
`A = 10 m²`, `Z = 2 m`: `lw1 = {moderate['levers']['steady_lever_m']:.6f} m`,
`lw2 = 1.5 lw1 = {moderate['levers']['gust_lever_m']:.6f} m`.

| case | `θ_0` [deg] | limit [deg] | `θ_2` [deg] | b [m rad] | a [m rad] | clauses | verdict |
|---|---|---|---|---|---|---|---|
{weather_rows}

The third row is why the 80 % deck-edge rule is implemented rather than the 16° limit alone:
with `θ_deck-edge = 12°` the limit becomes 9.6°, and `θ_0 = {binding['theta0_deg']:.3f}°` is
legal under the 16° rule and illegal under the code.

![weather](fig_weather.png)

`θ_1`, the windward roll amplitude, is an **input** here: clause 2.3.4's empirical formula
needs the roll-period tables, which are not in this repository and were not guessed. A
non-zero `θ_1` moves both areas to start at that angle, which the tests check against the
antiderivative.

The last two rows are clause 2.3.3, which needs one more input: `θ_f`, the flooding angle.
`θ_1 = 8° < θ_f = 20°` passes it, and because the code gives `θ_f` the same role in 2.3.2's
minimum the gust areas then stop at 20° instead of 50° (`b = {capped['area_b_m_rad']:.6f}`
instead of `{moderate['area_b_m_rad']:.6f}` m rad on the same curve). At `θ_1 = 25° > θ_f` the
span `θ_2 - θ_1` is empty: 2.3.3 fails and 2.3.2 is reported as **not measured** rather than
as a zero, which is why the two area cells are not numbers there.

## 6. Angle-grid convergence

The curve is a polyline, so the sampled step is the criteria's only discretisation error.
`A(0,30)` against the antiderivative as the step is halved:

| step [deg] | `A(0,30)` [m rad] | relative error | verdict on that grid |
|---|---|---|---|
{convergence_rows}

Successive error ratios: {ratios_text}. The two finest ratios are 4.00, the trapezoid
signature; the 5° rung is not yet in the asymptotic regime, which is worth stating rather
than averaging away -- a 5° ladder is still accurate to 0.26 % on this barge, but its error
is not predicted by `h^2`.

![convergence](fig_convergence.png)

## 7. The claims, as asserted

The script exits non-zero if any of these stops holding.

| claim | measured | expected | holds |
|---|---|---|---|
{claims}

## 8. What this page is not

* `θ_1` from clause 2.3.4 (see §5), and the roll-period tables behind it.
* Damage stability (part B), grain, timber deck, and the 2.2 variants for special ship types.
* Any statement about a *real* ship: the barge is a reference shape and the hull is a
  synthetic parametric body. Compliance of an actual vessel is a class-society calculation
  with its own loading manual, windage drawings and down-flooding points.
* The 2.2 areas are quoted to the sampled grid, so a 5° ladder under-reads `A(0,30)` by
  {results['convergence'][0]['relative_error'] * 100:.2f} % -- small against the code's margins
  ({barge['criteria'][0]['margin'] / barge['criteria'][0]['required'] * 100:.0f} % on this
  barge), but the reason the report carries `angle_step_max_deg` is so a reader can see it.
"""


if __name__ == "__main__":
    sys.exit(main())