#!/usr/bin/env python3
"""Cross-check the canonical ``GZ`` integrator against an independent station integrator.

:mod:`tensorcad.gz` is the canonical righting-arm route and ``examples/gz`` already holds it
against the two closed forms that exist. What that example *cannot* catch is a defect shared
by a closed form and its own implementation -- so this one puts a **third code path** on the
same meshes: ``examples/righting/oracle.py`` integrates heeled *sections* (Green's theorem on
clipped station polygons, Simpson along the hull, bisection on the tilted waterplane), while
``tensorcad.gz`` integrates the *closed surface* (divergence theorem, Newton on the
waterline). Different error sources, same answer if both are right.

The claims, each one checked and one of them able to fail:

1. **Box barge inside the wall-sided window** -- the oracle reproduces
   ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` to round-off, because it clips a rectangle exactly
   and a prism's station areas are constant along ``x``.
2. **... and stops reproducing it past the window** (negative control on claim 1): the same
   comparison above the deck-edge / bilge limit must show a visible deviation, otherwise
   claim 1 could be passing for the wrong reason.
3. **Deck-edge immersion onset** from the section grid lands within one sampling step of
   ``tensorcad.gz.box_wall_sided_limit_deg``.
4. **Hull: the two independent routes agree** over 0-90 degrees, to a fraction of the arm.
5. **Hull: the two routes float the same displacement** (upright bookkeeping, independent of
   the heel sweep).
6. **Station-grid convergence** -- refining ``station_intervals`` narrows the gap to the
   converged value, and the coarse grid is measurably worse than the reported one.
7. **Negative control on the whole comparison** -- shifting ``gz``'s displacement target
   moves the arm, so agreement in claim 4 is not agreement of two constants.

Writes ``results.json``, ``report.md`` and ``righting_arms.png`` into ``--output`` (default
``docs/righting``) and exits non-zero if any claim stops holding, so CI can run it as a check
rather than as a demo.

Usage::

    python examples/righting/run_example.py [--output docs/righting] [--no-plot]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
for path in (REPO_ROOT, Path(__file__).resolve().parent):
    if str(path) not in sys.path:  # allow running from a checkout without installing
        sys.path.insert(0, str(path))

from oracle import (
    DEFAULT_STATION_INTERVALS,
    righting_arm_curve as oracle_curve,
)
from tensorcad.geometry import Mesh
from tensorcad.gz import (
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    heeled_hydrostatics,
    righting_arm_curve,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull

#: Box barge: ``L=10``, ``B=2``, ``T=1`` keeps the wall-sided window up to 45 degrees.
BARGE = {"length": 10.0, "beam": 2.0, "draft": 1.0, "freeboard": 1.0, "kg": 0.8}
BARGE_ANGLE_STEP = 5.0
#: 0-85 degrees: ``tan`` in the closed form is undefined at 90, and a barge is on its side there.
BARGE_ANGLES = tuple(0.0 + BARGE_ANGLE_STEP * index for index in range(18))

#: The hull of ``examples/gz``, with the freeboard band a generated hull needs to be heeled.
HULL = {"length": 10.0, "beam": 2.0, "draft": 1.0, "nx": 41, "nz": 13, "kg": 0.5}
HULL_FREEBOARD = 1.0
HULL_ANGLE_STEP = 5.0
HULL_ANGLES = tuple(0.0 + HULL_ANGLE_STEP * index for index in range(19))

#: Station grids of the convergence study, coarse to fine.
CONVERGENCE_INTERVALS = (8, 16, 32, 64, 128)
CONVERGENCE_HEEL_DEG = 40.0

#: Tolerances, each tied to the error it is meant to catch, not to what was measured.
ROUND_OFF_TOLERANCE = 1e-12      # a barge is a polygon: exact clip, constant sections
VISIBLE_DEVIATION = 1e-3         # past the wall-sided window the closed form must drift
ROUTE_AGREEMENT = 5e-4           # two integrators, different discretisations, same body
VOLUME_AGREEMENT = 1e-3          # upright displacement, closed surface vs station sum
CONVERGENCE_GAP = 5e-4           # |GZ(64) - GZ(128)| on a 1 m-draft hull
NEGATIVE_CONTROL_SHIFT = 2e-4    # a 2 % displacement error must move GZ by more than this


def barge_mesh(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Watertight rectangular barge: ``draft`` below ``z = 0``, deck at ``z = freeboard``."""
    half = beam / 2.0
    vertices = (
        (0.0, -half, -draft),
        (length, -half, -draft),
        (length, half, -draft),
        (0.0, half, -draft),
        (0.0, -half, freeboard),
        (length, -half, freeboard),
        (length, half, freeboard),
        (0.0, half, freeboard),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return Mesh(vertices=vertices, faces=faces)


def _relative(deviation: float, scale: float) -> float:
    return deviation / scale if scale else math.inf


def build_results() -> tuple[dict, list[str]]:
    """Compute every claim, returning the record and the list of failures."""
    results: dict = {"barge": {}, "hull": {}, "convergence": {}, "controls": {}, "claims": []}
    failures: list[str] = []

    def claim(name: str, value: float, limit: float, unit: str) -> None:
        passed = bool(math.isfinite(value) and value <= limit)
        results["claims"].append(
            {"claim": name, "value": value, "limit": limit, "unit": unit, "passed": passed}
        )
        if not passed:
            failures.append(f"{name}: {value:.6e} {unit} > {limit:.1e} {unit}")

    # ------------------------------------------------------------------ 1/2. barge
    barge = barge_mesh(BARGE["length"], BARGE["beam"], BARGE["draft"], BARGE["freeboard"])
    depth = BARGE["draft"] + BARGE["freeboard"]
    limit = box_wall_sided_limit_deg(BARGE["beam"], BARGE["draft"], depth)
    inside = tuple(angle for angle in BARGE_ANGLES if angle <= limit)
    barge_curve = oracle_curve(
        barge, kg=BARGE["kg"], heel_angles_deg=BARGE_ANGLES, station_intervals=8
    )
    kb = BARGE["draft"] / 2.0
    bm = BARGE["beam"] ** 2 / (12.0 * BARGE["draft"])
    gm = kb + bm - BARGE["kg"]
    # The guarded helper refuses to evaluate above its window, which is right for a caller; the
    # *formula* has to keep producing numbers here, because claim 2 asks how far it drifts, so
    # the unguarded reference is used for the sweep and the guarded helper is checked against it
    # inside the window.
    reference = {angle: wall_sided_righting_arm(angle, gm=gm, bm=bm) for angle in BARGE_ANGLES}
    guarded = box_righting_arms_analytic(
        BARGE["beam"], BARGE["draft"], BARGE["kg"], inside, depth=depth
    )
    guarded_gap = max(
        abs(left - reference[angle]) for angle, left in zip(inside, guarded, strict=True)
    )
    if guarded_gap > ROUND_OFF_TOLERANCE:
        failures.append(
            f"box_righting_arms_analytic disagrees with wall_sided_righting_arm by "
            f"{guarded_gap:.3e} m"
        )
    rows = []
    worst_inside = 0.0
    worst_outside = 0.0
    for point in barge_curve.points:
        expected = reference[point.heel_deg]
        deviation = abs(point.gz - expected)
        if point.heel_deg <= limit:
            worst_inside = max(worst_inside, deviation)
        else:
            worst_outside = max(worst_outside, deviation)
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "gz_m": point.gz,
                "gz_wall_sided_m": expected,
                "abs_deviation_m": deviation,
                "submerged_volume_m3": point.submerged_volume,
                "y_b_m": point.transverse_centre_buoyancy,
                "z_b_m": point.vertical_centre_buoyancy,
                "deck_immersion_fraction": point.deck_immersion_fraction,
            }
        )
    results["barge"] = {
        "parameters": dict(BARGE, depth_m=depth),
        "wall_sided_limit_deg": limit,
        "gm_from_geometry_m": gm,
        "bm_m": bm,
        "kb_m": kb,
        "displacement_volume_m3": barge_curve.displacement_volume,
        "displacement_volume_closed_form_m3": BARGE["length"] * BARGE["beam"] * BARGE["draft"],
        "curve": rows,
    }
    claim("barge GZ vs wall-sided closed form, inside the window", worst_inside, ROUND_OFF_TOLERANCE, "m")
    # Claim 2 is the other direction: the comparison has to be able to see the window end.
    # The nearest out-of-window angle is the honest one to quote (further out, ``tan^2`` makes
    # the closed form meaningless rather than merely wrong), so both are kept.
    first_outside = next((row for row in rows if row["heel_deg"] > limit), None)
    if first_outside is None:  # pragma: no cover - the sampled angles always leave the window
        failures.append("no sampled heel lies above the wall-sided limit")
        first_outside = {"heel_deg": float("nan"), "abs_deviation_m": 0.0}
    results["barge"]["first_outside_window"] = {
        "heel_deg": first_outside["heel_deg"],
        "abs_deviation_m": first_outside["abs_deviation_m"],
    }
    results["claims"].append(
        {
            "claim": "closed form wrong just above its window (negative control)",
            "value": first_outside["abs_deviation_m"],
            "limit": VISIBLE_DEVIATION,
            "unit": "m (must exceed)",
            "passed": first_outside["abs_deviation_m"] >= VISIBLE_DEVIATION,
        }
    )
    if first_outside["abs_deviation_m"] < VISIBLE_DEVIATION:
        failures.append(
            f"the closed form is only {first_outside['abs_deviation_m']:.3e} m out just above "
            "its window, so claim 1 is vacuous"
        )

    # ------------------------------------------------- 3. deck-edge immersion onset
    onset = next(
        (point.heel_deg for point in barge_curve.points if point.deck_immersion_fraction >= 0.5),
        float("nan"),
    )
    onset_error = abs(onset - limit) if math.isfinite(onset) else math.inf
    results["barge"]["deck_immersion_onset_deg"] = onset
    results["barge"]["deck_immersion_onset_error_deg"] = onset_error
    claim("deck-edge immersion onset vs wall-sided limit", onset_error, BARGE_ANGLE_STEP, "deg")

    # ------------------------------------------------------------- 4/5/7. the hull
    hull = generate_hull(
        HullParameters(
            length=HULL["length"],
            beam=HULL["beam"],
            draft=HULL["draft"],
            nx=HULL["nx"],
            nz=HULL["nz"],
        )
    ).mesh
    banded = with_freeboard(hull, HULL_FREEBOARD)
    canonical = righting_arm_curve(banded, HULL["kg"], list(HULL_ANGLES), waterline_z=0.0)
    section_curve = oracle_curve(
        banded,
        kg=HULL["kg"],
        heel_angles_deg=HULL_ANGLES,
        station_intervals=DEFAULT_STATION_INTERVALS,
    )
    worst_route = 0.0
    hull_rows = []
    for angle, canonical_gz, heeled, section_point in zip(
        canonical.angles_deg,
        canonical.gz,
        canonical.points,
        section_curve.points,
        strict=True,
    ):
        difference = section_point.gz - canonical_gz
        worst_route = max(worst_route, abs(difference))
        hull_rows.append(
            {
                "heel_deg": angle,
                "gz_canonical_m": canonical_gz,
                "gz_station_m": section_point.gz,
                "abs_difference_m": abs(difference),
                "canonical_volume_residual_m3": heeled.displacement_volume
                - canonical.displacement_volume,
                "station_volume_residual_m3": section_point.submerged_volume
                - section_curve.displacement_volume,
            }
        )
    upright = heeled_hydrostatics(banded, 0.0, waterline_z=0.0)
    volume_gap = _relative(
        abs(section_curve.displacement_volume - upright.displacement_volume),
        upright.displacement_volume,
    )
    results["hull"] = {
        "parameters": dict(HULL),
        "freeboard_m": HULL_FREEBOARD,
        "station_intervals": DEFAULT_STATION_INTERVALS,
        "gz_max_canonical_m": canonical.max_gz,
        "gz_max_station_m": section_curve.gz_max,
        "angle_gz_max_canonical_deg": canonical.angle_of_max_gz_deg,
        "angle_gz_max_station_deg": section_curve.angle_gz_max_deg,
        "displacement_canonical_m3": canonical.displacement_volume,
        "displacement_station_m3": section_curve.displacement_volume,
        "displacement_closed_surface_m3": upright.displacement_volume,
        "max_abs_route_difference_m": worst_route,
        "max_relative_route_difference": _relative(worst_route, canonical.max_gz),
        "max_volume_error_canonical": canonical.max_volume_error,
        "max_volume_error_station": section_curve.max_volume_error,
        "curve": hull_rows,
    }
    claim("hull: independent routes agree (relative to max GZ)", _relative(worst_route, canonical.max_gz), ROUTE_AGREEMENT, "rel")
    claim("hull: both routes float the same upright displacement", volume_gap, VOLUME_AGREEMENT, "rel")

    # ------------------------------------------------------------ 7. negative control
    shifted = righting_arm_curve(
        banded,
        HULL["kg"],
        list(HULL_ANGLES),
        waterline_z=0.0,
        displacement_volume=canonical.displacement_volume * 1.02,
    )
    control_shift = max(np - n for np, n in zip(shifted.gz, canonical.gz, strict=True))
    results["controls"] = {
        "description": "canonical curve with its displacement target scaled by 1.02",
        "displacement_scale": 1.02,
        "gz_shift_at_20deg_m": shifted.gz[list(HULL_ANGLES).index(20.0)]
        - canonical.gz[list(HULL_ANGLES).index(20.0)],
        "max_gz_shift_m": control_shift,
    }
    results["claims"].append(
        {
            "claim": "2 % displacement error moves GZ (negative control)",
            "value": control_shift,
            "limit": NEGATIVE_CONTROL_SHIFT,
            "unit": "m (must exceed)",
            "passed": control_shift >= NEGATIVE_CONTROL_SHIFT,
        }
    )
    if control_shift < NEGATIVE_CONTROL_SHIFT:
        failures.append(
            f"a 2 % displacement error shifts GZ by only {control_shift:.3e} m, so the "
            "agreement test is insensitive"
        )

    # -------------------------------------------------------- 6. station convergence
    heel = CONVERGENCE_HEEL_DEG
    sweep = []
    for intervals in CONVERGENCE_INTERVALS:
        curve = oracle_curve(
            banded,
            kg=HULL["kg"],
            heel_angles_deg=(0.0, heel),
            station_intervals=intervals,
        )
        sweep.append({"station_intervals": intervals, "gz_m": curve.gz[1]})
    converged = sweep[-1]["gz_m"]
    reported = next(
        row["gz_station_m"]
        for row in hull_rows
        if row["heel_deg"] == heel
    )
    for row in sweep:
        row["abs_error_m"] = abs(row["gz_m"] - converged)
        row["abs_error_relative"] = _relative(row["abs_error_m"], converged)
    coarse_error = sweep[0]["abs_error_m"]
    reported_error = abs(reported - converged)
    ratio = coarse_error / reported_error if reported_error else math.inf
    results["convergence"] = {
        "heel_deg": heel,
        "intervals": list(CONVERGENCE_INTERVALS),
        "sweep": sweep,
        "converged_gz_m": converged,
        "reported_station_intervals": DEFAULT_STATION_INTERVALS,
        "reported_abs_error_m": reported_error,
        "coarsest_abs_error_m": coarse_error,
        "coarse_over_reported_error": ratio,
        "canonical_gz_m": canonical.gz[list(HULL_ANGLES).index(heel)],
        "canonical_minus_converged_m": canonical.gz[list(HULL_ANGLES).index(heel)] - converged,
    }
    claim(f"station grid {DEFAULT_STATION_INTERVALS} vs {CONVERGENCE_INTERVALS[-1]} at {heel:.0f} deg", reported_error, CONVERGENCE_GAP, "m")
    results["claims"].append(
        {
            "claim": "the coarse station grid is measurably worse (convergence is real)",
            "value": ratio,
            "limit": 4.0,
            "unit": "x (must exceed)",
            "passed": ratio >= 4.0,
        }
    )
    if ratio < 4.0:
        failures.append(
            f"the {CONVERGENCE_INTERVALS[0]}-interval error is only {ratio:.2f}x the "
            f"{DEFAULT_STATION_INTERVALS}-interval error, so convergence is not demonstrated"
        )

    results["checks"] = {"failures": failures, "passed": not failures}
    return results, failures


def _plot(results: dict, path: Path) -> bool:
    """Three panels: the barge against its closed form, the hull cross-route, convergence."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - optional dependency
        print("matplotlib not installed; skipping plot")
        return False

    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))

    rows = results["barge"]["curve"]
    heels = [row["heel_deg"] for row in rows]
    axes[0].plot(
        heels, [row["gz_wall_sided_m"] for row in rows], "k-", lw=1.2,
        label=r"closed form $\sin\phi\,(GM+\frac{1}{2}BM\tan^2\phi)$",
    )
    axes[0].plot(
        heels, [row["gz_m"] for row in rows], "o", ms=4, mfc="none",
        label="station integrator (oracle)",
    )
    axes[0].axvline(
        results["barge"]["wall_sided_limit_deg"], color="0.4", ls=":", lw=1.0,
        label=f"wall-sided limit {results['barge']['wall_sided_limit_deg']:.0f} deg",
    )
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].set_ylim(0.0, 1.25 * max(row["gz_m"] for row in rows))
    axes[0].set_title("Box barge 10x2x2 m, KG 0.8 m")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)

    rows = results["hull"]["curve"]
    heels = [row["heel_deg"] for row in rows]
    axes[1].plot(
        heels, [row["gz_canonical_m"] for row in rows], "-", lw=1.4,
        label="canonical (tensorcad.gz, divergence theorem)",
    )
    axes[1].plot(
        heels, [row["gz_station_m"] for row in rows], "o", ms=4, mfc="none",
        label="station integrator (PR #52 method)",
    )
    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].set_xlabel("heel angle [deg]")
    axes[1].set_ylabel("GZ [m]")
    axes[1].set_title(
        f"Hull 41x13, freeboard 1.0 m: max |diff| "
        f"{results['hull']['max_abs_route_difference_m'] * 1000:.3f} mm"
    )
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)

    sweep = results["convergence"]["sweep"]
    axes[2].loglog(
        [row["station_intervals"] for row in sweep],
        [max(row["abs_error_m"], 1e-16) for row in sweep],
        "o-",
        label=f"|GZ - GZ(128)| at {results['convergence']['heel_deg']:.0f} deg",
    )
    axes[2].axhline(
        abs(results["convergence"]["canonical_minus_converged_m"]), color="0.4", ls="--", lw=1.0,
        label="|canonical - GZ(128)|",
    )
    axes[2].set_xlabel("station intervals")
    axes[2].set_ylabel("error in GZ [m]")
    axes[2].set_title("Station-grid convergence")
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.3, which="both")

    figure.tight_layout()
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return True


def _report(results: dict) -> list[str]:
    barge = results["barge"]
    hull = results["hull"]
    convergence = results["convergence"]
    controls = results["controls"]
    inside = [row for row in barge["curve"] if row["heel_deg"] <= barge["wall_sided_limit_deg"]]
    outside = [row for row in barge["curve"] if row["heel_deg"] > barge["wall_sided_limit_deg"]]
    lines = [
        "# Righting arms: two integrators on one body",
        "",
        "Generated by `examples/righting/run_example.py`; every number below is recomputed on",
        "each run and the run exits non-zero if a claim breaks. The canonical route is",
        "`tensorcad.gz` (closed-surface integration, Newton on the waterline); the oracle is",
        "`examples/righting/oracle.py` (heeled sections clipped and integrated by Green's",
        "theorem, Simpson along the hull, bisection on the tilted waterplane).",
        "",
        "## 1. Box barge against the wall-sided closed form",
        "",
        f"- barge {barge['parameters']['length']:.0f} x {barge['parameters']['beam']:.0f} x"
        f" {barge['parameters']['depth_m']:.0f} m, draft {barge['parameters']['draft']:.1f} m,"
        f" KG {barge['parameters']['kg']:.1f} m; `KB = {barge['kb_m']:.3f} m`,",
        f"  `BM = {barge['bm_m']:.6f} m`, `GM = {barge['gm_from_geometry_m']:.6f} m`",
        f"- displacement {barge['displacement_volume_m3']:.6f} m3, closed form"
        f" {barge['displacement_volume_closed_form_m3']:.6f} m3",
        f"- worst |GZ - closed form| inside the window (0-{barge['wall_sided_limit_deg']:.0f} deg):"
        f" **{max(row['abs_deviation_m'] for row in inside):.3e} m**"
        f" over {len(inside)} angles",
        f"- just above the window ({barge['first_outside_window']['heel_deg']:.0f} deg) the closed"
        f" form is already **{barge['first_outside_window']['abs_deviation_m'] * 1000:.2f} mm** out"
        f" (and {max(row['abs_deviation_m'] for row in outside):.3e} m by"
        f" {max(outside, key=lambda row: row['abs_deviation_m'])['heel_deg']:.0f} deg, where"
        " `tan^2` makes it meaningless rather than merely wrong), so the agreement above is",
        "  not vacuous",
        f"- deck-edge immersion onset {barge['deck_immersion_onset_deg']:.0f} deg vs wall-sided limit"
        f" {barge['wall_sided_limit_deg']:.2f} deg",
        "",
        "## 2. Hull: the two independent routes",
        "",
        f"- hull {hull['parameters']['length']:.0f} x {hull['parameters']['beam']:.1f} x"
        f" {hull['parameters']['draft']:.1f} m, {hull['parameters']['nx']}x{hull['parameters']['nz']} mesh,"
        f" freeboard {hull['freeboard_m']:.1f} m, KG {hull['parameters']['kg']:.2f} m",
        f"- upright displacement: canonical {hull['displacement_canonical_m3']:.6f} m3,"
        f" station {hull['displacement_station_m3']:.6f} m3,"
        f" closed surface {hull['displacement_closed_surface_m3']:.6f} m3",
        f"- max |GZ_canonical - GZ_station| = **{hull['max_abs_route_difference_m'] * 1000:.3f} mm**"
        f" ({hull['max_relative_route_difference'] * 100:.4f} % of the largest arm"
        f" {hull['gz_max_canonical_m']:.4f} m)",
        f"- max GZ {hull['gz_max_canonical_m']:.4f} m at {hull['angle_gz_max_canonical_deg']:.1f} deg"
        f" (canonical) vs {hull['gz_max_station_m']:.4f} m at"
        f" {hull['angle_gz_max_station_deg']:.1f} deg (station)",
        f"- displacement residual along the curve: canonical"
        f" {hull['max_volume_error_canonical']:.2e} relative, station"
        f" {hull['max_volume_error_station']:.2e} relative",
        "",
        "## 3. Station-grid convergence",
        "",
        f"`GZ` at {convergence['heel_deg']:.0f} deg against the finest grid"
        f" ({convergence['intervals'][-1]} intervals):",
        "",
        "| station intervals | GZ [m] | |GZ - GZ(128)| [m] |",
        "|---|---|---|",
    ]
    for row in convergence["sweep"]:
        lines.append(
            f"| {row['station_intervals']} | {row['gz_m']:.9f} | {row['abs_error_m']:.3e} |"
        )
    lines += [
        "",
        f"The reported grid ({DEFAULT_STATION_INTERVALS} intervals) is"
        f" {convergence['reported_abs_error_m']:.2e} m from the finest, while the"
        f" {convergence['intervals'][0]}-interval grid is"
        f" {convergence['coarse_over_reported_error']:.1f}x worse -- the quadrature is the",
        "error source, and it is under control at the reported setting.",
        f"Canonical route at the same angle: {convergence['canonical_gz_m']:.9f} m, i.e."
        f" {convergence['canonical_minus_converged_m']:.2e} m from the finest station grid.",
        "",
        "## 4. Negative control",
        "",
        f"Scaling the canonical route's displacement target by"
        f" {controls['displacement_scale']:.2f} moves `GZ(20 deg)` by",
        f"{controls['gz_shift_at_20deg_m'] * 1000:.2f} mm and the curve by up to"
        f" {controls['max_gz_shift_m'] * 1000:.2f} mm, so the agreement in section 2 is a",
        "statement about two integrations and not about two constants.",
        "",
        "## Verdict",
        "",
    ]
    passed = results["checks"]["passed"]
    lines.append(
        "**All claims hold.**" if passed else "**FAILED:** " + "; ".join(results["checks"]["failures"])
    )
    lines.append("")
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "docs" / "righting",
        help="directory for results.json, report.md and righting_arms.png",
    )
    parser.add_argument("--no-plot", action="store_true", help="skip the figure")
    options = parser.parse_args(argv)
    output: Path = options.output
    output.mkdir(parents=True, exist_ok=True)

    results, failures = build_results()
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    lines = _report(results)
    (output / "report.md").write_text("\n".join(lines))
    if not options.no_plot:  # pragma: no branch - plotted by default
        _plot(results, output / "righting_arms.png")

    print("\n".join(lines))
    for failure in failures:
        print(f"FAIL: {failure}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
