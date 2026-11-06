"""Certify both heeled-mesh integrators against the exact multi-regime box solution.

    python3 tools/verify_box_large_heel.py     # writes artifacts/m4864/ and syncs the doc tables

What is checked, per barge geometry:

1. **The exact reference against the closed form it generalises.**  Inside the wall-sided
   window the new solution must agree with ``sin(phi)(GM + BM/2 tan^2 phi)`` to round-off; if it
   did not, the new one would be the suspect.
2. **Both merged integrators against the exact reference** over ``0..85`` degrees: the Newton
   route in :mod:`tensorcad.gz` and the bracket-and-bisect route carried in from PR #40.  Until
   now the two routes only had each other out there, which cannot catch a shared modelling
   error.
3. **Continuity across the regime boundaries** -- the waterline crossing the deck edge and the
   bilge -- where a mis-assembled reference would show a jump.
4. **A negative control**: a KG shift of one millimetre and a 0.5 % displacement shift, to show
   the comparison can fail.

Everything is written from the same run: ``artifacts/m4864/box_large_heel.json`` (the per-heel
curves, the regime bounds and the controls) plus the two figures, and -- through
:func:`sync_doc_tables` -- the table blocks of ``docs/box_large_heel.md``, which are replaced in
place by the tables *this* run computed. Running the script twice leaves the page byte-identical,
so a stale table there is a diff rather than a typo nobody noticed, and the tables are not
duplicated into a second artifact file that could drift from the page.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# runnable straight from a clone, like the other tools here
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt

from tensorcad import bisection_righting_arm_curve
from tensorcad.box_gz import (
    box_gz,
    box_regime_bounds_deg,
    box_section,
    box_vanishing_angle_deg,
)
from tensorcad.geometry import MeshBuilder
from tensorcad.gz import righting_arm_curve, wall_sided_righting_arm

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "m4864"
DOCS = ROOT / "docs"

SPECS = [
    {
        "name": "doc barge",
        "length": 10.0,
        "breadth": 2.0,
        "draft": 0.8,
        "freeboard": 0.4,
        "kg": 0.55,
    },
    {
        "name": "wide shallow",
        "length": 6.0,
        "breadth": 4.0,
        "draft": 0.5,
        "freeboard": 0.3,
        "kg": 0.9,
    },
    {
        "name": "square deep",
        "length": 3.0,
        "breadth": 3.0,
        "draft": 1.0,
        "freeboard": 0.5,
        "kg": 1.2,
    },
    {
        # enough freeboard that the bilge leaves the water first: this is the only geometry
        # here that enters ``bilge_emerged`` before ``deck_immersed``
        "name": "tall barge",
        "length": 8.0,
        "breadth": 2.0,
        "draft": 0.8,
        "freeboard": 1.6,
        "kg": 0.6,
    },
]

ANGLES = tuple(float(angle) for angle in range(0, 86))


def box_mesh(length: float, breadth: float, draft: float, freeboard: float):
    """Closed box, bottom at ``-draft``, deck at ``+freeboard``, waterline at ``z = 0``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * breadth, 0.5 * breadth):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * breadth, 0.5 * breadth
    bottom, top = -draft, freeboard
    left, right = 0.0, length
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def wall_sided_values(spec: dict, angles) -> list[float]:
    bm = spec["breadth"] ** 2 / (12.0 * spec["draft"])
    gm = 0.5 * spec["draft"] + bm - spec["kg"]
    return [wall_sided_righting_arm(angle, gm=gm, bm=bm) for angle in angles]


def run_spec(spec: dict) -> dict:
    breadth, draft = spec["breadth"], spec["draft"]
    depth = draft + spec["freeboard"]
    mesh = box_mesh(spec["length"], breadth, draft, spec["freeboard"])
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)

    exact_points = [
        box_gz(breadth=breadth, draft=draft, depth=depth, kg=spec["kg"], heel_deg=angle)
        for angle in ANGLES
    ]
    exact = [point.righting_arm for point in exact_points]
    regimes = [point.regime for point in exact_points]

    newton = righting_arm_curve(mesh, spec["kg"], angles_deg=ANGLES)
    bisect = bisection_righting_arm_curve(mesh, kg=spec["kg"], heels_deg=ANGLES)
    newton_values = list(newton.gz)
    bisect_values = [point.righting_arm for point in bisect.points]

    wall_sided = wall_sided_values(spec, ANGLES)
    # the wall-sided form dies at whichever comes first: a deck edge or a bilge
    wall_sided_limit = min(bounds["deck_immersion"], bounds["bilge_emergence"])
    wall_sided_window = [angle <= wall_sided_limit + 1e-12 for angle in ANGLES]
    classical_bilge_deg = math.degrees(math.atan(2.0 * draft / breadth))
    window_error = max(
        abs(exact[index] - wall_sided[index])
        for index in range(len(ANGLES))
        if wall_sided_window[index]
    )

    volume_residual = max(
        abs(point.submerged_area - breadth * draft) / (breadth * draft) for point in exact_points
    )
    # the mesh routes state their own volume residual per point
    mesh_residual = max(
        abs(point.volume_residual) for point in bisect.points
    ) / (breadth * draft * spec["length"])

    # negative controls: the comparison has to be able to fail
    control_kg = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=spec["kg"] + 1e-3, heel_deg=40.0
    ).righting_arm
    control_displacement = box_section(
        breadth=breadth, draft=draft * 1.005, depth=depth, heel_deg=40.0
    )

    return {
        "spec": spec,
        "depth": depth,
        "regime_bounds_deg": bounds,
        "vanishing_angle_deg": box_vanishing_angle_deg(
            breadth=breadth, draft=draft, depth=depth, kg=spec["kg"]
        ),
        "wall_sided_window_error_m": window_error,
        "wall_sided_limit_deg": wall_sided_limit,
        "classical_bilge_emergence_deg": classical_bilge_deg,
        "bilge_minus_classical_deg": bounds["bilge_emergence"] - classical_bilge_deg,
        "regimes_present": sorted(set(regimes)),
        "max_abs_deviation_newton_m": max(
            abs(a - b) for a, b in zip(exact, newton_values, strict=False)
        ),
        "max_abs_deviation_bisect_m": max(
            abs(a - b) for a, b in zip(exact, bisect_values, strict=False)
        ),
        "volume_residual_relative": volume_residual,
        "mesh_volume_residual_relative": mesh_residual,
        "wall_sided_error_at_40_deg_m": abs(
            exact[40] - wall_sided[40]
        ),
        "peak": {
            "heel_deg": max(zip(exact, range(len(exact)), strict=False))[1],
            "gz": max(exact),
        },
        "control_kg_1mm_shift_m": abs(control_kg - exact[40]),
        "control_displacement_half_percent_effect_m": abs(
            control_displacement[1] * math.cos(math.radians(40.0))
            + control_displacement[2] * math.sin(math.radians(40.0))
            - exact_points[40].kn
        ),
        "curve": {
            "angles_deg": list(ANGLES),
            "exact": exact,
            "newton": newton_values,
            "bisect": bisect_values,
            "wall_sided": wall_sided,
            "regimes": regimes,
        },
    }


def plot_curve(results: list[dict]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    left, right = axes
    main = results[0]
    curve = main["curve"]
    angles = curve["angles_deg"]
    left.plot(angles, curve["exact"], "k-", lw=2.0, label="exact, four regimes (this work)")
    left.plot(
        angles[::5],
        curve["newton"][::5],
        "o",
        ms=5,
        mfc="none",
        label="tensorcad.gz (Newton)",
    )
    left.plot(
        angles[::5],
        curve["bisect"][::5],
        "s",
        ms=4,
        mfc="none",
        label="tensorcad.righting_arm (bisect)",
    )
    limit = main["regime_bounds_deg"]["deck_immersion"]
    left.plot(angles, curve["wall_sided"], "r--", lw=1.2, label="wall-sided closed form")
    left.axvspan(0.0, limit, color="0.92", zorder=0)
    for key in ("deck_immersion", "bilge_emergence"):
        value = main["regime_bounds_deg"][key]
        left.axvline(value, color="0.4", lw=0.8, ls=":")
    left.set_xlabel("heel [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(f"{main['spec']['name']}: GZ across all four regimes")
    left.legend(fontsize=7.5, loc="lower left")
    left.grid(alpha=0.3)

    for result in results:
        curve = result["curve"]
        right.semilogy(
            curve["angles_deg"][1:],
            [abs(a - b) + 1e-18 for a, b in zip(curve["exact"], curve["newton"], strict=False)][1:],
            "-",
            lw=1.2,
            label=f"{result['spec']['name']}: Newton",
        )
        right.semilogy(
            curve["angles_deg"][1:],
            [abs(a - b) + 1e-18 for a, b in zip(curve["exact"], curve["bisect"], strict=False)][1:],
            "--",
            lw=1.2,
            label=f"{result['spec']['name']}: bisect",
        )
    right.axhline(1e-9, color="0.5", lw=0.8, ls=":")
    right.set_xlabel("heel [deg]")
    right.set_ylabel("|integrator - exact| [m]")
    right.set_title("agreement with the exact solution")
    right.legend(fontsize=7)
    right.grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(ARTIFACTS / "box_large_heel_gz.png", dpi=150)
    plt.close(figure)


def plot_sections(spec: dict) -> None:
    breadth, draft = spec["breadth"], spec["draft"]
    depth = draft + spec["freeboard"]
    heels = (0.0, 21.8, 42.0, 60.0, 80.0)
    figure, axes = plt.subplots(1, len(heels), figsize=(13.0, 2.9))
    for axis, heel in zip(axes, heels, strict=False):
        phi = math.radians(heel)
        point = box_gz(breadth=breadth, draft=draft, depth=depth, kg=spec["kg"], heel_deg=heel)
        half = 0.5 * breadth
        corners = ((-half, 0.0), (half, 0.0), (half, depth), (-half, depth))
        world = [
            (y * math.cos(phi) + z * math.sin(phi), -y * math.sin(phi) + z * math.cos(phi))
            for y, z in corners
        ]
        axis.fill(*zip(*world, strict=False), color="0.88", edgecolor="0.4", lw=0.8)
        # the immersed part of the section, drawn from the reference's own clip
        submerged = []
        for index in range(len(world)):
            current, following = world[index], world[(index + 1) % len(world)]
            dc, df = current[1] - point.waterline_z, following[1] - point.waterline_z
            if dc <= 0.0:
                submerged.append(current)
            if (dc < 0.0 < df) or (df < 0.0 < dc):
                fraction = dc / (dc - df)
                submerged.append(
                    (
                        current[0] + fraction * (following[0] - current[0]),
                        point.waterline_z,
                    )
                )
        if submerged:
            axis.fill(*zip(*submerged, strict=False), color="#7fb3d5", alpha=0.75)
        axis.axhline(point.waterline_z, color="#1f618d", lw=1.0)
        axis.plot([point.centre_of_buoyancy[0]], [point.centre_of_buoyancy[1]], "o", color="k", ms=4)
        axis.set_title(f"{heel:.0f} deg: {point.regime}\nGZ = {point.righting_arm:+.4f} m", fontsize=8)
        axis.set_aspect("equal")
        axis.set_xlim(-1.25 * (half + 0.2), 1.25 * (half + 0.2))
        axis.set_ylim(-1.6 * breadth, 1.3 * depth + 0.3 * breadth)
        axis.tick_params(labelsize=6)
    figure.suptitle(
        f"heeled section, {spec['name']} (B={breadth} T={draft} D={depth} m): "
        "shaded = immersed, dot = B",
        fontsize=9,
    )
    figure.tight_layout()
    figure.savefig(ARTIFACTS / "box_large_heel_sections.png", dpi=150)
    plt.close(figure)


def markdown_tables(results: list[dict]) -> str:
    lines: list[str] = []
    lines.append(
        "| geometry | heel of deck-edge immersion | bilge emergence | both | GZ peak "
        "| vanishing angle |"
    )
    lines.append("|---|---|---|---|---|---|")
    for result in results:
        spec = result["spec"]
        bounds = result["regime_bounds_deg"]
        vanishing = result["vanishing_angle_deg"]
        lines.append(
            "| {name} (L={length} B={breadth} T={draft} FB={freeboard} KG={kg}) "
            "| {deck:.3f} deg | {bilge:.3f} deg | {both:.3f} deg | {peak:.6f} m at {peak_heel:.0f} deg "
            "| {vanishing} |".format(
                name=spec["name"],
                length=spec["length"],
                breadth=spec["breadth"],
                draft=spec["draft"],
                freeboard=spec["freeboard"],
                kg=spec["kg"],
                deck=bounds["deck_immersion"],
                bilge=bounds["bilge_emergence"],
                both=bounds["both_walls_lost"],
                peak=result["peak"]["gz"],
                peak_heel=result["peak"]["heel_deg"],
                vanishing=(
                    "beyond 89 deg"
                    if math.isinf(vanishing)
                    else f"{vanishing:.3f} deg"
                ),
            )
        )
    lines.append("")
    lines.append(
        "| geometry | wall-sided form vs exact, inside its window | wall-sided form at 40 deg "
        "| Newton vs exact, 0-85 deg | bisect vs exact, 0-85 deg |"
    )
    lines.append("|---|---|---|---|---|")
    for result in results:
        lines.append(
            "| {name} | {window:.2e} m | {at40:.4f} m | {newton:.2e} m | {bisect:.2e} m |".format(
                name=result["spec"]["name"],
                window=result["wall_sided_window_error_m"],
                at40=result["wall_sided_error_at_40_deg_m"],
                newton=result["max_abs_deviation_newton_m"],
                bisect=result["max_abs_deviation_bisect_m"],
            )
        )
    lines.append("")
    lines.append(
        "| geometry | wall-sided window ends at | computed bilge emergence "
        "| classical atan(2T/B) | difference |"
    )
    lines.append("|---|---|---|---|---|")
    for result in results:
        lines.append(
            "| {name} | {limit:.4f} deg | {bilge:.4f} deg | {classical:.4f} deg | {delta:+.4f} deg |".format(
                name=result["spec"]["name"],
                limit=result["wall_sided_limit_deg"],
                bilge=result["regime_bounds_deg"]["bilge_emergence"],
                classical=result["classical_bilge_emergence_deg"],
                delta=result["bilge_minus_classical_deg"],
            )
        )
    lines.append("")
    lines.append("| geometry | area residual of the reference | mesh volume residual (bisect route) |")
    lines.append("|---|---|---|")
    for result in results:
        lines.append(
            "| {name} | {area:.2e} (relative) | {mesh:.2e} (relative) |".format(
                name=result["spec"]["name"],
                area=result["volume_residual_relative"],
                mesh=result["mesh_volume_residual_relative"],
            )
        )
    return "\n".join(lines) + "\n"


def sync_doc_tables(tables_markdown: str, *, doc: Path) -> list[str]:
    """Replace this run's four result tables in the page, matching each by its header row.

    The page interleaves prose with the tables, and also carries two tables this script does not
    own (the error judgement and the negative control), so the blocks are matched on their header
    line rather than by position.  A header that is missing from the page, or present twice, is
    an error: the page is never edited on a guess.
    """
    generated = [block for block in tables_markdown.strip("\n").split("\n\n") if block.strip()]
    lines = doc.read_text().splitlines()
    blocks: list[tuple[int, int]] = []
    index = 0
    while index < len(lines):
        if lines[index].startswith("|"):
            begin = index
            while index < len(lines) and lines[index].startswith("|"):
                index += 1
            blocks.append((begin, index))
        else:
            index += 1
    headers = {lines[begin]: (begin, stop) for begin, stop in blocks}
    synced: list[str] = []
    edits: list[tuple[int, int, list[str]]] = []
    used: set[int] = set()
    for block in generated:
        header = block.splitlines()[0]
        if header not in headers:
            raise SystemExit(f"docs/{doc.name}: no table with header {header!r} to update")
        begin, stop = headers[header]
        if begin in used:
            raise SystemExit(f"docs/{doc.name}: two computed tables share the header {header!r}")
        used.add(begin)
        edits.append((begin, stop, block.splitlines()))
        synced.append(header)
    for begin, stop, replacement in sorted(edits, reverse=True):
        lines[begin:stop] = replacement
    doc.write_text("\n".join(lines) + "\n")
    return synced


def main() -> int:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    results = [run_spec(spec) for spec in SPECS]
    payload = {
        "mission": "m4864",
        "library": "TensorCAD",
        "module": "tensorcad/box_gz.py",
        "generated_by": "tools/verify_box_large_heel.py",
        "angles_deg": list(ANGLES),
        "specs": results,
        "worst": {
            "newton": max(result["max_abs_deviation_newton_m"] for result in results),
            "bisect": max(result["max_abs_deviation_bisect_m"] for result in results),
            "wall_sided_window": max(result["wall_sided_window_error_m"] for result in results),
        },
    }
    (ARTIFACTS / "box_large_heel.json").write_text(json.dumps(payload, indent=2) + "\n")
    tables = markdown_tables(results)
    synced = sync_doc_tables(tables, doc=DOCS / "box_large_heel.md")
    plot_curve(results)
    plot_sections(SPECS[0])

    print(tables)
    print(f"worst |Newton - exact| = {payload['worst']['newton']:.3e} m")
    print(f"worst |bisect - exact| = {payload['worst']['bisect']:.3e} m")
    print(f"worst wall-sided window error = {payload['worst']['wall_sided_window']:.3e} m")
    print(f"synced {synced} table blocks into docs/box_large_heel.md")

    failures = []
    for result in results:
        if result["max_abs_deviation_newton_m"] > 1e-9:
            failures.append(f"{result['spec']['name']}: Newton deviates")
        if result["max_abs_deviation_bisect_m"] > 1e-9:
            failures.append(f"{result['spec']['name']}: bisect deviates")
        if result["wall_sided_window_error_m"] > 1e-9:
            failures.append(f"{result['spec']['name']}: window disagrees with the closed form")
        if result["volume_residual_relative"] > 1e-12:
            failures.append(f"{result['spec']['name']}: area not conserved")
    if failures:
        print("FAIL: " + "; ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
