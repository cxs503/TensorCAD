"""A waterline plane that contains mesh edges: the freeboard seam, measured once.

Run:  python examples/freeboard_seam/run_example.py --output examples/freeboard_seam/outputs

A hull is generated only up to its waterline and then given freeboard by
:func:`tensorcad.gz.with_freeboard`, because a GZ curve needs the waterline to stay on the
hull as it heels. The band is vertical, so the waterplane at the design waterline is
*identical* with and without it -- but the seam where band and hull meet is a whole ring of
edges lying **in** the cutting plane, and a plane that contains an edge is met by both of
that edge's faces. Every boundary segment was therefore emitted twice, the contour was
traced twice, and the two coincident loops were read as one loop nested inside the other
and subtracted: the Wigley hull reported a waterplane area of ``-13.33125 m^2`` (exactly
``-2 x 6.665625``) and refused to report moments at all.

The same seam also exposed a second, quieter defect: ``initial_stability`` integrated the
*whole closed mesh*, which is the displacement only for a hull capped at its waterline. The
banded Wigley handed back ``4.108733 m^3`` as its displacement instead of ``2.775608 m^3``
-- total body volume -- and a ``GM`` 28 % too large.

This example pins the repaired behaviour against closed forms, from one mesh:

1. the section through the banded mesh equals the section through the capped hull
   (area, ``I_T``, ``I_L``) to 1e-14 relative -- a waterplane must not depend on whether
   freeboard is modelled;
2. ``GM`` measured on the banded mesh equals ``GM`` measured on the capped hull to 1e-9,
   and both approach the analytic ``KB + I_T/nabla - KG`` as the grid refines;
3. the slope of the GZ curve at 2 degrees equals ``GM`` from the same banded mesh to
   0.1 %, widening to the expected ``BM/2 tan^2 phi`` wall-sided term at larger heel;
4. every heeled point floats at the upright displacement (< 1e-10 relative).

The Wigley closed forms used as the reference (``n = m = 2``, ``L``, ``B``, ``T``)::

    nabla = 4 L B T / 9
    KB    = 5 T / 8                     above the keel
    A_wp  = 2 L B / 3
    I_T   = 4 B^3 L / 105               about the centreline
    I_L   = B L^3 / 5 - A_wp (L/2)^2    about the transverse axis through the centroid

Writes ``results.json``, ``report.md`` and ``plot.png`` into the output directory and exits
non-zero if a claim stops holding, so CI can run it as a check rather than a demo.
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
from pathlib import Path

from tensorcad.gz import with_freeboard
from tensorcad.righting_arm import righting_arm_curve
from tensorcad.stability import initial_stability, waterplane_moments
from tensorcad.wigley import WigleyParameters, generate_wigley

LENGTH, BEAM, DRAFT = 10.0, 1.0, 0.625
FREEBOARD, KG, WATERLINE = 0.2, 0.3, 0.0
GRIDS = ((41, 13), (81, 21), (161, 41))
SLOPE_HEELS = (2.0, 5.0, 10.0)


def analytic_reference() -> dict[str, float]:
    """Wigley ``n = m = 2`` closed forms, derived in the module docstring."""
    volume = 4.0 * LENGTH * BEAM * DRAFT / 9.0
    kb = 5.0 * DRAFT / 8.0
    waterplane = 2.0 * LENGTH * BEAM / 3.0
    transverse = 4.0 * BEAM**3 * LENGTH / 105.0
    longitudinal = BEAM * LENGTH**3 / 5.0 - waterplane * (0.5 * LENGTH) ** 2
    return {
        "displacement_volume_m3": volume,
        "kb_m": kb,
        "waterplane_area_m2": waterplane,
        "transverse_second_moment_m4": transverse,
        "longitudinal_second_moment_m4": longitudinal,
        "bm_transverse_m": transverse / volume,
        "gm_transverse_m": kb + transverse / volume - KG,
    }


class Claims:
    """Collect ``(name, value, limit, detail)`` and fail the run if any is out of limit."""

    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(self, name: str, value: float, limit: float, detail: str) -> None:
        self.rows.append(
            {
                "claim": name,
                "value": value,
                "limit": limit,
                "unit": "relative error" if limit < 1.0 else "m^3",
                "holds": bool(value <= limit),
                "detail": detail,
            }
        )

    @property
    def failures(self) -> list[dict[str, object]]:
        return [row for row in self.rows if not row["holds"]]


def relative(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


def measure(nx: int, nz: int) -> dict[str, object]:
    """Capped hull and its freeboarded twin, all quantities from each mesh."""
    parameters = WigleyParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)
    model = generate_wigley(parameters)
    capped_mesh = model.mesh
    banded_mesh = with_freeboard(capped_mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)

    row: dict[str, object] = {"nx": nx, "nz": nz, "cells": nx * nz}
    capped_section = waterplane_moments(capped_mesh, WATERLINE)
    banded_section = waterplane_moments(banded_mesh, WATERLINE)
    row["capped"] = {
        "area": capped_section.area,
        "transverse": capped_section.transverse_second_moment,
        "longitudinal": capped_section.longitudinal_second_moment,
        "loops": capped_section.loops,
    }
    row["banded"] = {
        "area": banded_section.area,
        "transverse": banded_section.transverse_second_moment,
        "longitudinal": banded_section.longitudinal_second_moment,
        "loops": banded_section.loops,
    }

    capped_stability = initial_stability(capped_mesh, kg=KG, waterline_z=WATERLINE)
    banded_stability = initial_stability(banded_mesh, kg=KG, waterline_z=WATERLINE)
    row["capped"].update(
        {
            "volume": capped_stability.displacement_volume,
            "kb": capped_stability.kb,
            "bm": capped_stability.bm_transverse,
            "gm": capped_stability.gm_transverse,
        }
    )
    row["banded"].update(
        {
            "volume": banded_stability.displacement_volume,
            "kb": banded_stability.kb,
            "bm": banded_stability.bm_transverse,
            "gm": banded_stability.gm_transverse,
            "basis": banded_stability.basis,
        }
    )

    slopes = {}
    for heel in SLOPE_HEELS:
        curve = righting_arm_curve(
            banded_mesh, kg=KG, waterline_z=WATERLINE, heels_deg=(0.0, heel)
        )
        gz = curve.points[1].righting_arm
        slopes[str(heel)] = {
            "gz": gz,
            "slope": gz / math.sin(math.radians(heel)),
            "volume_error": abs(curve.points[1].displacement_volume - banded_stability.displacement_volume)
            / banded_stability.displacement_volume,
        }
    row["slopes"] = slopes
    return row


def make_plot(rows: list[dict[str, object]], reference: dict[str, float], path: Path) -> None:
    """Two panels: the GZ curve against ``GM sin(phi)``, and the mesh error budget."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    finest = rows[-1]
    parameters = WigleyParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=int(finest["nx"]), nz=int(finest["nz"]))
    banded = with_freeboard(generate_wigley(parameters).mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
    heels = [0.0, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5, 15.0]
    curve = righting_arm_curve(banded, kg=KG, waterline_z=WATERLINE, heels_deg=tuple(heels))
    gm = float(finest["banded"]["gm"])

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    radians = [math.radians(angle) for angle in heels]
    left.plot(heels, [point.righting_arm for point in curve.points], "o-", label="GZ from heeled volume integral")
    left.plot(heels, [gm * math.sin(value) for value in radians], "--", label=f"GM sin(phi), GM={gm:.6f} m")
    left.plot(heels, [float(finest["slopes"]["2.0"]["slope"]) * math.sin(value) for value in radians], ":", label="slope at 2 deg extrapolated")
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(f"Wigley + {FREEBOARD:g} m freeboard, {finest['nx']}x{finest['nz']}")
    left.legend(loc="upper left", fontsize=8)
    left.grid(alpha=0.3)

    cells = [int(row["cells"]) for row in rows]
    area_error = [relative(float(row["banded"]["area"]), reference["waterplane_area_m2"]) for row in rows]
    gm_error = [relative(float(row["banded"]["gm"]), reference["gm_transverse_m"]) for row in rows]
    seam_error = [
        relative(float(row["banded"]["area"]), float(row["capped"]["area"])) for row in rows
    ]
    right.loglog(cells, [1e2 * value for value in area_error], "s-", label="waterplane area vs analytic (x1e2)")
    right.loglog(cells, [1e2 * value for value in gm_error], "o-", label="GM vs analytic (x1e2)")
    right.loglog(cells, [max(value, 1e-17) for value in seam_error], "^-", label="banded vs capped (seam)")
    right.set_xlabel("cells nx*nz")
    right.set_ylabel("relative error")
    right.set_title("error budget under refinement")
    right.legend(fontsize=8)
    right.grid(alpha=0.3, which="both")

    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def write_report(rows, reference, claims, path: Path) -> None:
    lines = [
        "# The freeboard seam: a waterline plane that contains mesh edges",
        "",
        "Generated by `examples/freeboard_seam/run_example.py`; every number below is read",
        "from `results.json` in this directory.",
        "",
        "## Analytic reference (Wigley n = m = 2, L=10, B=1, T=0.625)",
        "",
        "| quantity | closed form | value |",
        "|---|---|---|",
        f"| displacement | 4LBT/9 | {reference['displacement_volume_m3']:.6f} m^3 |",
        f"| KB | 5T/8 | {reference['kb_m']:.6f} m |",
        f"| waterplane area | 2LB/3 | {reference['waterplane_area_m2']:.6f} m^2 |",
        f"| I_T | 4B^3L/105 | {reference['transverse_second_moment_m4']:.6f} m^4 |",
        f"| I_L (centroid) | BL^3/5 - A(L/2)^2 | {reference['longitudinal_second_moment_m4']:.6f} m^4 |",
        f"| BM_T | I_T/nabla | {reference['bm_transverse_m']:.6f} m |",
        f"| GM (KG=0.3) | KB+BM-KG | {reference['gm_transverse_m']:.6f} m |",
        "",
        "## Section through the banded mesh vs the capped hull",
        "",
        "| grid | capped A [m^2] | banded A [m^2] | rel diff | banded loops | banded I_T [m^4] |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {nx}x{nz} | {ca:.9f} | {ba:.9f} | {d:.2e} | {loops} | {it:.9f} |".format(
                nx=row["nx"],
                nz=row["nz"],
                ca=float(row["capped"]["area"]),
                ba=float(row["banded"]["area"]),
                d=relative(float(row["banded"]["area"]), float(row["capped"]["area"])),
                loops=int(row["banded"]["loops"]),
                it=float(row["banded"]["transverse"]),
            )
        )
    lines += [
        "",
        "## Metacentric height, banded mesh only",
        "",
        "| grid | nabla [m^3] | KB [m] | BM_T [m] | GM [m] | GM vs analytic |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        gm = float(row["banded"]["gm"])
        lines.append(
            "| {nx}x{nz} | {v:.6f} | {kb:.6f} | {bm:.6f} | {gm:.6f} | {e:+.3f}% |".format(
                nx=row["nx"],
                nz=row["nz"],
                v=float(row["banded"]["volume"]),
                kb=float(row["banded"]["kb"]),
                bm=float(row["banded"]["bm"]),
                gm=gm,
                e=100.0 * (gm - reference["gm_transverse_m"]) / reference["gm_transverse_m"],
            )
        )
    finest = rows[-1]
    lines += [
        "",
        "## Slope of the GZ curve against GM from the same mesh",
        "",
        f"finest grid {finest['nx']}x{finest['nz']}, GM = {float(finest['banded']['gm']):.6f} m",
        "",
        "| heel [deg] | GZ [m] | GZ/sin(phi) [m] | error vs GM | volume residual |",
        "|---|---|---|---|---|",
    ]
    for heel, values in finest["slopes"].items():
        slope = float(values["slope"])
        gm = float(finest["banded"]["gm"])
        lines.append(
            f"| {heel} | {float(values['gz']):.9f} | {slope:.9f} | "
            f"{100.0 * (slope - gm) / gm:+.4f}% | {float(values['volume_error']):.2e} |"
        )
    lines += [
        "",
        "## Claims",
        "",
        "| claim | value | limit | holds |",
        "|---|---|---|---|",
    ]
    for row in claims.rows:
        lines.append(
            f"| {row['claim']} | {float(row['value']):.3e} | {float(row['limit']):.3e} | "
            f"{'PASS' if row['holds'] else 'FAIL'} |"
        )
    lines += ["", f"Overall: {'PASS' if not claims.failures else 'FAIL'}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=Path("examples/freeboard_seam/outputs"))
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    reference = analytic_reference()
    rows = [measure(nx, nz) for nx, nz in GRIDS]

    claims = Claims()
    for row in rows:
        tag = f"{row['nx']}x{row['nz']}"
        claims.add(
            f"seam: banded area == capped area [{tag}]",
            relative(float(row["banded"]["area"]), float(row["capped"]["area"])),
            1e-14,
            "a vertical band cannot change the waterplane",
        )
        claims.add(
            f"seam: banded I_T == capped I_T [{tag}]",
            relative(float(row["banded"]["transverse"]), float(row["capped"]["transverse"])),
            1e-14,
            "second moment measured once, not twice",
        )
        claims.add(
            f"immersed body: banded nabla == capped nabla [{tag}]",
            relative(float(row["banded"]["volume"]), float(row["capped"]["volume"])),
            1e-9,
            "freeboard volume is not displacement",
        )
        claims.add(
            f"immersed body: banded GM == capped GM [{tag}]",
            relative(float(row["banded"]["gm"]), float(row["capped"]["gm"])),
            1e-9,
            "same loading, same metacentre",
        )
    for row in rows:
        tag = f"{row['nx']}x{row['nz']}"
        for heel, values in row["slopes"].items():
            gm = float(row["banded"]["gm"])
            limit = 0.001 * (float(heel) / 2.0) ** 2 + 0.0002
            claims.add(
                f"GZ slope vs GM at {heel} deg [{tag}]",
                relative(float(values["slope"]), gm),
                limit,
                "wall-sided BM/2 tan^2 term widens the tolerance with heel",
            )
            claims.add(
                f"displacement conserved at {heel} deg [{tag}]",
                float(values["volume_error"]),
                1e-10,
                "every heeled point floats at the upright displacement",
            )

    finest = rows[-1]
    claims.add(
        "convergence: GM reaches the closed form on the finest grid",
        relative(float(finest["banded"]["gm"]), reference["gm_transverse_m"]),
        0.001,
        f"finest grid {finest['nx']}x{finest['nz']}",
    )

    results = {
        "schema": "tensorcad.freeboard_seam/1",
        "inputs": {
            "length_m": LENGTH,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "freeboard_m": FREEBOARD,
            "kg_m": KG,
            "waterline_z_m": WATERLINE,
            "grids": [list(pair) for pair in GRIDS],
        },
        "analytic": reference,
        "grids": rows,
        "claims": claims.rows,
        "all_claims_hold": not claims.failures,
    }
    (arguments.output / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_report(rows, reference, claims, arguments.output / "report.md")
    try:
        make_plot(rows, reference, arguments.output / "plot.png")
    except ImportError:  # pragma: no cover - matplotlib is an examples-only dependency
        (arguments.output / "plot.png").unlink(missing_ok=True)
        print("matplotlib unavailable: plot.png skipped")

    for row in claims.rows:
        print(f"{'PASS' if row['holds'] else 'FAIL'}  {row['claim']}: {float(row['value']):.3e} (limit {float(row['limit']):.3e})")
    for row in rows:
        gm = float(row["banded"]["gm"])
        print(
            f"grid {row['nx']:>3}x{row['nz']:<3} area={float(row['banded']['area']):.6f} m^2 "
            f"I_T={float(row['banded']['transverse']):.6f} m^4 nabla={float(row['banded']['volume']):.6f} m^3 "
            f"KB={float(row['banded']['kb']):.6f} m GM={gm:.6f} m "
            f"(analytic {reference['gm_transverse_m']:.6f}, "
            f"{100.0 * (gm - reference['gm_transverse_m']) / reference['gm_transverse_m']:+.3f}%)"
        )
    print(f"all claims hold: {not claims.failures}  -> {arguments.output}")
    return 0 if not claims.failures else 1


if __name__ == "__main__":
    raise SystemExit(main())