#!/usr/bin/env python3
"""Verified benchmark: cross-tool displacement volume (XVAL-GEOM-VOLCROSS).

The decision this benchmark supports: **may the platform's displacement volume be
certified by its own integrator alone?**  The answer this campaign asserts is no --
the JBC ``-18.3 %`` episode is the standing counterexample -- and the acceptance
line copied from the platform's volume gate is 0.3 %.

What is measured, and against what
----------------------------------
Two sources of truth are used, and they are not the same kind of thing:

1. **Live rows (computed here)** -- the Wigley reference hull, whose displacement
   has the closed form ``V = L B T n m / ((n+1)(m+1)) = 4 L B T / 9``.  The closed
   triangle solid is rebuilt on a four-level ladder, and each level is measured by
   trimesh and numpy-stl (independent engines) and by a tessellation-free Simpson
   section quadrature.  The analytic value is the reference, so a deviation here is
   attributable: the engine rows must decay onto it at the measured order, and the
   quadrature row must hit it exactly (Simpson is exact for this polynomial field).
2. **Committed rows (re-validated, not transcribed)** -- the multi-hull evidence
   bundle in ``artifacts/xval_geom_volcross/results.json`` (Wigley, a
   KCS-dimensioned NURBS hull at full scale, and an asymmetric-transom integrator
   battery).  The KCS row's reference is the published displacement
   (52 030 m3, KRISO/SIMMAN); its caveat -- the form is calibrated to the published
   Cb, so the row's independent content is the three-route agreement -- is carried
   through into the table rather than dropped.

Gates (predeclared, module constants in ``tensorcad.volume_xval``):

* ``G1_engine_agreement`` -- the two independent engines, reading the same
  triangles, agree to :data:`tensorcad.volume_xval.ENGINE_AGREEMENT_TOLERANCE`;
* ``G2_reference_at_finest_wigley`` -- finest rung of the live ladder against the
  analytic value inside :data:`tensorcad.volume_xval.ACCEPTANCE_TOLERANCE`;
* ``G3_quadrature_is_exact_for_the_reference_hull`` -- the tessellation-free route
  reproduces the analytic Wigley volume to 1e-9 relative;
* ``G4_measured_order`` -- the observed order of the engine deficit falls inside
  ``ORDER_WINDOW`` (a tessellation deficit, not an unexplained bias);
* ``G5_committed_rows_hold`` -- every committed multi-hull row still passes its own
  reference band and still shows engine agreement.

Run from anywhere::

    python benchmarks/verified/volume_cross_engine/run_benchmark.py

Exit code 0 means every gate passed; a failing gate still writes its artifacts,
because an honest FAIL is evidence.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

CASEDIR = Path(__file__).resolve().parent
ROOT = CASEDIR.parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad import volume_xval as xv  # noqa: E402
from tensorcad import wigley as wg  # noqa: E402

LADDER = (21, 41, 81, 161)
COMMITTED = ROOT / "artifacts" / "xval_geom_volcross" / "results.json"


def _levels(stations: int) -> int:
    return max(3, stations // 4 + 1)


def live_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Wigley against its closed form, four rungs, three independent routes."""
    parameters = wg.reference_parameters()
    analytic = wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)

    def breadth(x: float, z: float) -> float:
        return wg.half_breadth(parameters, x, z)

    sectional = xv.volume_sectional(
        breadth, length=parameters.length, draft=parameters.draft, stations=401, levels=801
    )
    rows: list[dict[str, Any]] = []
    for stations in LADDER:
        levels = _levels(stations)
        solid = xv.closed_solid_from_breadth(
            breadth, length=parameters.length, draft=parameters.draft, stations=stations, levels=levels
        )
        topology = xv.topology_report(solid)
        engine_a = xv.volume_trimesh(solid)
        engine_b = xv.volume_numpy_stl(solid)
        rows.append(
            {
                "stations": stations,
                "levels": levels,
                "h_m": parameters.length / stations,
                "triangles": topology["faces"],
                "watertight": topology["watertight"],
                "closed_genus": topology["closed_genus"],
                "engine_a_trimesh_m3": engine_a,
                "engine_b_numpy_stl_m3": engine_b,
                "reference_analytic_m3": analytic,
                "engine_a_vs_reference": (engine_a - analytic) / analytic,
                "engine_b_vs_reference": (engine_b - analytic) / analytic,
                "engine_spread": abs(engine_a - engine_b) / abs(engine_a),
            }
        )
    deficits = [(row["h_m"], abs(analytic - row["engine_a_trimesh_m3"])) for row in rows]
    attribution = xv.classify_routes(
        divergence=xv.volume_divergence(
            xv.closed_solid_from_breadth(
                breadth,
                length=parameters.length,
                draft=parameters.draft,
                stations=LADDER[-1],
                levels=_levels(LADDER[-1]),
            )
        ),
        engine_a=rows[-1]["engine_a_trimesh_m3"],
        engine_b=rows[-1]["engine_b_numpy_stl_m3"],
        pipeline=sectional,
        deficits=deficits,
    )
    context = {
        "hull": "Wigley (n=m=2, ITTC reference proportions)",
        "principal_particulars": {
            "length_m": parameters.length,
            "beam_m": parameters.beam,
            "draft_m": parameters.draft,
        },
        "reference_volume_m3": analytic,
        "reference_kind": "analytic (exact)",
        "reference_source": "V = 4 L B T / 9 evaluated by tensorcad.wigley.analytic_volume",
        "sectional_volume_m3": sectional,
        "sectional_vs_reference": (sectional - analytic) / analytic,
        "attribution": attribution,
        "deficits": [{"h_m": h, "deficit_m3": d} for h, d in deficits],
    }
    return rows, context


def committed_rows() -> tuple[list[dict[str, Any]], list[str]]:
    """Re-validate the committed multi-hull bundle instead of quoting it."""
    if not COMMITTED.is_file():
        return [], [f"{COMMITTED.relative_to(ROOT)} is missing: run tools/volume_xval_report.py"]
    payload = json.loads(COMMITTED.read_text())
    rows: list[dict[str, Any]] = []
    problems: list[str] = []
    for case in payload.get("cases", []):
        for row in case.get("rows", []):
            engine_a = row.get("engine_a_trimesh_m3")
            engine_b = row.get("engine_b_numpy_stl_m3")
            reference = row.get("reference_volume_m3")
            if engine_a is None or engine_b is None:
                problems.append(f"{case.get('key')}: a row has no engine volume")
                continue
            spread = abs(engine_a - engine_b) / abs(engine_a)
            rows.append(
                {
                    "case": case.get("key"),
                    "name": case.get("name"),
                    "stations": row.get("stations"),
                    "levels": row.get("levels"),
                    "engine_a_trimesh_m3": engine_a,
                    "engine_b_numpy_stl_m3": engine_b,
                    "engine_spread": spread,
                    "reference_volume_m3": reference,
                    "reference_kind": row.get("reference_kind") or case.get("reference_kind"),
                    "engine_a_vs_reference": row.get("engine_a_vs_reference"),
                    "pipeline_volume_m3": row.get("pipeline_volume_m3"),
                    "pipeline_vs_reference": row.get("pipeline_vs_reference"),
                    "caveat": case.get("caveat", ""),
                }
            )
            if spread > xv.ENGINE_AGREEMENT_TOLERANCE:
                problems.append(f"{case.get('key')} @ {row.get('stations')}: engine spread {spread:.2e}")
    # The acceptance band belongs to the *finest* rung of each case: the coarse rungs
    # are the tessellation deficit itself, and demanding 0.3 % of a 21-station loft
    # would be demanding that refinement not matter.  What every rung must show is
    # engine agreement (checked above); what the finest must show is the band.
    for case_key in {row["case"] for row in rows}:
        case_rows = [
            row for row in rows if row["case"] == case_key and row.get("engine_a_vs_reference") is not None
        ]
        if not case_rows:
            continue
        finest = max(case_rows, key=lambda row: (row.get("stations") or 0, row.get("levels") or 0))
        if abs(finest["engine_a_vs_reference"]) > xv.ACCEPTANCE_TOLERANCE:
            problems.append(
                f"{case_key} @ finest ({finest['stations']}x{finest['levels']}): engine A vs "
                f"reference {100 * finest['engine_a_vs_reference']:+.3f} % outside the band"
            )
    return rows, problems


def figures(rows: list[dict[str, Any]], context: dict[str, Any], outdir: Path) -> list[dict[str, str]]:
    outdir.mkdir(parents=True, exist_ok=True)
    produced: list[dict[str, str]] = []

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    h = [row["h_m"] for row in rows]
    a = [abs(row["engine_a_vs_reference"]) * 100 for row in rows]
    b = [abs(row["engine_b_vs_reference"]) * 100 for row in rows]
    axes[0].loglog(h, a, "o-", label="trimesh (engine A)")
    axes[0].loglog(h, b, "s--", label="numpy-stl Mirtich (engine B)")
    axes[0].axhline(100 * xv.ACCEPTANCE_TOLERANCE, color="k", ls=":", label="0.3 % acceptance")
    axes[0].set_xlabel("h = L / stations [m]")
    axes[0].set_ylabel("|V_engine - V_analytic| [%]")
    axes[0].invert_xaxis()
    axes[0].set_title("Wigley: inscribed deficit vs analytic volume")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, which="both", alpha=0.3)

    spreads = [row["engine_spread"] for row in rows]
    axes[1].semilogy(h, spreads, "o-", label="measured engine spread")
    axes[1].axhline(xv.ENGINE_AGREEMENT_TOLERANCE, color="k", ls=":", label="1e-6 float32 budget")
    axes[1].set_xlabel("h = L / stations [m]")
    axes[1].set_ylabel("|A - B| / A")
    axes[1].invert_xaxis()
    order = context["attribution"]["observed_order"]
    axes[1].set_title(f"engines agree; observed deficit order {order:.2f}" if order else "engine agreement")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, which="both", alpha=0.3)
    path = outdir / "volume_cross_engine.png"
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    produced.append(
        {
            "name": path.name,
            "path": str(path.relative_to(ROOT)),
            "caption": "Left: the two independent engines' deficit against the analytic Wigley "
            "volume on a four-rung ladder, with the 0.3 % acceptance line; right: the engines' "
            "spread over the same ladder against the 1e-6 binary-STL float32 budget.",
        }
    )
    return produced


def main() -> int:
    started = time.perf_counter()
    rows, context = live_rows()
    committed, problems = committed_rows()
    gates = [
        {
            "name": "G1_engine_agreement",
            "value": max(row["engine_spread"] for row in rows),
            "limit": xv.ENGINE_AGREEMENT_TOLERANCE,
            "passed": max(row["engine_spread"] for row in rows) <= xv.ENGINE_AGREEMENT_TOLERANCE,
        },
        {
            "name": "G2_reference_at_finest_wigley",
            "value": abs(rows[-1]["engine_a_vs_reference"]),
            "limit": xv.ACCEPTANCE_TOLERANCE,
            "passed": abs(rows[-1]["engine_a_vs_reference"]) <= xv.ACCEPTANCE_TOLERANCE,
        },
        {
            "name": "G3_quadrature_is_exact_for_the_reference_hull",
            "value": abs(context["sectional_vs_reference"]),
            "limit": 1e-9,
            "passed": abs(context["sectional_vs_reference"]) <= 1e-9,
        },
        {
            "name": "G4_measured_order",
            "value": context["attribution"]["observed_order"],
            "limit": list(xv.ORDER_WINDOW),
            "passed": context["attribution"]["observed_order"] is not None
            and xv.ORDER_WINDOW[0] <= context["attribution"]["observed_order"] <= xv.ORDER_WINDOW[1],
        },
        {
            "name": "G5_committed_rows_hold",
            "value": problems,
            "limit": [],
            "passed": not problems,
        },
    ]
    figures_ = figures(rows, context, CASEDIR / "figures")
    payload = {
        "schema": "tensorcad.benchmark.volume_cross_engine/1.0",
        "feature": "XVAL-GEOM-VOLCROSS",
        "repro_cmd": "python benchmarks/verified/volume_cross_engine/run_benchmark.py",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "trimesh": _version("trimesh"),
            "numpy_stl": _version("stl"),
        },
        "tolerances": {
            "acceptance": xv.ACCEPTANCE_TOLERANCE,
            "engine_agreement": xv.ENGINE_AGREEMENT_TOLERANCE,
            "order_window": list(xv.ORDER_WINDOW),
        },
        "live": {"context": context, "rows": rows},
        "committed": {"source": str(COMMITTED.relative_to(ROOT)), "rows": committed, "problems": problems},
        "gates": gates,
        "passed": all(gate["passed"] for gate in gates),
        "figures": figures_,
        "seconds": time.perf_counter() - started,
    }
    CASEDIR.mkdir(parents=True, exist_ok=True)
    (CASEDIR / "result.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    (CASEDIR / "gates_verified.json").write_text(
        json.dumps({"feature": payload["feature"], "gates": gates, "passed": payload["passed"]}, indent=1),
        encoding="utf-8",
    )
    for gate in gates:
        state = "PASS" if gate["passed"] else "FAIL"
        print(f"[{state}] {gate['name']}: {gate['value']} (limit {gate['limit']})")
    print(f"live ladder (Wigley, analytic {context['reference_volume_m3']:.6f} m3):")
    for row in rows:
        print(
            f"  {row['stations']:4d} x {row['levels']:3d}  A {row['engine_a_trimesh_m3']:.6f} "
            f"({100 * row['engine_a_vs_reference']:+.4f} %)  B {row['engine_b_numpy_stl_m3']:.6f} "
            f"({100 * row['engine_b_vs_reference']:+.4f} %)  spread {row['engine_spread']:.2e}"
        )
    print(f"committed rows re-validated: {len(committed)}, problems {len(problems)}")
    return 0 if payload["passed"] else 1


def _version(module: str) -> str:
    try:
        return __import__(module).__version__
    except Exception:  # pragma: no cover - environment dependent
        return "absent"


if __name__ == "__main__":
    raise SystemExit(main())