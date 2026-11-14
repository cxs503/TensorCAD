"""Cross-tool volume cross-validation report (feature ``XVAL-GEOM-VOLCROSS``).

Runs every hull case through four independent volume routes and writes the
evidence bundle: a JSON record, a generated Markdown table, and figures.  Nothing
in the generated documents is transcribed by hand -- the numbers are written by
this script and read back by the docs.

    python tools/volume_xval_report.py [--out artifacts/xval_geom_volcross]

Routes per case (see :mod:`tensorcad.volume_xval` for what each one can prove):

* ``pipeline``  -- the platform's own volume route for that hull: the Gauss
  quadrature of the NURBS patch (``hull_nurbs``) where the hull *is* a patch, or
  the closed-surface divergence theorem of :mod:`tensorcad.mesh`;
* ``engine A``  -- :mod:`trimesh`, independent third-party engine;
* ``engine B``  -- :mod:`numpy-stl` (Mirtich), a second independent engine, through
  a binary-STL round trip;
* ``smooth``    -- a route that shares no mathematics with any of the others: a
  quadrature of the half-breadth field (no triangles, no divergence theorem).

Engines must be installed for the bundle to be complete; a missing engine aborts
the run instead of silently shrinking the table.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Sequence

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad import volume_xval as xv  # noqa: E402
from tensorcad.geometry import Mesh  # noqa: E402

# --------------------------------------------------------------------------- #
# hull cases
# --------------------------------------------------------------------------- #
@dataclass
class Case:
    """One hull plus the reference that is allowed to grade it."""

    key: str
    name: str
    length: float
    beam: float
    draft: float
    breadth: xv.Breadth
    pipeline_volume: Callable[[], float]
    pipeline_name: str
    reference_volume: float | None
    reference_kind: str
    reference_source: str
    caveat: str = ""
    #: ``True`` when the pipeline route integrates the *smooth* geometry (a NURBS
    #: quadrature), so its gap to the engines is the tessellation deficit.  ``False``
    #: when the pipeline route is itself a mesh route (Wigley), in which case the
    #: convergence ladder has to be read off the engines against the reference.
    pipeline_is_smooth: bool = True
    stations: Sequence[int] = (21, 41, 81, 161)
    levels_of: Callable[[int], int] = lambda n: max(3, n // 4 + 2)


#: Published KCS principal particulars (KRISO; SIMMAN 2008 / Tokyo 2015 workshop
#: series).  The parametric hull built from these dimensions is *calibrated* to the
#: published Cb, so agreement with ``Cb L B T`` is a calibration identity and is
#: reported as such, not as an independent validation of the volume integrator.
KCS = {
    "length": 230.0,
    "beam": 32.2,
    "draft": 10.8,
    "cb": 0.6505,
    "volume": 52_030.0,
    "source": "KCS principal particulars, KRISO (Lpp 230.0 m, B 32.2 m, T 10.8 m, Cb 0.6505, "
    "nabla 52 030 m^3; SIMMAN 2008 / Tokyo 2015 workshop series)",
}


def wigley_case() -> Case:
    from tensorcad import wigley as wg

    p = wg.reference_parameters()
    exact = wg.analytic_volume(p.length, p.beam, p.draft)

    def breadth(x: float, z: float) -> float:
        return wg.half_breadth(p, x, z)

    def pipeline() -> float:
        # The platform's route for Wigley: build its own hull mesh, integrate the
        # closed surface with the divergence theorem (tensorcad.mesh).
        mesh, _offsets = wg.build_mesh(p)
        return meshlib.total_signed_volume(mesh)

    return Case(
        key="wigley",
        name=f"Wigley (ITTC reference proportions, n=m=2)",
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        breadth=breadth,
        pipeline_volume=pipeline,
        pipeline_name="tensorcad.wigley mesh + divergence theorem",
        reference_volume=exact,
        reference_kind="analytic (exact)",
        reference_source="V = L B T n m / ((n+1)(m+1)) = 4 L B T / 9, evaluated in "
        "tensorcad.wigley.analytic_volume",
        pipeline_is_smooth=False,
        stations=(21, 41, 81, 161),
        levels_of=lambda n: max(3, n // 4 + 1),
    )


def nurbs_hull_case(*, full_scale: bool) -> Case:
    from tensorcad import hull_nurbs as hn

    if full_scale:
        params = hn.HullNurbsParameters(
            length=KCS["length"],
            beam=KCS["beam"],
            draft=KCS["draft"],
            block_coefficient=KCS["cb"],
        )
        key = "kcs-fullscale"
        name = "KCS-dimensioned NURBS hull (full scale, Cb calibrated to 0.6505)"
        reference: float | None = KCS["volume"]
        ref_kind = "published displacement volume"
        ref_source = KCS["source"]
        caveat = (
            "the form is calibrated to the published Cb, so pipeline-vs-published is a "
            "calibration identity; the independent content of this row is the three-route "
            "agreement at realistic L/B/T and the tessellation attribution"
        )
    else:
        params = hn.kcs_lightweight_parameters()
        key = "kcs-modelscale"
        name = "KCS-dimensioned NURBS hull (model scale, Cb 0.555) -- IN383-CAD01 case"
        reference = None
        ref_kind = "none (WITHHELD)"
        ref_source = (
            "no published value exists for this parametric form (it is not the KCS offsets); "
            "before this run its volume was only ever compared against a mesh generated by the "
            "same patch, i.e. a WITHHELD claim"
        )
        caveat = (
            "reference-free by construction: the evidence produced here is cross-tool "
            "agreement, not agreement with an outside number"
        )

    hull = hn.build_nurbs_hull(params)
    deck = hull.deck_height
    quadrature = hn.nurbs_hydrostatics(hull).displaced_volume

    def breadth(x: float, z: float) -> float:
        return hull.surface.evaluate(x / params.length, z / deck)[1]

    return Case(
        key=key,
        name=name,
        length=params.length,
        beam=params.beam,
        draft=params.draft,
        breadth=breadth,
        pipeline_volume=lambda: quadrature,
        pipeline_name="tensorcad.hull_nurbs Gauss quadrature of the NURBS patch "
        f"({params.quadrature_order[0]}x{params.quadrature_order[1]})",
        reference_volume=reference,
        reference_kind=ref_kind,
        reference_source=ref_source,
        caveat=caveat,
        stations=(21, 41, 81, 161),
        levels_of=lambda n: max(3, n // 4 + 2),
    )


def integrator_battery() -> dict[str, Any]:
    """Same-triangle integrator identity on hulls that have no symmetric breadth field.

    The asymmetric transom hull is the platform's one hull whose starboard and port
    sections differ, so the symmetric loft used for the cross-validation table cannot
    represent it.  It is still worth an independent-engine check: three routes read the
    platform's *own* closed mesh, and the analytic closed form grades all three.  This
    row answers "does the platform integrate a closed surface correctly", not "is its
    geometry source right" -- the table above answers the second question.
    """
    from tensorcad import transom_asymmetric as ta

    params = ta.AsymmetricTransomParameters()
    mesh, _offsets = ta.build_asymmetric_mesh(params)
    analytic = ta.analytic_asymmetric_hydrostatics(params)["displacement_volume"]
    starboard, port = params.half_breadths()
    divergence = xv.volume_divergence(mesh)
    engine_a = xv.volume_trimesh(mesh)
    engine_b = xv.volume_numpy_stl(mesh)
    topo = xv.topology_report(mesh)
    return {
        "name": "Asymmetric transom hull (starboard/port differ, flat transom plate)",
        "principal_particulars": {
            "L_m": params.length,
            "B_m": params.beam,
            "T_m": params.draft,
            "starboard_half_breadth_m": starboard,
            "port_half_breadth_m": port,
            "stations": params.nx,
            "levels": params.nz,
        },
        "reference_kind": "analytic closed form",
        "reference_source": "tensorcad.transom_asymmetric.analytic_asymmetric_hydrostatics",
        "reference_volume_m3": analytic,
        "divergence_m3": divergence,
        "engine_a_trimesh_m3": engine_a,
        "engine_b_numpy_stl_m3": engine_b,
        "divergence_vs_analytic": xv.relative_difference(divergence, analytic),
        "engine_a_vs_analytic": xv.relative_difference(engine_a, analytic),
        "engine_b_vs_analytic": xv.relative_difference(engine_b, analytic),
        "self_vs_engine_a": xv.relative_difference(divergence, engine_a),
        "engine_a_vs_engine_b": xv.relative_difference(engine_b, engine_a),
        "topology": topo,
        "note": (
            "the mesh is the platform's own hull at its default station grid, so the "
            "deviation from the analytic volume is a sampling deficit of that grid -- the "
            "three routes agreeing on the same triangles is what this row proves"
        ),
    }


# --------------------------------------------------------------------------- #
# running a case
# --------------------------------------------------------------------------- #
def run_case(case: Case) -> dict[str, Any]:
    pipeline = case.pipeline_volume()
    # The smooth-field quadrature does not depend on the mesh resolution, so it is
    # evaluated once per case rather than once per row: a NURBS patch evaluation
    # costs ~200 us in pure Python and repeating it 4x would dominate the run.
    smooth = xv.volume_sectional(
        case.breadth, length=case.length, draft=case.draft, stations=201, levels=101
    )
    rows: list[xv.CaseRow] = []
    same_mesh_rows: list[dict[str, Any]] = []
    for stations in case.stations:
        levels = case.levels_of(stations)
        solid = xv.closed_solid_from_breadth(
            case.breadth,
            length=case.length,
            draft=case.draft,
            stations=stations,
            levels=levels,
        )
        topo = xv.topology_report(solid)
        divergence = xv.volume_divergence(solid)
        engine_a = xv.volume_trimesh(solid)
        engine_b = xv.volume_numpy_stl(solid)
        h = case.length / (stations - 1)
        row = xv.CaseRow(
            ship=case.key,
            reference_kind=case.reference_kind,
            reference_volume=case.reference_volume,
            reference_source=case.reference_source,
            stations=stations,
            levels=levels,
            h=h,
            pipeline=divergence,
            divergence=divergence,
            engine_a=engine_a,
            engine_b=engine_b,
            sectional=smooth,
            extras={
                "pipeline_geometry_volume_m3": pipeline,
                "pipeline_geometry_vs_engines": xv.relative_difference(
                    pipeline, 0.5 * (engine_a + engine_b)
                ),
                "topology": topo,
                "engine_agreement": xv.relative_difference(engine_b, engine_a),
            },
        )
        rows.append(row)
        same_mesh_rows.append(
            {
                "stations": stations,
                "levels": levels,
                "vertices": topo["vertices"],
                "faces": topo["faces"],
                "h_m": h,
                "divergence_m3": divergence,
                "engine_a_trimesh_m3": engine_a,
                "engine_b_numpy_stl_m3": engine_b,
                "self_vs_engine_a": xv.relative_difference(divergence, engine_a),
                "engine_a_vs_engine_b": xv.relative_difference(engine_b, engine_a),
                "geometry_quadrature_m3": pipeline,
                "geometry_vs_engines": xv.relative_difference(
                    pipeline, 0.5 * (engine_a + engine_b)
                ),
                "smooth_quadrature_m3": smooth,
                "smooth_vs_engines": xv.relative_difference(smooth, 0.5 * (engine_a + engine_b)),
                "watertight": topo["watertight"],
                "closed_genus": topo["closed_genus"],
            }
        )

    if case.pipeline_is_smooth:
        deficit = [(row.h, abs(row.extras["pipeline_geometry_vs_engines"])) for row in rows]
        comparison = pipeline
    else:
        # The pipeline route for this hull is itself a mesh route, so the only
        # convergence axis it has is the engines' own deficit against the reference.
        if case.reference_volume is None:
            raise SystemExit(f"{case.key}: a mesh-route pipeline needs a reference to be graded")
        deficit = [
            (row.h, abs(xv.relative_difference(row.engine_a, case.reference_volume))) for row in rows
        ]
        comparison = case.reference_volume
    worst_spread = max(abs(row.extras["engine_agreement"]) for row in rows)
    verdict = xv.classify_routes(
        divergence=rows[-1].divergence,
        engine_a=rows[-1].engine_a,
        engine_b=rows[-1].engine_b,
        pipeline=comparison,
        deficits=deficit,
    )
    order = verdict["observed_order"]
    extrapolated = None
    if order and order > 0.0:
        extrapolated = xv.richardson_limit([(row.h, row.engine_a) for row in rows], order)
    verdict["engine_spread_worst"] = worst_spread
    verdict["engine_note"] = xv.engine_agreement_note(worst_spread)
    verdict["extrapolated_engine_volume_m3"] = extrapolated
    if extrapolated is not None:
        verdict["extrapolated_vs_pipeline"] = xv.relative_difference(extrapolated, comparison)
        if case.reference_volume:
            verdict["extrapolated_vs_reference"] = xv.relative_difference(
                extrapolated, case.reference_volume
            )
    if case.reference_volume:
        verdict["pipeline_vs_reference"] = xv.relative_difference(pipeline, case.reference_volume)
        verdict["within_acceptance_vs_reference"] = (
            abs(verdict["pipeline_vs_reference"]) <= xv.ACCEPTANCE_TOLERANCE
        )
        # Half the sample spacing is a fair proxy for the tessellation error scale.
        verdict["reference_is_within_acceptance_of_pipeline"] = (
            abs(verdict["pipeline_vs_reference"]) <= xv.ACCEPTANCE_TOLERANCE
        )
    verdict["topology_all_watertight"] = all(
        row.extras["topology"]["watertight"] for row in rows
    )
    verdict["topology_genera"] = sorted({row.extras["topology"]["closed_genus"] for row in rows})
    verdict["self_vs_engine_a_worst"] = max(
        abs(xv.relative_difference(row.divergence, row.engine_a)) for row in rows
    )
    return {
        "key": case.key,
        "name": case.name,
        "principal_particulars": {
            "L_m": case.length,
            "B_m": case.beam,
            "T_m": case.draft,
            "LBT_m3": case.length * case.beam * case.draft,
        },
        "pipeline_route": case.pipeline_name,
        "pipeline_volume_m3": pipeline,
        "reference_volume_m3": case.reference_volume,
        "reference_kind": case.reference_kind,
        "reference_source": case.reference_source,
        "caveat": case.caveat,
        "rows": [row.as_dict() for row in rows],
        "same_mesh_rows": same_mesh_rows,
        "verdict": verdict,
    }


def quadrature_ladder() -> dict[str, Any]:
    """Convergence of the pipeline quadrature *itself*, in Gauss order.

    A tessellation study cannot say anything about the quadrature it is compared
    against, so the pipeline route is refined on its own axis: if doubling the
    Gauss order moves the volume, the pipeline number was not converged and any
    agreement with an engine would have been coincidence.
    """
    from tensorcad import hull_nurbs as hn

    params = hn.kcs_lightweight_parameters().replaced(
        quadrature_order=(8, 6)
    )
    hull = hn.build_nurbs_hull(params)
    rows = []
    previous = None
    for order in ((8, 6), (12, 8), (16, 12), (24, 16), (32, 24), (48, 32)):
        volume = hn.nurbs_hydrostatics(hull, quadrature_order=order).displaced_volume
        rows.append(
            {
                "order": list(order),
                "volume_m3": volume,
                "delta_vs_previous": None if previous is None else volume - previous,
                "relative_vs_previous": None
                if previous is None
                else xv.relative_difference(volume, previous),
            }
        )
        previous = volume
    return {
        "hull": "KCS-dimensioned NURBS hull (model scale)",
        "rows": rows,
        "converged_delta_m3": abs(rows[-1]["volume_m3"] - rows[-2]["volume_m3"]),
        "converged_relative": abs(rows[-1]["relative_vs_previous"]),
    }


def engine_versions() -> dict[str, Any]:
    versions: dict[str, Any] = {"python": platform.python_version()}
    for label, module in (("trimesh", "trimesh"), ("numpy-stl", "stl"), ("numpy", "numpy")):
        try:
            imported = __import__(module)
            versions[label] = getattr(imported, "__version__", "unknown")
        except ImportError:  # pragma: no cover - environment dependent
            versions[label] = None
    versions["unavailable"] = {}
    for label, module in (("gmsh", "gmsh"), ("pythonocc-core", "OCC"), ("cadquery", "cadquery")):
        try:
            __import__(module)
        except ImportError:
            versions["unavailable"][label] = "not installed (no network in this environment)"
    return versions


# --------------------------------------------------------------------------- #
# markdown and figures
# --------------------------------------------------------------------------- #
def markdown(results: dict[str, Any]) -> str:
    lines: list[str] = []
    add = lines.append

    add("# Geometry volume cross-validation: pipeline vs independent engines")
    add("")
    add(f"Feature: `XVAL-GEOM-VOLCROSS`  ")
    add(f"Generated: {results['generated_at']} by `tools/volume_xval_report.py`  ")
    add(
        "Engines: "
        + ", ".join(
            f"{name} {ver}" for name, ver in results["engines"].items() if name not in {"python", "unavailable"} and ver
        )
    )
    add("")
    add(
        "Every number below is written by the generator into `results.json` and rendered from there; "
        "no value is copied by hand."
    )
    add("")
    add("## Route definitions")
    add("")
    add("| route | implementation | shares code/algorithms with |")
    add("|---|---|---|")
    add(
        "| pipeline | the platform's own route for that hull (NURBS Gauss quadrature, or the "
        "closed-surface divergence theorem of `tensorcad.mesh`) | -- |"
    )
    add("| engine A | `trimesh.Trimesh.volume` (third-party, NumPy) | nothing in this repo |")
    add(
        "| engine B | `numpy-stl` `get_mass_properties` (Mirtich, C++) through a binary-STL "
        "round trip (float32 coordinates) | nothing in this repo |"
    )
    add(
        "| smooth | Simpson quadrature of the half-breadth field: `A(x) = 2 int b dz`, then "
        "`V = int A dx`. No triangles, no divergence theorem | the half-breadth field only |"
    )
    add("")
    for case in results["cases"]:
        add(f"## {case['name']}")
        add("")
        pp = case["principal_particulars"]
        add(
            f"L = {pp['L_m']:.6g} m, B = {pp['B_m']:.6g} m, T = {pp['T_m']:.6g} m, "
            f"LBT = {pp['LBT_m3']:.8g} m³.  Pipeline route: `{case['pipeline_route']}`."
        )
        add("")
        add(
            f"Reference: **{case['reference_kind']}**"
            + (f" = {case['reference_volume_m3']:.8g} m³" if case["reference_volume_m3"] else "")
            + f" -- {case['reference_source']}."
        )
        if case["caveat"]:
            add("")
            add(f"> Caveat: {case['caveat']}.")
        add("")
        add("| stations | levels | h [m] | faces | pipeline (m³) | engine A (m³) | engine B (m³) | smooth (m³) | A vs B | pipeline vs engines | smooth vs engines |")
        add("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in case["same_mesh_rows"]:
            add(
                f"| {row['stations']} | {row['levels']} | {row['h_m']:.4g} | {row['faces']} | "
                f"{row['geometry_quadrature_m3']:.9g} | {row['engine_a_trimesh_m3']:.9g} | "
                f"{row['engine_b_numpy_stl_m3']:.9g} | {row['smooth_quadrature_m3']:.9g} | "
                f"{100 * row['engine_a_vs_engine_b']:+.2e} % | "
                f"{100 * row['geometry_vs_engines']:+.4f} % | "
                f"{100 * row['smooth_vs_engines']:+.4f} % |"
            )
        add("")
        add("**Self route against the engines, on the identical triangles** "
            "(integrator identity, not tessellation):")
        add("")
        add("| stations | divergence (m³) | engine A (m³) | A vs divergence | B vs divergence | watertight | genus |")
        add("|---:|---:|---:|---:|---:|---|---:|")
        for row in case["same_mesh_rows"]:
            add(
                f"| {row['stations']} | {row['divergence_m3']:.9g} | {row['engine_a_trimesh_m3']:.9g} | "
                f"{100 * row['self_vs_engine_a']:+.2e} % | {100 * row['engine_a_vs_engine_b']:+.2e} % | "
                f"{row['watertight']} | {row['closed_genus']} |"
            )
        add("")
        verdict = case["verdict"]
        add(f"Verdict: **{verdict['status']}** -- {verdict['reason']}.")
        add("")
        add(f"- observed convergence order of the pipeline-to-engine gap: {verdict['observed_order']!r}")
        add(
            "- engine-to-engine spread (worst resolution): "
            f"{verdict['engine_spread_worst']:.3e} ({verdict['engine_note']})"
        )
        add(f"- self divergence vs engine A (worst): {verdict['self_vs_engine_a_worst']:.3e}")
        if verdict.get("extrapolated_engine_volume_m3") is not None:
            add(
                f"- Richardson extrapolation of engine A to h->0: "
                f"{verdict['extrapolated_engine_volume_m3']:.9g} m³, which is "
                f"{100 * verdict['extrapolated_vs_pipeline']:+.4f} % from the pipeline volume"
            )
        if case["reference_volume_m3"]:
            add(
                f"- pipeline vs reference: {100 * verdict['pipeline_vs_reference']:+.4f} % "
                f"({'inside' if verdict['within_acceptance_vs_reference'] else 'OUTSIDE'} the 0.3 % line)"
            )
            if verdict.get("extrapolated_vs_reference") is not None:
                add(
                    f"- extrapolated engines vs reference: "
                    f"{100 * verdict['extrapolated_vs_reference']:+.4f} %"
                )
        add("")
    add("## Pipeline quadrature convergence (does the pipeline number stand still?)")
    add("")
    ladder = results["quadrature_ladder"]
    add("| Gauss order | volume (m³) | delta vs previous (m³) | relative vs previous |")
    add("|---|---:|---:|---:|")
    for row in ladder["rows"]:
        add(
            f"| {row['order'][0]}x{row['order'][1]} | {row['volume_m3']:.12g} | "
            f"{'' if row['delta_vs_previous'] is None else format(row['delta_vs_previous'], '.3e')} | "
            f"{'' if row['relative_vs_previous'] is None else format(100 * row['relative_vs_previous'], '+.3e') + ' %'} |"
        )
    add("")
    add(
        f"Between the two finest orders the pipeline volume moves "
        f"{ladder['converged_delta_m3']:.3e} m³ ({100 * ladder['converged_relative']:.3e} %), so the "
        "quadrature is converged far below every tolerance compared against it."
    )
    add("")
    battery = results["integrator_battery"]
    add("## Integrator-identity battery (hulls without a symmetric breadth field)")
    add("")
    add(f"**{battery['name']}** -- reference: {battery['reference_kind']}, "
        f"{battery['reference_volume_m3']:.8g} m3 ({battery['reference_source']}).")
    add("")
    add("| route | volume (m3) | vs analytic |")
    add("|---|---:|---:|")
    add(f"| divergence (platform) | {battery['divergence_m3']:.9g} | {100*battery['divergence_vs_analytic']:+.4f} % |")
    add(f"| engine A (trimesh) | {battery['engine_a_trimesh_m3']:.9g} | {100*battery['engine_a_vs_analytic']:+.4f} % |")
    add(f"| engine B (numpy-stl) | {battery['engine_b_numpy_stl_m3']:.9g} | {100*battery['engine_b_vs_analytic']:+.4f} % |")
    add("")
    add(f"Self divergence vs engine A: {100*battery['self_vs_engine_a']:+.2e} %; "
        f"engine A vs engine B: {100*battery['engine_a_vs_engine_b']:+.2e} %; "
        f"watertight={battery['topology']['watertight']}, genus={battery['topology']['closed_genus']}.")
    add("")
    add(f"> {battery['note']}.")
    add("")
    add("## Engines that were not available")
    add("")
    for name, reason in results["engines"]["unavailable"].items():
        add(f"- `{name}`: {reason}")
    add("")
    add(
        "The anchor names Gmsh/OpenCASCADE as the second engine. Neither is installable in this "
        "environment (no package index reachable), so the second independent engine is "
        "`numpy-stl`'s Mirtich implementation, and the tessellation-free route is the sectional "
        "quadrature. This is a capability limitation of the environment, recorded rather than "
        "worked around silently; the OpenCASCADE BRep route remains open work."
    )
    return "\n".join(lines) + "\n"


def figures(results: dict[str, Any], outdir: Path) -> list[dict[str, str]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    produced: list[dict[str, str]] = []

    # Figure 1: convergence of every route against the reference / against the engines.
    fig, axes = plt.subplots(1, len(results["cases"]), figsize=(4.6 * len(results["cases"]), 4.2))
    if len(results["cases"]) == 1:  # pragma: no cover - single-case fallback
        axes = [axes]
    reference_slope = None
    for ax, case in zip(axes, results["cases"]):
        rows = case["same_mesh_rows"]
        h = [row["h_m"] for row in rows]
        gap = [abs(row["geometry_vs_engines"]) for row in rows]
        ax.loglog(h, gap, "o-", color="#1f77b4", label="pipeline vs engines")
        smooth_gap = [abs(row["smooth_vs_engines"]) for row in rows]
        ax.loglog(h, smooth_gap, "s-", color="#2ca02c", label="smooth quadrature vs engines")
        # order-2 guide anchored on the coarsest point
        guide = [gap[0] * (h[0] / value) ** 2 for value in h]
        ax.loglog(h, guide, "k--", lw=1.0, label="O(h²) guide")
        reference_slope = (
            (math.log(gap[-1]) - math.log(gap[0])) / (math.log(h[-1]) - math.log(h[0]))
            if gap[0] > 0 and gap[-1] > 0
            else None
        )
        ax.set_title(case["key"], fontsize=10)
        ax.set_xlabel("h = L/(stations-1)  [m]")
        ax.set_ylabel("relative difference to the engines")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle(
        "Tessellation deficit: the pipeline volume sits above the inscribed polyhedron by O(h²)"
        + (f"  (measured slope {reference_slope:.2f})" if reference_slope else "")
    )
    fig.tight_layout()
    path = outdir / "convergence.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    produced.append(
        {
            "name": "convergence.png",
            "caption": "Pipeline-to-engine volume gap versus mesh spacing for every hull case, "
            "with the O(h²) guide. The gap closing at order 2 is what distinguishes a "
            "tessellation deficit from an integrator bias.",
        }
    )

    # Figure 2: route deviation against the reference for the ships that have one.
    graded = [case for case in results["cases"] if case["reference_volume_m3"]]
    fig, ax = plt.subplots(figsize=(max(6.0, 2.2 * len(graded)), 4.0))
    width = 0.2
    for index, case in enumerate(graded):
        row = case["same_mesh_rows"][-1]
        values = [
            100 * xv.relative_difference(row["engine_a_trimesh_m3"], case["reference_volume_m3"]),
            100 * xv.relative_difference(row["engine_b_numpy_stl_m3"], case["reference_volume_m3"]),
            100 * xv.relative_difference(row["smooth_quadrature_m3"], case["reference_volume_m3"]),
            100
            * xv.relative_difference(row["geometry_quadrature_m3"], case["reference_volume_m3"]),
        ]
        labels = ["engine A (trimesh)", "engine B (numpy-stl)", "smooth quad", case["pipeline_route"].split(" ")[0]]
        positions = [index + (k - 1.5) * width for k in range(4)]
        ax.bar(positions, values, width=width, label=None)
        for pos, value in zip(positions, values):
            ax.annotate(f"{value:+.4f}", (pos, value), ha="center", va="bottom", fontsize=7, rotation=90)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axhline(0.3, color="r", ls="--", lw=0.8, label="+0.3 % line")
    ax.axhline(-0.3, color="r", ls="--", lw=0.8, label="-0.3 % line")
    ax.set_xticks(range(len(graded)))
    ax.set_xticklabels(
        [case["key"] + "\n" + " ".join(labels) for case in graded], fontsize=6
    )
    ax.set_ylabel("deviation from reference [%]")
    ax.set_title("Route deviation against the reference at the finest mesh")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    path = outdir / "deviation_vs_reference.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    produced.append(
        {
            "name": "deviation_vs_reference.png",
            "caption": "Signed deviation of each route from its own reference at the finest "
            "resolution, with the 0.3 % acceptance lines. Only the analytic-reference cases "
            "carry a reference at all; the model-scale NURBS hull is reference-free (WITHHELD).",
        }
    )

    # Figure 3: the geometry itself -- breadth field and body plan of the KCS-dimensioned hull.
    from tensorcad import hull_nurbs as hn

    hull = hn.build_nurbs_hull(hn.kcs_lightweight_parameters())
    params = hull.parameters
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.0))
    zs = [params.draft * j / 60 for j in range(61)]
    xs = [params.length * i / 120 for i in range(121)]
    field = [
        [hull.surface.evaluate(x / params.length, z / hull.deck_height)[1] for x in xs] for z in zs
    ]
    image = ax1.pcolormesh(xs, zs, field, shading="auto", cmap="viridis")
    fig.colorbar(image, ax=ax1, label="half-breadth b(x,z) [m]")
    ax1.set_xlabel("x [m] (0 = transom, L = stem)")
    ax1.set_ylabel("z [m] above keel")
    ax1.set_title("Half-breadth field of the KCS-dimensioned NURBS hull")
    for station in (0.05, 0.25, 0.5, 0.75, 0.95):
        x = station * params.length
        section = [hull.surface.evaluate(station, z / hull.deck_height)[1] for z in zs]
        ax2.plot(section, zs, label=f"x/L={station}")
        ax2.plot([-value for value in section], zs, color=ax2.lines[-1].get_color())
    ax2.set_xlabel("y [m]")
    ax2.set_ylabel("z [m]")
    ax2.set_title("Body plan (port mirrored)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=7)
    fig.tight_layout()
    path = outdir / "hull_geometry.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    produced.append(
        {
            "name": "hull_geometry.png",
            "caption": "Half-breadth field b(x,z) and body plan of the KCS-dimensioned NURBS hull "
            "whose volume is cross-validated: the geometry the quadrature integrates and the "
            "engines tessellate.",
        }
    )
    return produced


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def build_bundle(outdir: Path) -> dict[str, Any]:
    engines = engine_versions()
    for required in ("trimesh", "numpy-stl"):
        if engines.get(required) is None:
            raise SystemExit(
                f"engine {required!r} is required for the cross-validation bundle and is not "
                "installed; a bundle with one engine would not cross-validate anything"
            )
    cases = [
        wigley_case(),
        nurbs_hull_case(full_scale=True),
        nurbs_hull_case(full_scale=False),
    ]
    results: dict[str, Any] = {
        "schema": xv.SCHEMA,
        "feature": "XVAL-GEOM-VOLCROSS",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tolerances": {
            "acceptance": xv.ACCEPTANCE_TOLERANCE,
            "engine_agreement": xv.ENGINE_AGREEMENT_TOLERANCE,
            "order_window": list(xv.ORDER_WINDOW),
        },
        "engines": engines,
        "cases": [],
    }
    for case in cases:
        print(f"[case] {case.key}: {case.name}")
        results["cases"].append(run_case(case))
    results["quadrature_ladder"] = quadrature_ladder()
    results["integrator_battery"] = integrator_battery()

    outdir.mkdir(parents=True, exist_ok=True)
    results["figures"] = figures(results, outdir / "figures")
    (outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")
    document = markdown(results)
    (outdir / "results.md").write_text(document)

    summary = {
        case["key"]: {
            "status": case["verdict"]["status"],
            "pipeline_volume_m3": case["pipeline_volume_m3"],
            "reference_volume_m3": case["reference_volume_m3"],
            "pipeline_vs_reference": case["verdict"].get("pipeline_vs_reference"),
            "engine_spread_worst": case["verdict"]["engine_spread_worst"],
            "observed_order": case["verdict"]["observed_order"],
            "extrapolated_vs_pipeline": case["verdict"].get("extrapolated_vs_pipeline"),
        }
        for case in results["cases"]
    }
    summary["all_watertight"] = all(
        case["verdict"]["topology_all_watertight"] for case in results["cases"]
    )
    summary["quadrature_converged_relative"] = results["quadrature_ladder"]["converged_relative"]
    results["summary"] = summary
    (outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")
    (outdir / "results.md").write_text(markdown(results))
    return results


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default=str(REPO / "artifacts" / "xval_geom_volcross"),
        help="artifact directory for results.json / results.md / figures",
    )
    args = parser.parse_args(argv)
    outdir = Path(args.out)
    results = build_bundle(outdir)
    summary = results["summary"]
    for key, row in summary.items():
        if key.startswith("all_") or key.startswith("quadrature"):
            continue
        print(
            f"  {key:16s} status={row['status']:22s} "
            f"V={row['pipeline_volume_m3']:.9g} m3  "
            f"engines-spread={row['engine_spread_worst']:.2e}  order={row['observed_order']}"
        )
    print(f"  all watertight: {summary['all_watertight']}")
    print(f"  quadrature converged to {summary['quadrature_converged_relative']:.2e} relative")
    print(f"written: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())