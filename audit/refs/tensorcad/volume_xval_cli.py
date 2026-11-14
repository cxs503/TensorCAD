"""Command-line entry for the cross-tool displacement-volume cross-check.

Why a CLI and not only a library call
-------------------------------------
The value of :mod:`tensorcad.volume_xval` is that it produces a *second, independent*
number in the place where the first one is produced.  A capability that only exists
inside a script and a test cannot do that: the operator who suspects a volume has
to be able to ask for the second opinion from the product surface, on a hull file
they already have, and get an exit code a pipeline can read.

Two cases, both offline:

``wigley``
    the classical reference hull whose displacement has a closed form
    (``V = 4 L B T / 9``).  Builds the closed triangle solid from the half-breadth
    field, then measures it with every route: the platform's divergence sum,
    trimesh, numpy-stl (Mirtich through a binary-STL round trip) and a
    tessellation-free Simpson quadrature.  The analytic value is the reference, so
    this case separates "integrator wrong" from "tessellation coarse" on a hull
    where the truth is known.

``mesh FILE.stl``
    a hull mesh an operator already has.  No breadth field, so no sectional route
    and no analytic reference: the three triangle routes are compared with each
    other, which is exactly what localises a shared-code bias.  ``--reference-volume``
    adds the published displacement when the operator has one.

Exit codes: ``0`` the route comparison passed its gates, ``1`` a gate failed (the
report is still written -- an honest FAIL is evidence), ``2`` the input was refused.

The numbers printed here are the numbers in the JSON; nothing is transcribed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

from . import volume_xval as xv

SCHEMA = "tensorcad.volume_xval_cli/1.0"


def _wigley_case(
    *, stations: int, levels: int, sectional: tuple[int, int], ladder: tuple[int, ...]
) -> dict[str, Any]:
    from . import mesh as meshlib
    from . import wigley as wg

    parameters = wg.reference_parameters()
    mesh, _offsets = wg.build_mesh(parameters)

    def breadth(x: float, z: float) -> float:
        return wg.half_breadth(parameters, x, z)

    solid = xv.closed_solid_from_breadth(
        breadth,
        length=parameters.length,
        draft=parameters.draft,
        stations=stations,
        levels=levels,
    )
    pipeline = meshlib.total_signed_volume(mesh)
    reference = wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)
    sectional_volume = xv.volume_sectional(
        breadth,
        length=parameters.length,
        draft=parameters.draft,
        stations=sectional[0],
        levels=sectional[1],
    )
    return {
        "deficits": _deficits(
            breadth,
            length=parameters.length,
            draft=parameters.draft,
            smooth_total=sectional_volume,
            ladder=ladder,
        ),
        "case": "wigley",
        "ship": "Wigley (n=m=2)",
        "length_m": parameters.length,
        "beam_m": parameters.beam,
        "draft_m": parameters.draft,
        "stations": stations,
        "levels": levels,
        "topology": xv.topology_report(solid),
        "routes": {
            "pipeline": pipeline,
            "pipeline_name": "tensorcad.wigley mesh + divergence theorem",
            "divergence": xv.volume_divergence(solid),
            "engine_a_trimesh": xv.volume_trimesh(solid),
            "engine_b_numpy_stl": xv.volume_numpy_stl(solid),
            "sectional": sectional_volume,
        },
        "reference": {
            "volume_m3": reference,
            "kind": "analytic (exact)",
            "source": "V = L B T n m / ((n+1)(m+1)) = 4 L B T / 9 "
            "(tensorcad.wigley.analytic_volume)",
        },
    }


def _deficits(
    breadth,
    *,
    length: float,
    draft: float,
    smooth_total: float,
    ladder: tuple[int, ...],
) -> list[tuple[float, float]]:
    """Measure the (h, |smooth - inscribed|) points the attribution needs.

    A gap without an order is not a diagnosis, so the ladder is measured, not
    assumed: the closed solid is rebuilt at each station count, averaged over the
    two independent engines (they agree to ~1e-8, so the average is the polyhedron
    value), and compared with the tessellation-free smooth-field integral.
    """
    points: list[tuple[float, float]] = []
    for stations in ladder:
        levels = max(3, stations // 4 + 1)
        solid = xv.closed_solid_from_breadth(
            breadth, length=length, draft=draft, stations=stations, levels=levels
        )
        engines = 0.5 * (xv.volume_trimesh(solid) + xv.volume_numpy_stl(solid))
        points.append((length / stations, abs(smooth_total - engines)))
    return points


def _mesh_case(
    path: Path, *, reference_volume: float | None
) -> dict[str, Any]:
    from . import mesh as meshlib
    from .fleet_gate import read_stl

    # An STL stores every triangle's three corners separately, so the file as read is
    # a triangle soup: no edge is shared and every engine would (correctly) refuse to
    # call it a closed surface.  The weld is part of the exchange path a solver sees,
    # so it is done here, at a stated tolerance, and reported.
    mesh = meshlib.weld(read_stl(path), tolerance=1e-6)
    routes = {
        "pipeline": None,
        "pipeline_name": "not applicable: the file *is* the geometry",
        "divergence": xv.volume_divergence(mesh),
        "engine_a_trimesh": xv.volume_trimesh(mesh),
        "engine_b_numpy_stl": xv.volume_numpy_stl(mesh),
        "sectional": None,
    }
    reference = (
        None
        if reference_volume is None
        else {
            "volume_m3": float(reference_volume),
            "kind": "published displacement (operator-supplied)",
            "source": "--reference-volume",
        }
    )
    return {
        "case": "mesh",
        "ship": path.name,
        "path": str(path),
        "exchange": {"weld_tolerance": 1e-6, "note": "STL triangle soup welded before measuring"},
        "topology": xv.topology_report(mesh),
        "routes": routes,
        "reference": reference,
    }


def _verdict(record: dict[str, Any]) -> dict[str, Any]:
    """Classify the routes and state which gates the row passes."""
    routes = record["routes"]
    reference = record.get("reference")
    # The "pipeline" number is what the platform would report for this hull.  It has to
    # be the *smooth*-geometry claim for the deficit to mean anything: the Simpson
    # section quadrature integrates the half-breadth field and shares no mathematics
    # with the engines, so it is used when the case has one; otherwise the platform's
    # own divergence sum is the number under test.
    if routes.get("sectional") is not None:
        pipeline, pipeline_route = routes["sectional"], "Simpson section quadrature (smooth field)"
    else:
        pipeline, pipeline_route = routes["divergence"], "tensorcad divergence theorem"
    classification = xv.classify_routes(
        divergence=routes["divergence"],
        engine_a=routes["engine_a_trimesh"],
        engine_b=routes["engine_b_numpy_stl"],
        pipeline=pipeline,
        deficits=record.get("deficits") or (),
    )
    engine_spread = classification["engine_spread"]
    gates: list[dict[str, Any]] = [
        {
            "name": "engines_agree",
            "value": engine_spread,
            "limit": xv.ENGINE_AGREEMENT_TOLERANCE,
            "passed": engine_spread <= xv.ENGINE_AGREEMENT_TOLERANCE,
            "detail": "independent engines reading the same triangles",
        }
    ]
    reference_error = None
    if reference and reference.get("volume_m3"):
        measured = routes["engine_a_trimesh"]
        reference_error = (measured - reference["volume_m3"]) / reference["volume_m3"]
        gates.append(
            {
                "name": "reference_error",
                "value": reference_error,
                "limit": xv.ACCEPTANCE_TOLERANCE,
                "passed": abs(reference_error) <= xv.ACCEPTANCE_TOLERANCE,
                "detail": f"engine A trimesh vs {reference['kind']}",
            }
        )
    gates.append(
        {
            "name": "no_systematic_bias",
            "value": classification["status"],
            "limit": "not 'systematic-bias' / 'engine-disagreement'",
            "passed": classification["status"] in ("agrees", "tessellation-deficit"),
            "detail": classification["reason"],
        }
    )
    return {
        "classification": classification,
        "pipeline_route": pipeline_route,
        "engine_note": xv.engine_agreement_note(engine_spread),
        "reference_error": reference_error,
        "gates": gates,
        "passed": all(gate["passed"] for gate in gates),
    }


def _table(record: dict[str, Any], verdict: dict[str, Any]) -> str:
    routes = record["routes"]
    reference = record.get("reference")
    lines = [
        f"case: {record['case']}  ship: {record['ship']}",
        f"triangles: {record['topology']['faces']}  watertight: {record['topology']['watertight']}"
        f"  closed genus: {record['topology']['closed_genus']}",
        "",
        f"{'route':38s} {'volume [m3]':>16s} {'vs engine A':>13s}",
    ]
    engine_a = routes["engine_a_trimesh"]
    for key, label in (
        ("pipeline", routes["pipeline_name"]),
        ("divergence", "tensorcad divergence theorem"),
        ("engine_a_trimesh", "trimesh (independent engine A)"),
        ("engine_b_numpy_stl", "numpy-stl Mirtich (engine B, via STL)"),
        ("sectional", "Simpson section quadrature (no triangles)"),
    ):
        value = routes.get(key)
        if value is None:
            lines.append(f"{label:38s} {'withheld':>16s} {'-':>13s}")
            continue
        lines.append(f"{label:38s} {value:16.6f} {(value - engine_a) / engine_a * 100:+12.3e}%")
    if reference and reference.get("volume_m3"):
        lines.append("")
        lines.append(
            f"reference ({reference['kind']}): {reference['volume_m3']:.6f} m3; "
            f"engine A error {verdict['reference_error'] * 100:+.4f}%"
        )
    lines.append("")
    lines.append(f"verdict: {verdict['classification']['status']} -- {verdict['classification']['reason']}")
    lines.append(f"engines: {verdict['engine_note']}")
    lines.append(f"pipeline route under test: {verdict['pipeline_route']}")
    for gate in verdict["gates"]:
        state = "PASS" if gate["passed"] else "FAIL"
        lines.append(f"  [{state}] {gate['name']}: {gate['value']} (limit {gate['limit']})")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-hull volume-xval",
        description="Cross-check a displacement volume against engines that share no code.",
    )
    sub = parser.add_subparsers(dest="command")
    wigley = sub.add_parser("wigley", help="reference hull with a closed-form volume")
    wigley.add_argument("--stations", type=int, default=81, help="loft stations (default 81)")
    wigley.add_argument("--levels", type=int, default=17, help="loft waterline levels (default 17)")
    wigley.add_argument(
        "--sectional-stations", type=int, default=201, help="Simpson stations (default 201)"
    )
    wigley.add_argument(
        "--sectional-levels", type=int, default=401, help="Simpson levels (default 401)"
    )
    wigley.add_argument(
        "--ladder",
        type=int,
        nargs="+",
        default=[21, 41, 81],
        help="station counts of the refinement ladder used to attribute the gap",
    )
    mesh = sub.add_parser("mesh", help="a hull STL the operator already has")
    mesh.add_argument("path", type=Path)
    mesh.add_argument("--reference-volume", type=float, default=None, help="published displacement")
    for target in (wigley, mesh):
        target.add_argument("--out", type=Path, default=None, help="write the JSON report here")
        target.add_argument("--json", action="store_true", help="print JSON instead of a table")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    try:
        if args.command == "wigley":
            record = _wigley_case(
                stations=args.stations,
                levels=args.levels,
                sectional=(args.sectional_stations, args.sectional_levels),
                ladder=tuple(args.ladder),
            )
        else:
            record = _mesh_case(args.path, reference_volume=args.reference_volume)
        verdict = _verdict(record)
    except Exception as exc:  # refused input is a reportable state, not a traceback
        message = f"{type(exc).__name__}: {exc}" if not isinstance(exc, ValueError) else str(exc)
        print(f"tensorcad-hull volume-xval: refused: {message}", file=sys.stderr)
        if getattr(args, "json", False):
            print(json.dumps({"error": message}, sort_keys=True))
        return 2

    payload = {
        "schema": SCHEMA,
        "case": record,
        "verdict": verdict,
        "tolerances": {
            "acceptance": xv.ACCEPTANCE_TOLERANCE,
            "engine_agreement": xv.ENGINE_AGREEMENT_TOLERANCE,
            "order_window": list(xv.ORDER_WINDOW),
        },
    }
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=1, sort_keys=True))
    else:
        print(_table(record, verdict))
        if args.out is not None:
            print(f"\njson: {args.out}")
    return 0 if verdict["passed"] else 1


if __name__ == "__main__":  # pragma: no cover - exercised through the console script
    raise SystemExit(main())