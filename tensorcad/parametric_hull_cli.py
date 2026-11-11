"""Command-line entry point for :mod:`tensorcad.parametric_hull`.

This is the executable face of the parametric hull generator: one command goes from a
hull definition (Wigley family, DARPA SUBOFF, or a tabulated offset file) to a lofted
NURBS surface, a triangulated surface mesh, hydrostatics and mesh-quality numbers,
optionally written out as a gmsh 2.2 file for the FEM/CFD handoff.

Three definitions are supported::

    python -m tensorcad.parametric_hull_cli wigley --length 10 --beam 1 --draft 0.625
    python -m tensorcad.parametric_hull_cli suboff --station-fractions
    python -m tensorcad.parametric_hull_cli offsets sections.json --length 10 \\
        --gmsh hull.msh --json hull.json

Every number printed is measured on the mesh the command just built; nothing is
transcribed.  The JSON document written by ``--json`` is the audit record: hull
parameters, section count, surface degrees, mesh counts, hydrostatics, mesh quality and
the watertight/manifold verdict of the exported file.

Exit codes: ``0`` success, ``1`` the mesh failed the quality/closure gate (the JSON is
still written), ``2`` usage/definition error.

Only the standard library is used, exactly like the core: the entry point must run in a
constrained environment with no CAD kernel installed.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

from .parametric_hull import (
    CrossSection,
    ParametricHull,
    ParametricHullError,
    build_parametric_hull,
    cross_section_from_offsets,
    depth_exponent_for_block_coefficient,
    mesh_density_for,
    suboff_cross_sections,
    suboff_station_fractions,
    wigley_analytic_properties,
    wigley_family_sections,
)

__all__ = ["build_hull_from_args", "hull_summary", "main"]


# --------------------------------------------------------------------------- #
# definition -> cross-sections
# --------------------------------------------------------------------------- #
def _sections_wigley(args: argparse.Namespace) -> tuple[tuple[CrossSection, ...], dict[str, Any]]:
    depth_exponent = 2.0
    if args.block_coefficient is not None:
        depth_exponent = depth_exponent_for_block_coefficient(
            args.block_coefficient, fore_exponent=args.fore_exponent
        )
    sections = wigley_family_sections(
        length=args.length,
        beam=args.beam,
        draft=args.draft,
        fore_exponent=args.fore_exponent,
        depth_exponent=depth_exponent,
        n_stations=args.stations,
        n_points=args.points,
    )
    definition = {
        "family": "wigley",
        "length": args.length,
        "beam": args.beam,
        "draft": args.draft,
        "fore_exponent": args.fore_exponent,
        "depth_exponent": depth_exponent,
        "requested_block_coefficient": args.block_coefficient,
        "stations": len(sections),
        "points_per_section": args.points,
    }
    return sections, definition


def _sections_suboff(args: argparse.Namespace) -> tuple[tuple[CrossSection, ...], dict[str, Any]]:
    from .suboff import SUBOFF_LENGTH_M

    length = args.length if args.length is not None else SUBOFF_LENGTH_M
    fractions = suboff_station_fractions(args.stations) if args.station_fractions else None
    sections = suboff_cross_sections(
        length_m=length, n_stations=args.stations, n_points=args.points,
        station_fractions=fractions,
    )
    definition = {
        "family": "suboff",
        "length": length,
        "stations": len(sections),
        "points_per_section": args.points,
        "clustered_stations": bool(args.station_fractions),
    }
    return sections, definition


def _sections_offsets(args: argparse.Namespace) -> tuple[tuple[CrossSection, ...], dict[str, Any]]:
    """Read a tabulated hull definition.

    Accepted JSON forms (all lengths in metres, ``y >= 0`` half-breadths)::

        {"length": 10.0, "stations": [{"x": 0.0, "half_breadths": [...], "heights": [...]}, ...]}
        {"length": 10.0, "stations": [{"x_fraction": 0.5, "points": [[y, z], ...]}, ...]}

    ``x`` may be a station fraction in ``[0, 1]`` (``x_fraction``) or an absolute
    longitudinal position (``x``); the two are not mixed inside one file.
    """
    path = pathlib.Path(args.table)
    if not path.is_file():
        raise ParametricHullError(f"offset table {path} does not exist")
    payload = json.loads(path.read_text())
    length = float(payload.get("length", args.length))
    stations = payload.get("stations")
    if not stations:
        raise ParametricHullError(f"{path} has no 'stations' list")
    sections = []
    for index, station in enumerate(stations):
        if "x_fraction" in station:
            x = float(station["x_fraction"]) * length
        elif "x" in station:
            x = float(station["x"])
        else:
            raise ParametricHullError(f"station {index} has neither 'x' nor 'x_fraction'")
        if "points" in station:
            points = [(float(y), float(z)) for y, z in station["points"]]
            if len(points) < 2:
                raise ParametricHullError(f"station {index} needs at least 2 points")
            sections.append(CrossSection(x=x, points=tuple(points)))
        else:
            sections.append(
                cross_section_from_offsets(x, station["half_breadths"], station["heights"])
            )
    sections.sort(key=lambda s: s.x)
    definition = {
        "family": "offset_table",
        "table": str(path),
        "length": length,
        "stations": len(sections),
    }
    return tuple(sections), definition


_BUILDERS = {"wigley": _sections_wigley, "suboff": _sections_suboff, "offsets": _sections_offsets}


# --------------------------------------------------------------------------- #
# build + summary
# --------------------------------------------------------------------------- #
def build_hull_from_args(args: argparse.Namespace) -> tuple[ParametricHull, dict[str, Any]]:
    """Build the hull described by parsed CLI arguments plus its definition record."""
    sections, definition = _BUILDERS[args.family](args)
    if args.mesh_n_u is not None and args.mesh_n_v is not None:
        n_u, n_v = args.mesh_n_u, args.mesh_n_v
    else:
        # Default density: keep triangles near-isotropic (dx ~ 2.5 dz), which is the
        # rule the acceptance criterion for mesh quality is built on.
        n_u, n_v = mesh_density_for(definition["length"], _profile_extent(sections))
        if args.mesh_n_u is not None:
            n_u = args.mesh_n_u
        if args.mesh_n_v is not None:
            n_v = args.mesh_n_v
    # A waterplane lid needs a planar top edge: a canoe body's sections end at the
    # waterline, a body of revolution's do not.  For the SUBOFF hull the lid is
    # therefore geometrically impossible, and the record says so instead of failing
    # with a message about a flag the caller never touched.
    close_waterplane = not args.open_waterplane
    if args.family == "suboff" and close_waterplane:
        close_waterplane = False
        definition["waterplane_lid"] = (
            "omitted: body of revolution, section tops are not coplanar (shell, not a body)"
        )
    hull = build_parametric_hull(
        sections,
        mesh_n_u=n_u,
        mesh_n_v=n_v,
        waterline_z=args.waterline,
        close_waterplane=close_waterplane,
    )
    return hull, definition


def _profile_extent(sections: tuple[CrossSection, ...]) -> float:
    """Vertical extent of the sections (draft for a canoe body)."""
    low = min(s.min_z for s in sections)
    high = max(s.max_height for s in sections)
    extent = high - low
    return extent if extent > 0.0 else 1.0


def hull_summary(hull: ParametricHull, definition: dict[str, Any]) -> dict[str, Any]:
    """Machine-readable record of one build: geometry, hydrostatics, quality."""
    hydro = hull.hydrostatics.to_dict()
    quality = dict(hull.quality)
    summary = {
        "definition": definition,
        "parameters": dict(hull.parameters),
        "sections": len(hull.sections),
        "surface": {
            "degree": list(hull.surface.degree),
            "control_points": [
                len(hull.surface.control_net),
                len(hull.surface.control_net[0]),
            ],
        },
        "mesh": {
            "n_u": hull.mesh_data.n_u,
            "n_v": hull.mesh_data.n_v,
            "vertices": hull.mesh_data.n_vertices,
            "faces": hull.mesh_data.n_faces,
            "hull_faces": hull.mesh_data.hull_face_count,
            "waterplane_closed": hull.mesh_data.waterplane_z is not None,
        },
        "hydrostatics": hydro,
        "quality": quality,
    }
    return summary


def _gate_failures(summary: dict[str, Any]) -> list[str]:
    """Checks the command insists on before it returns success."""
    failures = []
    quality = summary["quality"]
    jacobian = quality.get("min_scaled_jacobian")
    aspect = quality.get("max_aspect_ratio")
    if jacobian is not None and jacobian <= 0.0:
        failures.append(f"min_scaled_jacobian {jacobian:.4g} <= 0 (inverted triangle)")
    if aspect is not None and aspect < 1.0:
        failures.append(f"max_aspect_ratio {aspect:.4g} < 1 (degenerate)")
    volume = summary["hydrostatics"].get("displaced_volume_m3")
    if volume is not None and volume <= 0.0:
        failures.append(f"displacement volume {volume:.6g} <= 0 (surface faces inward)")
    return failures


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def _write_gmsh(hull: ParametricHull, path: pathlib.Path) -> dict[str, Any]:
    """Hand the built surface to the FEM adapter and write gmsh 2.2 ASCII."""
    from .adapters.fem import parametric_hull_surface, write_gmsh

    surface = parametric_hull_surface(hull)
    written = write_gmsh(surface, path)
    return {
        "path": str(written),
        "watertight": surface.watertight,
        "edge_manifold": surface.edge_manifold,
        "euler_characteristic": surface.euler_characteristic,
        "volume_m3": surface.volume_m3,
        "area_m2": surface.area_m2,
        "digest": surface.digest,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _add_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mesh-n-u", type=int, default=None, help="mesh columns along the hull")
    parser.add_argument("--mesh-n-v", type=int, default=None, help="mesh rows across the profile")
    parser.add_argument("--waterline", type=float, default=0.0, help="waterline z (default 0)")
    parser.add_argument(
        "--open-waterplane",
        action="store_true",
        help="leave the waterplane open (no lid). The FEM handoff then refuses the mesh.",
    )


def _add_output_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--json", type=pathlib.Path, default=None, help="write the record here")
    parser.add_argument("--gmsh", type=pathlib.Path, default=None, help="write gmsh 2.2 ASCII here")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-parametric-hull",
        description="parametric hull surface generation, surface mesh and hydrostatics",
    )
    subparsers = parser.add_subparsers(dest="family", required=True)

    wigley = subparsers.add_parser("wigley", help="generalised Wigley hull (p = 2 default)")
    wigley.add_argument("--length", type=float, default=10.0)
    wigley.add_argument("--beam", type=float, default=1.0)
    wigley.add_argument("--draft", type=float, default=0.625)
    wigley.add_argument("--stations", type=int, default=11)
    wigley.add_argument("--points", type=int, default=12)
    wigley.add_argument("--fore-exponent", type=float, default=2.0)
    wigley.add_argument(
        "--block-coefficient",
        type=float,
        default=None,
        help="target Cb; the depth exponent is solved for it",
    )
    _add_mesh_arguments(wigley)
    _add_output_arguments(wigley)

    suboff = subparsers.add_parser("suboff", help="DARPA SUBOFF bare hull")
    suboff.add_argument("--length", type=float, default=None, help="scale length (default 4.356 m)")
    suboff.add_argument("--stations", type=int, default=21)
    suboff.add_argument("--points", type=int, default=20)
    suboff.add_argument(
        "--station-fractions",
        action="store_true",
        help="cosine-clustered stations that resolve the elliptical stern cap",
    )
    _add_mesh_arguments(suboff)
    _add_output_arguments(suboff)

    offsets = subparsers.add_parser("offsets", help="tabulated hull from a JSON offset file")
    offsets.add_argument("table", type=str)
    offsets.add_argument("--length", type=float, default=10.0, help="used if the table omits it")
    _add_mesh_arguments(offsets)
    _add_output_arguments(offsets)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        hull, definition = build_hull_from_args(args)
    except (ParametricHullError, ValueError, KeyError) as error:
        print(f"tensorcad-parametric-hull: {error}", file=sys.stderr)
        return 2

    summary = hull_summary(hull, definition)
    if args.family == "wigley":
        # The closed form is the reference the measured numbers are judged against; it
        # is carried in the record so a consumer never has to recompute it by hand.
        summary["analytic_reference"] = wigley_analytic_properties(
            definition["length"],
            definition["beam"],
            definition["draft"],
            fore_exponent=definition["fore_exponent"],
            depth_exponent=definition["depth_exponent"],
        )
        reference = summary["analytic_reference"]
        measured = summary["hydrostatics"]
        if reference.get("displaced_volume"):
            summary["analytic_comparison"] = {
                "volume_rel_error": (
                    measured["displaced_volume_m3"] - reference["displaced_volume"]
                ) / reference["displaced_volume"],
                "waterplane_area_rel_error": (
                    measured["waterplane_area_m2"] - reference["waterplane_area"]
                ) / reference["waterplane_area"],
            }
    if args.gmsh is not None:
        args.gmsh.parent.mkdir(parents=True, exist_ok=True)
        try:
            summary["gmsh"] = _write_gmsh(hull, args.gmsh)
        except (ParametricHullError, ValueError) as error:
            # The FEM handoff refuses an open shell or an inward-facing surface; that is
            # a definition error of this invocation, not a traceback for the caller.
            print(f"tensorcad-parametric-hull: {error}", file=sys.stderr)
            return 2

    summary["gate_failures"] = _gate_failures(summary)
    text = json.dumps(summary, indent=2, sort_keys=True)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n")
    print(text)
    return 0 if not summary["gate_failures"] else 1


if __name__ == "__main__":  # pragma: no cover - exercised through main() as well
    raise SystemExit(main())
