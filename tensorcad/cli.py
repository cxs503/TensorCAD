"""Command-line entry point for the TensorCAD hull geometry model.

The CLI is the agent-facing surface: it takes parameters, prints measured
results (not just "done"), and exits non-zero with the full validation list when
a model is unusable. ``--json`` emits the complete reproducibility record for a
machine consumer.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from .export import export_hull
from .hull import HullParameters, generate_hull

SUMMARY_KEYS = (
    ("displacement_volume_m3", "displacement volume"),
    ("displacement_mass_kg", "displacement mass"),
    ("wetted_surface_area_m2", "wetted surface"),
    ("waterplane_area_m2", "waterplane area"),
    ("block_coefficient", "block coeff Cb"),
    ("prismatic_coefficient", "prismatic coeff Cp"),
    ("longitudinal_centre_buoyancy_m", "LCB from x=0"),
    ("vertical_centre_buoyancy_m", "VCB"),
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-hull",
        description="Generate a closed, validated parametric hull surface with metadata.",
    )
    parser.add_argument("--length", type=float, default=10.0, help="length in metres (required to be > 0)")
    parser.add_argument("--beam", type=float, default=2.0, help="maximum moulded breadth in metres")
    parser.add_argument("--draft", type=float, default=0.8, help="keel-to-waterline draft in metres")
    parser.add_argument("--nx", type=int, default=41,
                        help="longitudinal stations (an odd value lands on midship)")
    parser.add_argument("--nz", type=int, default=13, help="vertical levels from keel to waterline")
    parser.add_argument("--density", type=float, default=1025.0,
                        help="water density in kg/m^3 used for displacement mass (1025 seawater, 1000 fresh)")
    parser.add_argument("--output", default="tensorcad-output", help="directory for the exported artifacts")
    parser.add_argument("--no-export", action="store_true", help="validate and report only, write no files")
    parser.add_argument("--no-convergence", action="store_true",
                        help="skip the refined-grid discretisation estimate (faster, less informative)")
    parser.add_argument("--json", action="store_true",
                        help="print the full metadata record instead of a summary")
    return parser


def _summary(model: Any) -> dict[str, Any]:
    metadata = model.metadata()
    hydro = metadata.get("hydrostatics", {})
    lines: dict[str, Any] = {
        "parameters": metadata["parameters"],
        "mesh": {
            "vertex_count": metadata["mesh"]["vertex_count"],
            "triangle_count": metadata["mesh"]["triangle_count"],
            "bounding_box_size": metadata["mesh"]["bounding_box"]["size"],
        },
        "validation": {
            "passed": metadata["validation"]["passed"],
            "errors": metadata["validation"]["errors"],
            "warnings": metadata["validation"]["warnings"],
        },
        "hydrostatics": {key: hydro.get(key) for key, _ in SUMMARY_KEYS},
    }
    return lines


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    parameters = HullParameters(args.length, args.beam, args.draft, args.nx, args.nz, args.density)

    def fail(stage: str, exc: Exception) -> int:
        """One error contract: humans on stderr, machines on stdout, exit 2."""
        message = f"{type(exc).__name__}: {exc}" if not isinstance(exc, ValueError) else str(exc)
        print(f"tensorcad-hull: {stage} failed: {message}", file=sys.stderr)
        if args.json:
            print(json.dumps({"error": message, "stage": stage}, sort_keys=True))
        return 2

    try:
        model = generate_hull(parameters)
    except ValueError as exc:
        print(f"tensorcad-hull: invalid model: {exc}", file=sys.stderr)
        if args.json:
            print(json.dumps({"error": str(exc), "stage": "validate"}, sort_keys=True))
        return 2

    # One metadata document, computed once: the convergence estimate regenerates
    # a refined mesh, so paying for it twice (file + stdout) is real work.
    try:
        metadata = model.metadata(with_discretisation=not args.no_convergence)
    except Exception as exc:
        return fail("metadata", exc)

    artifacts: dict[str, str] = {}
    if not args.no_export:
        try:
            artifacts = export_hull(model, args.output, metadata=metadata)
        except Exception as exc:
            return fail("export", exc)

    if args.json:
        payload: dict[str, Any] = dict(metadata)
        payload["artifacts"] = artifacts
        try:
            print(json.dumps(payload, indent=2, sort_keys=True))
        except (TypeError, ValueError) as exc:
            return fail("serialise", exc)
    else:
        report = _summary(model)
        report["artifacts"] = artifacts
        print(json.dumps(report, indent=2))
    return 0 if model.validation.passed else 1


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main(sys.argv[1:]))
