"""Command-line entry point for the first TensorCAD geometry example."""
from __future__ import annotations

import argparse
import json

from .hull import HullParameters, export_hull, generate_hull


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a reproducible TensorCAD parametric hull example.")
    parser.add_argument("--length", type=float, default=10.0, help="length in metres")
    parser.add_argument("--beam", type=float, default=2.0, help="maximum beam in metres")
    parser.add_argument("--draft", type=float, default=0.8, help="draft in metres")
    parser.add_argument("--nx", type=int, default=41, help="longitudinal sample count")
    parser.add_argument("--nz", type=int, default=13, help="vertical sample count")
    parser.add_argument("--output", default="tensorcad-output", help="output directory")
    args = parser.parse_args()

    model = generate_hull(HullParameters(args.length, args.beam, args.draft, args.nx, args.nz))
    artifacts = export_hull(model, args.output)
    print(json.dumps({"artifacts": artifacts, "metadata": model.metadata()}, indent=2))


if __name__ == "__main__":
    main()
