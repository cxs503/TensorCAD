"""CLI for the IGES/STEP → watertight mesh import pipeline (GAP-CAD-T1-IGES-PIPE).

Usage::

    tensorcad-iges-pipe input.igs [--output mesh.stl] [--tolerance 1e-4]
                                  [--u-samples 64] [--v-samples 32]
                                  [--report report.json] [--quiet]

The pipeline reads an IGES or STEP file, tessellates trimmed B-spline faces,
repairs the mesh (weld, drop degenerate/duplicate, fill small holes, re-orient),
and writes a watertight triangle mesh.  Every step is measured and the gate
verdicts are printed and optionally saved as JSON.

Exit codes
----------
0  All gates passed (watertight, manifold, single shell, positive volume).
1  One or more gates failed (report still written).
2  File could not be read or parsed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .geometry import Mesh
from .geomquality.trim_mesh import import_cad


def _parse_planes(specs: list[str], default_tolerance: float | None) -> list[tuple]:
    """Parse ``--closure-plane NAME=OX,OY,OZ,NX,NY,NZ[,TOL]`` into audit entries."""
    planes: list[tuple] = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"--closure-plane expects NAME=OX,OY,OZ,NX,NY,NZ[,TOL], got {spec!r}")
        name, raw = spec.split("=", 1)
        numbers = [piece.strip() for piece in raw.split(",")]
        if len(numbers) not in (6, 7):
            raise ValueError(f"--closure-plane {name!r} needs 6 or 7 numbers, got {len(numbers)}")
        try:
            values = [float(piece) for piece in numbers]
        except ValueError as exc:
            raise ValueError(f"--closure-plane {name!r}: {exc}") from exc
        tolerance = values[6] if len(values) == 7 else (
            default_tolerance if default_tolerance else 1e-6)
        planes.append((name, (tuple(values[0:3]), tuple(values[3:6])), tolerance))
    return planes


def _gate_summary(report: dict) -> str:
    lines = []
    for name, verdict in report.get("gates", {}).items():
        status = "PASS" if verdict.get("passed") else "FAIL"
        detail = verdict.get("detail", "")
        value = verdict.get("value", "")
        limit = verdict.get("limit")
        limit_str = f" (limit {limit})" if limit is not None else ""
        lines.append(f"  [{status}] {name}: {value}{limit_str} — {detail}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="tensorcad-iges-pipe",
        description="IGES/STEP → watertight mesh import pipeline",
    )
    parser.add_argument("input", help="IGES or STEP file path")
    parser.add_argument("--output", "-o", help="Output mesh file (STL/OBJ/VTK)")
    parser.add_argument("--tolerance", type=float, default=1e-4,
                        help="Weld tolerance in model units (default 1e-4)")
    parser.add_argument("--u-samples", type=int, default=64,
                        help="U-direction sample density for untrimmed patches")
    parser.add_argument("--v-samples", type=int, default=32,
                        help="V-direction sample density for untrimmed patches")
    parser.add_argument("--seal-openings", action="store_true",
                        help="Cap boundary loops that lie on a plane (the file's own cut "
                             "planes); leave defect gaps open and reported")
    parser.add_argument("--seal-planarity", type=float, default=None,
                        help="Absolute planarity tolerance in model units for the seal "
                             "stage (default: 1e-3 of the model extent)")
    parser.add_argument("--seam-tolerance", type=float, default=None,
                        metavar="DIST",
                        help="Seam coincidence budget in model units. Neighbouring patches "
                             "exported within DIST of each other are stitched, moving the "
                             "surface by at most DIST. Omit to measure the open boundary "
                             "without changing the mesh.")
    parser.add_argument("--symmetry-closure", action="store_true",
                        help="Detect a symmetry plane carrying the open boundary (the "
                             "publisher's half-hull centre plane), mirror the mesh across "
                             "it and weld the seam, producing a full watertight hull")
    parser.add_argument("--symmetry-min-coverage", type=float, default=0.5,
                        metavar="FRACTION",
                        help="Fraction of the open-boundary length the plane must carry "
                             "before it is accepted as the symmetry cut (default 0.5)")
    parser.add_argument("--closure-tolerance", type=float, default=None,
                        metavar="DIST",
                        help="Classify the open boundary within DIST (seam / patch "
                             "overlap / declared plane / unexplained) and record the "
                             "length in metres per cause. Measurement only: the mesh is "
                             "not changed.")
    parser.add_argument("--closure-plane", action="append", default=[],
                        metavar="NAME=OX,OY,OZ,NX,NY,NZ[,TOL]",
                        help="Declare an opening plane (e.g. a half-hull symmetry cut). "
                             "Boundary chords lying on it are reported as plane:<NAME>. "
                             "TOL defaults to the closure tolerance. Repeatable.")
    parser.add_argument("--report", "-r", help="Save gate report as JSON")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress stdout output")
    args = parser.parse_args(argv)
    try:
        closure_planes = _parse_planes(args.closure_plane, args.closure_tolerance)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    source = Path(args.input)
    if not source.exists():
        print(f"Error: file not found: {source}", file=sys.stderr)
        return 2

    try:
        result = import_cad(
            source,
            u_samples=args.u_samples,
            v_samples=args.v_samples,
            weld_tolerance=args.tolerance,
            seal_openings=args.seal_openings,
            seal_planarity=args.seal_planarity,
            seam_tolerance=args.seam_tolerance,
            symmetry_closure=args.symmetry_closure,
            symmetry_min_coverage=args.symmetry_min_coverage,
            closure_tolerance=args.closure_tolerance,
            closure_planes=closure_planes,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    mesh = result.mesh
    report = result.report

    # Write output mesh if requested
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        suffix = out_path.suffix.lower()
        if suffix == ".stl":
            _write_stl(mesh, out_path)
        elif suffix == ".obj":
            _write_obj(mesh, out_path)
        elif suffix == ".vtk":
            _write_vtk(mesh, out_path)
        else:
            # Default to STL
            _write_stl(mesh, out_path.with_suffix(".stl"))

    # Write report if requested
    if args.report:
        report_path = Path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(report, indent=2, default=str) + "\n",
            encoding="utf-8",
        )

    # Print summary
    if not args.quiet:
        after = report.get("after_repair", {})
        gates = report.get("gates", {})
        all_passed = all(v.get("passed") for v in gates.values())

        print(f"Source: {source.name}")
        print(f"Vertices: {len(mesh.vertices)}")
        print(f"Triangles: {len(mesh.faces)}")
        print(f"Volume: {after.get('signed_volume', 'N/A'):.4f} {mesh.units}³")
        print(f"Watertight: {after.get('watertight', False)}")
        print(f"Components: {after.get('components', 'N/A')} "
              f"(closed: {after.get('closed_components', 'N/A')})")
        print()
        print("Gate verdicts:")
        print(_gate_summary(report))
        closure = report.get("open_boundary") or {}
        if closure:
            family = closure.get("family", {})
            print()
            print(f"Open boundary by cause (tolerance {closure['tolerance']:g}, "
                  f"{closure['free_edges']} chords, {closure['open_length']:.2f} total):")
            for name in ("seam", "overlap", "declared_plane", "unexplained"):
                entry = family.get(name, {})
                print(f"  {name:<15} {entry.get('edges', 0):6d} chords "
                      f"{entry.get('length', 0.0):10.2f}")
        print()
        if all_passed:
            print("Result: ALL GATES PASSED")
        else:
            failed = [n for n, v in gates.items() if not v.get("passed")]
            print(f"Result: {len(failed)} GATE(S) FAILED: {', '.join(failed)}")

    # Exit code: 0 if all gates passed, 1 otherwise
    all_passed = all(v.get("passed") for v in report.get("gates", {}).values())
    return 0 if all_passed else 1


def _write_stl(mesh: Mesh, path: Path) -> None:
    """Write binary STL."""
    import struct
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)  # header
        f.write(struct.pack("<I", len(mesh.faces)))
        for a, b, c in mesh.faces:
            va, vb, vc = mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
            # Normal (zero = auto-computed by most readers)
            f.write(struct.pack("<12f", 0, 0, 0,
                                *va, *vb, *vc))
            f.write(b"\x00\x00")  # attribute byte count


def _write_obj(mesh: Mesh, path: Path) -> None:
    """Write OBJ."""
    lines = []
    for x, y, z in mesh.vertices:
        lines.append(f"v {x:.8g} {y:.8g} {z:.8g}")
    for a, b, c in mesh.faces:
        lines.append(f"f {a+1} {b+1} {c+1}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_vtk(mesh: Mesh, path: Path) -> None:
    """Write legacy VTK."""
    lines = ["# vtk DataFile Version 3.0", "TensorCAD mesh", "ASCII",
             "DATASET POLYDATA", f"POINTS {len(mesh.vertices)} float"]
    for x, y, z in mesh.vertices:
        lines.append(f"{x:.8g} {y:.8g} {z:.8g}")
    lines.append(f"POLYGONS {len(mesh.faces)} {4 * len(mesh.faces)}")
    for a, b, c in mesh.faces:
        lines.append(f"3 {a} {b} {c}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
