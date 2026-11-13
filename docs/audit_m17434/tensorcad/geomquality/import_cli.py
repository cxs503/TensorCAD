"""Command-line entry point for the CAD import pipeline.

    python -m tensorcad.geomquality.import_cli hull.igs --out-dir runs/import

For every input file the CLI writes

    <out-dir>/<stem>.report.json   the full audit trail (source hash, per-face statistics,
                                   before/after topology, gate verdicts)
    <out-dir>/<stem>.stl           the repaired, watertight-if-possible triangle mesh
    <out-dir>/manifest.json        one aggregate record: route, units, gate verdicts

and prints one line per file.  The exit status is 0 only when every gate passed (or when
``--allow-open`` says an open surface model -- a propeller, a lifting surface -- is an
expected outcome).  A refusal is printed with the measured reason; the CLI never turns a
failed gate into a success.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

from .iges_trim import IgesError
from .trim_mesh import import_cad, write_mesh

__all__ = ["main", "run"]


def _summarise(report: dict[str, Any]) -> dict[str, Any]:
    gates = report.get("gates", {})
    after = report.get("after_repair", {})
    before = report.get("before_repair", {})
    return {
        "path": report.get("source", {}).get("path"),
        "sha256": report.get("source", {}).get("sha256"),
        "route": report.get("route", "trimmed"),
        "units": report.get("source", {}).get("units"),
        "scale_applied": report.get("source", {}).get("scale_applied"),
        "faces_in_file": report.get("import", {}).get("faces_in_file"),
        "faces_tessellated": report.get("import", {}).get("faces_tessellated"),
        "faces_refused": report.get("import", {}).get("faces_refused"),
        "triangles": after.get("triangles"),
        "free_edges_before": before.get("free_edges"),
        "free_edges_after": after.get("free_edges"),
        "watertight": after.get("watertight"),
        "signed_volume": after.get("signed_volume"),
        "gates": {name: verdict.get("passed") for name, verdict in gates.items()},
        "gates_passed": all(verdict.get("passed") for verdict in gates.values()) if gates else False,
        "route_notes": report.get("route_notes", []),
    }


def run(inputs: Sequence[str | Path], out_dir: str | Path, *, tolerance: float | None = None,
        scale: float | None = None, repair: bool = True, fill_holes: bool = True,
        max_hole_edges: int = 32, write_stl: bool = True,
        allow_patch_fallback: bool = True) -> dict[str, Any]:
    """Import every input and write the per-file reports, meshes and the manifest."""
    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, Any]] = []
    for item in inputs:
        source = Path(item)
        try:
            result = import_cad(
                source, tolerance=tolerance, scale=scale, repair=repair,
                fill_holes=fill_holes, max_hole_edges=max_hole_edges,
                allow_patch_fallback=allow_patch_fallback,
            )
            report_path = target / f"{source.stem}.report.json"
            report_path.write_text(json.dumps(result.report, indent=2, default=str))
            if write_stl:
                write_mesh(result.mesh, target / f"{source.stem}.stl")
            summary = _summarise(result.report)
            summary["report"] = str(report_path)
            entries.append(summary)
        except (IgesError, ValueError, OSError) as exc:
            entries.append({"path": str(source), "error": str(exc), "gates_passed": False})
    manifest = {
        "tool": "tensorcad.geomquality.import_cli",
        "files": entries,
        "all_gates_passed": all(entry.get("gates_passed") for entry in entries) if entries else False,
    }
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", help="IGES or STEP files")
    parser.add_argument("--out-dir", default="runs/cad_import", help="where reports and meshes go")
    parser.add_argument("--tolerance", type=float, default=None,
                        help="chord tolerance for boundary sampling, in file units")
    parser.add_argument("--scale", type=float, default=None,
                        help="metres per file unit; defaults to the file's declared units")
    parser.add_argument("--no-repair", action="store_true", help="weld only; skip hole filling")
    parser.add_argument("--no-fill-holes", action="store_true", help="repair but leave holes open")
    parser.add_argument("--max-hole-edges", type=int, default=32,
                        help="largest boundary loop, in edges, that repair may close")
    parser.add_argument("--no-patch-fallback", action="store_true",
                        help="refuse a file the trimmed reader cannot read instead of falling back")
    parser.add_argument("--no-stl", action="store_true", help="skip mesh export")
    parser.add_argument("--allow-open", action="store_true",
                        help="treat an open surface model (propeller, lifting surface) as expected")
    args = parser.parse_args(argv)
    manifest = run(
        args.inputs, args.out_dir, tolerance=args.tolerance, scale=args.scale,
        repair=not args.no_repair, fill_holes=not args.no_fill_holes,
        max_hole_edges=args.max_hole_edges, write_stl=not args.no_stl,
        allow_patch_fallback=not args.no_patch_fallback,
    )
    for entry in manifest["files"]:
        if "error" in entry:
            print(f"FAIL {entry['path']}: {entry['error']}")
            continue
        gates = " ".join(
            f"{name}={'PASS' if value else 'FAIL'}" for name, value in entry["gates"].items()
        )
        print(
            f"{'PASS' if entry['gates_passed'] else 'FAIL'} {entry['path']} "
            f"route={entry['route']} faces={entry['faces_tessellated']}/{entry['faces_in_file']} "
            f"triangles={entry['triangles']} free_edges={entry['free_edges_before']}->"
            f"{entry['free_edges_after']} volume={entry['signed_volume']} {gates}"
        )
    if args.allow_open:
        return 0 if all("error" not in entry for entry in manifest["files"]) else 1
    return 0 if manifest["all_gates_passed"] else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())