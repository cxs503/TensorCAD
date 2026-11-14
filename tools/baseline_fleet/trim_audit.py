#!/usr/bin/env python3
"""Audit the trim geometry of a hull archive: entity identity, curve sanity, trim loops.

Why this tool exists
--------------------
Two ledger reasons for the baseline-fleet anchor were built on claims that the bytes
do not support:

* "the file carries N type-126 *surfaces* and this reader decodes type 128 only" --
  type 126 is the rational B-spline **curve**, and the reader decodes it;
* "curve_uv is field 5 of a type 142" -- stated as a field offset with no measurement.

This tool replaces both with output: for every archive it reports

1. ``census`` -- what entity types are present, and how many basis surfaces are bound
   to a trimmed patch (``IgesFile.surface_census``);
2. ``curves`` -- a sanity verdict per type-126 record: a knot vector must be
   non-decreasing and all coordinates must lie inside the file's model-space box;
   a decoder that indexes the tail in the wrong order fails both at once
   (it reads coordinates as knots);
3. ``trim_pointer_audit`` -- per trimmed patch, which candidate curve-on-surface
   pointer lands inside the basis surface's ``(u, v)`` domain.  A parameter-space loop
   must be inside; a model-space curve cannot be.  Where neither candidate wins, the
   patch is reported undecided instead of defaulted.

Usage
-----
    python tools/baseline_fleet/trim_audit.py FILE.igs [FILE.igs ...] --out out.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.iges import BSplineCurve, IgesError, read_iges  # noqa: E402


def curve_sanity(model, samples: int = 12) -> dict:
    """Monotonicity and model-space containment over every type-126 record."""
    surfaces = model.surfaces()
    boxes = [surface.bounding_box() for surface in surfaces]
    if boxes:
        lo = [min(box[0][axis] for box in boxes) for axis in range(3)]
        hi = [max(box[1][axis] for box in boxes) for axis in range(3)]
    else:
        lo = hi = None
    entries = model.of_type(126)
    monotone = decoded = outside = 0
    failures: list[str] = []
    for entry in entries:
        try:
            curve = BSplineCurve.from_parameters(model.parameters(entry.sequence))
        except (IgesError, ValueError, IndexError) as exc:
            failures.append(f"DE {entry.sequence}: {type(exc).__name__}: {exc}")
            continue
        decoded += 1
        if all(a <= b for a, b in zip(curve.knots, curve.knots[1:])):
            monotone += 1
        elif len(failures) < 5:
            failures.append(
                f"DE {entry.sequence}: knot vector not non-decreasing "
                f"(first three {[round(k, 4) for k in curve.knots[:3]]})"
            )
        if lo is not None:
            span = [hi[axis] - lo[axis] for axis in range(3)]
            tolerance = 0.5 * max(span) if max(span) > 0 else 1.0
            for point in curve.sample(samples):
                if any(
                    point[axis] < lo[axis] - tolerance or point[axis] > hi[axis] + tolerance
                    for axis in range(3)
                ):
                    outside += 1
                    break
    return {
        "type126_count": len(entries),
        "decoded": decoded,
        "knots_non_decreasing": monotone,
        "outside_model_box": outside,
        "model_box": None if lo is None else [[round(lo[a], 3), round(hi[a], 3)] for a in range(3)],
        "failures": failures,
    }


def audit(path: str) -> dict:
    model = read_iges(path)
    return {
        "file": Path(path).name,
        "unit_flag": model.unit_flag,
        "metres_per_unit": model.metres_per_unit,
        "census": model.surface_census(),
        "curves": curve_sanity(model),
        "trim_pointer_audit": model.trim_pointer_audit(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="IGES archives to audit")
    parser.add_argument("--out", default=None, help="write the report here as JSON")
    args = parser.parse_args(argv)

    report = {"schema": "tensorcad.trim_audit/1.0", "files": []}
    for path in args.files:
        row = audit(path)
        report["files"].append(row)
        census = row["census"]
        curves = row["curves"]
        print(
            f"{row['file']}: 128={census['surfaces_type128']} "
            f"126={census['curves_type126']} trimmed={census['basis_surfaces_trimmed']} "
            f"untrimmed={census['basis_surfaces_untrimmed']} "
            f"| curves decoded {curves['decoded']}/{curves['type126_count']} "
            f"non-decreasing {curves['knots_non_decreasing']} outside box {curves['outside_model_box']}"
        )
        audit_row = row["trim_pointer_audit"]
        if audit_row["patch_count"]:
            print(
                f"    trim patches {audit_row['patch_count']}: CRV wins "
                f"{audit_row['decided_field_crv']}, BPTR wins {audit_row['decided_field_bptr']}, "
                f"undecided {audit_row['undecided']}"
            )
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())