"""Diagnose the authoritative baseline-fleet IGES sources, byte by byte.

The fleet acceptance (``GAP-CAD-T1-BASELINE-FLEET``) grades a hull on
*provenance + topology + displaced volume*.  Every time a hull failed that gate
the interesting question was "why", and the answers that were written down by
hand kept drifting: "the reader ignores entity 124", "the file is at model
scale", "the z origin does not match a convention".  This script replaces those
narratives with measurements that a reviewer can re-run:

* what entity types the file actually contains (so a "missing transform matrix"
  claim can be checked instead of believed);
* the declared units and the reader's metres-per-unit;
* per patch: control-net size, chosen layout, whether the weight test -- not the
  smoothness score -- decided it, and whether the choice is ambiguous;
* the coordinate-class ranges of the whole control-point cloud, i.e. what the
  file's three axes really span, compared with the publisher's L / B / T / D;
* trailing values per type-128 record that no layout consumes (writer noise).

Nothing here is a verdict.  It prints facts plus the residual of each principal
dimension against the published figure, which is what tells a frame bug apart
from a scale bug.

Usage::

    python tools/baseline_fleet/diagnose_sources.py \
        --raw artifacts/baseline_fleet/raw \
        --json benchmarks/verified/baseline_fleet/geometry_diagnostics.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tensorcad import iges  # noqa: E402
from tensorcad.fleet_baseline import specs_from_json  # noqa: E402

#: Principal particulars used only as *diagnostic references* (they come from the
#: same publisher tables that ``fleet_specs.json`` carries, with page SHA-256).
#: ``depth`` is an optional extra that only KVLCC2 publishes, which is what makes
#: the axis-range test below decisive for that hull.
REFERENCES = {
    "KCS": {"l_pp": 230.0, "beam": 32.2, "draft": 10.8, "depth": None},
    "KVLCC2": {"l_pp": 320.0, "beam": 58.0, "draft": 20.8, "depth": 30.0},
    "JBC": {"l_pp": 280.0, "beam": 45.0, "draft": 16.5, "depth": None},
    "ONRT": {"l_pp": 154.0, "beam": 18.78, "draft": 5.494, "depth": None},
}


def _bbox(values: list[float]) -> dict[str, float]:
    return {"min": min(values), "max": max(values), "span": max(values) - min(values)}


def diagnose(name: str, path: Path) -> dict[str, object]:
    parsed = iges.read_iges(str(path))
    scale = parsed.metres_per_unit
    patches = []
    class_values: list[list[float]] = [[], [], []]
    trailing: list[int] = []
    for entry in parsed.of_type(iges.ENTITY_BSPLINE_SURFACE):
        params = parsed.parameters(entry.sequence)
        k1, k2, m1, m2 = (int(params[i]) for i in (1, 2, 3, 4))
        count_u, count_v = k1 + 1, k2 + 1
        points = count_u * count_v
        rational = int(params[7]) == 0
        start = 10 + (k1 + m1 + 2) + (k2 + m2 + 2)
        body = params[start:]
        surface = iges.BSplineSurface.from_parameters(params)
        weights = list(surface.weights)
        patches.append(
            {
                "de": entry.sequence,
                "net": [count_u, count_v],
                "rational": rational,
                "layout": surface.layout,
                "layout_ambiguous": surface.layout_ambiguous,
                "decided_by_weight_test": iges._uniform_weights(weights),
                "weight_count": len(set(round(w, 12) for w in weights)),
            }
        )
        trailing.append(len(body) - len(surface.control) * 3 - len(weights))
        for k in range(3):
            class_values[k].extend(point[k] * scale for point in surface.control)

    ranges = [_bbox(values) for values in class_values]
    reference = REFERENCES.get(name, {})
    # The three axes are reported in file order; the test below asks whether *one*
    # of them matches the published length, beam or depth, which is what an axis
    # permutation defect looks like from the outside.
    matches = {}
    for axis, box in enumerate(ranges):
        for quantity in ("l_pp", "beam", "depth"):
            target = reference.get(quantity)
            if not target:
                continue
            matches[f"axis{axis}_vs_{quantity}"] = round(100.0 * (box["span"] - target) / target, 2)
    kinds: dict[str, int] = {}
    for patch in patches:
        kinds[patch["layout"]] = kinds.get(patch["layout"], 0) + 1
    return {
        "file": str(path),
        "sha256": parsed.__dict__.get("sha256"),
        "units": parsed.unit_name or parsed.units,
        "unit_flag": parsed.unit_flag,
        "metres_per_unit": scale,
        "entity_counts": parsed.stats(),
        "transformation_matrices": len(parsed.of_type(124)),
        "patch_count": len(patches),
        "layout_histogram": kinds,
        "ambiguous_patches": sum(1 for p in patches if p["layout_ambiguous"]),
        "patches_decided_by_weight_test": sum(1 for p in patches if p["decided_by_weight_test"]),
        "trailing_values_max": max(trailing) if trailing else 0,
        "control_cloud_axes_m": ranges,
        "diagnostic_residuals_percent": matches,
        "patches": patches,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", default="artifacts/baseline_fleet/raw")
    parser.add_argument("--specs", default="benchmarks/verified/baseline_fleet/fleet_specs.json")
    parser.add_argument("--json", default="")
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args(argv)

    raw = Path(args.raw)
    if not raw.is_absolute():
        raw = REPO / raw
    specs_path = Path(args.specs)
    if not specs_path.is_absolute():
        specs_path = REPO / specs_path
    report: dict[str, object] = {"schema": "tensorcad.baseline_fleet.diagnostics/1.0", "hulls": {}}
    for spec in specs_from_json(specs_path):
        local = raw / spec.name / Path(spec.iges_path).name
        if not local.exists():
            report["hulls"][spec.name] = {
                "file": str(local),
                "status": "NOT_FETCHED",
                "remedy": "tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet",
            }
            continue
        report["hulls"][spec.name] = diagnose(spec.name, local)

    if args.json:
        out = Path(args.json)
        if not out.is_absolute():
            out = REPO / out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, sort_keys=False))
        print(f"wrote {out}")
    if args.summary:
        for name, entry in report["hulls"].items():  # type: ignore[union-attr]
            if entry.get("status") == "NOT_FETCHED":
                print(f"{name:8s} NOT_FETCHED")
                continue
            print(
                f"{name:8s} units={entry['units']}(x{entry['metres_per_unit']:g}) "
                f"patches={entry['patch_count']} layouts={entry['layout_histogram']} "
                f"ambig={entry['ambiguous_patches']} xform124={entry['transformation_matrices']}"
            )
            for axis, box in enumerate(entry["control_cloud_axes_m"]):
                print(f"         axis{axis}: {box['min']:+.3f} .. {box['max']:+.3f} m (span {box['span']:.3f})")
            print(f"         residuals {entry['diagnostic_residuals_percent']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())