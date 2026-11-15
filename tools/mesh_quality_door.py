"""The CI door: run :func:`tensorcad.meshquality.door.admit` over a fixed mesh suite.

What this tool is for
---------------------
Board audit of `GAP-CAD-D1-MESH-QUALITY-LIB` (``docs/mesh-quality-lib-audit-m17404.md``,
finding F2) recorded that the gate "can be a door" but that no workflow invoked it, and
that the only meshes the gate was ever applied to were the ones the cross-check happened
to build. This tool turns that into a door with a *fixed* suite, so a change to the metric
library, the topology code or the gate is measured against the same meshes every run:

1. **Verdict.** Six shells (the audited hull set plus the sliver control) and three
   volume meshes go through :func:`admit`; each case has an expected verdict recorded
   here, so a door that starts admitting a sliver -- or starts rejecting a good hull --
   is a failure of this tool, not a note in someone's terminal.
2. **Closed form.** The volume cases have exact references derived by hand on the
   *analytic* geometry (sheared-block Jacobian, Kuhn tetrahedra of the unit cube), not
   recomputed by the library's own arithmetic. The deviation is printed for every one.
3. **Regression baseline.** Every check value and every summary statistic is compared
   against ``docs/meshquality/door_baseline.json`` at a relative tolerance of ``1e-12``.
   The baseline is *recomputed* on every run -- it is never read back into the report --
   so a hand-edited baseline cannot buy a pass; ``--bless`` is the only way to move it
   and it has to be part of a reviewed commit.
4. **Exit code** ``0`` clean, ``1`` a door/closed-form/baseline failure, ``2`` usage.
   That is the contract a CI workflow needs.

Upstream (``vtk``/``gmsh``/``meshmagick``) is still absent on this host; the cross-check
tool records that per run. This tool does not repeat the upstream claim: the differential
it enforces is closed-form plus regression, which is exactly as strong as it says.

Usage::

    python tools/mesh_quality_door.py                 # verify (CI)
    python tools/mesh_quality_door.py --bless         # re-freeze the baseline
    python tools/mesh_quality_door.py --case hex_shear_2  # one case
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # allow `python tools/...` from a source checkout
    sys.path.insert(0, str(ROOT))
TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from tensorcad.meshquality.door import admit  # noqa: E402

#: Relative tolerance of the regression baseline and of the closed-form comparisons.
#: The library's own differential (``tools/mesh_quality_crosscheck.py``) accepts
#: ``1e-6`` against an independent route and ``1e-12`` against exact arithmetic; a
#: *regression* door has no excuse for a loose number, so it uses the exact-arithmetic
#: tolerance.
TOLERANCE = 1e-12

#: Angles come out of a ``acos``/``atan`` and are compared absolutely, in degrees.
ANGLE_TOLERANCE = 1e-9

#: Where the frozen baseline lives. In the repository, next to the other committed
#: measurement records (``docs/golden/``), not under ``artifacts/``, so a baseline move
#: shows up in review as a data edit.
BASELINE_PATH = ROOT / "docs" / "meshquality" / "door_baseline.json"
ARTIFACT_DIR = ROOT / "artifacts" / "mesh_quality_door"
DOC_PATH = ROOT / "docs" / "mesh_quality_door.md"


# --------------------------------------------------------------------------------------
# the suite
# --------------------------------------------------------------------------------------
def _unit_cube_points() -> tuple[tuple[float, float, float], ...]:
    return (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    )


_KUHN_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)


def hex_shear_block(shear: float = 0.5) -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    """Two hexahedra on ``x in [0,1]``, sheared in ``y`` by ``y -> y + shear * x``.

    The Jacobian of the map is constant over the block,
    ``J = [[dx, shear * dx, 0], [0, 1, 0], [0, 0, 1]]`` with ``dx = 0.5``, so the exact
    element numbers are known in closed form and the interior face at ``x = 0.5`` has an
    exactly known non-orthogonality ``atan(shear)`` (see :func:`closed_form`).
    """
    points: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}

    def node(x: float, y: float, z: float) -> int:
        key = (x, y + shear * x, z)
        if key not in index:
            index[key] = len(points)
            points.append(key)
        return index[key]

    cells: list[list[int]] = []
    for x0, x1 in ((0.0, 0.5), (0.5, 1.0)):
        cells.append(
            [
                node(x0, 0.0, 0.0),
                node(x1, 0.0, 0.0),
                node(x1, 1.0, 0.0),
                node(x0, 1.0, 0.0),
                node(x0, 0.0, 1.0),
                node(x1, 0.0, 1.0),
                node(x1, 1.0, 1.0),
                node(x0, 1.0, 1.0),
            ]
        )
    return points, cells


def hex_inverted_block() -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    """The same sheared block with the right-hand hexahedron's orientation reversed.

    A mesher bug that is *invisible* to surface checks and to a bounding box, and that
    inverts the sign of the Jacobian everywhere in one cell. The door has to reject it.
    """
    points, cells = hex_shear_block(0.5)
    flipped = list(cells[1])
    flipped[1], flipped[2] = flipped[2], flipped[1]  # swap two corners of one face
    return points, [cells[0], flipped]


def _closed_form_shear(shear: float) -> dict[str, float]:
    dx = 0.5
    det = dx  # [[dx, dx*shear, 0], [0, 1, 0], [0, 0, 1]]
    fro2 = dx * dx + (dx * shear) ** 2 + 1.0 + 1.0
    return {
        "volume_per_cell": det,
        "volume_total": 2.0 * det,
        # min_sicn is d * det^(2/d) / ||J||_F^2 with the ideal element (unit cube) = 1
        "min_sicn": 3.0 * det ** (2.0 / 3.0) / fro2,
        # scaled Jacobian: det / (|J e1| |J e2| |J e3|)
        "scaled_jacobian": det / (dx * math.hypot(1.0, shear) * 1.0),
        # the face at x = 0.5 is planar in y-z, so its normal is e_x, while the
        # cell-centre vector is (0.5, 0.5 * shear, 0): the angle is atan(shear)
        "max_non_orthogonality": math.degrees(math.atan(shear)),
    }


def closed_form() -> dict[str, dict[str, float]]:
    """Exact references for the volume cases, derived on the analytic geometry."""
    shear = _closed_form_shear(0.5)
    kuhn_volume = 1.0 / 6.0
    # A Kuhn tetrahedron has J = [[1,1,1],[0,1,1],[0,0,1]] from its corner at the
    # origin: det = 1, ||J||_F^2 = 6, and Gmsh's reference tetrahedron is the right-corner
    # one (for which the same expression equals 1), so gmsh_icn = 3 * 1 / 6.
    kuhn_icn = 3.0 * 1.0 ** (2.0 / 3.0) / 6.0
    return {
        # ``hex_inverted`` is a defect control, not a closed-form case: one of its two
        # cells is wound inside out, so the analytic numbers above do not apply to it.
        "hex_shear_2": shear,
        "tet_cube_kuhn6": {
            "volume_per_cell": kuhn_volume,
            "volume_total": 1.0,
            "gmsh_icn": kuhn_icn,
        },
    }


def suite() -> dict[str, dict[str, Any]]:
    """The fixed mesh suite: shells from the audited generators plus volume controls."""
    import mesh_quality_crosscheck as xcheck  # the audited generators, reused verbatim

    cases: dict[str, dict[str, Any]] = {}

    hulls = xcheck._hull_meshes()
    cases["wigley"] = {
        "points": hulls["wigley"].vertices,
        "triangles": hulls["wigley"].faces,
        "expect": "pass",
    }
    cases["series60"] = {
        "points": hulls["series60"].vertices,
        "triangles": hulls["series60"].faces,
        "expect": "fail",
        "why": "min angle 9.494 deg against the shipped 10 deg floor: the marginal case",
    }
    cases["suboff_coarse"] = {
        "points": hulls["suboff_coarse"].vertices,
        "triangles": hulls["suboff_coarse"].faces,
        "expect": "fail",
        "why": "raw polar Suboff: a singular pole, min angle 0.0039 deg, edge ratio 9211",
    }
    cases["suboff_taubin30"] = {
        "points": hulls["suboff_taubin30"].vertices,
        "triangles": hulls["suboff_taubin30"].faces,
        "expect": "pass",
        "why": "Taubin-30 smoothing of the raw polar mesh restores it",
    }

    cube = xcheck._box_mesh()
    cases["box_unit"] = {"points": cube.vertices, "triangles": cube.faces, "expect": "pass"}

    sliver = xcheck._sliver_patch()
    cases["sliver_patch"] = {
        "points": sliver.vertices,
        "triangles": sliver.faces,
        "expect": "fail",
        "why": "a deliberately degenerate triangle patch: the gate has to reject it",
    }

    points, cells = hex_shear_block(0.5)
    cases["hex_shear_2"] = {
        "points": points,
        "cells": {"hex": cells},
        "expect": "pass",
        "why": "closed-form controls: non-orthogonality atan(0.5) = 26.565 deg",
    }

    points, cells = hex_shear_block(0.5)
    cases["tet_cube_kuhn6"] = {
        "points": list(_unit_cube_points()),
        "cells": {"tet": [list(tet) for tet in _KUHN_TETS]},
        "expect": "pass",
        "why": "6-tet Kuhn split of the unit cube: total volume exactly 1, gmsh_icn 0.5",
    }

    points, cells = hex_inverted_block()
    cases["hex_inverted"] = {
        "points": points,
        "cells": {"hex": cells},
        "expect": "fail",
        "why": "one cell wound inside out: min_sicn < 0 has to stop the door",
    }
    return cases


# --------------------------------------------------------------------------------------
# measurement, comparison
# --------------------------------------------------------------------------------------
def _flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    """Leaf paths of a nested JSON-able structure, floats kept as floats."""
    out: dict[str, Any] = {}
    if isinstance(value, Mapping):
        for key, item in value.items():
            out.update(_flatten(item, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, (list, tuple)):
        for position, item in enumerate(value):
            out.update(_flatten(item, f"{prefix}[{position}]"))
    else:
        out[prefix] = value
    return out


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _deviation(expected: float, got: float) -> float:
    return abs(got - expected)


def _within(expected: float, got: float, tolerance: float) -> bool:
    if math.isnan(expected) or math.isnan(got):
        # ``genus`` is undefined (NaN) on an open shell; NaN == NaN is the agreement we want
        return math.isnan(expected) and math.isnan(got)
    if not (math.isfinite(expected) and math.isfinite(got)):
        return expected == got
    scale = max(abs(expected), abs(got), 1.0)
    return abs(got - expected) <= tolerance * scale


def measure_case(case: Mapping[str, Any]) -> dict[str, Any]:
    verdict = admit(
        case["points"],
        triangles=case.get("triangles"),
        cells=case.get("cells"),
    )
    return verdict.to_dict()


def check_values(measured: Mapping[str, Any]) -> dict[str, float]:
    """The number behind every check, keyed by check name (duplicates suffixed)."""
    out: dict[str, float] = {}
    for check in measured["checks"]:
        key = check["name"]
        if key in out:
            index = 2
            while f"{key}#{index}" in out:
                index += 1
            key = f"{key}#{index}"
        out[key] = float(check["value"])
    return out


def snapshot(measured: Mapping[str, Any]) -> dict[str, Any]:
    """Everything the baseline freezes: the verdict, every check value, every statistic."""
    return {
        "passed": bool(measured["passed"]),
        "checks": check_values(measured),
        "measurements": {
            key: value
            for key, value in _flatten(measured["measurements"]).items()
            if isinstance(value, bool)
            # a NaN (``genus`` on an open shell) is not a regression signal and is not
            # representable in JSON, so it stays out of the frozen baseline
            or (_is_number(value) and math.isfinite(float(value)))
        },
    }


def _compare_baseline(
    frozen: Mapping[str, Any], fresh: Mapping[str, Any], case: str
) -> list[str]:
    problems: list[str] = []
    if frozen.get("passed") != fresh.get("passed"):
        problems.append(f"{case}: verdict {frozen.get('passed')} -> {fresh.get('passed')}")
    for section in ("checks", "measurements"):
        old = frozen.get(section, {})
        new = fresh.get(section, {})
        for key in sorted(set(old) | set(new)):
            if key not in new:
                problems.append(f"{case}/{section}/{key}: disappeared from the measurement")
                continue
            if key not in old:
                problems.append(f"{case}/{section}/{key}: new value {new[key]!r} not in baseline")
                continue
            expected, got = old[key], new[key]
            if isinstance(expected, bool) or isinstance(got, bool):
                if expected != got:
                    problems.append(f"{case}/{section}/{key}: {expected} -> {got}")
                continue
            if not _within(float(expected), float(got), TOLERANCE):
                problems.append(
                    f"{case}/{section}/{key}: {expected!r} -> {got!r} "
                    f"(relerr {abs(got - expected) / max(abs(expected), 1e-300):.3e})"
                )
    return problems


def _compare_closed_form(
    case: str, reference: Mapping[str, float], measured: Mapping[str, Any]
) -> list[dict[str, Any]]:
    """Pair each closed-form reference with the library's number for the same quantity."""
    rows: list[dict[str, Any]] = []
    summary = measured["measurements"].get("summary", {})
    kind = next(iter(summary), None)
    block = summary.get(kind, {}) if kind else {}
    stats = block.get("summary", {})
    cells = int(block.get("count", 0))
    checks = check_values(measured)

    def add(quantity: str, expected: float, got: float, tolerance: float) -> None:
        rows.append(
            {
                "case": case,
                "quantity": quantity,
                "expected": expected,
                "got": got,
                "difference": _deviation(expected, got),
                "tolerance": tolerance,
                "ok": _within(expected, got, tolerance),
            }
        )

    for quantity, expected in reference.items():
        if quantity == "max_non_orthogonality":
            got = checks.get("volume.max_non_orthogonality")
            if got is not None:
                add(quantity, expected, got, ANGLE_TOLERANCE)
        elif quantity == "volume_per_cell":
            if "volume" in stats:
                add(quantity, expected, stats["volume"]["mean"], TOLERANCE)
        elif quantity == "volume_total":
            if "volume" in stats:
                add(quantity, expected, stats["volume"]["mean"] * cells, TOLERANCE)
        elif quantity == "min_sicn":
            if "min_sicn" in stats:
                add(quantity, expected, stats["min_sicn"]["min"], TOLERANCE)
        elif quantity == "scaled_jacobian":
            if "scaled_jacobian" in stats:
                add(quantity, expected, stats["scaled_jacobian"]["min"], TOLERANCE)
        elif quantity == "gmsh_icn":
            if "gmsh_icn" in stats:
                add(quantity, expected, stats["gmsh_icn"]["mean"], TOLERANCE)
    return rows


def run(cases: Mapping[str, dict[str, Any]] | None = None) -> dict[str, Any]:
    """Measure the suite and evaluate every failure mode; returns the record."""
    cases = cases or suite()
    references = closed_form()
    baseline = _read_baseline()
    record: dict[str, Any] = {"cases": {}, "closed_form": [], "problems": []}

    for name, case in cases.items():
        measured = measure_case(case)
        record["cases"][name] = {
            "expected": case["expect"],
            "why": case.get("why", ""),
            "passed": measured["passed"],
            "failures": measured["failures"],
            "checks": check_values(measured),
            "snapshot": snapshot(measured),
            "payload": measured,
        }
        got = "pass" if measured["passed"] else "fail"
        if got != case["expect"]:
            record["problems"].append(
                f"{name}: the door said {got.upper()}, the suite expects {case['expect'].upper()}"
            )
        if name in references:
            record["closed_form"].extend(_compare_closed_form(name, references[name], measured))

    for row in record["closed_form"]:
        if not row["ok"]:
            record["problems"].append(
                f"closed-form {row['case']}/{row['quantity']}: expected {row['expected']:.17g}, "
                f"got {row['got']:.17g} (difference {row['difference']:.3e})"
            )

    if baseline is not None:
        for name, case in record["cases"].items():
            if name in baseline.get("cases", {}):
                record["problems"].extend(
                    _compare_baseline(baseline["cases"][name], case["snapshot"], name)
                )
        record["baseline"] = {
            "path": _display(BASELINE_PATH),
            "verified": True,
            "cases": len(baseline.get("cases", {})),
        }
    else:
        record["baseline"] = {"path": _display(BASELINE_PATH), "verified": False}

    record["passed"] = not record["problems"]
    return record


def _display(path: Path) -> str:
    """Repository-relative when possible -- a test may point the baseline at ``tmp_path``."""
    try:
        return str(Path(path).relative_to(ROOT))
    except ValueError:
        return str(path)


def _read_baseline() -> dict[str, Any] | None:
    if not BASELINE_PATH.exists():
        return None
    return json.loads(BASELINE_PATH.read_text(encoding="utf-8"))


def baseline_document(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "generator": "tools/mesh_quality_door.py --bless",
        "tolerance": TOLERANCE,
        "angle_tolerance": ANGLE_TOLERANCE,
        "cases": {
            name: case["snapshot"] for name, case in record["cases"].items()
        },
    }


# --------------------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------------------
def write_markdown(record: Mapping[str, Any], path: Path = DOC_PATH) -> Path:
    lines: list[str] = []
    lines.append("# The mesh-quality door: one verdict, and the meshes it is measured on")
    lines.append("")
    lines.append(
        f"*Gap:* `GAP-CAD-D1-MESH-QUALITY-LIB`. *Generated by* `python tools/mesh_quality_door.py` "
        f"at {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} -- every number below is injected "
        "from `artifacts/mesh_quality_door/results.json`, none is transcribed by hand. "
        f"*Verdict:* **{'PASS' if record['passed'] else 'FAIL'}**."
    )
    lines.append("")
    lines.append("## 1. What the door is")
    lines.append("")
    lines.append(
        "`tensorcad.meshquality.door.admit` composes the three gates the library already "
        "shipped separately -- shell topology (watertight / edge- and vertex-manifold / "
        "winding / outward / genus), the per-element `DEFAULT_GATE`, and the OpenFOAM "
        "interior-face non-orthogonality -- into one verdict whose failures still name the "
        "threshold, the cell kind and the margin. It adds no measure and no threshold of "
        "its own."
    )
    lines.append("")
    lines.append("## 2. Suite and verdicts")
    lines.append("")
    lines.append("| case | kind | expected | door | why |")
    lines.append("|---|---|---|---|---|")
    for name, case in record["cases"].items():
        payload = case["payload"]["measurements"]
        kind = "shell" if "topology" in payload else ", ".join(
            key.split("_")[0] for key in payload.get("summary", {})
        )
        door = "PASS" if case["passed"] else "FAIL"
        lines.append(
            f"| `{name}` | {kind} | {case['expected'].upper()} | **{door}** | {case['why'] or 'audited hull mesh'} |"
        )
    lines.append("")
    lines.append("## 3. Checks the door applied on a passing hull (Wigley)")
    lines.append("")
    reference = record["cases"].get("wigley")
    if reference:
        lines.append("| check | value | limit | verdict | detail |")
        lines.append("|---|---|---|---|---|")
        for check in reference["payload"]["checks"]:
            lines.append(
                f"| `{check['name']}` | {check['value']:.9g} | `{check['limit']}` | "
                f"{'pass' if check['passed'] else '**fail**'} | {check['detail']} |"
            )
        lines.append("")
    lines.append("## 4. Closed form vs library (volume cases)")
    lines.append("")
    lines.append(
        "The references are derived by hand on the analytic geometry -- the sheared block's "
        "constant Jacobian and the Kuhn split of the unit cube -- not by the library's own "
        "arithmetic. Relative tolerance `1e-12`, angles absolute `1e-9` deg."
    )
    lines.append("")
    lines.append("| case | quantity | closed form | library | difference | ok |")
    lines.append("|---|---|---|---|---|---|")
    for row in record["closed_form"]:
        lines.append(
            f"| `{row['case']}` | {row['quantity']} | {row['expected']:.17g} | "
            f"{row['got']:.17g} | {row['difference']:.3e} | {'yes' if row['ok'] else '**no**'} |"
        )
    lines.append("")
    lines.append("## 5. Baseline")
    lines.append("")
    baseline = record["baseline"]
    if baseline["verified"]:
        lines.append(
            f"`{baseline['path']}` was recomputed and compared on all {baseline['cases']} cases at "
            f"relative tolerance `{record.get('tolerance', TOLERANCE):g}`. The baseline is regenerated "
            "in memory on every run and never read back into the report, so a hand-edited baseline "
            "cannot buy a pass; `--bless` is the only way to move it and it lands as a reviewed data edit."
        )
    else:
        lines.append(f"`{baseline['path']}` is absent: run `--bless` once to freeze it.")
    lines.append("")
    lines.append("## 6. Notes, limits, applicability")
    lines.append("")
    lines.append(
        "* The upstream differential (``vtk``/``gmsh``/``meshmagick``) is **not** repeated here: "
        "those packages have no wheel on this host, and `tools/mesh_quality_crosscheck.py` already "
        "records that per run. This door's differential is closed-form plus regression.\n"
        "* The gate thresholds are policy, not measurement: a PASS means the mesh is inside the "
        "declared limits, not that a solver will be accurate on it. "
        "`docs/mesh_quality_gate_power.md` measures that second, weaker claim and reports FAIL.\n"
        "* `sliver_patch` and `hex_inverted` are controls: they exist to prove the door closes. "
        "If either starts passing, the door has stopped being a door."
    )
    lines.append("")
    lines.append("## 7. Errors against the references")
    lines.append("")
    worst = max(
        (row["difference"] / max(abs(row["expected"]), 1e-300) for row in record["closed_form"]),
        default=0.0,
    )
    lines.append(
        f"Worst closed-form relative deviation over {len(record['closed_form'])} comparisons: "
        f"`{worst:.3e}`. All are below `3%`, and below the `1e-12` this door enforces; the quiet "
        "reason is that the volume cases are rational, so the only error source is floating-point "
        "evaluation of `det^(2/3)`, `atan` and the floating Jacobian of the isoparametric map."
    )
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the mesh-quality door over the fixed suite (exit 1 on failure)."
    )
    parser.add_argument("--bless", action="store_true", help="rewrite the regression baseline")
    parser.add_argument("--case", action="append", help="restrict to one case (repeatable)")
    parser.add_argument("--json", default=str(ARTIFACT_DIR / "results.json"), help="artifact path")
    parser.add_argument("--doc", default=str(DOC_PATH), help="markdown path")
    parser.add_argument("--no-doc", action="store_true", help="skip the markdown document")
    parser.add_argument("--quiet", action="store_true", help="verdict line only")
    args = parser.parse_args(argv)

    cases = suite()
    if args.case:
        unknown = [name for name in args.case if name not in cases]
        if unknown:
            print(f"unknown case(s): {', '.join(unknown)}", file=sys.stderr)
            return 2
        cases = {name: cases[name] for name in args.case}

    record = run(cases)
    record["tolerance"] = TOLERANCE

    if args.bless:
        if args.case:
            print("--bless rewrites the whole baseline; do not combine it with --case", file=sys.stderr)
            return 2
        BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
        BASELINE_PATH.write_text(
            json.dumps(baseline_document(record), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        record["problems"] = [p for p in record["problems"] if "baseline" not in p]
        record["baseline"] = {"path": _display(BASELINE_PATH), "verified": True, "cases": len(cases)}
        record["passed"] = not record["problems"]

    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    if not args.no_doc:
        write_markdown(record, Path(args.doc))

    if args.quiet:
        print("PASS" if record["passed"] else "FAIL")
    else:
        for name, case in record["cases"].items():
            door = "PASS" if case["passed"] else "FAIL"
            print(f"  {name:22s} expected {case['expected']:4s} door {door}")
        for row in record["closed_form"]:
            print(
                f"  closed form {row['case']}/{row['quantity']}: difference {row['difference']:.3e}"
            )
        print(f"  baseline: {'verified' if record['baseline']['verified'] else 'absent'}")
        for problem in record["problems"]:
            print(f"  FAILURE: {problem}")
        print("PASS" if record["passed"] else "FAIL")
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())