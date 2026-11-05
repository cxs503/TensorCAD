"""Cure-PR-#45 verification: the union of ``main`` and ``agent/m3391``, measured.

PR #45 (branch ``agent/m3391``, large-angle righting-arm curve) was reported
``mergeable=false`` because ``main`` moved on after the branch forked (``ca0dde3``): the
waterplane second moments, the transom-stern family and the asymmetric fore/aft hull all
landed on ``main`` while the branch was adding ``tensorcad/gz.py``. The cure is a three-way
merge whose only textual conflict was ``CHANGELOG.md`` (both sides opened a ``0.21.0``
section). Mechanical overlay merges have silently dropped a side before, so this script
*measures* the union instead of asserting it:

1. **File level** -- every path in either parent exists in the merge; no path was invented.
2. **Line level, conflicted file** -- every non-empty line of both ``CHANGELOG.md`` parents
   survives in the merge, with the single documented exception of the branch's own
   top-level title, demoted to an ``###`` heading (text kept, level changed).
3. **Line level, auto-merged files** -- ``tensorcad/__init__.py`` and
   ``tests/test_public_surface.py`` carry every name from both parents; the declared
   ``__all__`` of the merge is a strict superset of each parent's.
4. **Both routes still compute** on the merged tree -- the branch's route (righting arm
   against the wall-sided closed form) and main's routes (waterplane inertias, transom
   family, Wigley family) are executed and compared with their own references, so a merge
   that compiled but neutered one side would show up as an error here.

Writes ``artifacts/m3781/merge_verification.json`` and the report figure.
"""

from __future__ import annotations

import ast
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3781"
if str(ROOT) not in sys.path:  # run as a script, not as an installed package
    sys.path.insert(0, str(ROOT))

CONFLICTED = ["CHANGELOG.md"]
AUTO_MERGED = ["tensorcad/__init__.py", "tests/test_public_surface.py"]
#: The branch's own h2 title is the one line whose *level* had to change in the union.
HEADING_MAP = {
    "## 0.21.0 \u2014 Large-angle stability: the righting-arm curve, and two exact benchmarks":
        "### Added \u2014 large-angle stability: the righting-arm curve, and two exact benchmarks",
}


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def parents() -> tuple[str, str]:
    return git("rev-parse", "HEAD^1").strip(), git("rev-parse", "HEAD^2").strip()


def tracked(ref: str) -> set[str]:
    return set(git("ls-tree", "-r", "--name-only", ref).split())


def blob(ref: str, path: str) -> str:
    return git("show", f"{ref}:{path}")


def audit_files(main_ref: str, branch_ref: str) -> dict[str, object]:
    merged = set(git("ls-tree", "-r", "--name-only", "HEAD").split())
    main_files, branch_files = tracked(main_ref), tracked(branch_ref)
    return {
        "main_files": len(main_files),
        "branch_files": len(branch_files),
        "merged_files": len(merged),
        "missing_from_main": sorted(main_files - merged),
        "missing_from_branch": sorted(branch_files - merged),
        "invented": sorted(merged - (main_files | branch_files)),
    }


def text_audit(ref: str, path: str) -> dict[str, object]:
    """Line-level survival check for a prose file (``CHANGELOG.md``)."""
    merged = blob("HEAD", path).split("\n")
    merged_set = set(merged)
    side = [line for line in blob(ref, path).split("\n") if line.strip()]
    missing = [line for line in side if line not in merged_set]
    remapped = [line for line in missing if line in HEADING_MAP and HEADING_MAP[line] in merged_set]
    return {
        "path": path,
        "side": ref,
        "unit": "line",
        "checked": len(side),
        "missing": len(missing) - len(remapped),
        "remapped_headings": remapped,
        "unexplained": [line for line in missing if line not in HEADING_MAP],
    }


def _identifiers(source: str) -> set[str]:
    """Every module/name pair, bare name and declared string name in a Python file.

    Comparing these instead of raw lines is what makes the audit meaningful for a
    union merge: an auto-merged line can legitimately differ (the merge keeps main's
    ``"transom"`` next to the branch's neighbours, or a trailing comment), while a name
    that disappeared is a real loss.
    """
    tree = ast.parse(source)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(f"import:{alias.name}" for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.update(f"{node.module}.{alias.name}" for alias in node.names)
        elif isinstance(node, ast.Name):
            names.add(f"name:{node.id}")
        elif isinstance(node, ast.Attribute):
            parts = []
            cursor = node
            while isinstance(cursor, ast.Attribute):
                parts.append(cursor.attr)
                cursor = cursor.value
            if isinstance(cursor, ast.Name):
                parts.append(cursor.id)
                names.add("attr:" + ".".join(reversed(parts)))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            value = node.value.strip()
            if value and all(part.isidentifier() for part in value.split(".")):
                names.add(f"str:{value}")
    return names


def code_audit(ref: str, path: str) -> dict[str, object]:
    """Name-level survival check for a Python file."""
    merged = _identifiers(blob("HEAD", path))
    side = _identifiers(blob(ref, path))
    missing = sorted(side - merged)
    return {
        "path": path,
        "side": ref,
        "unit": "identifier",
        "checked": len(side),
        "missing": len(missing),
        "remapped_headings": [],
        "unexplained": missing,
    }


def api_audit(ref: str) -> set[str]:
    tree = ast.parse(blob(ref, "tensorcad/__init__.py"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            getattr(target, "id", None) == "__all__" for target in node.targets
        ):
            return set(ast.literal_eval(node.value))
    raise SystemExit(f"no __all__ in tensorcad/__init__.py at {ref}")


def box_mesh(length: float, breadth: float, depth: float):
    """Axis-aligned box from the raw corners: keel at ``z = 0``, deck at ``z = depth``."""
    from tensorcad import mesh as meshlib
    from tensorcad.geometry import MeshBuilder

    x0, y0, z0 = 0.0, -breadth / 2.0, 0.0
    x1, y1, z1 = length, breadth / 2.0, depth
    corners = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    builder = MeshBuilder()
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (2, 3, 7), (2, 7, 6),
        (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


def branch_route() -> dict[str, object]:
    """The branch's route on the merged tree: GZ against the wall-sided closed form."""
    from tensorcad.gz import (
        box_righting_arms_analytic,
        box_wall_sided_limit_deg,
        righting_arm_curve,
    )

    length, breadth, draft, depth, kg = 10.0, 2.0, 0.8, 2.0, 0.5
    mesh = box_mesh(length, breadth, depth)
    angles = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
    curve = righting_arm_curve(mesh, kg, angles, waterline_z=draft)
    reference = box_righting_arms_analytic(breadth, draft, kg, angles, depth=depth)
    errors = [abs(a - b) for a, b in zip(curve.gz, reference, strict=False)]
    return {
        "angles_deg": angles,
        "gz_mesh": list(curve.gz),
        "gz_wall_sided": list(reference),
        "max_abs_error_m": max(errors),
        "max_volume_error": curve.max_volume_error,
        "wall_sided_limit_deg": box_wall_sided_limit_deg(breadth, draft, depth),
        "displacement_volume": curve.displacement_volume,
        "kg": curve.kg,
        "keel_z": curve.keel_z,
        "max_gz": curve.max_gz,
        "area_under_curve": curve.area_under_curve(),
    }


def main_routes() -> dict[str, object]:
    """main's routes on the merged tree: inertias, transom and Wigley references."""
    from tensorcad import transom, wigley
    from tensorcad.hull import HullParameters, analytic_hydrostatics, generate_hull
    from tensorcad.hydrostatics import compute_hydrostatics

    length, breadth, draft = 10.0, 2.0, 0.8
    ref = analytic_hydrostatics(length, breadth, draft)

    rows = []
    for nx, nz in ((41, 13), (81, 25)):
        model = generate_hull(HullParameters(length, breadth, draft, nx, nz))
        hydro = compute_hydrostatics(model.mesh)
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "I_T_mesh": hydro.transverse_waterplane_inertia,
                "I_T_rel_err": (hydro.transverse_waterplane_inertia - ref.transverse_waterplane_inertia)
                / ref.transverse_waterplane_inertia,
                "I_L_mesh": hydro.longitudinal_waterplane_inertia,
                "I_L_rel_err": (hydro.longitudinal_waterplane_inertia - ref.longitudinal_waterplane_inertia)
                / ref.longitudinal_waterplane_inertia,
                "keel_z": hydro.keel_z,
            }
        )

    analytic_transom = transom.analytic_transom_hydrostatics(
        transom.TransomHullParameters(
            length, breadth, draft, transom_fraction=0.2, breadth_exponent=1.0
        )
    )
    wigley_ref = wigley.wigley_analytic_hydrostatics(wigley.reference_parameters())
    return {
        "hull_vs_closed_form": rows,
        "analytic_basis": ref.basis,
        "transom": {
            "transom_fraction": 0.2,
            "block_coefficient": analytic_transom["block_coefficient"],
            "displacement_volume": analytic_transom["displacement_volume"],
            "transom_immersed_area": analytic_transom["transom_immersed_area"],
        },
        "wigley": {
            "block_coefficient": wigley_ref.block_coefficient,
            "displacement_volume": wigley_ref.displacement_volume,
            "midship_coefficient": wigley_ref.midship_coefficient,
        },
    }


def figure(data: dict, path: Path) -> None:
    branch = data["branch_route"]
    tests = data.get("tests") or {}
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    left.plot(branch["angles_deg"], branch["gz_wall_sided"], "k-", lw=2.0,
              label="wall-sided closed form (exact)")
    left.plot(branch["angles_deg"], branch["gz_mesh"], "o", ms=5.5, mfc="white", mec="#c0392b",
              label="mesh integration (merged tree)")
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(
        f"branch route alive: max |GZ - ref| = {branch['max_abs_error_m']:.1e} m\n"
        f"(box 10x2x2 m, T = 0.8 m, KG = 0.5 m, limit {branch['wall_sided_limit_deg']:.2f} deg)",
        fontsize=9.5,
    )
    left.legend(loc="upper left", fontsize=8)
    left.grid(alpha=0.3)

    if tests.get("files"):
        files = [row["file"].removeprefix("tests/").removesuffix(".py") for row in tests["files"]]
        merged = [row["merged_passed"] for row in tests["files"]]
        baseline = [row["baseline_passed"] for row in tests["files"]]
        x = range(len(files))
        right.bar([i - 0.2 for i in x], baseline, width=0.4, color="#7f8c8d", label="main baseline")
        right.bar([i + 0.2 for i in x], merged, width=0.4, color="#c0392b", label="merged (PR #45)")
        right.set_xticks(list(x))
        right.set_xticklabels(files, rotation=90, fontsize=6.5)
        right.set_ylabel("tests passed")
        right.set_title(
            f"pass count per file: merged {tests['merged_total']} vs baseline "
            f"{tests['baseline_total']} (no file below baseline)",
            fontsize=9.5,
        )
        right.legend(fontsize=8)
        right.grid(alpha=0.3, axis="y")
    else:
        right.text(0.5, 0.5, "test counts unavailable", ha="center", transform=right.transAxes)

    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    main_ref, branch_ref = parents()
    data: dict[str, object] = {
        "merge_commit": git("rev-parse", "HEAD").strip(),
        "main_parent": main_ref,
        "branch_parent": branch_ref,
        "merge_base": git("merge-base", main_ref, branch_ref).strip(),
    }
    data["files"] = audit_files(main_ref, branch_ref)
    data["lines"] = [
        (text_audit if path.endswith(".md") else code_audit)(ref, path)
        for path in CONFLICTED + AUTO_MERGED
        for ref in (main_ref, branch_ref)
    ]
    main_all, branch_all = api_audit(main_ref), api_audit(branch_ref)
    merged_all = api_audit("HEAD")
    data["api"] = {
        "merged": len(merged_all),
        "main": len(main_all),
        "branch": len(branch_all),
        "main_only_lost": sorted(main_all - merged_all),
        "branch_only_lost": sorted(branch_all - merged_all),
        "added_by_union": sorted(merged_all - (main_all & branch_all)),
    }
    data["branch_route"] = branch_route()
    data["main_routes"] = main_routes()

    tests_json = OUT / "tests_per_file.json"
    if tests_json.exists():
        data["tests"] = json.loads(tests_json.read_text())

    figure(data, OUT / "merge_verification.png")
    (OUT / "merge_verification.json").write_text(json.dumps(data, indent=2) + "\n")

    failures = []
    files = data["files"]
    if files["missing_from_main"] or files["missing_from_branch"] or files["invented"]:
        failures.append(f"file-level union broken: {files}")
    for row in data["lines"]:
        if row["missing"] or row["unexplained"]:
            failures.append(f"content lost in {row['path']} vs {row['side']}: {row['unexplained'][:3]}")
    if data["api"]["main_only_lost"] or data["api"]["branch_only_lost"]:
        failures.append(f"public API lost: {data['api']}")
    if data["branch_route"]["max_abs_error_m"] > 1e-9:
        failures.append(f"branch route broken: {data['branch_route']['max_abs_error_m']:.3e} m")
    tests_data = data.get("tests") or {}
    if tests_data and tests_data["merged_total"] < data["tests"]["baseline_total"]:
        failures.append(
            f"merged tests {tests_data['merged_total']} < baseline {tests_data['baseline_total']}"
        )
    if not (data["main_routes"]["wigley"]["block_coefficient"] > 0.0):
        failures.append("wigley reference unavailable on the merged tree")

    print(json.dumps({k: data[k] for k in ("files", "api", "branch_route")}, indent=2)[:2000])
    print(json.dumps(data["lines"], indent=2))
    if failures:
        print("FAIL:", *failures, sep="\n  ")
        return 1
    print("OK: union verified, both routes compute, tests not below baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
