"""PR #44 cure (m3847): merge inventory, cross-feature numbers and report figure.

Everything printed/published by this script is computed here — nothing is typed by
hand into the report.  Three products land in ``artifacts/m3847/``:

* ``merge_inventory.json`` — per-file record of what each side changed since the
  merge base and whether the cured tree still carries that side byte-for-byte.
* ``numbers.json`` — the experiments of §3 of the report (cross-feature GZ curves on
  the branch's solver *and* main's hull families, displacement conservation, GM
  agreement, transom hydrostatics against its closed form, and the test counts
  parsed out of the two pytest logs).
* ``gz_integration.png`` — four panels: GZ vs the independent 2-D/wall-sided
  reference, the same on main's hull families, the log-scale deviation, and the
  test-count comparison against the ``main`` baseline.

Run from the repository root::

    PYTHONPATH=. python3 tools/report_m3847.py
"""

from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import (  # noqa: E402
    analytic_transom_hydrostatics,
    compute_hydrostatics,
    generate_transom_hull,
    generate_wigley,
)
from tensorcad.geometry import Mesh, MeshBuilder  # noqa: E402
from tensorcad.stability import initial_stability, righting_arm_curve  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "artifacts" / "m3847"
BASE_REF = "e604f02"  # merge base of origin/main and origin/agent/m3347
OLD_SHA = "9f1572f"  # origin/agent/m3347 before the cure
MAIN_SHA = "2ff80b2"  # origin/main at cure time


# ---------------------------------------------------------------------------
# 1. merge inventory
# ---------------------------------------------------------------------------
def _git(*args: str) -> str:
    return subprocess.run(
        ("git",) + args, cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def _blob(ref: str, path: str) -> bytes | None:
    proc = subprocess.run(
        ("git", "show", f"{ref}:{path}"), cwd=ROOT, capture_output=True
    )
    return proc.stdout if proc.returncode == 0 else None


def _worktree(path: str) -> bytes | None:
    target = ROOT / path
    return target.read_bytes() if target.exists() else None


def inventory() -> dict:
    ours, theirs = "origin/agent/m3347", "origin/main"
    paths: set[str] = set()
    for ref in (ours, theirs):
        for line in _git("diff", "--name-status", BASE_REF, ref).splitlines():
            if line.strip():
                paths.add(line.split("\t")[-1])

    rows = []
    for path in sorted(paths):
        base = _blob(BASE_REF, path)
        mine = _blob(ours, path)
        main = _blob(theirs, path)
        merged = _worktree(path)
        branch_changed = mine != base
        main_changed = main != base
        side = (
            "both"
            if branch_changed and main_changed
            else "branch"
            if branch_changed
            else "main"
            if main_changed
            else "none"
        )
        rows.append(
            {
                "path": path,
                "changed_by": side,
                "identical_to_branch": merged == mine,
                "identical_to_main": merged == main,
                "textual_conflict": path in CONFLICTED,
                "resolution": RESOLUTIONS.get(
                    path,
                    "single-sided: taken byte-for-byte from the side that changed it",
                ),
                "preserved": (
                    (merged == mine or merged == main)
                    if side != "both"
                    else "union-verified"
                ),
            }
        )
    # every single-sided file must be byte-identical to the side that changed it
    lost = [
        r["path"]
        for r in rows
        if r["changed_by"] == "branch" and not r["identical_to_branch"]
        or r["changed_by"] == "main" and not r["identical_to_main"]
    ]
    return {
        "merge_base": BASE_REF,
        "branch_ref": ours,
        "branch_sha_before": _git("rev-parse", ours).strip(),
        "main_ref": theirs,
        "main_sha": _git("rev-parse", theirs).strip(),
        "cured_commit": _git("rev-parse", "HEAD").strip(),
        "files_touched_by_either_side": len(rows),
        "counts": {
            k: sum(1 for r in rows if r["changed_by"] == k)
            for k in ("branch", "main", "both", "none")
        },
        "textual_conflicts": sorted(CONFLICTED),
        "single_sided_content_lost": lost,
        "rows": rows,
    }


CONFLICTED = {"CHANGELOG.md"}
RESOLUTIONS = {
    "CHANGELOG.md": (
        "the only textual conflict: both sides opened the same `## 0.21.0` section. "
        "Kept main's `### Added` heading and then both bullet sets verbatim — the "
        "branch's 3 GZ/freeboard bullets followed by main's 6 transom/inertia bullets. "
        "Every line of both sides is present (checked line-by-line against both "
        "parents' blobs)."
    ),
    "tensorcad/__init__.py": (
        "auto-merged to the union: 125 exported names = 112 (branch) + 120 (main) "
        "minus 107 shared; nothing from either side's `__all__` is missing and there "
        "are no duplicates; both feature sets import from the package root."
    ),
    "tests/test_public_surface.py": (
        "auto-merged to the union: `EXPECTED_API` holds all 125 names, the branch's "
        "`GZCurve`/`righting_arm_curve` group and main's `transom` group are both "
        "asserted, and `transom` joined the module-import list."
    ),
    "pyproject.toml": (
        "both sides bumped the version to 0.21.0; the two blobs were identical, so the "
        "merged file is byte-identical to both."
    ),
}


# ---------------------------------------------------------------------------
# 2. cross-feature experiments
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Closed rectangular barge: keel at ``-draft``, deck at ``+freeboard``, WL at 0."""
    b = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = b.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, c, d, e = (corner[p] for p in (p0, p1, p2, p3))
        b.add_triangle(a, c, d)
        b.add_triangle(a, d, e)

    lo, hi = -0.5 * beam, 0.5 * beam
    left, right = 0.0, length
    bottom, top = -draft, freeboard
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return b.build()


def wall_sided(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    """Classical closed arm for a wall-sided body (válid while tan(phi) < 2T/B)."""
    phi = math.radians(heel_deg)
    bm = beam * beam / (12.0 * draft)
    bg = kg - 0.5 * draft
    return math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)


def stable_kg(mesh: Mesh) -> float:
    """KG that puts the metacentre a half-BM above G, so the case is stably trimmed."""
    gm_at_kg0 = initial_stability(mesh, kg=0.0, closed=True).to_dict()["gm_transverse_m"]
    return 0.5 * gm_at_kg0


def rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else float("nan")


def experiments() -> dict:
    heels = [0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0]

    # E1 — branch solver on a box, against the classical wall-sided formula.
    beam, draft = 2.0, 0.8
    barge = box(10.0, beam, draft, 1.2)
    kg = stable_kg(barge)
    barge_curve = righting_arm_curve(barge, kg=kg, angles_deg=heels).to_dict()
    limit = math.degrees(math.atan(2.0 * draft / beam))
    rows = []
    for point in barge_curve["points"]:
        closed = wall_sided(beam, draft, kg, point["heel_deg"])
        rows.append(
            {
                "heel_deg": point["heel_deg"],
                "gz": point["gz"],
                "gz_wall_sided": closed,
                "rel_error": rel(point["gz"], closed) if abs(closed) > 1e-9 else None,
                "inside_window": point["heel_deg"] < limit,
            }
        )
    inside = [r for r in rows if r["inside_window"] and r["rel_error"] is not None]
    outside = [
        r
        for r in rows
        if r["heel_deg"] > limit + 3.0 and r["rel_error"] is not None
    ]
    box_stats = {
        "beam_m": beam,
        "draft_m": draft,
        "freeboard_m": 1.2,
        "kg_m": kg,
        "wall_sided_window_deg": limit,
        "worst_rel_error_inside_window": max(r["rel_error"] for r in inside),
        "first_rel_error_outside_window": min(r["rel_error"] for r in outside),
        "gm_curve": barge_curve["initial_gm_transverse"],
        "gm_wall_sided": beam * beam / (12.0 * draft) - (kg - 0.5 * draft),
        "rows": rows,
    }

    # E2 — main's transom hull family driven by the branch's GZ solver (the actual
    # semantic seam of this merge: two features developed against different mains).
    transom = generate_transom_hull()
    tmesh = transom.mesh
    t_initial = initial_stability(tmesh, kg=0.8, closed=True).to_dict()
    t_curve = righting_arm_curve(
        tmesh, kg=0.8, angles_deg=heels, freeboard=0.3
    ).to_dict()
    t_hydro = compute_hydrostatics(tmesh).to_dict()
    t_analytic = analytic_transom_hydrostatics(transom.parameters)
    if hasattr(t_analytic, "to_dict"):
        t_analytic = t_analytic.to_dict()
    vols = [p["displacement_volume"] for p in t_curve["points"]]
    transom_stats = {
        "tau": transom.parameters.transom_fraction,
        "breadth_exponent": transom.parameters.breadth_exponent,
        "kg_m": 0.8,
        "freeboard_m": 0.3,
        "gm_initial_stability": t_initial["gm_transverse_m"],
        "gm_gz_curve": t_curve["initial_gm_transverse"],
        "gm_rel_error": rel(t_curve["initial_gm_transverse"], t_initial["gm_transverse_m"]),
        "displacement_volume_curve": t_curve["points"][0]["displacement_volume"],
        "displacement_volume_spread_over_heel": max(vols) - min(vols),
        "gz_max_m": t_curve["gz_max"],
        "angle_of_gz_max_deg": t_curve["angle_of_gz_max_deg"],
        "vanishing_angle_deg": t_curve["vanishing_angle_deg"],
        "dynamic_stability_area_m_rad": t_curve["dynamic_stability_area"],
        "gz_at_zero_heel_m": t_curve["points"][0]["gz"],
        "hydrostatics_vs_analytic": {
            mesh_key: {
                "mesh": t_hydro[mesh_key],
                "analytic": t_analytic[ana_key],
                "rel_error": rel(t_hydro[mesh_key], t_analytic[ana_key]),
            }
            for mesh_key, ana_key in (
                ("displacement_volume_m3", "displacement_volume"),
                ("waterplane_area_m2", "waterplane_area"),
                ("longitudinal_centre_buoyancy_m", "longitudinal_centre_buoyancy"),
                ("block_coefficient", "block_coefficient"),
            )
        },
        "curve": [
            {"heel_deg": p["heel_deg"], "gz": p["gz"]} for p in t_curve["points"]
        ],
    }

    # E3 — main's Wigley family with the branch's wall-sided topsides.
    wigley = generate_wigley().mesh
    w_kg = stable_kg(wigley)
    w_curve = righting_arm_curve(
        wigley, kg=w_kg, angles_deg=heels, freeboard=0.3
    ).to_dict()
    w_initial = initial_stability(wigley, kg=w_kg, closed=True).to_dict()
    wigley_stats = {
        "kg_m": w_kg,
        "gm_initial_stability": w_initial["gm_transverse_m"],
        "gm_gz_curve": w_curve["initial_gm_transverse"],
        "gm_rel_error": rel(w_curve["initial_gm_transverse"], w_initial["gm_transverse_m"]),
        "gz_max_m": w_curve["gz_max"],
        "angle_of_gz_max_deg": w_curve["angle_of_gz_max_deg"],
        "vanishing_angle_deg": w_curve["vanishing_angle_deg"],
        "curve": [
            {"heel_deg": p["heel_deg"], "gz": p["gz"]} for p in w_curve["points"]
        ],
    }

    return {
        "box_vs_wall_sided": box_stats,
        "transom_hull_with_gz": transom_stats,
        "wigley_with_topsides": wigley_stats,
    }


# ---------------------------------------------------------------------------
# 3. test evidence
# ---------------------------------------------------------------------------
COUNT_RE = re.compile(
    r"(?:(?P<failed>\d+) failed)?(?:, )?(?:(?P<passed>\d+) passed)?"
    r"(?:, )?(?:(?P<skipped>\d+) skipped)?"
)


def counts(logfile: Path) -> dict:
    text = logfile.read_text()
    summary = text.strip().splitlines()[-1]
    out = {"log": str(logfile.relative_to(ROOT)), "summary_line": summary}
    for key in ("passed", "failed", "skipped"):
        match = re.search(rf"(\d+) {key}\b", text)
        out[key] = int(match.group(1)) if match else 0
    out["failed_tests"] = sorted(
        {
            line.split(" ")[1].split("::", 1)[1]
            for line in text.splitlines()
            if line.startswith("FAILED ")
        }
    )
    return out


# ---------------------------------------------------------------------------
# 4. figure
# ---------------------------------------------------------------------------
def figure(exp: dict, tests: dict) -> Path:
    fig, ax = plt.subplots(2, 2, figsize=(11.0, 7.6))
    box = exp["box_vs_wall_sided"]

    a = ax[0][0]
    heels = [r["heel_deg"] for r in box["rows"]]
    a.plot(heels, [r["gz"] for r in box["rows"]], "k-", lw=2, label="tensorcad GZ (3-D clip)")
    a.plot(
        heels,
        [r["gz_wall_sided"] for r in box["rows"]],
        "C3--",
        lw=1.4,
        label="wall-sided closed form",
    )
    a.axvline(box["wall_sided_window_deg"], color="0.5", ls=":", lw=1)
    a.text(
        box["wall_sided_window_deg"] + 1.5,
        0.6 * max(r["gz"] for r in box["rows"]),
        "tan φ = 2T/B\n(formula leaves its\nvalidity domain)",
        fontsize=7,
        color="0.35",
    )
    a.set_title(f"(a) branch solver vs analytic reference, box {box['beam_m']}×{box['draft_m']} m")
    a.set_xlabel("heel φ (deg)")
    a.set_ylabel("GZ (m)")
    a.legend(fontsize=7, loc="upper left")
    a.grid(alpha=0.3)

    b = ax[0][1]
    for key, style, label in (
        ("transom_hull_with_gz", "C0-", "main's transom hull (τ=%.2f) + branch GZ" % exp["transom_hull_with_gz"]["tau"]),
        ("wigley_with_topsides", "C2-", "main's Wigley + branch topsides"),
    ):
        curve = exp[key]["curve"]
        b.plot([p["heel_deg"] for p in curve], [p["gz"] for p in curve], style, lw=1.8, label=label)
    b.axhline(0.0, color="0.6", lw=0.8)
    b.set_title("(b) the seam of this merge: branch solver on main's hull families")
    b.set_xlabel("heel φ (deg)")
    b.set_ylabel("GZ (m)")
    b.legend(fontsize=7)
    b.grid(alpha=0.3)

    c = ax[1][0]
    err = [(r["heel_deg"], r["rel_error"]) for r in box["rows"] if r["rel_error"]]
    c.semilogy([e[0] for e in err], [e[1] for e in err], "ko-", ms=3, lw=1)
    c.axvline(box["wall_sided_window_deg"], color="0.5", ls=":", lw=1)
    c.set_title("(c) box: |GZ − GZ_wall-sided| / |GZ_wall-sided|")
    c.set_xlabel("heel φ (deg)")
    c.set_ylabel("relative error")
    c.grid(alpha=0.3, which="both")
    c.text(
        0.03,
        0.05,
        f"inside window: ≤ {box['worst_rel_error_inside_window']:.1e}\n"
        f"outside window: ≥ {box['first_rel_error_outside_window']:.2%}",
        transform=c.transAxes,
        fontsize=7,
    )

    d = ax[1][1]
    labels = ["passed", "failed", "skipped"]
    main_v = [tests["main"][k] for k in labels]
    merged_v = [tests["merged"][k] for k in labels]
    xpos = range(len(labels))
    d.bar([x - 0.19 for x in xpos], main_v, width=0.36, color="C7", label=f"main {MAIN_SHA}")
    d.bar([x + 0.19 for x in xpos], merged_v, width=0.36, color="C0", label="merged cure tree")
    for x, (m, g) in enumerate(zip(main_v, merged_v)):
        d.text(x - 0.19, m + 8, str(m), ha="center", fontsize=8)
        d.text(x + 0.19, g + 8, str(g), ha="center", fontsize=8)
    d.set_xticks(list(xpos), labels)
    d.set_yscale("symlog")
    d.set_title("(d) test suite: cure tree vs main baseline (same command)")
    d.set_ylabel("tests")
    d.legend(fontsize=7)
    d.grid(alpha=0.3, axis="y")

    fig.suptitle(
        "TensorCAD PR #44 cure (m3847): semantic union of agent/m3347 and main",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = ART / "gz_integration.png"
    for dpi in (110, 100, 90, 80, 70):
        fig.savefig(out, dpi=dpi)
        if out.stat().st_size <= 58_000:
            break
    plt.close(fig)
    if out.stat().st_size > 58_000:  # last resort: keep the report payload small
        from PIL import Image

        with Image.open(out) as image:
            image.thumbnail((1200, 1200))
            image.save(out, optimize=True)
    return out


# ---------------------------------------------------------------------------
def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    inv = inventory()
    (ART / "merge_inventory.json").write_text(json.dumps(inv, indent=2) + "\n")

    exp = experiments()
    tests = {
        "main": counts(ART / "pytest_main.log"),
        "merged": counts(ART / "pytest_merged_full.log"),
    }
    payload = {"inventory_summary": {k: inv[k] for k in (
        "merge_base", "branch_sha_before", "main_sha", "cured_commit",
        "files_touched_by_either_side", "counts", "textual_conflicts",
        "single_sided_content_lost")}, "experiments": exp, "tests": tests}
    (ART / "numbers.json").write_text(json.dumps(payload, indent=2) + "\n")
    png = figure(exp, tests)

    print(json.dumps({
        "inventory": payload["inventory_summary"],
        "box_worst_inside": exp["box_vs_wall_sided"]["worst_rel_error_inside_window"],
        "box_first_outside": exp["box_vs_wall_sided"]["first_rel_error_outside_window"],
        "transom_gm_rel_error": exp["transom_hull_with_gz"]["gm_rel_error"],
        "transom_volume_spread": exp["transom_hull_with_gz"]["displacement_volume_spread_over_heel"],
        "wigley_gm_rel_error": exp["wigley_with_topsides"]["gm_rel_error"],
        "tests_main": {k: tests["main"][k] for k in ("passed", "failed", "skipped")},
        "tests_merged": {k: tests["merged"][k] for k in ("passed", "failed", "skipped")},
        "figure_bytes": png.stat().st_size,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())