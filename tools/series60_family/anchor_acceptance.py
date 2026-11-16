#!/usr/bin/env python3
"""Executable acceptance verdict for ``GAP-CAD-T2-PARAM-SERIES``.

The ledger entry for this anchor carries three acceptance clauses::

    Series60 五档 Cb 生成 + 静水力对照公开值；gate-power 突变测试沿用 >=80%；文档

Until now those three clauses were only checkable by reading prose: the five-Cb
numbers lived in a report, the gate-power record sat outside the tree, and the
"published comparison" was a *closure* check (the published particulars are the
inputs) whose machine-readable role field could be -- and once was -- read as an
external verification.  A board audit (m17308, findings F1/F2/F3/F4) called
exactly that out: ``external anchor = WITHHELD``, no registered case in
``tensor_benchmark_cases``, artefacts that disagreed.

This tool turns the anchor's acceptance into something a machine can adjudicate
and a reader cannot over-read:

* every clause reports ``pass`` / ``fail`` / ``not_run`` / ``withheld``;
* every comparison carries an evidence **grade** drawn from a closed set, and a
  ``closure`` comparison is never allowed to report itself as external;
* a missing gate-power record is ``not_run`` -- *unknown*, never "passed";
* the external shape anchor (Todd's tabulated station offsets) is declared
  ``WITHHELD`` in the artifact itself instead of only in a docstring.

Outputs (relative to the repository root)::

    artifacts/series60_anchor_acceptance/anchor_acceptance.json
    artifacts/series60_anchor_acceptance/anchor_acceptance.csv
    artifacts/series60_anchor_acceptance/figures/anchor_acceptance.png

Usage::

    python3 tools/series60_family/anchor_acceptance.py
    python3 tools/series60_family/anchor_acceptance.py --mutations /path/to/gate_power.json

Nothing is typed by hand: the published particulars are read from
``docs/data/series60_published_parents.json`` (which carries its own provenance)
and every model number is computed here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

DATA = ROOT / "docs" / "data" / "series60_published_parents.json"
DEFAULT_MUTATIONS = ROOT / "artifacts" / "series60_gate_power" / "gate_power_series60.json"
OUT = ROOT / "artifacts" / "series60_anchor_acceptance"

#: Principal dimensions of the systematic-series model the ladder is quoted for.
L, B, T = 10.0, 1.5, 0.6
#: Grid the anchor's numbers are quoted on (161 stations x 41 levels, as shipped).
GRID = (161, 41)

#: Predeclared bars.  These are fixed *before* the numbers are inspected, and the
#: artifact records them so a later reader can see they were not moved.
TOLERANCE_PUBLISHED = 0.03      # model vs published particular
#: Length, beam and draft are honoured by construction, so a 0.1 % gate is honest.
TOLERANCE_DIMENSION = 0.001
#: ``Cb`` is *integrated off the loft*, so it carries the grid's first-order
#: deficit (the same 2 % the module already declares as its discretisation bar).
#: Holding it to 0.1 % would fail every converged hull for the wrong reason.
TOLERANCE_REALISED_CB = 0.02
GATE_POWER_MINIMUM = 0.80       # anchor's mutation-kill bar

ANCHOR = "GAP-CAD-T2-PARAM-SERIES"

#: Closed evidence-grade vocabulary.  ``closure`` = the published value is an
#: *input* to the shape solve, so agreement measures algebra, not fidelity.
#: ``external`` = the reference is independent of the model's inputs.
GRADES = ("external", "closure", "withheld")

#: Clause statuses.  ``not_run`` is deliberately distinct from ``pass``.
#: ``unverified`` = a record exists but cannot be read as evidence of a real run
#: (hand-transcribed, or an injection that could not have changed behaviour).
STATUSES = ("pass", "fail", "not_run", "unverified", "withheld")

#: Which published particular plays which role.  Recorded per quantity so the
#: closure/external mistake the audit found cannot be re-introduced silently.
QUANTITY_GRADES = {
    "block_coefficient": "closure",        # requested as Cb, recovered by integration
    "prismatic_coefficient": "closure",    # read off the published Cp ladder
    "midship_coefficient": "closure",      # Cm = Cb / Cp follows from that ladder
    "waterplane_coefficient": "closure",   # solved for, from the published Cwp
    "lcb_over_length": "closure",          # solved for, from the published LCB
}

REQUIRED_DOCS = (
    "docs/series60.md",
    "docs/series60_verification.md",
    "docs/series60_fidelity.md",
    "docs/series60_bulb_transom.md",
    "docs/series60_gate_power.md",
)


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_published(path: Path = DATA) -> dict:
    return json.loads(path.read_text())


def relative_error(model: float, reference: float) -> float:
    return model / reference - 1.0


# --------------------------------------------------------------------- clause 1
def geometry_precondition(model, request: dict, tolerance: float = TOLERANCE_DIMENSION,
                          cb_tolerance: float = TOLERANCE_REALISED_CB) -> dict:
    """Realised L/B/T and Cb against the request -- the geometry gate.

    The platform's rule is that a hull's Cb/L/B must be validated *before* any
    hydrostatic or resistance number is quoted off it; a loft that misses its
    moulded breadth silently moves every downstream coefficient.
    """
    hydro = model.hydrostatics
    box_x = [v[0] for v in model.mesh.vertices]
    box_y = [v[1] for v in model.mesh.vertices]
    box_z = [v[2] for v in model.mesh.vertices]
    length = max(box_x) - min(box_x)
    beam = 2.0 * max(abs(y) for y in box_y)
    draft = max(box_z) - min(box_z)
    realised = {
        "length": length,
        "beam": beam,
        "draft": draft,
        # Cb realised from the integrated volume, not from the requested parameter.
        "block_coefficient": hydro.displacement_volume / (length * beam * draft),
    }
    requested = {
        "length": request["length"],
        "beam": request["beam"],
        "draft": request["draft"],
        "block_coefficient": request["cb"],
    }
    rows = []
    for key in ("length", "beam", "draft", "block_coefficient"):
        err = relative_error(realised[key], requested[key])
        bar = cb_tolerance if key == "block_coefficient" else tolerance
        rows.append({"quantity": key, "realised": realised[key],
                     "requested": requested[key], "relative": err,
                     "bar": bar,
                     "within_bar": abs(err) <= bar,
                     "grade": "closure"})
    return {
        "id": "geometry_precondition_cb_l_b",
        "statement": "realised length / beam / draft / Cb match the request before "
                     "any coefficient is quoted",
        "tolerance": tolerance,
        "cb_tolerance": cb_tolerance,
        "rows": rows,
        "passed": all(row["within_bar"] for row in rows),
    }


# --------------------------------------------------------------------- clause 2
def grade_gradient(published: dict, nx: int = GRID[0], nz: int = GRID[1]) -> dict:
    """Five Cb parents generated; coefficients against the published particulars."""
    parents = []
    for parent in published["parents"]:
        cb = parent["cb"]
        request = dict(length=L, beam=B, draft=T, cb=cb, nx=nx, nz=nz,
                       cwp=parent["cwp"], lcb=parent["lcb_over_length"])
        model = generate_series60(Series60Parameters(**request))
        hydro = model.hydrostatics
        closed = analytic_series60_hydrostatics(Series60Parameters(**request))
        model_values = {
            "block_coefficient": hydro.block_coefficient,
            "prismatic_coefficient": hydro.prismatic_coefficient,
            "midship_coefficient": hydro.maximum_section_area / (B * T),
            "waterplane_coefficient": hydro.waterplane_area / (L * B),
            "lcb_over_length": hydro.longitudinal_centre_buoyancy / L,
        }
        published_values = {
            "block_coefficient": cb,
            "prismatic_coefficient": parent["cp"],
            "midship_coefficient": parent["cm"],
            "waterplane_coefficient": parent["cwp"],
            "lcb_over_length": parent["lcb_over_length"],
        }
        comparisons = {}
        for key, value in model_values.items():
            err = relative_error(value, published_values[key])
            comparisons[key] = {
                "model": value,
                "published": published_values[key],
                "relative": err,
                "grade": QUANTITY_GRADES[key],
                "within_bar": abs(err) <= TOLERANCE_PUBLISHED,
            }
        parents.append({
            "cb": cb,
            "comparisons": comparisons,
            "displacement_volume": hydro.displacement_volume,
            "vcb_over_draft": (hydro.vertical_centre_buoyancy + T) / T,
            "closed_form_displacement_volume": closed["displacement_volume"],
            "mesh_vs_closed_form_volume": relative_error(
                hydro.displacement_volume, closed["displacement_volume"]),
            "validation_passed": bool(model.validation.passed),
            "validation_errors": list(model.validation.errors),
            "triangles": len(model.mesh.faces),
            "geometry": geometry_precondition(model, request),
        })
    worst = max(
        (abs(c["relative"]), parent["cb"], key)
        for parent in parents
        for key, c in parent["comparisons"].items()
    )
    return {
        "id": "five_cb_generation_and_published_hydrostatics",
        "statement": "Cb 0.60/0.65/0.70/0.75/0.80 generate a valid closed hull and the "
                     "five particulars agree with the published ladder within 3%",
        "tolerance": TOLERANCE_PUBLISHED,
        "grid": {"nx": nx, "nz": nz},
        "parents": parents,
        "worst_relative": worst[0],
        "worst_quantity": worst[2],
        "worst_cb": worst[1],
        "grade": "closure",
        "grade_note": (
            "All five particulars are inputs to, or arithmetic consequences of, the "
            "published Cp ladder (_PARENT_FORMS): Cb/Cwp/LCB are solved for and Cp/Cm "
            "are read from the ladder. Agreement therefore measures the algebra, not "
            "fidelity to the real ship. This clause is a closure check and must never "
            "be scored as an external verification."
        ),
        "passed": all(
            parent["validation_passed"]
            and parent["geometry"]["passed"]
            and all(c["within_bar"] for c in parent["comparisons"].values())
            for parent in parents
        ),
    }


#: What a kill claim must carry before it can be read as evidence: the injection
#: that was made (``mutant``) and the measured exit status of the run (``run``).
#: A record without these is a summary of claims, not a mutation run.
PROVENANCE_MUTANT_KEYS = ("id", "kind", "file", "old", "new", "claim")


def provenance_defects(entry: dict) -> list[str]:
    """Why one record entry cannot be counted as a mutation kill.

    The failure this guards against is concrete: a hand-written record with
    ``{"status": "ok", "run": {"killed": true}}`` per mutant used to be scored as
    an *external* kill rate, sha256-sealed into the verdict.  Nothing was broken,
    nothing was run, and the acceptance clause read as passed.
    """
    defects: list[str] = []
    mutant = entry.get("mutant")
    if not isinstance(mutant, dict):
        return ["no 'mutant' block: the record does not say what was broken"]
    defects += [f"mutant.{key} missing" for key in PROVENANCE_MUTANT_KEYS
                if not mutant.get(key)]
    if mutant.get("old") and mutant.get("new") and mutant["old"] == mutant["new"]:
        defects.append("mutant.old == mutant.new: a no-op injection cannot change behaviour")
    run = entry.get("run")
    if not isinstance(run, dict):
        return defects + ["no 'run' block: no measured exit status"]
    returncode = run.get("returncode")
    killed = run.get("killed")
    if not isinstance(returncode, int):
        defects.append("run.returncode missing: a kill must be a measured exit status")
    if not isinstance(killed, bool):
        defects.append("run.killed is not a boolean")
    elif isinstance(returncode, int) and killed != (returncode != 0):
        defects.append(f"run.killed={killed} contradicts run.returncode={returncode}")
    if not run.get("summary"):
        defects.append("run.summary empty: no pytest summary line")
    if not isinstance(run.get("failed"), list):
        defects.append("run.failed missing: which test killed it is not recorded")
    wall = run.get("wall_seconds")
    if not (isinstance(wall, (int, float)) and wall > 0):
        defects.append("run.wall_seconds missing or zero")
    return defects


# --------------------------------------------------------------------- clause 3
def gate_power_clause(record_path: Path, minimum: float = GATE_POWER_MINIMUM) -> dict:
    """Mutation kill rate from a machine-readable record, or ``not_run``."""
    clause = {
        "id": "gate_power_mutation_kill_rate",
        "statement": f"Series 60 mutation kill rate >= {minimum:.0%}",
        "minimum": minimum,
        "record": None,
        "grade": "external",
    }
    if not record_path.exists():
        clause.update({
            "status": "not_run",
            "passed": None,
            "kill_rate": None,
            "detail": (f"no mutation record at {record_path}; the kill rate is UNKNOWN, "
                       "which is not the same as passed -- regenerate it with "
                       "tools/hull_gate_power/mutate.py --only <S60-*>"),
        })
        return clause
    record = json.loads(record_path.read_text())
    killed = survived = 0
    survivors: list[str] = []
    unverifiable: dict[str, list[str]] = {}
    for mutant_id, entry in sorted(record.items()):
        if entry.get("status") != "ok":
            continue
        defects = provenance_defects(entry)
        if defects:
            unverifiable[mutant_id] = defects
            continue
        if entry["run"]["killed"]:
            killed += 1
        else:
            survived += 1
            survivors.append(mutant_id)
    total = killed + survived
    rate = killed / total if total else None
    clause.update({
        "record": str(record_path.relative_to(ROOT)) if record_path.is_relative_to(ROOT)
                  else str(record_path),
        "record_sha256": sha256_of(record_path),
        "killed": killed,
        "survived": survived,
        "total": total,
        "survivors": survivors,
        "kill_rate": rate,
        "unverifiable": unverifiable,
        "passed": bool(rate is not None and rate >= minimum and not unverifiable),
        "status": ("unverified" if unverifiable
                   else "pass" if (rate is not None and rate >= minimum)
                   else "fail"),
    })
    if unverifiable:
        first = sorted(unverifiable)[0]
        clause["detail"] = (
            f"{len(unverifiable)} of {len(unverifiable) + total} entries carry no "
            f"injection/exit-status provenance (first: {first}: "
            f"{'; '.join(unverifiable[first])}). A record that does not say what was "
            "broken and what the run returned is a summary of claims, not a mutation "
            "run: regenerate it with tools/hull_gate_power/mutate.py --only <S60-*>."
        )
    return clause


# --------------------------------------------------------------------- clause 4
def documentation_clause(paths=REQUIRED_DOCS) -> dict:
    rows = []
    for rel in paths:
        path = ROOT / rel
        rows.append({"path": rel, "exists": path.exists(),
                     "bytes": path.stat().st_size if path.exists() else 0})
    passed = all(row["exists"] and row["bytes"] > 0 for row in rows)
    return {
        "id": "documentation",
        "statement": "the family has a reference page, verification page, fidelity page, "
                     "bulb/transom page and gate-power page",
        "rows": rows,
        "passed": passed,
        "status": "pass" if passed else "fail",
    }


# --------------------------------------------------------------------- clause 5
def external_anchor_clause() -> dict:
    """The honest gap, stated in the artifact rather than only in prose."""
    return {
        "id": "external_shape_anchor_todd_offsets",
        "statement": "station offsets reproduce Todd's tabulated Series 60 offsets",
        "status": "withheld",
        "passed": None,
        "grade": "withheld",
        "detail": (
            "No Todd (1951) offset table is available in this sandbox and no offset table "
            "was ever landed for the family: tests/test_series60_stations_published_comparison.py "
            "compares the mesh against the module's own analytic station law, which is a "
            "discretisation check, not a shape validation. The family is an analytic stand-in "
            "with the right integral particulars; it does not reproduce Todd's ordinates, and "
            "no resistance claim should be built on it as if it did. Recorded as WITHHELD, "
            "not as passed."
        ),
    }


# --------------------------------------------------------------------- assembly
def build_report(published: dict, mutations: Path) -> dict:
    geometry = grade_gradient(published)
    parents_geometry = [p["geometry"] for p in geometry["parents"]]
    geometry_clause = {
        "id": "geometry_precondition_cb_l_b",
        "statement": "realised length / beam / draft / Cb match the request before any "
                     "coefficient is quoted",
        "tolerance": TOLERANCE_DIMENSION,
        "parents": parents_geometry,
        "passed": geometry["passed"] and all(g["passed"] for g in parents_geometry),
    }
    geometry_clause["status"] = "pass" if geometry_clause["passed"] else "fail"

    clauses = [geometry_clause, geometry, gate_power_clause(mutations),
               documentation_clause(), external_anchor_clause()]
    for clause in clauses:
        if "status" not in clause:
            clause["status"] = "pass" if clause["passed"] else "fail"

    graded = [c for c in clauses
              if c["status"] not in ("withheld", "not_run", "unverified")]
    failed = [c["id"] for c in graded if c["status"] == "fail"]
    unrun = [c["id"] for c in clauses if c["status"] == "not_run"]
    unverifiable = [c["id"] for c in clauses if c["status"] == "unverified"]
    withheld = [c["id"] for c in clauses if c["status"] == "withheld"]
    if failed:
        verdict = "fail"
    elif unrun or unverifiable:
        verdict = "incomplete"
    else:
        verdict = "pass"

    # Self-consistency: a closure comparison may never be reported as external.
    overstated = [
        (clause["id"], key)
        for clause in clauses
        for key, value in _graded_items(clause)
        if value == "external" and clause["id"] != "gate_power_mutation_kill_rate"
    ]

    return {
        "schema": "tensorcad.anchor_acceptance/1",
        "anchor": ANCHOR,
        "mission": "m22471",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "series60_module_sha256": sha256_of(ROOT / "tensorcad" / "series60.py"),
            "published_data": str(DATA.relative_to(ROOT)),
            "published_data_sha256": sha256_of(DATA),
            "tool_sha256": sha256_of(Path(__file__).resolve()),
            "repro_cmd": "python3 tools/series60_family/anchor_acceptance.py",
            "test_cmd": "python3 -m pytest tests/test_series60_anchor_acceptance.py -q",
        },
        "model": {"length": L, "beam": B, "draft": T, "grid": {"nx": GRID[0], "nz": GRID[1]}},
        "bars": {
            "published_relative": TOLERANCE_PUBLISHED,
            "dimension_relative": TOLERANCE_DIMENSION,
            "gate_power_minimum": GATE_POWER_MINIMUM,
            "platform_accuracy_rule": "a key quantity off by more than 15% is unacceptable; "
                                      "this family holds itself to 3%",
        },
        "grades": {key: value for key, value in QUANTITY_GRADES.items()},
        "clauses": clauses,
        "verdict": verdict,
        "failed_clauses": failed,
        "clauses_not_run": unrun,
        "clauses_unverified": unverifiable,
        "clauses_withheld": withheld,
        "closure_where_external_should_be": overstated,
        "overall_evidence": (
            "The five-Cb gradient and its published comparison are REAL arithmetic on a "
            "REAL closure; the shape is an analytic stand-in, so the external anchor is "
            "WITHHELD. The verdict is about the code doing what it claims, not about the "
            "hull matching Todd's lines."
        ),
    }


def _graded_items(clause: dict):
    """Yield ``(quantity, grade)`` pairs declared by a clause."""
    if "comparisons" in clause:
        for key, value in clause["comparisons"].items():
            yield key, value["grade"]
    for parent in clause.get("parents", []):
        for key, value in parent.get("comparisons", {}).items():
            yield key, value["grade"]
        for row in parent.get("rows", []):
            yield row["quantity"], row.get("grade", "closure")
    if "grade" in clause and clause["grade"] in GRADES:
        yield clause["id"], clause["grade"]


# ------------------------------------------------------------------------ write
def write_csv(report: dict, path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cb", "quantity", "model", "published", "relative",
                         "grade", "within_bar"])
        for clause in report["clauses"]:
            if clause["id"] != "five_cb_generation_and_published_hydrostatics":
                continue
            for parent in clause["parents"]:
                for key, comp in sorted(parent["comparisons"].items()):
                    writer.writerow([parent["cb"], key, repr(comp["model"]),
                                     repr(comp["published"]), repr(comp["relative"]),
                                     comp["grade"], comp["within_bar"]])


def figure(report: dict, path: Path) -> None:
    gradient = next(c for c in report["clauses"]
                    if c["id"] == "five_cb_generation_and_published_hydrostatics")
    quantities = sorted(QUANTITY_GRADES)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    colors = plt.cm.tab10.colors
    for index, key in enumerate(quantities):
        xs = [p["cb"] for p in gradient["parents"]]
        ys = [100.0 * p["comparisons"][key]["relative"] for p in gradient["parents"]]
        axes[0].plot(xs, ys, marker="o", ms=4, lw=1.4, color=colors[index],
                     label=f"{key} ({QUANTITY_GRADES[key]})")
    axes[0].axhline(0.0, color="k", lw=0.8)
    axes[0].axhline(100 * TOLERANCE_PUBLISHED, color="crimson", ls="--", lw=1.0,
                    label=f"+{100 * TOLERANCE_PUBLISHED:.0f}% bar")
    axes[0].axhline(-100 * TOLERANCE_PUBLISHED, color="crimson", ls="--", lw=1.0)
    axes[0].set_xlabel("Cb (Todd parent form)")
    axes[0].set_ylabel("model / published - 1  [%]")
    axes[0].set_title("Five-Cb gradient vs published particulars\n(all five are closure grades)")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)

    labels, values, tints = [], [], []
    for clause in report["clauses"]:
        labels.append(clause["id"].replace("_", "\n", 1))
        status = clause["status"]
        values.append({"pass": 1.0, "fail": 0.0, "not_run": 0.5,
                       "withheld": 0.25}[status])
        tints.append({"pass": "#2a9d8f", "fail": "#e63946", "not_run": "#f4a261",
                      "withheld": "#8d99ae"}[status])
    axes[1].barh(labels, values, color=tints)
    axes[1].set_xlim(0, 1.0)
    axes[1].set_xlabel("status score (1 = pass, 0.5 = not run, 0.25 = withheld, 0 = fail)")
    axes[1].set_title(f"Anchor clause verdict: {report['verdict'].upper()}")
    for tick, clause in zip(axes[1].get_yticks(), report["clauses"]):
        axes[1].text(0.02, tick, clause["status"], va="center", fontsize=8, color="white")
    axes[1].grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", default=str(DEFAULT_MUTATIONS))
    parser.add_argument("--out", default=str(OUT))
    args = parser.parse_args()

    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    (out / "figures").mkdir(parents=True, exist_ok=True)

    published = load_published()
    report = build_report(published, Path(args.mutations))

    (out / "anchor_acceptance.json").write_text(json.dumps(report, indent=2))
    write_csv(report, out / "anchor_acceptance.csv")
    figure(report, out / "figures" / "anchor_acceptance.png")

    for clause in report["clauses"]:
        print(f"[{clause['status']:>8}] {clause['id']}", flush=True)
    print(f"[verdict] {report['verdict']}", flush=True)
    if report["clauses_not_run"]:
        print(f"[not_run] {', '.join(report['clauses_not_run'])} -> UNKNOWN, not passed",
              flush=True)
    print(f"[written] {out}", flush=True)
    # A failed clause is a real failure and must fail the process; an unrun clause
    # is unknown, which is also not a pass.
    return 0 if report["verdict"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())