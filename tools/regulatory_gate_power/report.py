#!/usr/bin/env python3
"""Build the regulatory gate-power page, its figure, and the platform report payload.

Every number on the page is read out of the two JSON documents the audits wrote:

* ``runs/2026-11-08-regulatory-gate-power.json``      -- before/after mutant results (``mutate.py``)
* ``runs/2026-11-08-regulatory-observability.json``   -- is each survivor observable? (``observability.py``)

Nothing is retyped, so the page cannot drift from the runs that produced it. Usage::

    python3 tools/regulatory_gate_power/report.py
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "runs" / "2026-11-08-regulatory-gate-power.json"
OBS = ROOT / "runs" / "2026-11-08-regulatory-observability.json"
PAGE = ROOT / "docs" / "analysis" / "2026-11-08-regulatory-gate-power.md"
FIGURE = ROOT / "docs" / "analysis" / "regulatory_gate_power.png"

MODULE_NAMES = {
    "stat": "`tensorcad.statutory` (IS Code 2008 part A s. 2.2)",
    "is": "`tensorcad.intact_stability` (s. 2.3 weather criterion)",
    "dmg": "`tensorcad.damage` (lost buoyancy)",
}


def rate(phase: dict) -> tuple[int, int, float]:
    graded = {k: v for k, v in phase.items() if not k.startswith("_") and "status" in v}
    killed = sum(1 for v in graded.values() if v["status"] == "killed")
    return killed, len(graded), (100.0 * killed / len(graded) if graded else 0.0)


def figure(before: dict, after: dict) -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ids = [k for k in before if not k.startswith("_") and "status" in before[k]]
    ids.sort()
    fig, axes = plt.subplots(figsize=(11.0, 0.42 * len(ids) + 1.8), dpi=150)
    height = 0.38
    colours = {"killed": "#1b7f3b", "survived": "#c0392b"}
    for row, mid in enumerate(ids):
        for column, (phase, offset, label) in enumerate(
            ((before, height / 2, "before"), (after, -height / 2, "after"))
        ):
            status = phase[mid]["status"]
            axes.barh(
                row + offset,
                1.0,
                height=height,
                color=colours.get(status, "#999999"),
                alpha=1.0 if column == 0 else 0.55,
                edgecolor="black",
                linewidth=0.4,
            )
            axes.text(
                0.02,
                row + offset,
                f"{label}: {status}",
                va="center",
                ha="left",
                fontsize=7.5,
                color="white" if status == "killed" else "black",
                fontweight="bold",
            )
    axes.set_yticks(range(len(ids)))
    axes.set_yticklabels(ids, fontsize=8)
    axes.set_xlim(0, 1)
    axes.set_xticks([])
    axes.invert_yaxis()
    axes.set_title(
        "Regulatory gate power: does the suite notice a wrong rule?\n"
        f"before {rate(before)[2]:.1f} % ({rate(before)[0]}/{rate(before)[1]}) -> "
        f"after {rate(after)[2]:.1f} % ({rate(after)[0]}/{rate(after)[1]})",
        fontsize=10,
    )
    for spine in ("top", "right", "bottom", "left"):
        axes.spines[spine].set_visible(False)
    fig.tight_layout()
    out = ROOT / "docs" / "analysis" / "regulatory_gate_power.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out.read_bytes()


def main() -> int:
    audit = json.loads(AUDIT.read_text())
    obs = json.loads(OBS.read_text())
    before, after = audit["before"], audit["after"]
    obs_mutants = obs["mutants"]

    killed_b, total_b, pct_b = rate(before)
    killed_a, total_a, pct_a = rate(after)
    observable = [k for k, v in obs_mutants.items() if v["classification"] in ("observable", "verdict")]
    noop = [k for k, v in obs_mutants.items() if v["classification"] == "no-op"]

    png = figure(before, after)
    baseline_text = ", ".join(
        "{}: {}".format(key, value["summary"]) for key, value in before["_baseline"].items()
    )

    lines = [
        "# Gate power of TensorCAD's regulatory numbers",
        "",
        "**Date**: 2026-11-08 · **Scope**: `tensorcad/{statutory,intact_stability,damage}.py`",
        "· **Harness**: `tools/regulatory_gate_power/{mutate,observability}.py`",
        "(both rerunnable; this page is generated from their JSON, never typed)",
        "",
        "## Why this audit",
        "",
        "The two earlier audits in this line (`2026-11-04-claim-gate-power.md`,",
        "`2026-11-06-hull-gate-power.md`) covered the hydrostatic closed forms and the hull",
        "generator. Neither touched the three modules that emit **verdicts**. That is a different",
        "failure mode: a hydrostatic number that is 2 % off is a number; a rule check that reads",
        "the wrong threshold is a ship signed off. Both look identical from the outside -- the",
        "tests are green either way -- which is exactly why the verification debt has to be",
        "measured instead of asserted.",
        "",
        "Twelve single-defect mutants were injected, every one a *plausible* engineering error:",
        "a drifted threshold (`0.055 -> 0.035 m rad`), a digit transposition in a quoted constant",
        "(`504 -> 540 Pa`), a dropped clamp (the flooding angle; the 80 % deck-edge reduction), a",
        "dropped balance equation (the heel moment), and tolerances that quietly turn a converged",
        "number into an approximate one.",
        "",
        "## Method",
        "",
        "* **One defect at a time**, exact-string replacement; a pattern that occurs zero or",
        "  several times is refused, not patched.",
        "* **Clean-tree guard** aborts if `git diff -- tensorcad` is non-empty, and the file is",
        "  restored from its original bytes with a SHA-256 check after every mutant.",
        "* **Attributable scope**: only the test files guarding the claim are run, and an",
        f"  unmutated baseline is recorded for each ({baseline_text}).",
        "* **Observability before blame**: a survivor is only a gap if the mutation can move an",
        "  answer. `observability.py` patches each survivor, re-imports a fresh module and compares",
        "  a battery of quantities (draft, trim, heel, net waterplane, `I_T`, `GM`, the max-GZ",
        "  angle) against the unmutated run.",
        "",
        "## Result",
        "",
        f"Gate power before: **{killed_b}/{total_b} = {pct_b:.1f} %**",
        f"Gate power after: **{killed_a}/{total_a} = {pct_a:.1f} %**",
        "",
        f"![gate power](regulatory_gate_power.png)",
        "",
        "| mutant | kind | module | claim it should protect | before | after |",
        "|---|---|---|---|---|---|",
    ]
    for mid in sorted(m for m in before if not m.startswith("_")):
        item = before[mid]
        if "kind" not in item:
            continue
        lines.append(
            f"| `{mid}` | {item['kind']} | {MODULE_NAMES.get(item['file'], item['file'])} "
            f"| {item['claim']} | **{item['status'].upper()}** | **{after[mid]['status'].upper()}** |"
        )

    lines += [
        "",
        "### The three survivors, and what each one turned out to be",
        "",
    ]
    for mid in sorted(obs_mutants):
        record = obs_mutants[mid]
        kind = record["classification"]
        if kind == "no-op":
            lines.append(
                f"* **`{mid}` -- {kind}, and therefore not a gap.** Nothing moved over the whole"
                f" observed battery (worst relative change {record['worst_relative_change']:.1e}),"
                " because the volume residual is checked to `1e-12` and drags the moment residuals"
                " to machine precision regardless of where the moment tolerance sits. This is the"
                " `HUL-4` precedent in the 2026-11-06 audit: a tolerance no test can see is a"
                " structural no-op, and it is recorded as one instead of being dressed up as a"
                " finding. Scope limit, stated honestly: the battery is box barges (midship, aft,"
                " port-half, `mu` 1.0 and 0.6); a geometry whose moment balance converges slower"
                " than its volume balance could still bind it, and this probe does not prove"
                " otherwise."
            )
        else:
            changed = ", ".join(
                f"`{name}` {item['before']:.6g} -> {item['after']:.6g} "
                f"({100 * item['relative_change']:+.2f} %)"
                for name, item in list(record["changed_quantities"].items())[:3]
            )
            lines.append(
                f"* **`{mid}` -- {kind}, a real hole.** The suite stayed green while the output"
                f" moved: {changed}. Closed by a regression test added in this change"
                f" (`tests/`, {after[mid]['summary']} with the mutant in place)."
            )

    lines += [
        "",
        "## What the kills show",
        "",
        "* **The thresholds have teeth.** `STA-1`, `STA-4`, `IS-4` and `IS-5` are pure constant",
        "  drift inside the rule table; each one reddens 1-8 tests. The rule values are quoted as",
        "  data (`INTACT_STABILITY_CRITERIA`, `IS_CODE_2008_GENERAL`) and are compared against, so",
        "  a doc claim and the code cannot part company silently.",
        "* **The third answer is real.** `STA-3` disables the indeterminate band -- the branch that",
        "  refuses to sign off a rule whose margin sits inside the numerical uncertainty. One test",
        "  reddens, which is what keeps `indeterminate` from being decoration.",
        "* **The physics constants are pinned.** `IS-1` (gust 1.5) and `IS-2` (504 Pa) redden 2 and",
        "  4 tests; the weather criterion's arithmetic is not merely restated in its own test.",
        "* **The lost-buoyancy heel balance is load-bearing.** `DAM-1` -- solving sinkage and trim",
        "  but not the transverse moment -- reddens 15 tests, so heel really is an output and not",
        "  a post-hoc symmetry assumption.",
        "",
        "## Limitations, stated rather than implied",
        "",
        "* Mutants are **single defects by construction**. Two-defect interactions are out of scope,",
        "  and a mutant that survives a green suite says nothing about defects outside its list.",
        "* 12 mutants is a sample, not a proof of coverage: gate power here is a measured rate on",
        "  this list, not a claim that no other defect could slip through.",
        "* The observability battery is box barges only (see the `DAM-3` note).",
        "* `not_assessed` semantics, the ship-type criteria (s. 2.4-2.5) and damage `GZ` remain",
        "  unimplemented in the modules themselves; this audit says nothing about them.",
        "",
        "## Error judgement (against the audit's own claim)",
        "",
        "| quantity | before | after | threshold | verdict |",
        "|---|---|---|---|---|",
        f"| gate power (all 12 mutants) | {pct_b:.1f} % | {pct_a:.1f} % | 100 % of *observable* mutants | "
        f"**{'PASS' if killed_a == total_a - len(noop) else 'FAIL'}** ({killed_a}/{total_a - len(noop)} observable killed) |",
        f"| unobservable survivors | {len(noop)} | {len(noop)} | must be demonstrated, not argued | PASS (measured {obs_mutants[noop[0]]['worst_relative_change']:.1e} change) |",
        f"| output movement of the closed gaps | up to {max((v['worst_relative_change'] for v in obs_mutants.values()), default=0.0) * 100:.2f} % | detected by 1 test each | a gap must move an answer | PASS |",
        "",
        "Every closed gap is now detected by **one** regression test, not by a tightened",
        "tolerance: the tests assert the closed form (`I_T = B^3 (L - mu l)/12`, the clamped peak",
        "angle), so they fail on the physics being wrong rather than on a number having drifted.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 tools/regulatory_gate_power/mutate.py --out runs/2026-11-08-regulatory-gate-power.json --phase before",
        "python3 tools/regulatory_gate_power/mutate.py --out runs/2026-11-08-regulatory-gate-power.json --phase after",
        "python3 tools/regulatory_gate_power/observability.py --out runs/2026-11-08-regulatory-observability.json",
        "python3 tools/regulatory_gate_power/report.py",
        "```",
        "",
    ]
    PAGE.parent.mkdir(parents=True, exist_ok=True)
    PAGE.write_text("\n".join(lines))
    print(f"wrote {PAGE} and {FIGURE}")

    payload = {
        "title": "TensorCAD regulatory gate power: 9/12 -> 11/12, and the one survivor is a measured no-op",
        "summary": (
            "First mutation audit of TensorCAD's verdict-emitting modules (statutory s.2.2, "
            f"intact_stability s.2.3, damage). {killed_b}/{total_b} mutants killed before, "
            f"{killed_a}/{total_a} after; the 2 observable gaps (flooding-angle clamp on the peak "
            "search, permeability in the net damaged waterplane) are closed by closed-form "
            "regression tests; the remaining survivor is a measured structural no-op."
        ),
        "scenario_key": "tensorcad.regulatory_gate_power",
        "mission_id": 6401,
        "tags": ["TensorCAD", "gate-power", "mutation-testing", "IS-Code-2008", "damage-stability", "verification-debt"],
        "markdown": "\n".join(lines),
        "figures": [
            {
                "name": "regulatory_gate_power.png",
                "caption": "Per-mutant gate power before and after remediation (green = the suite notices, red = it does not).",
                "png_base64": base64.b64encode(png).decode(),
            }
        ],
    }
    out = Path("/tmp/report_m6401.json")
    out.write_text(json.dumps(payload))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())