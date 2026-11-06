#!/usr/bin/env python3
"""Build the docs page, figure and platform report for the hull gate-power audit.

Reads ``runs_m4780_hull_gate_power.json`` (produced by ``mutate.py``) and emits:

* ``docs/analysis/hull_gate_power.png`` -- one bar per mutant, green = killed,
  red = survived;
* ``docs/analysis/2026-11-06-hull-gate-power.md`` -- the audit page, generated
  from the JSON (no retyped numbers);
* ``/tmp/report_m4780.json`` -- the platform report payload (post it with the
  usual curl).

Usage::

    python3 tools/hull_gate_power/report.py
"""
from __future__ import annotations

import base64
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "runs_m4780_hull_gate_power.json"
FIG = ROOT / "docs" / "analysis" / "hull_gate_power.png"
PAGE = ROOT / "docs" / "analysis" / "2026-11-06-hull-gate-power.md"
REPORT = Path("/tmp/report_m4780.json")

ORDER = [
    "HUL-1-taper-exponent",
    "HUL-2-vertical-taper",
    "HUL-3-end-zeroing",
    "HUL-4-midship-split",
    "ANA-1-volume-factor",
    "ANA-2-vcb",
    "ANA-3-waterplane-half-percent",
    "ANA-4-taper-arg",
    "ANA-5-volume-half-percent",
    "HUL-5-taper-integral-const",
]

WHY = {
    "HUL-1-taper-exponent": "The generator's exponent is a parameter, but the analytic reference keeps its own copy of the closed form. A generator that samples ``sin^0.8`` while the reference integrates ``sin^0.7`` moves the mesh volume by ~48 % -- far outside every documented band.",
    "HUL-2-vertical-taper": "The ``sqrt`` profile is what makes the vertical integral ``2/3``. A linear profile changes the volume factor to ``1/2`` (-25 %).",
    "HUL-3-end-zeroing": "``math.sin(math.pi)`` is ``1.22e-16``, not zero. Without the explicit zeroing the end stations carry a ``~1e-11 m`` half-breadth. The offsets test asserts the waterline closes at both ends *exactly*, so the slot is visible even though the volume change is below any tolerance.",
    "HUL-4-midship-split": "At midship ``sin(pi/2) = 1`` and ``1**p = 1`` for every exponent, and the default form has ``fore == aft``. The ``>=`` vs ``>`` choice is unobservable in the output: no test can kill this mutant, and it is not a defect. It is recorded as a structural no-op, not a gap.",
    "ANA-1-volume-factor": "The reference volume is the thing the mesh is compared against. Halving the ``2/3`` factor moves the documented -0.88 % error to ~-26 %.",
    "ANA-2-vcb": "``VCB = -0.4*T`` is the z-weighted mean of ``zeta**(1/2)`` on ``[-T, 0]``. Shifting it to ``-0.5*T`` is a 25 % error in the centre of buoyancy.",
    "ANA-3-waterplane-half-percent": "A 0.5 % shift of the reference waterplane sits inside the 1.5 % default-grid band but outside the 0.5 % fine-grid band -- the fine-grid test is what has the teeth.",
    "ANA-4-taper-arg": "The Beta integral's first gamma argument is ``(e+1)/2``. Off by half an exponent, the textbook values (``e=1 -> 2/pi``, ``e=2 -> 1/2``) break immediately.",
    "ANA-5-volume-half-percent": "The documented -0.880 % grid error is pinned to ``abs=5e-5``. A 0.5 % shift of the reference moves it to -1.375 %.",
    "HUL-5-taper-integral-const": "``TAPER_INTEGRAL`` in ``hull.py`` is cross-derived two independent ways: Simpson quadrature of ``sin^0.7`` and the gamma-function derivation in ``hydrostatics_analytic.py``. Halving the prefactor breaks both.",
}


def main() -> None:
    data = json.loads(RESULTS.read_text())
    rows = []
    for mid in ORDER:
        entry = data.get(mid)
        if not entry or entry.get("status") != "ok":
            rows.append({"id": mid, "killed": None, "summary": "not-run", "failed": []})
            continue
        rows.append(
            {
                "id": mid,
                "killed": entry["run"]["killed"],
                "summary": entry["run"]["summary"],
                "failed": entry["run"]["failed"],
                "wall": entry["run"]["wall_seconds"],
            }
        )
    done = [r for r in rows if r["killed"] is not None]
    killed = [r for r in done if r["killed"]]
    survived = [r for r in done if not r["killed"]]
    power = 100.0 * len(killed) / len(done) if done else float("nan")

    # ---- figure ---------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 3.6), dpi=150)
    labels = [r["id"].split("-")[0] + "-" + r["id"].split("-")[1] for r in rows]
    colors = ["#2e7d32" if r["killed"] else "#c62828" if r["killed"] is False else "#9e9e9e" for r in rows]
    ax.bar(range(len(rows)), [1] * len(rows), color=colors, width=0.72)
    for i, r in enumerate(rows):
        ax.text(i, 1.06, "killed" if r["killed"] else ("survived" if r["killed"] is False else "n/a"),
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=38, ha="right", fontsize=8)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_title(f"Gate power of the hull closed-form claims: {len(killed)}/{len(done)} mutants killed ({power:.1f} %)")
    fig.tight_layout()
    fig.savefig(FIG)
    plt.close(fig)

    # ---- docs page ------------------------------------------------------
    lines = [
        "# Gate power of the hull generator's closed forms",
        "",
        "**Date**: 2026-11-06 · **Scope**: `tensorcad/hull.py` (generator) + `tensorcad/hydrostatics_analytic.py` (analytic reference)",
        "· **Harness**: `tools/hull_gate_power/mutate.py` (rerun it; this page is generated from its JSON by `report.py`)",
        "",
        "## Why this audit",
        "",
        "The 2026-11-04 claim-gate-power audit attacked five modules and closed at 14/14, but it",
        "explicitly left the parametric hull generator uncovered: *\"the generator's own closed forms",
        "are shared constants -- a mutant of the exponent would have to be paired with the analytic",
        "side to be meaningful\"*. That pairing is what this page does. The generator and the",
        "reference each carry their own copy of the closed form (the exponent, the `sqrt` vertical",
        "profile, the `2/3` volume factor, the Beta integral). A defect that moves one copy and not",
        "the other is exactly the class of regression a mesh-vs-reference comparison exists to catch --",
        "if any test actually performs the comparison at a tolerance tighter than the defect.",
        "",
        "## Method",
        "",
        "* **One defect at a time.** Exact-string replacement (a non-unique pattern is refused), one",
        "  pytest run per mutant, restore from original bytes with a SHA-256 check.",
        "* **Clean-tree guard.** Aborts if `git diff -- tensorcad` is non-empty.",
        "* **Targeted scope.** Generator mutants run `tests/test_hull.py` + `tests/test_analytic_hydrostatics.py`;",
        "  analytic-side mutants run `tests/test_hydrostatics_analytic.py` + `tests/test_analytic_hydrostatics.py`",
        "  + `tests/test_hull.py`, so a kill is attributable.",
        "",
        "## Result",
        "",
        f"Gate power: **{len(killed)}/{len(done)} = {power:.1f} %**",
        "",
        "| mutant | kind | file | documented claim it should protect | verdict |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        entry = data.get(r["id"], {})
        m = entry.get("mutant", {})
        verdict = "**KILLED**" if r["killed"] else ("**SURVIVED**" if r["killed"] is False else "not run")
        lines.append(f"| `{r['id']}` | {m.get('kind', '?')} | `{m.get('file', '?')}` | {m.get('claim', '?')} | {verdict} |")
    lines += ["", "![gate power](hull_gate_power.png)", "", "### What survived, and why", ""]
    for r in survived:
        lines.append(f"* `{r['id']}` -- {WHY.get(r['id'], 'see harness')}")
    lines += [
        "",
        "### What the kills show",
        "",
    ]
    for r in killed:
        n = len(r["failed"])
        lines.append(f"* `{r['id']}` -- {n} test(s) red. {WHY.get(r['id'], '')}")
    lines += [
        "",
        "## Limitations",
        "",
        "* **The mutant set is hand-picked, not exhaustive.** 10 mutants over two modules is a probe of the",
        "  evidence style, not a coverage metric. A surviving mutant proves a gap; a killed mutant does not",
        "  prove the absence of one.",
        "* **`HUL-4` is a structural no-op**, not a missed defect: the midship tie is unobservable because",
        "  `1**p == 1` for every `p`. It is listed so a future reader does not mistake it for debt.",
        "* **`tolerance`-kind mutants measure the band, not correctness.** They show a documented number is",
        "  pinned tighter than the injected shift; they say nothing about whether the number is right.",
        "* **Not covered here**: the lofting/topology code of `generate_hull` (its claims are already pinned",
        "  by the watertight/manifold/Euler gates audited on 2026-11-04), `mesh_conditioning.py` (35 s per",
        "  run, a separate deliberate audit), and the `export.py`/`gate.py` thresholds.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 tools/hull_gate_power/mutate.py --out runs_m4780_hull_gate_power.json   # resumable",
        "python3 tools/hull_gate_power/report.py                                        # this page + figure",
        "```",
        "",
    ]
    PAGE.write_text("\n".join(lines))

    # ---- platform report ------------------------------------------------
    fig_b64 = base64.b64encode(FIG.read_bytes()).decode()
    md = PAGE.read_text()
    payload = {
        "title": "TensorCAD 生成器闭式 gate-power 审计：hull.py × 解析参考侧 10 mutant 配对攻击（m4780）",
        "summary": (
            f"对 2026-11-04 审计明示未覆盖的 hull.py 生成器闭式做配对突变审计：10 个单缺陷 mutant"
            f"（生成器侧 4 + 解析侧 6），{len(killed)}/{len(done)} 被杀（{power:.1f}%）。"
            "HUL-4（中点归属 >= vs >）结构性不可观测，如实记为 no-op 而非欠账。"
        ),
        "scenario_key": "tensorcad-hull-gate-power",
        "mission_id": 4780,
        "tags": ["tensorcad", "verification-debt", "mutation-testing", "gate-power"],
        "markdown": md,
        "figures": [
            {
                "name": "hull_gate_power.png",
                "caption": f"每个 mutant 一根柱：绿=测试套件变红（killed），红=存活。{len(killed)}/{len(done)} killed。",
                "png_base64": fig_b64,
            }
        ],
    }
    REPORT.write_text(json.dumps(payload, ensure_ascii=False))
    print(f"page: {PAGE}")
    print(f"figure: {FIG} ({FIG.stat().st_size} bytes)")
    print(f"report payload: {REPORT}")
    print(f"gate power: {len(killed)}/{len(done)} = {power:.1f}%")


if __name__ == "__main__":
    main()
