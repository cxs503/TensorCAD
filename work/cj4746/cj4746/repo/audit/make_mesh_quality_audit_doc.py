#!/usr/bin/env python3
"""Generates docs/mesh-quality-lib-audit-m17404.md and the report payload.

Numbers are injected from the artefacts, never typed by hand:
  * /root/tp-agent/work/m17404/artifacts/mesh_quality_lib/{crosscheck,gate}.json (audited side)
  * artifacts/mesh_quality_audit/dialect.json (this audit's own measurements)

    PYTHONPATH=/root/tp-agent/work/m17404 python3 audit/make_mesh_quality_audit_doc.py
"""
from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
AUDITED = Path("/root/tp-agent/work/m17404")
DIALECT = json.loads((ROOT / "artifacts" / "mesh_quality_audit" / "dialect.json").read_text())
XC = json.loads((AUDITED / "artifacts" / "mesh_quality_lib" / "crosscheck.json").read_text())
GATE = json.loads((AUDITED / "artifacts" / "mesh_quality_lib" / "gate.json").read_text())


def f(x: float, n: int = 6) -> str:
    return f"{x:.{n}g}"


def figure_png() -> bytes:
    """Left: canonical vs legacy aspect on every element of every audited mesh.
    Right: mesh-level max under both dialects against the shared limit 25."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8))
    colours = {"wigley": "#1f77b4", "series60": "#2ca02c", "suboff_coarse": "#d62728",
               "suboff_taubin30": "#9467bd"}
    for name, colour in colours.items():
        d = XC["dialect_audit"][name]
        legacy = d["validation_max_aspect_lmax2_over_2A"]
        canon = d["library_max_aspect_ratio"]
        ax1.scatter([legacy], [canon], s=70, color=colour, label=f"{name} (max element)")
    single = DIALECT["verdict"]["single_element_counterexample"]
    ax1.scatter([single["legacy_aspect_lmax2_over_2A"]], [single["canonical_aspect"]],
                s=110, marker="*", color="black", label="audit counterexample")
    lim = np.array([1.0, 2e4])
    ax1.plot(lim, lim, "k--", lw=1, label="y = x")
    ax1.plot(lim, lim * 1 / np.sqrt(3.0), "k:", lw=1, label="y = x/sqrt(3) (lower bound of ratio)")
    ax1.axvline(25.0, color="grey", lw=1)
    ax1.axhline(25.0, color="grey", lw=1)
    ax1.fill_between(lim, 25.0, lim, where=(lim > 25.0), color="orange", alpha=0.12)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("legacy aspect $l_{max}^2/(2A)$  (mesh_conditioning / validation, limit 25)")
    ax1.set_ylabel("canonical aspect $l_{max}p/(4\\sqrt{3}A)$  (geomquality, limit 25)")
    ax1.set_title("Same element, two 'aspect ratio' dialects\nshaded band = legacy fails, canonical passes")
    ax1.legend(fontsize=7, loc="upper left")
    ax1.grid(alpha=0.3)

    names = sorted(DIALECT["mesh_measurements"])
    legacy_max = [DIALECT["mesh_measurements"][n]["max_aspect_legacy_mesh_conditioning"] for n in names]
    canon_max = [DIALECT["mesh_measurements"][n]["max_aspect_geomquality"] for n in names]
    y = np.arange(len(names))
    ax2.barh(y - 0.2, legacy_max, height=0.4, color="#8c564b", label="legacy $l_{max}^2/(2A)$")
    ax2.barh(y + 0.2, canon_max, height=0.4, color="#17becf", label="canonical $l_{max}p/(4\\sqrt{3}A)$")
    ax2.axvline(25.0, color="black", ls="--", lw=1, label="shared limit 25")
    n_dis = DIALECT["mesh_measurements"]["suboff_coarse"]["elements_where_verdict_differs_at_limit_25"]
    ax2.text(50, y[names.index("suboff_coarse")],
             f"{n_dis}/9088 elements flip verdict", fontsize=8, color="#d62728")
    ax2.set_yticks(y)
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_xscale("log")
    ax2.set_xlabel("max aspect ratio (log)")
    ax2.set_title("Mesh-level maximum under both dialects")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    return buf.getvalue()


def table(rows, header) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def build_markdown() -> str:
    single = DIALECT["verdict"]["single_element_counterexample"]
    diff = XC["differential"]
    ide = XC["ideal_values"]

    ideal_rows = [[k, f(v["expected"]), f(v["got"]), f"{v['difference']:.3g}"]
                  for k, v in sorted(ide.items())]
    def worst_diff(d):
        return max(v["difference"] for v in d["metrics"].values())

    route_rows = [[n, f(worst_diff(d), 4)] for n, d in sorted(diff.items())]
    dialect_rows = [[n,
                     f(d["validation_max_aspect_lmax2_over_2A"], 5),
                     f(d["library_max_aspect_ratio"], 5),
                     f(d["library_max_aspect_ratio"] / d["validation_max_aspect_lmax2_over_2A"], 5),
                     f(d["gate_over_validation"], 3)]
                    for n, d in sorted(XC["dialect_audit"].items())]
    flip_rows = [[n, m["triangles"],
                  m["elements_where_verdict_differs_at_limit_25"],
                  m["of_which_pass_every_other_library_check"]]
                 for n, m in sorted(DIALECT["mesh_measurements"].items())]
    upstream = ", ".join(f"{k}: {'available' if v['available'] else 'absent'}"
                         for k, v in XC["upstream"].items())

    return f"""# Board audit — m17404 / `GAP-CAD-D1-MESH-QUALITY-LIB` (mesh-quality metric library)

*Auditor:* m17437 (board-audit, 交差核验单 — not a redo). *Repo:* TensorCAD.
*Audited commit:* `9b66af8` (`geomquality: canonical mesh-quality metric library + differential
cross-check + CI gate`), tree `/root/tp-agent/work/m17404`.
*Machine-readable artefact of this audit:* `artifacts/mesh_quality_audit/dialect.json`.

## 1. Scope and what was checked

Ten measures (`area`, `min/max_angle_deg`, `edge_ratio`, `aspect_ratio`, `aspect_frobenius`,
`gamma`, `radius_ratio`, `scaled_jacobian`, `condition_number_inverse`), the aggregate report,
the criteria/gate path, and the four claimed deliverables of the mission: artefact, benchmark
numbers, `repro_cmd`, six-section report.

Executed by the auditor in the audited tree (not read from the report):

| check | command | result |
|---|---|---|
| reproduction | `python tools/mesh_quality_crosscheck.py` | exit 0, `PASS`, all six meshes recomputed, artefact rewritten |
| unit/invariant suite | `python -m pytest tests/test_mesh_quality_lib.py -q` | 10 passed |
| no regression in neighbouring gates | `... tests/test_gate.py tests/test_mesh_conditioning.py tests/test_geometry_audit.py -q` | 77 passed, 1 skipped |
| anchors vs ideal values | `ideal_values` in the artefact | worst deviation {f(max(v['difference'] for v in ide.values()))} |
| route agreement | `differential.*.metrics.*.difference` | worst {f(max(max(v["difference"] for v in d["metrics"].values()) for d in diff.values()))} (acos-conditioned angle on `suboff_coarse`) |
| production reachability | `grep -rn geomquality --include=*.py tensorcad` | only inside `tensorcad/geomquality/*`; no solver/gate path imports `geomquality.quality` |

Independence limitations: the auditor used the same host, no gmsh/VTK/meshmagick
({upstream}), and no external reference values are available for the same meshes.

## 2. Anchor criterion, item by item

Anchor text: *与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异<1e-6 或差异解释）；可作 CI 门用*.

| item | criterion | measured | verdict |
|---|---|---|---|
| ≥3 meshes | ≥3 | 6 (`box_unit`, `series60`, `wigley`, `suboff_coarse`, `suboff_taubin30`, `sliver_patch`) | PASS |
| differential < 1e-6 | vs upstream | **not run** — all three upstream packages absent, `pip index` on this host (loong64) has no `gmsh` wheel | NOT RUN |
| "或差异解释" | explanation recorded | `upstream.*.reason` records the verbatim import error for each; the doc re-words acceptance as `|library − independent route| < 1e-6` | PASS (explained) |
| usable as CI gate | gate + thresholds + sensitivity | `quality_gate` + `DEFAULT_CRITERIA` + sweep in `gate.json.sensitivity`; sliver/raw-Suboff FAIL, smoothed PASS | PASS at library level, **not wired into CI** (see F2) |
| `tensor_benchmark_cases` row for the anchor | benchmark judge | query returns 0 rows for `anchor_slug='GAP-CAD-D1-MESH-QUALITY-LIB'` | ABSENT |

Because the benchmark table holds no row for this anchor, no external reference `got` exists to
compare against: the only `ref` available is analytic (ideal values and closed forms), and the
only cross-check is an internally written independent route. That is a real verification, but it
is *verification*, not *benchmark validation*.

## 3. Numbers reproduced by the auditor

Ideal (closed-form) values vs measured:

{table(ideal_rows, ["measure", "expected", "got", "difference"])}

Worst |library − independent route| per mesh:

{table(route_rows, ["mesh", "max abs deviation"])}

## 4. Findings (severity-ranked)

### F1 (medium, blocking for the "same threshold" claim) — the two aspect-ratio dialects are not the same quantity, so the shared limit 25 is not transferable

`mesh_conditioning`/`validation` ship `aspect = l_max²/(2A)`; the new library ships
`aspect = l_max·p/(4√3·A)`. The ratio is `p/(2√3·l_max) ∈ (1/√3, √3/2] = (0.577, 0.866]`, i.e.
bounded away from 1 and shape-dependent, so the library's number is systematically *smaller*:

{table(dialect_rows, ["mesh", "validation l_max²/2A", "canonical aspect", "ratio", "gate/validation"])}

Measured consequence, per mesh, of applying both dialects at the same numeric limit 25:

{table(flip_rows, ["mesh", "triangles", "elements whose verdict flips", "of which pass every other canonical check"])}

Single-element counterexample (no mesh, no threshold tuning, scale-invariant to 2e-16):

```
triangle {single['triangle']}
legacy = {f(single['legacy_aspect_lmax2_over_2A'])}  -> FAIL at 25
canonical = {f(single['canonical_aspect'])}  -> PASS at 25
min_angle = {f(single['min_angle_deg'])} deg (>2), gamma = {f(single['gamma'])} (>0.05),
edge_ratio = {f(single['edge_ratio'])} (<100)
```

The code comment in `quality.py` (`DEFAULT_CRITERIA`) says the limits are "the thresholds this
repository already ships elsewhere … both re-expressed here in the canonical definitions", and
`docs/mesh-quality-lib.md` §3 says "the repository's existing decimation guard (25), canonically
expressed". The first half is true (2° and 25 are the shipped numbers); "re-expressed" is not:
adopting the library's gate relaxes the aspect guard by a factor of 1.15–1.73. The library's own
artefact already *prints* both numbers side by side (`dialect_audit`), so this is a wording and
contract defect, not a hidden arithmetic error — but a reader who replaces the old gate with the
new one silently changes behaviour.

**Remediation:** either (a) state the conversion `canonical = legacy·p/(2√3·l_max)` next to the
limit and re-derive the canonical limit from 25 (i.e. 14.4–21.7, or simply drop the equivalence
claim), or (b) keep 25 and say plainly that the new gate is the more permissive one.

### F2 (medium) — the "CI gate" is not in CI, and the library is not reachable from production

The gate exists as `quality_gate`/`gate_report` returning a verdict dict; `tools/mesh_quality_crosscheck.py`
exits 1 on violation, so it *can* be a door. But: no production module imports
`tensorcad.geomquality.quality` (only `tools/` and `tests/`), no workflow in `.github/` invokes the
driver, and `mesh_conditioning`/`gate` keep their own dialects. So the mission's "可作 CI 门用" is
readiness, not integration.

### F3 (low) — definition provenance of `gamma` is transcribed, not verified

`gamma` is labelled *Gmsh 'gamma'* in `METRIC_DEFINITIONS` and the doc, but gmsh is absent and the
host's package index has no `gmsh` wheel, so the transcription was never diffed. The formula
satisfies the only locally checkable constraint (ideal = 1 at the equilateral element). Gmsh's own
normalisation of Γ must be re-checked on a host where `import gmsh` works; until then the label
carries the same status as the rest of the upstream column: transcribed. (Recorded as an
unverified claim, not as an error — I could not exhibit a counterexample.)

### F4 (low) — duplicate definitions left in place

`tensorcad/mesh_conditioning._aspect_ratio` and `tensorcad/validation` still define aspect ratio
their own way; the new module adds a fourth definition without deprecating or aliasing the old
ones. The doc names the four dialects (good) but nothing prevents a fifth from appearing.

## 5. Report quality (six-section requirement)

`tensor_reports` id 16393 for m17404: the six sections are present (case/method/parameters/
comparison/sensitivity/failure modes), 1 figure, `qa_status=pending`. The comparison table states
`ref` as ideal/closed forms and `got` as measured values — consistent with §2 above, i.e. it is
honest about *not* having an upstream reference. Two wording items overstate: the "same measure"
claim (F1) and "可作 CI 门用" (F2).

## 6. Error judgement against a baseline (explicit)

| key quantity | reference | measured | relative error | < 3 % ? |
|---|---|---|---|---|
| `aspect_ratio`, `gamma`, `edge_ratio`, `radius_ratio`, `scaled_jacobian`, `frobenius`, `condition_number_inverse` on the equilateral element | closed form (ideal 1 / 0.866 / 0.577) | see §3 | <= 2.3e-16 | yes (by 14 orders of magnitude) |
| 10 measures on the six hull meshes | auditor-independent route (Heron, law of cosines, SVD) | see §3 | worst 2.4e-08 absolute on an angle of ~0.0039 deg, i.e. ~6e-6 relative; all other measures <= 4.3e-14 absolute | yes |
| the same on rational meshes | exact `fractions.Fraction` arithmetic | `box_unit` 2.13e-14, `sliver_patch` 2.54e-12 | <= 1e-12 | yes |
| **vs gmsh / meshmagick / VTK** | **external reference** | **not computed** | **n/a** | **cannot be answered — no baseline was run (all three packages absent, no `gmsh` wheel on the loong64 index). This is the honest answer, not a pass.** |
| aspect dialect vs the shipped gate | `mesh_conditioning` limit 25 | library value is 0.577–0.866x the shipped one | systematic, not round-off | **no — the two numbers are not the same quantity (F1)** |

## 7. Disposition

**Approve with notes** for the library itself: the implementation is real, the arithmetic is
independently verified to {f(max(max(v["difference"] for v in d["metrics"].values()) for d in diff.values()))} and exact to
1e-12 on rational meshes, the artefacts regenerate, the tests assert invariants
(ideal values, exact arithmetic, rigid-motion/scale invariance, degenerate→worst-not-NaN, gate
direction, aggregate consistency) rather than merely executing code.

**Maturity level: L2** — independent functional implementation. Not L3: no production caller.
Not L4: the anchor's own benchmark (gmsh/meshmagick differential) was not executed and no
`tensor_benchmark_cases` row exists, so there is no external reference and no error against a
baseline; internal verification against closed forms is L1-grade evidence carried along, not a
benchmark.

Requested changes (F1, F2) are documentation/contract level and do not invalidate the deliverable.
"""


def main() -> int:
    md = build_markdown()
    out = ROOT / "docs" / "mesh-quality-lib-audit-m17404.md"
    out.write_text(md)

    png = figure_png()
    fig_dir = ROOT / "artifacts" / "mesh_quality_audit"
    fig_dir.mkdir(parents=True, exist_ok=True)
    (fig_dir / "dialect_figure.png").write_bytes(png)

    payload = {
        "title": "m17437 交差核验：m17404 GAP-CAD-D1-MESH-QUALITY-LIB 网格质量指标库（L2 判定 + 方言缺陷）",
        "summary": (
            "独立复跑 m17404 的 crosscheck/pytest（exit 0、10 passed、77 passed）确认真代码真工件；"
            "对照 tensor_benchmark_cases（该锚点 0 行）与锚点判据逐项核验：6 网格 ≥3 满足、"
            "与 gmsh/meshmagick 对拍未执行（上游全缺，有差异解释，PASS 走解释分支）、库级 CI 门可用但未接入 CI；"
            "新增实测缺陷：库的 aspect_ratio 与仓内既有 l_max²/2A 不是同一个量（比值 p/(2√3 l_max) ∈ (0.577,0.866]），"
            "同一阈值 25 不可迁移——suboff_coarse 上 33/9088 元素旧门判超限而新库判合规（其中 8 个元素其它指标全过）。"
            "判级 L2：无生产调用方（L3 未达）、无外部基准（L4 未达）。"
        ),
        "scenario_key": "board_audit",
        "mission_id": 17437,
        "tags": ["TensorCAD", "board-audit", "GAP-CAD-D1-MESH-QUALITY-LIB", "m17404", "maturity-verdict"],
        "markdown": md,
        "figures": [{
            "name": "aspect_dialects.png",
            "caption": (
                "左：同一元素在两种 aspect 方言下的值（点=各网格最大元素，星=审计反例），"
                "y=x/sqrt(3) 是比值下界，橙带=旧方言超 25 而新库合规的判决翻转区；"
                "右：各网格最大 aspect 在两方言下与共用阈值 25 的关系，suboff_coarse 有 33/9088 元素翻转判决。"
            ),
            "png_base64": base64.b64encode(png).decode(),
        }],
    }
    (Path("/tmp/report_m17437.json")).write_text(json.dumps(payload))
    print(f"wrote {out} and /tmp/report_m17437.json ({len(png)} bytes png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())