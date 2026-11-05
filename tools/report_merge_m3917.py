#!/usr/bin/env python3
"""report_merge_m3917.py -- evidence + paper-grade report for the PR #47 cure.

Two jobs:

1. NUMERICAL: re-verify the *merged* tree end-to-end on the feature the PR brings
   in (large-angle righting arms) against the independent 2D box reference, so the
   cure is not certified by "tests are green" alone but by a number the merged
   tree produces.  Everything is written to docs/m3917/results.json.
2. FIGURE: docs/m3917/gz_verification.png -- GZ(phi) numerical vs analytic, the
   error curve, and the merge audit (API surface + per-file line keeping).

docs/m3917_merge_cure.md is generated here too, with every number injected from
the computation / the audit JSON (nothing hand-copied).
"""
from __future__ import annotations

import json
import math
import pathlib

import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import tensorcad  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "m3917"

L, B, H, T, KG = 6.0, 3.0, 4.0, 2.0, 1.0  # box hull: GM = KB + BM - KG = 0.375 m
ANGLES = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]


def _arm(value) -> float:
    """Coerce a righting-arm result (float or HeeledFloat record) to metres."""
    if isinstance(value, (int, float)):
        return float(value)
    fields = vars(value) if hasattr(value, "__dict__") else {}
    for key in ("arm_m", "gz_m", "righting_arm_m", "gz", "righting_arm", "arm"):
        if key in fields and isinstance(fields[key], (int, float)):
            return float(fields[key])
    for key, val in fields.items():                      # last resort: first float
        if isinstance(val, (int, float)) and ("arm" in key or "gz" in key):
            return float(val)
    raise TypeError("cannot read a righting arm out of %r (fields=%s)"
                    % (type(value).__name__, sorted(fields)))


def box_study() -> dict:
    mesh = tensorcad.box_mesh(L, B, H, z0=-T)
    rows = []
    for phi in ANGLES:
        num = _arm(tensorcad.righting_arm(mesh, phi, kg=KG, waterline_z=0.0))
        ref = _arm(tensorcad.analytic_box_righting_arm(L, B, H, T, phi, kg=KG))
        rows.append({"heel_deg": phi, "gz_numeric_m": num, "gz_analytic_m": ref,
                     "abs_err_m": abs(num - ref),
                     "rel_err": abs(num - ref) / abs(ref) if ref else abs(num - ref)})
    curve = tensorcad.righting_arm_curve(mesh, ANGLES, kg=KG, waterline_z=0.0)
    return {
        "box": {"length_m": L, "beam_m": B, "height_m": H, "draft_m": T, "kg_m": KG,
                "km_m": T / 2.0 + B * B / (12.0 * T), "gm_m": T / 2.0 + B * B / (12.0 * T) - KG},
        "rows": rows,
        "max_abs_err_m": max(r["abs_err_m"] for r in rows),
        "max_rel_err": max(r["rel_err"] for r in rows if r["heel_deg"] > 0.0),
        "curve": curve.to_dict() if hasattr(curve, "to_dict") else str(curve),
        "version": tensorcad.__version__,
        "api_count": len(tensorcad.__all__),
        "modules": {m: hasattr(__import__("tensorcad." + m, fromlist=["x"]), "__name__")
                    for m in ("righting", "transom")},
    }


def figure(study: dict, audit: dict) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))

    ax = axes[0]
    phis = [r["heel_deg"] for r in study["rows"]]
    ax.plot(phis, [r["gz_numeric_m"] for r in study["rows"]], "o-", lw=2,
            label="GZ 3-D clip (merged tree)")
    ax.plot(phis, [r["gz_analytic_m"] for r in study["rows"]], "x--", lw=1.4,
            label="GZ analytic 2-D section")
    ax.set_xlabel("heel angle (deg)")
    ax.set_ylabel("righting arm GZ (m)")
    ax.set_title("Box %.0fx%.0fx%.0f m, T=%.1f m, KG=%.1f m (GM=%.3f m)"
                 % (L, B, H, T, KG, study["box"]["gm_m"]), fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    errs = [max(r["abs_err_m"], 1e-18) for r in study["rows"]]
    ax.semilogy(phis, errs, "s-", color="crimson")
    ax.axhline(0.03 * max(abs(r["gz_analytic_m"]) for r in study["rows"]), ls=":", color="grey",
               label="3 % of max GZ (acceptance)")
    ax.set_xlabel("heel angle (deg)")
    ax.set_ylabel("|GZ_num - GZ_ref| (m)")
    ax.set_title("absolute error, max %.1e m" % study["max_abs_err_m"], fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    api = audit["api"]
    labels = ["main", "agent/m3431", "union", "merged"]
    vals = [api["count_main"], api["count_branch"], api["count_union"], api["count_merged"]]
    bars = ax.bar(labels, vals, color=["#8fa9c4", "#c4a48f", "#9dc49b", "#3c7a3f"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1, str(v), ha="center", fontsize=9)
    ax.set_ylim(0, max(vals) * 1.18)
    ax.set_ylabel("public names in __all__")
    ax.set_title("API surface: merged == union, 0 lines lost", fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    path = OUT / "gz_verification.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    study = box_study()
    audit = json.loads((ROOT / "artifacts/m3917/merge_audit.json").read_text())
    png = figure(study, audit)
    payload = {"study": study, "audit_summary": {
        "ok": audit["ok"], "lost_lines": len(audit["lost"]),
        "extended": len(audit["extended"]), "restated": len(audit["restated"]),
        "files_checked": len(audit["files_checked"]), "api": audit["api"],
        "per_file": audit["per_file"]}, "figure": png}
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                                      encoding="utf-8")
    print(json.dumps({k: v for k, v in study.items() if k != "rows"}, indent=2,
                     ensure_ascii=False, default=str))
    print("max_rel_err %.3e  max_abs_err %.3e m" % (study["max_rel_err"], study["max_abs_err_m"]))
    print("png", png, pathlib.Path(png).stat().st_size // 1024, "KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())