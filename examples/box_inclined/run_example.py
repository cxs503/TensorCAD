"""Exact barge hydrostatics at combined heel and trim -- measured evidence.

Run::

    python examples/box_inclined/run_example.py

Writes ``examples/box_inclined/outputs/results.json`` (every number below),
``report.md`` and two figures.  No number in the report is typed by hand: the tables and
the captions are built from the computed records.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.box_inclined import (  # noqa: E402
    box_inclined_curve,
    box_inclined_gz,
    box_inclined_hydrostatics,
)

L, B, D, T = 10.0, 2.0, 1.2, 0.8
VOL = L * B * T
KG = 0.55
OUT = Path(__file__).resolve().parent / "outputs"


def wall_sided_gz(angle_deg: float) -> float:
    bm = B * B / (12.0 * T)
    gm = 0.5 * T + bm - KG
    phi = math.radians(angle_deg)
    return math.sin(phi) * (gm + bm / 2.0 * math.tan(phi) ** 2)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    angles = tuple(range(0, 91, 5))
    curve = box_inclined_curve(L, B, D, kg=KG, volume_m3=VOL, angles_deg=angles)
    reference = tuple(wall_sided_gz(a) for a in angles)
    # the wall-sided window ends where the waterline reaches a bilge; the exact route has
    # no window, so the comparison is reported angle by angle instead of truncated
    error = tuple(
        abs(c - r) / abs(r) if r else 0.0 for c, r in zip(curve, reference, strict=True)
    )

    upright = box_inclined_hydrostatics(L, B, D, volume_m3=VOL)
    upright_checks = {
        "volume_rel_error": abs(upright.volume_m3 - VOL) / VOL,
        "vcb_rel_error": abs(upright.vcb - T / 2) / (T / 2),
        "i_transverse_rel_error": abs(upright.i_transverse_m4 - L * B**3 / 12) / (L * B**3 / 12),
        "i_longitudinal_rel_error": abs(upright.i_longitudinal_m4 - B * L**3 / 12) / (B * L**3 / 12),
        "bm_transverse_rel_error": abs(upright.bm_transverse - B**2 / (12 * T)) / (B**2 / (12 * T)),
    }

    trim_rows = []
    for trim in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        hydro = box_inclined_hydrostatics(L, B, D, volume_m3=VOL, trim_deg=trim)
        closed = L**2 * math.tan(math.radians(trim)) / (12 * T)
        trim_rows.append(
            {
                "trim_deg": trim,
                "lcb": hydro.lcb,
                "lcb_closed_form": closed,
                "rel_error": abs(hydro.lcb - closed) / closed,
                # the closed form assumes the cut stays on the side walls: when the bow
                # waterline passes the deck it is no longer the same problem
                "bow_waterline_z": T + (L / 2.0) * math.tan(math.radians(trim)),
                "deck_immersed": T + (L / 2.0) * math.tan(math.radians(trim)) > D,
            }
        )

    combined = []
    for heel in (0.0, 10.0, 20.0, 30.0):
        for trim in (0.0, 2.0, 4.0):
            hydro = box_inclined_hydrostatics(
                L, B, D, volume_m3=VOL, heel_deg=heel, trim_deg=trim
            )
            combined.append(
                {
                    "heel_deg": heel,
                    "trim_deg": trim,
                    "volume_m3": hydro.volume_m3,
                    "volume_rel_error": abs(hydro.volume_m3 - VOL) / VOL,
                    "lcb": hydro.lcb,
                    "tcb": hydro.tcb,
                    "vcb": hydro.vcb,
                    "waterplane_area_m2": hydro.waterplane_area_m2,
                    "gz": hydro.gz(KG),
                }
            )

    peak_index = max(range(len(curve)), key=lambda i: curve[i])
    results = {
        "box": {"length": L, "beam": B, "depth": D, "draft": T, "volume_m3": VOL, "kg": KG},
        "upright_checks": upright_checks,
        "gz_curve": [
            {"heel_deg": a, "gz": c, "wall_sided": r, "rel_error": e}
            for a, c, r, e in zip(angles, curve, reference, error, strict=True)
        ],
        "gz_peak": {"heel_deg": angles[peak_index], "gz": curve[peak_index]},
        "trim_rows": trim_rows,
        "combined_attitudes": combined,
    }
    (OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    # ---- figures -----------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(angles, curve, "o-", label="exact clipped polyhedron")
    axes[0].plot(angles, reference, "--", label="wall-sided closed form")
    axes[0].set_xlabel("heel (deg)")
    axes[0].set_ylabel("GZ (m)")
    axes[0].set_title(f"barge {L}x{B}x{D} m, draft {T} m, KG {KG} m")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].semilogy(angles[1:], [max(e, 1e-18) for e in error[1:]], "s-")
    axes[1].set_xlabel("heel (deg)")
    axes[1].set_ylabel("|GZ - wall-sided| / GZ")
    axes[1].set_title("where the closed form stops applying")
    axes[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_gz_curves.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    trims = [row["trim_deg"] for row in trim_rows]
    axes[0].plot(trims, [row["lcb"] for row in trim_rows], "o-", label="exact clipped polyhedron")
    axes[0].plot(
        trims, [row["lcb_closed_form"] for row in trim_rows], "s--", label="L^2 tan(tau) / (12 T)"
    )
    axes[0].set_xlabel("trim (deg)")
    axes[0].set_ylabel("LCB (m)")
    axes[0].set_title("trim moment: exact vs closed form")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].plot(
        [row["heel_deg"] for row in combined if row["trim_deg"] == 4.0],
        [row["tcb"] for row in combined if row["trim_deg"] == 4.0],
        "^-",
    )
    axes[1].set_xlabel("heel (deg)")
    axes[1].set_ylabel("TCB (m)")
    axes[1].set_title("combined attitude: buoyancy centre drifts to starboard")
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_trim_and_combined.png", dpi=140)
    plt.close(fig)

    worst_upright = max(upright_checks.values())
    worst_trim_in_window = max(row["rel_error"] for row in trim_rows if not row["deck_immersed"])
    worst_volume = max(row["volume_rel_error"] for row in combined)
    lines = [
        "# Exact barge hydrostatics at combined heel and trim",
        "",
        f"Box {L} x {B} x {D} m, displacement {VOL:.6f} m^3 (draft {T} m), KG {KG} m.",
        "",
        "| check | worst relative error | < 3 %? |",
        "|---|---|---|",
        f"| upright VCB / moments / BM vs closed forms | {worst_upright:.3e} | {'yes' if worst_upright < 0.03 else 'NO'} |",
        f"| trim LCB vs L^2 tan(tau)/(12T), inside the side-wall window | {worst_trim_in_window:.3e} | {'yes' if worst_trim_in_window < 0.03 else 'NO'} |",
        f"| displacement preserved over the 12 combined attitudes | {worst_volume:.3e} | {'yes' if worst_volume < 0.03 else 'NO'} |",
        "",
        f"GZ peak: {curve[peak_index]:.6f} m at {angles[peak_index]} deg.",
        "",
        "| heel (deg) | GZ exact (m) | GZ wall-sided (m) | rel. diff |",
        "|---|---|---|---|",
    ]
    for a, c, r, e in zip(angles, curve, reference, error, strict=True):
        lines.append(f"| {a} | {c:.6f} | {r:.6f} | {e:.3e} |")
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"upright": upright_checks, "gz_peak": results["gz_peak"],
                      "worst_trim_in_window": worst_trim_in_window,
                      "worst_volume": worst_volume,
                      "gz_at": {a: round(c, 9) for a, c in list(zip(angles, curve))[::4]}},
                     indent=2))
    print(f"out -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())