"""Asymmetric transom family: convergence, the skew sweep, and figures.

Run from a checkout::

    python3 examples/transom_asymmetric/run_example.py

Writes ``artifacts/transom_asymmetric/`` (results.json, manifest.json, figures) and
``docs/transom_asymmetric/tables.md``. Every number in ``docs/transom_asymmetric.md``
comes from one of those files; nothing is typed by hand. The script exits non-zero
if any claim it makes stops holding.
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # run from a clone without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import hashlib
import json
import math
from pathlib import Path

from tensorcad.transom_asymmetric import (
    AsymmetricTransomParameters,
    analytic_asymmetric_hydrostatics,
    generate_asymmetric_transom_hull,
    simpson_reference,
)

GRIDS = ((41, 13), (81, 25), (161, 41), (321, 81))
SKEWS = (-0.6, -0.3, 0.0, 0.3, 0.6)
REFERENCE = dict(length=10.0, beam=2.0, draft=0.8, transom_fraction=0.02)


def convergence_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for nx, nz in GRIDS:
        model = generate_asymmetric_transom_hull(
            AsymmetricTransomParameters(**REFERENCE, skew=0.4, nx=nx, nz=nz)
        )
        errors = model.relative_errors()
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(model.mesh.faces),
                "volume_m3": model.hydrostatics.displacement_volume,
                "volume_error": errors["displacement_volume"],
                "waterplane_error": errors["waterplane_area"],
                "lcb_error": errors["longitudinal_centre_buoyancy"],
                "y_b_m": model.transverse_centre_buoyancy,
                "y_b_analytic_m": model.analytic["transverse_centre_buoyancy"],
                "y_b_error": errors["transverse_centre_buoyancy"],
                "inertia_measured_m4": model.hydrostatics.transverse_waterplane_inertia,
                "inertia_analytic_m4": model.analytic["transverse_waterplane_inertia"],
            }
        )
    return rows


def skew_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for skew in SKEWS:
        a = analytic_asymmetric_hydrostatics(AsymmetricTransomParameters(**REFERENCE, skew=skew))
        quad = simpson_reference(AsymmetricTransomParameters(**REFERENCE, skew=skew), panels=4000)
        rows.append(
            {
                "skew": skew,
                "starboard_half_breadth_m": a["starboard_half_breadth"],
                "port_half_breadth_m": a["port_half_breadth"],
                "volume_m3": a["displacement_volume"],
                "waterplane_area_m2": a["waterplane_area"],
                "lcb_m": a["longitudinal_centre_buoyancy"],
                "block_coefficient": a["block_coefficient"],
                "y_b_m": a["transverse_centre_buoyancy"],
                "y_f_m": a["transverse_centre_flotation"],
                "y_f_over_y_b": (
                    a["transverse_centre_flotation"] / a["transverse_centre_buoyancy"]
                    if a["transverse_centre_buoyancy"]
                    else float("nan")
                ),
                "inertia_centreline_m4": a["transverse_waterplane_inertia"],
                "inertia_centroid_m4": a["transverse_waterplane_inertia_centroid"],
                "volume_simpson_m3": quad["displacement_volume"],
                "y_b_simpson_m": quad["transverse_centre_buoyancy"],
                "volume_series_vs_simpson": quad["displacement_volume"]
                / a["displacement_volume"]
                - 1.0,
            }
        )
    return rows


def observed_order(rows: list[dict[str, float]], key: str) -> list[float]:
    orders: list[float] = []
    for coarse, fine in zip(rows, rows[1:], strict=False):
        ratio = (coarse["nx"] - 1) / (fine["nx"] - 1)
        e_coarse, e_fine = abs(coarse[key]), abs(fine[key])
        orders.append(math.log(e_coarse / e_fine) / math.log(ratio) if e_fine else float("nan"))
    return orders


def write_tables(path: Path, convergence: list[dict[str, float]], skew: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Generated tables (do not edit by hand)", ""]
    lines += ["## Mesh against the closed forms (skew = 0.4)", "",
              "| nx | nz | faces | V [m3] | V error | Aw error | LCB error | y_B [m] | y_B error |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in convergence:
        lines.append(
            f"| {r['nx']} | {r['nz']} | {r['faces']} | {r['volume_m3']:.6f} | "
            f"{100 * r['volume_error']:+.3f} % | {100 * r['waterplane_error']:+.3f} % | "
            f"{100 * r['lcb_error']:+.3f} % | {r['y_b_m']:.6f} | {100 * r['y_b_error']:+.3f} % |"
        )
    lines += ["", "## The skew sweep at the reference body", "",
              "| skew | B_s/2 [m] | B_p/2 [m] | V [m3] | Aw [m2] | Cb | y_B [m] | y_F [m] | y_F/y_B | I_T cl [m4] | I_T cg [m4] |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in skew:
        ratio = "&mdash;" if math.isnan(r["y_f_over_y_b"]) else f"{r['y_f_over_y_b']:.9f}"
        lines.append(
            f"| {r['skew']:+.2f} | {r['starboard_half_breadth_m']:.4f} | "
            f"{r['port_half_breadth_m']:.4f} | {r['volume_m3']:.6f} | "
            f"{r['waterplane_area_m2']:.6f} | {r['block_coefficient']:.6f} | "
            f"{r['y_b_m']:+.6f} | {r['y_f_m']:+.6f} | {ratio} | "
            f"{r['inertia_centreline_m4']:.6f} | {r['inertia_centroid_m4']:.6f} |"
        )
    path.write_text("\n".join(lines) + "\n")


def make_figures(output: Path, convergence: list[dict[str, float]], skew: list[dict[str, float]]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures: list[tuple[str, Path]] = []

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    spacing = [1.0 / (row["nx"] - 1) for row in convergence]
    axes[0].loglog(spacing, [abs(row["volume_error"]) for row in convergence], "o-", label="volume")
    axes[0].loglog(spacing, [abs(row["y_b_error"]) for row in convergence], "s-", label="y_B")
    axes[0].set_xlabel("longitudinal spacing / LOA")
    axes[0].set_ylabel("|relative error| vs closed form")
    axes[0].set_title("Asymmetric body converges (skew = 0.4)")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend()
    skews = [row["skew"] for row in skew]
    axes[1].plot(skews, [row["y_b_m"] for row in skew], "o-", label="closed form y_B")
    axes[1].plot(skews, [row["y_f_m"] for row in skew], "s--", label="closed form y_F")
    axes[1].plot(skews, [row["y_b_simpson_m"] for row in skew], "x", ms=10,
                 label="Simpson route y_B")
    axes[1].set_xlabel("breadth skew")
    axes[1].set_ylabel("transverse centre [m]")
    axes[1].set_title("Centres move linearly, and in a 4:3 ratio")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.tight_layout()
    path = output / "convergence_and_skew.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(("convergence_and_skew.png", path))

    body = AsymmetricTransomParameters(**REFERENCE, skew=0.5, nx=161, nz=41)
    model = generate_asymmetric_transom_hull(body)
    box = model.mesh.bounding_box()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for face in model.mesh.faces:
        a, b, c = (model.mesh.vertices[i] for i in face)
        if abs(a[2]) <= 1e-12 and abs(b[2]) <= 1e-12 and abs(c[2]) <= 1e-12:
            axes[0].plot([a[0], b[0], c[0], a[0]], [a[1], b[1], c[1], a[1]],
                         color="0.5", lw=0.3)
    axes[0].axhline(0.0, color="0.2", lw=0.8, ls=":")
    axes[0].set_xlabel("x [m] (transom at 0)")
    axes[0].set_ylabel("y [m]")
    axes[0].set_title(f"Waterline plan, skew = {body.skew:+.2f} (starboard up)")
    axes[0].grid(True, alpha=0.3)
    axes[0].set_aspect("equal")
    for label, station in (("transom", 0), ("midship", (body.nx - 1) // 2)):
        x = body.length * station / (body.nx - 1)
        points = [v for f in model.mesh.faces for v in (model.mesh.vertices[i] for i in f)
                  if abs(v[0] - x) <= 1e-12]
        points = sorted(set(points), key=lambda p: math.atan2(p[2], p[1]))
        axes[1].plot([p[1] for p in points], [p[2] for p in points], marker=".", ms=2,
                     lw=0.8, label=f"{label} (x = {x:.2f} m)")
    axes[1].axhline(0.0, color="tab:blue", lw=1.0, ls="--", label="design waterline")
    axes[1].set_xlabel("y [m]")
    axes[1].set_ylabel("z [m]")
    axes[1].set_title("Sections: the two sides differ, the span does not")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8)
    axes[1].set_aspect("equal")
    fig.tight_layout()
    path = output / "geometry.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(("geometry.png", path))
    return figures, box


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="artifacts/transom_asymmetric")
    args = parser.parse_args()
    output = _Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    convergence = convergence_rows()
    skew = skew_rows()
    orders = {
        key: observed_order(convergence, key)
        for key in ("volume_error", "y_b_error", "waterplane_error")
    }
    figures, box = make_figures(output, convergence, skew)
    write_tables(_REPO_ROOT / "docs" / "transom_asymmetric" / "tables.md", convergence, skew)

    # The claims the script refuses to publish a failure for.
    coarse = convergence[0]
    fine = convergence[-1]
    claims = {
        "volume_error_at_41x13_below_3pct": abs(coarse["volume_error"]) < 0.03,
        "volume_error_at_321x81_below_0.1pct": abs(fine["volume_error"]) < 1e-3,
        "y_b_error_at_321x81_below_0.1pct": abs(fine["y_b_error"]) < 1e-3,
        "y_b_converges_monotonically": abs(convergence[0]["y_b_error"])
        > abs(convergence[-1]["y_b_error"]),
        "volume_is_skew_free": max(abs(row["volume_m3"] - skew[0]["volume_m3"]) for row in skew)
        < 1e-12 * skew[0]["volume_m3"],
        "y_f_over_y_b_is_4_3": all(
            abs(row["y_f_over_y_b"] - 4.0 / 3.0) < 1e-9
            for row in skew
            if not math.isnan(row["y_f_over_y_b"])
        ),
        "series_and_simpson_agree_1e-6": max(
            abs(row["volume_series_vs_simpson"]) for row in skew
        )
        < 1e-6,
        "bounding_box_breadth_is_the_declared_beam": abs(
            (box.maximum[1] - box.minimum[1]) - 2.0
        )
        < 1e-3,
    }
    payload = {
        "case": {
            "length_m": REFERENCE["length"],
            "beam_m": REFERENCE["beam"],
            "draft_m": REFERENCE["draft"],
            "transom_fraction": REFERENCE["transom_fraction"],
            "grids": [list(g) for g in GRIDS],
            "skews": list(SKEWS),
        },
        "convergence": convergence,
        "skew_sweep": skew,
        "observed_orders": orders,
        "claims": claims,
        "figures": [name for name, _ in figures],
    }
    results = output / "results.json"
    results.write_text(json.dumps(payload, indent=2, sort_keys=True))
    blobs = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for _, p in figures}
    blobs[results.name] = hashlib.sha256(results.read_bytes()).hexdigest()
    tables = _REPO_ROOT / "docs" / "transom_asymmetric" / "tables.md"
    blobs["tables.md"] = hashlib.sha256(tables.read_bytes()).hexdigest()
    (output / "manifest.json").write_text(json.dumps({"sha256": blobs}, indent=2, sort_keys=True))
    print(json.dumps({"claims": claims, "orders": orders}, indent=2, sort_keys=True))
    for row in convergence:
        print(
            f"nx={row['nx']:4d} nz={row['nz']:3d} V={row['volume_m3']:.9f} "
            f"dV={100 * row['volume_error']:+.4f}% y_B={row['y_b_m']:+.9f} "
            f"dy_B={100 * row['y_b_error']:+.4f}%"
        )
    for row in skew:
        print(
            f"skew={row['skew']:+.2f} V={row['volume_m3']:.9f} Aw={row['waterplane_area_m2']:.9f} "
            f"y_B={row['y_b_m']:+.9f} y_F={row['y_f_m']:+.9f} "
            f"IT={row['inertia_centreline_m4']:.6f} ITcg={row['inertia_centroid_m4']:.6f}"
        )
    failures = [name for name, ok in claims.items() if not ok]
    if failures:
        print("CLAIM FAILURES:", failures, file=_sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())