#!/usr/bin/env python3
"""Close the two remaining gate holes of the baseline fleet: JBC's scale, ONRT's frame.

Why this tool exists
--------------------
``GAP-CAD-T1-BASELINE-FLEET`` asks for six standard hulls, watertight, with the
displacement reproduced to 1 % *at the publisher's own scale*.  Two archives were
``WITHHELD`` and neither refusal was about the shape:

* **JBC** ships in millimetres at a model scale the geometry page does not state, so the
  producer divided the *untrimmed* breadth by B = 45 m and got lambda = 39.7285.  The
  workshop report's own statement is lambda = 40 (280 m / 7.00 m), and the archive's
  half-breadth 562.5 mm x 40 = 45.000 m reproduces the published beam exactly, so the
  derived factor carried a 0.68 % length bias -- 2.06 % on a volume.
* **ONRT** ships as a *unit-Lpp* model (the bare-hull sheet spans 1.00614 in x) whose
  superstructure sits above z = 0.  The spec's 48.9355 is the workshop *model* scale of a
  3.147 m model, which a normalised archive does not use; worse, the sampling spacing was
  expressed in file units, so "8 samples per metre" meant *eight samples over the whole
  hull*.

Both are repaired as data, not as a fitted number: the scale comes from scale-free ratios
the archive itself carries (B/Lpp, T/Lpp) and the waterline is declared and then *checked*
against the published draft and beam before any volume is graded.

Usage
-----
    python tools/baseline_fleet/frame_and_scale_closure.py --ships JBC ONRT \
        --out runs/m21806_fleet_closure/evidence.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet_baseline import ShipSpec, build_ship, specs_from_json  # noqa: E402

#: Where each archive lives when it has not been ingested into the asset store.  The bytes
#: are not redistributed; their digests are pinned in the spec file.
LOCAL_ARCHIVES = {
    "KCS": "KCS_hull_SVA.igs",
    "KVLCC2": "kvlcc2.igs",
    "JBC": "jbc.igs",
    "ONRT": "onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
}


def _resolve(spec: ShipSpec, iges_dir: Path | None) -> str:
    if Path(spec.iges_path).is_file():
        return spec.iges_path
    if iges_dir is not None:
        candidate = iges_dir / LOCAL_ARCHIVES.get(spec.name, "")
        if candidate.is_file():
            return str(candidate)
    return spec.iges_path


def _compact(report: dict) -> dict:
    """The subset of the producer's report the closure table needs, copied verbatim."""
    keep_source = (
        "unit_gate",
        "unit_flag",
        "unit_name",
        "geometry_scale_to_full",
        "scale_basis",
        "scale_inferred",
        "frame_shift_z",
        "route",
        "sample_count",
        "sampling_spacing_m",
        "sampling_spacing_file_units",
        "metres_per_file_unit_for_spacing",
        "trim_sampling",
    )
    return {
        "verdict": report["verdict"],
        "reasons": report["reasons"],
        "checks": report["checks"],
        "volume_gate_basis": report["volume_gate_basis"],
        "source": {key: report["source"].get(key) for key in keep_source},
        "frame_check": report.get("frame_check"),
        "mesh": {
            key: report["mesh"][key]
            for key in (
                "volume_m3",
                "is_watertight",
                "is_edge_manifold",
                "boundary_edges",
                "zero_area_faces",
                "triangle_count",
                "extent_m",
            )
        },
        "derived": report["derived"],
        "field": {
            key: report["field"][key]
            for key in ("nx", "nz", "cover", "max_half_breadth_m", "levels_with_data")
            if key in report["field"]
        },
        "published": {
            "volume_m3": report["spec"]["published_volume_m3"],
            "cb": report["spec"]["published_cb"],
            "lpp": report["spec"]["length_pp"],
            "beam": report["spec"]["beam"],
            "draft": report["spec"]["draft"],
        },
    }


def _figure(hulls: dict, path: Path) -> str | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is optional
        return f"matplotlib unavailable: {exc}"

    names = list(hulls)
    figure, axes = plt.subplots(1, 1 + len(names), figsize=(4.6 * (1 + len(names)), 4.2))
    errors = [
        100.0 * hulls[name]["derived"]["published_volume_rel_error_percent"] / 100.0
        for name in names
    ]
    bars = axes[0].bar(names, errors, color="#2b6cb0")
    axes[0].axhline(1.0, color="#c53030", linestyle="--", linewidth=1.2, label="1 % tolerance")
    axes[0].axhline(-1.0, color="#c53030", linestyle="--", linewidth=1.2)
    axes[0].axhline(0.0, color="#4a5568", linewidth=0.8)
    axes[0].set_ylabel("V_mesh - V_published   [%]")
    axes[0].set_title("Displacement agreement after the fix")
    axes[0].legend(loc="lower right", fontsize=8)
    for bar, value in zip(bars, errors):
        axes[0].annotate(
            f"{value:+.2f}%",
            (bar.get_x() + bar.get_width() / 2.0, value),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    for axis, name in zip(axes[1:], names):
        cloud = hulls[name]["_cloud"]
        axis.scatter(cloud["x"], cloud["y"], c=cloud["z"], s=0.6, cmap="viridis")
        axis.set_aspect("equal")
        axis.set_xlabel("x [m]")
        axis.set_ylabel("y [m]")
        axis.set_title(f"{name}: reconstructed hull (colour = z)")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=140)
    return str(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specs", default="benchmarks/verified/baseline_fleet/fleet_specs.json")
    parser.add_argument("--iges-dir", default="artifacts/baseline_fleet/iges")
    parser.add_argument("--ships", nargs="*", default=None)
    parser.add_argument("--out", default="runs/m21806_fleet_closure/evidence.json")
    parser.add_argument("--figure", default="runs/m21806_fleet_closure/closure.png")
    parser.add_argument("--dx", type=float, default=0.5)
    parser.add_argument("--dz", type=float, default=0.1)
    parser.add_argument("--cloud-stride", type=int, default=40)
    args = parser.parse_args(argv)

    iges_dir = Path(args.iges_dir) if args.iges_dir else None
    specs = specs_from_json(args.specs)
    wanted = set(args.ships) if args.ships else None
    out: dict = {
        "schema": "tensorcad.fleet.frame_scale_closure/1.0",
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "resolution": {"dx_m": args.dx, "dz_m": args.dz},
        "hulls": {},
    }
    for spec in specs:
        if wanted is not None and spec.name not in wanted:
            continue
        spec.iges_path = _resolve(spec, iges_dir)
        if not Path(spec.iges_path).is_file():
            out["hulls"][spec.name] = {"verdict": "WITHHELD", "reasons": ["archive absent"]}
            continue
        mesh, report = build_ship(spec, dx=args.dx, dz=args.dz)
        block = _compact(report)
        vertices = mesh.vertices[:: max(args.cloud_stride, 1)]
        block["_cloud"] = {
            "x": [round(vertex[0], 2) for vertex in vertices],
            "y": [round(vertex[1], 2) for vertex in vertices],
            "z": [round(vertex[2], 2) for vertex in vertices],
        }
        out["hulls"][spec.name] = block

    clouds = {name: block.pop("_cloud") for name, block in out["hulls"].items()}
    figure_note = _figure(
        {name: {**block, "_cloud": clouds[name]} for name, block in out["hulls"].items()},
        Path(args.figure),
    )
    out["figure"] = {"path": figure_note, "caption": (
        "left: reconstructed displacement minus the publisher's value, per hull, against the "
        "1 % tolerance of the anchor; right: the reconstructed watertight hull below the design "
        "waterline (x-y projection, colour = z), the object the volume is integrated over."
    )}
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(out, indent=1))
    for name, block in out["hulls"].items():
        error = (block.get("derived") or {}).get("published_volume_rel_error_percent")
        volume = (block.get("mesh") or {}).get("volume_m3")
        print(
            f"{name:8s} {block['verdict']:9s} "
            f"V={volume if volume is None else round(volume, 1)} m3 "
            f"err={error if error is None else round(error, 4)} %"
        )
    return 0 if all(b["verdict"] == "VERIFIED" for b in out["hulls"].values()) else 1


if __name__ == "__main__":
    sys.exit(main())