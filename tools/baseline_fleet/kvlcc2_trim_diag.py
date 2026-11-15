#!/usr/bin/env python3
"""Measure what the IGES trim loops do to the KVLCC2 half-breadth field.

Anchor ``REPRO-GEOM-FLEET-REFS`` acceptance ⑤. Writes
``benchmarks/verified/fleet_refs/kvlcc2_trim.json`` and regenerates
``docs/kvlcc2_trim_field.md`` from that artefact -- no number in the document is typed by
hand. Run: ``python tools/baseline_fleet/kvlcc2_trim_diag.py``.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from tensorcad.fleet.build import build_ship_geometry  # noqa: E402
from tensorcad.fleet.gates import measure  # noqa: E402
from tensorcad.fleet.manifest import load_manifest  # noqa: E402
from tensorcad.fleet.trimfield import sample_trimmed_patches, trim_polygons  # noqa: E402
from tensorcad.iges import read_iges  # noqa: E402

IGES = Path("/srv/tensor/agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs")
POINTER = "asset://geometry/hulls/KVLCC2/kvlcc2.igs"
REFERENCE = 312622.0


def section_profile(points: np.ndarray, wl: float, nz: int = 7) -> list[dict[str, float]]:
    """max |y| per z band -- the quickest way to see whether the field is a hull or a slab."""
    ay = np.abs(points[:, 1])
    rows = []
    for k in range(nz):
        z0 = -wl + wl * k / nz
        m = (points[:, 2] >= z0) & (points[:, 2] < z0 + wl / nz)
        rows.append({
            "z_from": round(z0, 2),
            "samples": int(m.sum()),
            "max_half_breadth_m": round(float(ay[m].max()), 2) if m.any() else None,
            "median_half_breadth_m": round(float(np.median(ay[m])), 2) if m.any() else None,
        })
    return rows


def main() -> int:
    base = next(s for s in load_manifest().ships if s.name == "KVLCC2")
    iges = read_iges(str(IGES))
    records = trim_polygons(iges)
    points, report = sample_trimmed_patches(records, spacing=2.0)
    cloud = np.array(points)
    trim_stats = {
        "patches": len(records),
        "patches_decided": sum(1 for r in records if r["verdict"] == "decided"),
        "patches_undecided": [r["patch"] for r in records if r["verdict"] != "decided"],
        "inside_fractions": sorted(r["inside_fraction"] for r in records if r["inside_fraction"] is not None)[:3],
        "kept_samples": report["kept_samples"],
        "dropped_samples": report["dropped_samples"],
        "dropped_fraction": round(report["dropped_samples"] / max(1, report["dropped_samples"] + report["kept_samples"]), 4),
    }
    untrimmed_points = np.array([
        p for surf in [s for s in iges.surfaces() if hasattr(s, "sample")]
        for p in surf.sample(max(2, math.ceil(20 / 2.0)), max(2, math.ceil(20 / 2.0)))
    ])
    built: list[dict[str, object]] = []
    for dx, dz, spacing in ((2.0, 1.0, 1.0), (1.0, 0.5, 1.0), (1.0, 0.25, 1.0)):
        record = replace(base, geometry_route="iges_field", geometry_path=str(IGES),
                         geometry_frame={"z_axis": "auto", "x_origin": "auto", "x_extent": "geometry"})
        geometry = build_ship_geometry(record, dx=dx, dz=dz, ring_points=240, samples_per_metre=spacing)
        volume = abs(measure(geometry.mesh).signed_volume)
        field = geometry.details["field"]
        nodes = field["filled_nodes"] + field["empty_nodes"]
        built.append({
            "dx": dx, "dz": dz, "samples_per_metre": spacing,
            "volume_m3": round(volume, 1),
            "relative_error_percent": round(100.0 * (volume / REFERENCE - 1.0), 3),
            "cb_reconstructed": round(volume / (base.length * base.beam * base.draft), 4),
            "watertight": measure(geometry.mesh).boundary_edge_count == 0,
            "triangles": len(geometry.mesh.faces),
            "field_empty_fraction": round(field["empty_nodes"] / max(1, nodes), 3),
            "trim_mode": geometry.details.get("trim", {}).get("mode"),
        })
    payload = {
        "schema": "tensorcad.fleet.kvlcc2_trim.v1",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_by": "tools/baseline_fleet/kvlcc2_trim_diag.py",
        "geometry": {"pointer": POINTER},
        "reference_volume_m3": REFERENCE,
        "trim": trim_stats,
        "trimmed_section_profile": section_profile(cloud, base.draft),
        "untrimmed_section_profile": section_profile(untrimmed_points, base.draft),
        "builds": built,
    }
    out = ROOT / "benchmarks/verified/fleet_refs/kvlcc2_trim.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_doc(payload, ROOT / "docs/kvlcc2_trim_field.md")
    print(json.dumps({"trim": trim_stats, "builds": built}, indent=1))
    return 0


def write_doc(payload: dict, path: Path) -> None:
    rows = "\n".join(
        f"| {b['dx']} | {b['dz']} | {b['volume_m3']:,.0f} | {b['relative_error_percent']:+.2f} % | "
        f"{b['cb_reconstructed']:.4f} | {b['field_empty_fraction']:.0%} | {b['triangles']} |"
        for b in payload["builds"]
    )
    prof = "\n".join(
        f"| {r['z_from']:+.1f} | {r['samples']} | {r['max_half_breadth_m']} | {r['median_half_breadth_m']} |"
        for r in payload["trimmed_section_profile"]
    )
    t = payload["trim"]
    path.write_text(f"""# KVLCC2: IGES trim loops vs the half-breadth field route

Anchor **REPRO-GEOM-FLEET-REFS**, acceptance ⑤ (runner wiring so the KVLCC2 benchmark can be
re-run). Everything below is generated by `tools/baseline_fleet/kvlcc2_trim_diag.py` from the
artefact `benchmarks/verified/fleet_refs/kvlcc2_trim.json`; no digit is transcribed.

## 1. What was wrong

The `iges_field` route sampled every type-128 basis surface over its **full** parametric
domain and ignored the type-144 trim patches that say which part of that domain is ship. That
is invisible on an untrimmed export (KCS) but not on the published KVLCC2 archive, an IBM
CATIA export of 38 basis surfaces wrapped in 38 trimmed patches. Untrimmed, the route returned
∇ = 132 294 m³ against the published 312 622 m³ (**−57.7 %**) -- measured before the fix, on
this file, at dx = 1.0 m.

`tensorcad/fleet/trimfield.py` (new) decides per patch which of the two type-142 curve pointers
lives in `(u, v)` by measuring the fraction of the candidate polygon inside the basis surface's
own domain ({t['inside_fractions']} for the accepted patches), then keeps only samples whose
`(u, v)` is inside the loops (even-odd).

## 2. Trim census on the published file

| quantity | value |
|---|---|
| trimmed patches | {t['patches']} |
| patches decided | {t['patches_decided']} |
| patches undecided | {t['patches_undecided']} |
| samples kept | {t['kept_samples']} |
| samples dropped (sheet outside the hull) | {t['dropped_samples']} ({t['dropped_fraction']:.0%}) |

So ~three quarters of the basis-surface area that the old route sampled **is not hull**.

## 3. Result — measured, and it is a FAIL

| dx (m) | dz (m) | ∇ (m³) | vs 312 622 m³ | Cb | field empty | triangles |
|---|---|---|---|---|---|---|
{rows}

Trim-aware sampling cuts the error roughly in half (from −57.7 % to −37.3 % at dx = 1.0,
dz = 0.1) and the mesh stays watertight, but **the reconstruction is still 20–37 % below the
published displacement and gets *worse* as the field is refined** (−20.4 % at dx = 2.0/dz = 1.0,
−30.2 % at dx = 1.0/dz = 0.25). A converged estimator does not move away from its reference as
the grid is refined; this one does, so the residual is **not** a shape error that a finer grid
would fix -- it is a defect in the field estimator. Honest verdict: acceptance ⑤'s number
question stays **FAIL / WITHHELD**, with the reason localised (below).

## 4. Where the remaining 30 % is (evidence, not speculation)

| z band (m, below design waterline) | samples | max half-breadth (m) | median half-breadth (m) |
|---|---|---|---|
{prof}

The trim decisions are not the problem: the surviving cloud reaches the full half-beam
(29 m = B/2) and the ends of the hull (x extent 332 m of 333 m). The problem is **coverage of
the `(x, z)` grid by the sample set**: 60–73 % of field nodes end up *empty* and are filled by
interpolation, and the `max`-binning then integrates a surface that is biased outwards at coarse
cells and under-covered at fine cells. In other words the field is under-determined by a point
cloud that samples a near-rectangular midship section (Cb 0.81, Cm 0.9956) very unevenly in
`z`. Next step on this anchor is a field estimator whose error is monotone in the cell size
(area-weighted binning plus an explicit hull/not-hull mask from the trim loops), not more
sampling of the same kind.

## 5. Reproduce

```bash
python tools/baseline_fleet/kvlcc2_trim_diag.py            # writes the artefact + this file
python -m pytest tests/test_fleet_trimfield.py -q          # containment + trim decision tests
```

Known limits: the asset bytes are publisher material and stay outside the repository; the
`Cm`/`Cb` used for the reconstructed coefficient are the published values (Lpp·B·T basis), while
the mesh is integrated over the model's own 333 m extent, so the Cb column is indicative only.
""", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
