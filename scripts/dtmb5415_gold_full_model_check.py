"""The anchor's gold reference is a both-sides model: measure the halving on it.

``DTMB5415_Hull_exact.stl`` (the gold regression of GAP-CAD-T1-IGES-PIPE) carries both
sides of the hull, so ``classify_model`` reports ``symmetry_factor = 1`` and the pre-fix
``hull_characteristics`` read ``max |y|`` -- half the breadth -- at every station.  This
script measures that on the real gold file with the shipped code path and judges the
repaired reading against an **independent** estimator: the area of the section polylines
(``gold_loft.section_polyline``, i.e. the anchor's own gold tooling) closed at the
waterline and integrated with the shoelace formula -- no ``envelope_breadth`` anywhere.

Reproduce::

    python scripts/dtmb5415_gold_full_model_check.py     # -> runs/dtmb5415_full_model/
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "dtmb5415_full_model"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "benchmarks" / "verified" / "iges_import"))

from gold_loft import read_stl, section_polyline  # noqa: E402

from tensorcad.geomquality.hull_surface import (  # noqa: E402
    classify_model,
    hull_characteristics,
)

GOLD = Path("/root/ship-performance-platform/backend/data/geometry/ships/"
            "DTMB5415_Hull_exact.stl")
WL = 6.15  # design draft of the DTMB5415 model in the platform's geometry frame
STATIONS, LEVELS = 21, 24  # the pair that fits the local wall clock; disclosed in the doc


def shoelace_section_area(triangles: np.ndarray, x0: float, waterline: float) -> tuple[float, int]:
    """Immersed area of the section closed at the waterline: an independent estimator."""
    total = 0.0
    used = 0
    for chain in section_polyline(triangles, x0):
        if len(chain) < 2:
            continue
        pts = [(float(p[1]), float(p[2])) for p in chain if p[2] <= waterline + 1e-9]
        if len(pts) < 2:
            continue
        # close the open section along the waterline between its two ends
        pts = pts + [(pts[-1][0], waterline), (pts[0][0], waterline)]
        area = 0.0
        for (y0, z0), (y1, z1) in zip(pts, pts[1:]):
            area += y0 * z1 - y1 * z0
        total += abs(area) / 2.0
        used += 1
    return total, used


def main() -> int:
    triangles = read_stl(GOLD)
    flat = triangles.reshape(-1, 3)
    unique, inverse = np.unique(flat, axis=0, return_inverse=True)
    vertices = [tuple(map(float, point)) for point in unique]
    faces = [tuple(map(int, row)) for row in inverse.reshape(-1, 3)]
    xs = triangles[:, :, 0]
    x_lo, x_hi = float(xs.min()), float(xs.max())
    step = (x_hi - x_lo) / STATIONS
    grid = [x_lo + (x_hi - x_lo) * (i + 0.5) / STATIONS for i in range(STATIONS)]

    evidence: dict = {
        "anchor": "GAP-CAD-T1-IGES-PIPE",
        "input": {"path": str(GOLD), "triangles": int(len(faces)),
                  "vertices": int(len(vertices))},
        "waterline_m": WL,
        "stations": STATIONS,
        "levels": LEVELS,
        "model": classify_model(vertices),
        "repro": "python scripts/dtmb5415_gold_full_model_check.py",
    }

    runs = {}
    for flag, label in ((False, "before_max_abs_y"), (True, "after_full_span")):
        t0 = time.time()
        characteristic = hull_characteristics(vertices, faces, waterline=WL,
                                              stations=STATIONS, levels=LEVELS,
                                              dedup=False, full_breadth=flag)
        runs[label] = {
            "full_breadth": flag,
            "seconds": round(time.time() - t0, 3),
            "volume_m3": characteristic.volume,
            "beam_m": characteristic.beam,
            "lwl_m": characteristic.lwl,
            "cb": characteristic.cb,
            "cm": characteristic.cm,
            "s_wetted_naive_m2": characteristic.s_wetted_naive,
        }
    evidence["runs"] = runs
    evidence["halving"] = {
        "volume_ratio_after_over_before": (runs["after_full_span"]["volume_m3"]
                                           / runs["before_max_abs_y"]["volume_m3"]),
        "beam_ratio_after_over_before": (runs["after_full_span"]["beam_m"]
                                         / runs["before_max_abs_y"]["beam_m"]),
        "cb_ratio_after_over_before": (runs["after_full_span"]["cb"]
                                       / runs["before_max_abs_y"]["cb"]),
    }

    # independent estimator: shoelace on the section polylines
    t0 = time.time()
    areas = []
    for station in grid:
        area, _ = shoelace_section_area(triangles, station, WL)
        areas.append(area)
    reference_volume = sum(areas) * step
    evidence["independent_reference"] = {
        "method": "shoelace area of gold_loft section polylines closed at the waterline",
        "seconds": round(time.time() - t0, 3),
        "volume_m3": reference_volume,
        "sectional_areas_m2": areas,
        "stations_x_m": grid,
    }
    for label in runs:
        runs[label]["error_vs_independent_percent"] = (
            100.0 * (runs[label]["volume_m3"] - reference_volume) / reference_volume)
    evidence["verdict"] = {
        "before_error_percent": runs["before_max_abs_y"]["error_vs_independent_percent"],
        "after_error_percent": runs["after_full_span"]["error_vs_independent_percent"],
        "after_within_3pct": bool(
            abs(runs["after_full_span"]["error_vs_independent_percent"]) <= 3.0),
        "before_was_a_halving": bool(
            abs(evidence["halving"]["volume_ratio_after_over_before"] - 2.0) < 0.05),
    }
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2, default=str))
    print(json.dumps(evidence["verdict"], indent=2))
    print(json.dumps(evidence["halving"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
