#!/usr/bin/env python3
"""Before/after measurement of the JBC geometry bias (anchor GAP-CAD-T1-SOURCE-HUNT ③).

    python tools/baseline_fleet/jbc_scale_impact.py            # build both variants
    python tools/baseline_fleet/jbc_scale_impact.py --verify-only

Builds the JBC hull twice from the *same* IGES bytes through the production route
(:func:`tensorcad.fleet_baseline.build_ship`) and only differs in ``geometry_scale``:

* ``inferred`` -- what the manifest carried before the derived scale was adopted:
  ``geometry_scale = None`` makes the builder fit the hull's **x-extent** to the
  published ``L_PP``.  JBC's x-extent is the deck-level extreme (bow flare + stern
  overhang), 7.5428 m against 7.0478 m of waterline-parallel body, so the fit shrinks
  the ship: full-scale ``L_PP`` comes out ~6 % short while every volume gate still
  passes, because the same factor rescales the volume.
* ``derived`` -- the scale from ``artifacts/baseline_fleet/jbc_model_scale.json``
  (1:39.7285 from two independent published lengths).

That is the whole point of the anchor item: a length-fit left uncorrected is invisible
to a displacement gate.  The script writes
``artifacts/baseline_fleet/jbc_scale_impact.json``; nothing is typed by hand, every
number is read back from the two build reports and the record.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from tensorcad.fleet_baseline import ShipSpec, build_ship  # noqa: E402

RECORD = REPO / "artifacts" / "baseline_fleet" / "jbc_model_scale.json"
SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
OUT = REPO / "artifacts" / "baseline_fleet" / "jbc_scale_impact.json"

#: Coarse resolution: this tool measures *scale*, not convergence, and the bias under
#: test is 6 %, two orders of magnitude above the discretisation noise at this grid.
COARSE = {"dx": 0.5, "dz": 0.1, "ring_points": 240, "samples_per_metre": 8.0}


def locate_iges() -> Path:
    """First available copy of the JBC hull IGES.

    The publisher bytes are not redistributed, so the asset store is authoritative on
    a development host but absent on a cluster node; a checkout that ran
    ``jbc_model_scale.py`` keeps an extracted copy under ``artifacts/baseline_fleet/raw``
    (git-ignored).  ``TENSORCAD_JBC_IGES`` overrides both.
    """
    import os

    candidates = []
    override = os.environ.get("TENSORCAD_JBC_IGES")
    if override:
        candidates.append(Path(override))
    registry_path = REPO / "assets" / "geometry_registry.json"
    if registry_path.exists():
        registry = json.loads(registry_path.read_text())
        candidates.append(Path(registry["store_root"]) / "geometry" / "hulls" / "JBC" / "jbc.igs")
    candidates.append(REPO / "artifacts" / "baseline_fleet" / "raw" / "jbc.igs")
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit(
        "JBC hull IGES not found; run tools/baseline_fleet/jbc_model_scale.py first "
        "or set TENSORCAD_JBC_IGES"
    )


def jbc_spec(geometry_scale: float | None) -> ShipSpec:
    """JBC spec from the committed manifest, with the scale overridden for the test."""
    document = json.loads(SPECS.read_text())
    entry = next(ship for ship in document["ships"] if ship["name"] == "JBC")
    iges = locate_iges()
    return ShipSpec(
        name="JBC",
        draft=entry["draft"],
        length_pp=entry["length_pp"],
        beam=entry["beam"],
        published_volume_m3=entry["published_volume_m3"],
        published_cb=entry["published_cb"],
        published_lcb_percent=entry.get("published_lcb_percent"),
        geometry_scale=geometry_scale,
        scale_basis="impact test",
        published_source=entry.get("published_source", ""),
        source_url=entry.get("source_url", ""),
        publisher=entry.get("publisher", ""),
        iges_path=str(iges),
    )


def measure(report: dict) -> dict:
    """Pull the comparison quantities out of a build report.

    ``mesh.extent_m`` is *not* used as a length: it is clipped at the design waterline,
    so it carries the hull's largest submerged x-extent rather than ``L_PP`` or ``L_wl``.
    The hydrostatic block is the only honest source for a length.
    """
    hydro = report.get("derived", {})
    integrals = report.get("integrals", {})
    source = report.get("source", {})
    half_beam = integrals.get("half_breadth_max_m")
    return {
        "volume_m3": hydro.get("volume_m3"),
        "length_waterline_m": hydro.get("length_waterline_m"),
        "beam_waterline_m": None if half_beam is None else 2.0 * half_beam,
        "cb": hydro.get("cb"),
        "cwp": hydro.get("cwp"),
        "waterplane_area_m2": integrals.get("waterplane_area_m2"),
        "geometry_scale_to_full": source.get("geometry_scale_to_full"),
        "scale_inferred": source.get("scale_inferred"),
        "volume_error_percent": hydro.get("published_volume_rel_error_percent"),
        "extent_over_lpp": report.get("checks", {}).get("extent_over_lpp"),
        "lwl_vs_lpp_percent": report.get("checks", {}).get("lwl_vs_lpp_percent"),
        "verdict": report.get("verdict"),
        "reasons": report.get("reasons"),
    }


def run() -> dict:
    record = json.loads(RECORD.read_text())
    derived = float(record["scale_derived"]["value"])
    out: dict = {
        "schema": "tensorcad.fleet.jbc_scale_impact.v1",
        "anchor": "GAP-CAD-T1-SOURCE-HUNT",
        "hull": "JBC",
        "generated_by": "tools/baseline_fleet/jbc_scale_impact.py",
        "geometry": "geometry/hulls/JBC/jbc.igs (publisher CAD, model space, mm)",
        "derived_scale": derived,
        "record": str(RECORD.relative_to(REPO)),
        "published": {
            "length_waterline_m": 285.0,
            "beam_waterline_m": 45.0,
            "draft_m": 16.5,
            "volume_m3": 178369.9,
            "cb": 0.8580,
        },
        "variants": {},
    }
    for label, scale in (("inferred", None), ("derived", derived)):
        _, report = build_ship(jbc_spec(scale), **COARSE)
        out["variants"][label] = measure(report)
    # Errors against the published particulars, printed as the acceptance table.
    for label, got in out["variants"].items():
        errors = {}
        for key in ("length_waterline_m", "beam_waterline_m", "volume_m3", "cb"):
            reference = out["published"].get(key)
            value = got.get(key)
            errors[key] = (
                None if value is None or not reference
                else 100.0 * (value - reference) / reference
            )
        got["error_percent"] = errors
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {OUT}")
    for label, got in out["variants"].items():
        print(f"{label:9s} Lwl {got['length_waterline_m']:.2f} m "
              f"({got['error_percent']['length_waterline_m']:+.2f} %)  "
              f"B {got['beam_waterline_m']:.2f} m  Cb {got['cb']:.4f}  "
              f"nabla {got['volume_m3']:.0f} m3  scale {got['geometry_scale_to_full']:.5f}")
    return out


def verify() -> int:
    data = json.loads(OUT.read_text())
    problems = []
    grown = data["variants"]["derived"]
    shrunk = data["variants"]["inferred"]
    if shrunk["length_waterline_m"] >= grown["length_waterline_m"]:
        problems.append("the x-extent fit did not shrink the hull: bias model is wrong")
    error = grown["error_percent"]["length_waterline_m"]
    if abs(error) > 3.0:
        problems.append(f"derived-scale Lwl is {error:+.2f} % off the published value")
    beam_error = grown["error_percent"]["beam_waterline_m"]
    if abs(beam_error) > 3.0:
        problems.append(f"derived-scale B_wl is {beam_error:+.2f} % off the published value")
    bias = shrunk["error_percent"]["length_waterline_m"]
    if not -10.0 < bias < -3.0:
        problems.append(f"inferred-fit Lwl bias {bias:+.2f} % outside the 3-10 % band")
    for problem in problems:
        print(f"FAIL {problem}")
    if not problems:
        print(f"OK {OUT}: derived-scale Lpp error {error:+.2f} %, "
              f"x-extent fit {bias:+.2f} %")
    return 0 if not problems else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        return verify()
    run()
    return verify()


if __name__ == "__main__":
    raise SystemExit(main())