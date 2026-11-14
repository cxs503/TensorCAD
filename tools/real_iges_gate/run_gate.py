"""Run the CAD import pipeline on *publisher* IGES files and score it against published values.

The anchor's inputs used to be IGES files the repository wrote to itself with
``geomquality.write_iges`` (audit m17343, §3.2): a round trip, not an import.  This gate
takes the real specimens -- the geometry the publisher actually ships -- out of the asset
store (``asset://`` pointers, ``assets/geometry_registry.json``), runs the same
``import_cad`` the CLI uses, and compares the *mesh* against the principal dimensions and
displacement published for that hull in
``benchmarks/verified/baseline_fleet/fleet_specs.json``.

Nothing is hard-coded here: specimens, published numbers and tolerances all come from
those two JSON files, and the only measuring instrument is ``trimesh`` (an independent
library, not this package's own audit) plus the pipeline's own report.

    python tools/real_iges_gate/run_gate.py --ships KCS KVLCC2 KP505 \
        --out artifacts/real_iges_gate

Exit status is 0 only when every requested specimen was imported without a refused face
and every quantity that has a published reference agrees within ``--tolerance``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.trim_mesh import import_cad  # noqa: E402


def load_registry(root: Path) -> dict[str, Any]:
    return json.loads((root / "assets" / "geometry_registry.json").read_text())


def resolve(uri: str, registry: dict[str, Any], root: Path) -> Path:
    """``asset://geometry/hulls/KCS/...`` -> the store path (or a repo-relative path)."""
    if not uri.startswith("asset://"):
        return root / uri
    return Path(registry.get("store_root", "/srv/tensor/agent/assets")) / uri[len("asset://"):]


def measure(mesh_path: Path, published: dict[str, Any] | None) -> dict[str, Any]:
    """Independent measurement of the written STL: extents, area, displacement."""
    import trimesh

    mesh = trimesh.load(mesh_path, force="mesh")
    extents = [float(value) for value in mesh.extents]
    out: dict[str, Any] = {
        "triangles": int(len(mesh.faces)),
        "extents": extents,
        "area_m2": float(mesh.area),
        "watertight": bool(mesh.is_watertight),
    }
    if published:
        draft = float(published.get("draft", 0.0))
        low, high = mesh.bounds
        # The publisher ships x along the length, y to starboard, z up.  KCS/KVLCC2 are
        # *half* hulls (the shell is mirrored on the centre plane), so the beam is twice
        # the y extent; the code measures the mirror instead of assuming it.
        half = abs(float(low[1])) <= 1e-6 * max(1.0, float(high[1]))
        out["axes_assumed"] = {"length": "x", "beam": "y", "depth": "z", "half_hull": half}
        out["length_bbox"] = float(high[0] - low[0])
        out["beam_bbox"] = float(2.0 * (high[1] - low[1]) if half else high[1] - low[1])
        out["depth_bbox"] = float(high[2] - low[2])
        out["published_length_pp"] = published.get("length_pp")
        out["published_beam"] = published.get("beam")
        out["published_draft"] = draft
    volume = None
    if mesh.is_watertight:
        volume = abs(float(mesh.volume))
    else:
        # The KCS/KVLCC2 specimens are the shell *below the design waterline*, so the only
        # open boundary is the waterline/transom cut.  Close it (an independent library,
        # not this package's repair) and take the enclosed volume; say so in the record.
        capped = mesh.copy()
        trimesh.repair.fill_holes(capped)
        if capped.is_watertight:
            volume = abs(float(capped.volume))
            out["volume_source"] = "waterline cap (trimesh.repair.fill_holes)"
    if volume is None:
        out["displacement_m3"] = None
        out["displacement_note"] = ("the specimen is the shell below the design waterline and its "
                                    "waterline/transom cut is not planar-closable, so no enclosed "
                                    "volume is claimed (WITHHELD, not estimated)")
    if volume is not None:
        out["displacement_m3"] = volume
        if published and published.get("published_volume_m3"):
            reference = float(published["published_volume_m3"])
            out["published_volume_m3"] = reference
            out["volume_relative_error_pct"] = 100.0 * (volume - reference) / reference
        length = out.get("length_bbox")
        beam = out.get("beam_bbox")
        draft = out.get("published_draft")
        if published and length and beam and draft:
            out["cb_from_mesh"] = volume / (length * beam * draft)
            if published.get("published_cb"):
                out["published_cb"] = published["published_cb"]
                out["cb_relative_error_pct"] = 100.0 * (
                    out["cb_from_mesh"] - float(published["published_cb"])) / float(published["published_cb"])
    return out


def run_one(ship: str, spec: dict[str, Any], registry: dict[str, Any], root: Path,
            out_dir: Path, tolerance: float, kwargs: dict[str, Any]) -> dict[str, Any]:
    source = resolve(spec["iges_path"], registry, root)
    record: dict[str, Any] = {
        "ship": ship,
        "asset": spec["iges_path"],
        "path": str(source),
        "published_source": spec.get("published_source"),
        "publisher": spec.get("publisher"),
    }
    if not source.exists():
        record["status"] = "absent"
        record["note"] = "geometry not present in the asset store; nothing was measured"
        return record
    record["sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    record["bytes"] = source.stat().st_size
    record["sha256_matches_spec"] = record["sha256"] == spec.get("iges_sha256")
    started = time.time()
    result = import_cad(source, **kwargs)
    record["elapsed_seconds"] = round(time.time() - started, 3)
    report = result.report
    record["entity_census"] = report["import"]["entity_census"]
    record["faces_in_file"] = report["import"]["faces_in_file"]
    record["faces_tessellated"] = report["import"]["faces_tessellated"]
    record["faces_refused"] = report["import"]["faces_refused"]
    record["reader_refused"] = report["import"]["reader_refused"]
    record["max_inversion_residual"] = report["import"]["max_inversion_residual"]
    record["residual_tolerance"] = report["import"]["residual_tolerance"]
    record["chord_tolerance"] = report["import"]["chord_tolerance"]
    record["pipeline"] = {
        "triangles_before_repair": report["before_repair"]["triangles"],
        "free_edges_before_repair": report["before_repair"]["free_edges"],
        "triangles_after_repair": report["after_repair"]["triangles"],
        "free_edges_after_repair": report["after_repair"]["free_edges"],
        "components_after_repair": report["after_repair"]["components"],
        "watertight_after_repair": report["after_repair"]["watertight"],
    }
    record["pipeline_gates"] = {name: verdict["passed"] for name, verdict in report["gates"].items()}
    stl = out_dir / f"{ship}.stl"
    from tensorcad.geomquality.trim_mesh import write_mesh
    write_mesh(result.mesh, stl)
    record["stl"] = str(stl)
    record["stl_sha256"] = hashlib.sha256(stl.read_bytes()).hexdigest()
    measured = measure(stl, spec)
    record["mesh"] = measured
    errors = {}
    if "volume_relative_error_pct" in measured:
        errors["displacement"] = abs(measured["volume_relative_error_pct"])
    if "cb_relative_error_pct" in measured:
        errors["cb"] = abs(measured["cb_relative_error_pct"])
    for axis, published_key in (("length_bbox", "length_pp"), ("beam_bbox", "beam")):
        reference = spec.get(published_key)
        if reference and axis in measured:
            errors[axis] = 100.0 * abs(measured[axis] - float(reference)) / float(reference)
    record["relative_errors_pct"] = errors
    record["within_tolerance"] = {name: value <= tolerance for name, value in errors.items()}
    record["status"] = "measured"
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ships", nargs="+", default=["KCS"],
                        help="ship names in benchmarks/verified/baseline_fleet/fleet_specs.json")
    parser.add_argument("--out", default="artifacts/real_iges_gate")
    parser.add_argument("--tolerance", type=float, default=3.0,
                        help="acceptance in percent (the COG standard: >15%% is unacceptable)")
    parser.add_argument("--chord-tolerance", type=float, default=None)
    parser.add_argument("--repair/--no-repair", dest="repair", default=True)
    args = parser.parse_args(argv)

    registry = load_registry(ROOT)
    specs = json.loads((ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json").read_text())
    by_name = {ship["name"]: ship for ship in specs["ships"]}
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    kwargs: dict[str, Any] = {"repair": args.repair}
    if args.chord_tolerance is not None:
        kwargs["chord_tolerance"] = args.chord_tolerance

    records = []
    for ship in args.ships:
        spec = by_name.get(ship)
        if spec is None:
            records.append({"ship": ship, "status": "unknown_spec"})
            continue
        record = run_one(ship, spec, registry, ROOT, out_dir, args.tolerance, kwargs)
        records.append(record)
        if record["status"] == "measured":
            print(f"{ship}: faces {record['faces_tessellated']}/{record['faces_in_file']} "
                  f"tri {record['pipeline']['triangles_before_repair']} "
                  f"errors {json.dumps({k: round(v, 3) for k, v in record['relative_errors_pct'].items()})}")
        else:
            print(f"{ship}: {record['status']} -- {record.get('note', '')}")

    summary = {
        "tool": "tools/real_iges_gate/run_gate.py",
        "tolerance_pct": args.tolerance,
        "ships": records,
        "all_within_tolerance": all(
            record.get("status") == "measured" and all(record.get("within_tolerance", {}).values())
            and record.get("faces_refused") == 0
            for record in records
        ) if records else False,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"all_within_tolerance": summary["all_within_tolerance"],
                      "out": str(out_dir / "summary.json")}))
    return 0 if summary["all_within_tolerance"] else 1


if __name__ == "__main__":
    sys.exit(main())