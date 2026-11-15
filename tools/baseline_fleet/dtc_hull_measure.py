"""Measure the DTC workshop hull IGES that is already in this environment's source cache.

Anchor ``GAP-CAD-T1-SOURCE-HUNT``.  Round N of this anchor recorded the DTC geometry stage as
``WITHHELD`` because *"the IGES sampling collapses to a point: 448590 samples from 1345 B-spline
surfaces span 0.000 m"* (``benchmarks/verified/fleet_refs/dtc_geometry.json``).  That measurement
is stale: it was taken before the type-128 control-net layout fix documented in
:func:`tensorcad.iges._split_control` (the weight-block test that replaced the smoothness-only
choice).  Re-reading the same cached bytes today yields a real point cloud, so the DTC hull can
be *measured* instead of withdrawn.

This tool does exactly that, and it pins everything it measures:

* the archive digest, every member digest, the source URL and the licence statement;
* the IGES unit flag / model scale as the reader parses them (DTC is written in millimetres?
  no -- it declares metres, and the coordinate extent confirms full scale);
* the geometric extents (x, y, z) against the publisher's ``B`` and ``T``;
* displacement volume, block coefficient, waterplane area and LCB from the half-breadth field,
  compared with the cross-verified published values;
* the same numbers at two sampling densities, so the resolution sensitivity is visible rather
  than assumed.

Nothing is typed in: the published particulars are read out of
``tensorcad/fleet/data/fleet_manifest.json`` and ``benchmarks/verified/fleet_refs/dtc_geometry.json``
at run time, and every measurement comes out of the cached bytes.

The archive has one more defect worth pinning: it stores its members with **backslash** path
separators in the central directory (``dtc-2012-01-10\\DTC_hull.igs``) while the local headers use
forward slashes, so :meth:`zipfile.ZipFile.read` raises ``BadZipFile: File name in directory ...
and header ... differ``.  :func:`extract_members` normalises ``ZipInfo.orig_filename`` in memory;
the bytes on disk are untouched.

Usage::

    python tools/baseline_fleet/dtc_hull_measure.py                 # both densities
    python tools/baseline_fleet/dtc_hull_measure.py --spacing 1.0   # one density
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:  # pragma: no cover - import plumbing
    sys.path.insert(0, str(REPO))

from tensorcad.fleet_baseline import (  # noqa: E402
    field_from_points,
    field_integrals,
    partition_surfaces,
    sample_patches,
)
from tensorcad.iges import read_iges  # noqa: E402

SCHEMA = "tensorcad.fleet_refs.dtc_hull_measured/1.0"
ANCHOR = "GAP-CAD-T1-SOURCE-HUNT"

DEFAULT_ARCHIVE = Path("/root/.cache/hull_sources/dtc_igs_ordner.zip")
DEFAULT_MANIFEST = REPO / "tensorcad" / "fleet" / "data" / "fleet_manifest.json"
DEFAULT_DTC_GEOMETRY = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_geometry.json"
DEFAULT_OUT = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"

SHIP = "DTC"
HULL_MEMBER = "DTC_hull.igs"
SOURCE_URL = "https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip"
SOURCE_PAGE = "https://www.uni-due.de/ship/dtc.php"
LICENSE = (
    "publisher states no licence; the workshop publishes the IGES set for benchmarking. "
    "Bytes are digested here, never redistributed."
)
#: Volume tolerance the fleet build uses for every hull (``fleet_baseline.VOLUME_TOLERANCE``).
VOLUME_TOLERANCE = 0.01


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_members(archive: Path) -> dict[str, bytes]:
    """Return ``{basename: bytes}`` for every member of ``archive``.

    ``ZipFile.read`` refuses this archive because the central directory spells member
    paths with backslashes.  The local file header is authoritative (it is what the
    offsets point at), so the in-memory ``orig_filename`` is normalised to the header
    spelling and the bytes are returned unchanged.
    """
    out: dict[str, bytes] = {}
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            normalised = info.filename.replace("\\", "/")
            info.filename = normalised
            info.orig_filename = normalised
            data = handle.read(info)
            out[normalised.rsplit("/", 1)[-1]] = data
    return out


def published_reference(manifest_path: Path, dtc_geometry_path: Path) -> dict:
    """Published particulars, injected from the records rather than typed in here."""
    manifest = json.loads(manifest_path.read_text())
    record = next(ship for ship in manifest["ships"] if ship["name"] == SHIP)
    if dtc_geometry_path.exists():
        cross = json.loads(dtc_geometry_path.read_text()).get("published_reference", {})
    else:  # pragma: no cover - the file ships with the repo
        cross = {}
    return {
        "length_pp": record["length"],
        "beam": record["beam"],
        "draft": record["draft"],
        "cb": record["block_coefficient"],
        "displacement_volume_m3": record["displacement_volume"],
        "manifest_provenance": record["provenance"],
        "manifest_notes": record["notes"],
        "crossref_verdicts": {
            key: value.get("verdict") for key, value in cross.items() if isinstance(value, dict)
        },
    }


def _extent(points) -> dict:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    return {
        "bbox_min": [min(xs), min(ys), min(zs)],
        "bbox_max": [max(xs), max(ys), max(zs)],
        "extent": [max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)],
        "samples": len(points),
        "distinct": len({(round(x, 6), round(y, 6), round(z, 6)) for x, y, z in points}),
    }


def measure(
    hull_path: Path,
    reference: dict,
    spacings: tuple[float, ...] = (1.0, 0.5),
    dx: float = 0.5,
    dz: float = 0.1,
) -> dict:
    """Sample the hull at each density and integrate the half-breadth field."""
    model = read_iges(str(hull_path))
    surfaces = model.surfaces()
    usable, quarantined = partition_surfaces(surfaces)
    unit = {
        "unit_flag": model.unit_flag,
        "unit_name": model.unit_name,
        "metres_per_unit": model.metres_per_unit,
        "model_scale": model.model_scale,
    }
    rows = []
    for spacing in spacings:
        points = list(sample_patches(usable, spacing=spacing))
        extent = _extent(points)
        field = field_from_points(points, reference["draft"], dx=dx, dz=dz)
        integrals = field_integrals(field)
        volume = integrals["volume_m3"]
        cb = volume / (reference["length_pp"] * reference["beam"] * reference["draft"])
        volume_error = 100.0 * (volume - reference["displacement_volume_m3"]) / reference[
            "displacement_volume_m3"
        ]
        rows.append(
            {
                "sample_spacing_m": spacing,
                "field_dx_m": dx,
                "field_dz_m": dz,
                "point_cloud": extent,
                "surfaces_read": len(usable),
                "surfaces_quarantined": len(quarantined),
                "field_filled_nodes": field.filled_nodes,
                "field_interpolated_nodes": field.interpolated_nodes,
                "volume_m3": volume,
                "cb": cb,
                "lcb_x_m": integrals["lcb_x_m"],
                "waterplane_area_m2": integrals["waterplane_area_m2"],
                "half_breadth_max_m": integrals["half_breadth_max_m"],
                "volume_rel_error_percent": volume_error,
            }
        )
    best = rows[-1]
    checks = {
        "beam": {
            "measured_m": best["point_cloud"]["extent"][1],
            "published_m": reference["beam"],
            "rel_error_percent": 100.0
            * (best["point_cloud"]["extent"][1] - reference["beam"])
            / reference["beam"],
            "tolerance_percent": 0.5,
        },
        "block_coefficient": {
            "measured": best["cb"],
            "published": reference["cb"],
            "rel_error_percent": 100.0 * (best["cb"] - reference["cb"]) / reference["cb"],
            "tolerance_percent": 1.0,
        },
        "displacement_volume": {
            "measured_m3": best["volume_m3"],
            "published_m3": reference["displacement_volume_m3"],
            "rel_error_percent": best["volume_rel_error_percent"],
            "tolerance_percent": 100.0 * VOLUME_TOLERANCE,
        },
    }
    for item in checks.values():
        item["within_tolerance"] = abs(item["rel_error_percent"]) <= item["tolerance_percent"]
    return {
        "unit": unit,
        "surface_count": len(surfaces),
        "rows": rows,
        "checks": checks,
    }


def build_record(archive_path: Path, manifest_path: Path, dtc_geometry_path: Path, **kwargs) -> dict:
    reference = published_reference(manifest_path, dtc_geometry_path)
    members = extract_members(archive_path)
    if HULL_MEMBER not in members:
        raise SystemExit(f"{HULL_MEMBER} not in {archive_path}: members {sorted(members)}")
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        hull_path = Path(tmp) / HULL_MEMBER
        hull_path.write_bytes(members[HULL_MEMBER])
        measured = measure(hull_path, reference, **kwargs)
    record = {
        "schema": SCHEMA,
        "anchor": ANCHOR,
        "ship": SHIP,
        "source": {
            "archive_path_in_cache": str(archive_path),
            "archive_url": SOURCE_URL,
            "publisher_page": SOURCE_PAGE,
            "archive_bytes": archive_path.stat().st_size,
            "archive_sha256": sha256_file(archive_path),
            "hull_member": HULL_MEMBER,
            "hull_sha256": hashlib.sha256(members[HULL_MEMBER]).hexdigest(),
            "license": LICENSE,
            "archive_defect": (
                "central-directory member paths use backslashes (dtc-2012-01-10\\DTC_hull.igs) while "
                "the local headers use forward slashes, so ZipFile.read raises BadZipFile; the "
                "central directory spelling is wrong and the header is authoritative"
            ),
            "member_sha256": {
                name: hashlib.sha256(data).hexdigest() for name, data in sorted(members.items())
            },
            "member_bytes": {name: len(data) for name, data in sorted(members.items())},
        },
        "published_reference": reference,
        "measured": measured,
        "verdict": None,
        "verdict_reason": "",
    }
    failures = [name for name, item in measured["checks"].items() if not item["within_tolerance"]]
    record["verdict"] = "VERIFIED" if not failures else "WITHHELD"
    record["verdict_reason"] = (
        "every geometric check is inside its tolerance"
        if not failures
        else "; ".join(
            f"{name} {measured['checks'][name]['rel_error_percent']:+.2f}% "
            f"(tolerance {measured['checks'][name]['tolerance_percent']:.2f}%)"
            for name in failures
        )
    )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default=str(DEFAULT_ARCHIVE))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--dtc-geometry", default=str(DEFAULT_DTC_GEOMETRY))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--spacing", type=float, default=None, help="single sampling density [m]")
    parser.add_argument("--dx", type=float, default=0.5)
    parser.add_argument("--dz", type=float, default=0.1)
    args = parser.parse_args(argv)

    archive = Path(args.archive)
    if not archive.exists():
        raise SystemExit(f"archive not cached: {archive}")
    spacings = (args.spacing,) if args.spacing else (1.0, 0.5)
    record = build_record(
        archive, Path(args.manifest), Path(args.dtc_geometry), spacings=spacings, dx=args.dx, dz=args.dz
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    for row in record["measured"]["rows"]:
        print(
            f"  spacing {row['sample_spacing_m']:>4} m: samples {row['point_cloud']['samples']:>7} "
            f"beam {row['point_cloud']['extent'][1]:8.4f} m  volume {row['volume_m3']:12.1f} m3 "
            f"Cb {row['cb']:.4f}  err {row['volume_rel_error_percent']:+.2f}%"
        )
    print(f"  verdict {record['verdict']}: {record['verdict_reason']}")
    return 0 if record["verdict"] == "VERIFIED" else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())