"""Build the baseline fleet from the *publishers' own* IGES bytes and judge it.

The fleet records in :mod:`tensorcad.fleet.manifest` describe hulls through a
form-family loft (`geometry_route="series60_form"`) or a parametric section law, and
their provenance honestly says ``retrieval="offline_recollection"`` /
``verified=false``: a proxy is a shape that closes on the published particulars, not
the published shape.  This module is the other route.  It takes the official IGES
archive of each hull, builds the watertight mesh below the design waterline, and
returns a verdict whose every number can be recomputed from the bytes it names:

* the **digest of the geometry bytes** (so a re-fetch can be compared, and a source
  change cannot pass silently);
* the **digest of the publisher page** the particulars were read from, for the ships
  whose particulars are page-parsed rather than recalled;
* the topology facts (watertight, edge-manifold, boundary edges, Euler characteristic,
  zero-area faces), which are properties of the mesh and cannot be "close enough";
* the **volume against the published displacement volume**, which is the anchor's <1 %
  acceptance, and the immersed block coefficient against the published ``Cb`` as an
  independent shape check -- a mesh can hit the volume with the wrong distribution,
  but it cannot hit both volume and ``Cb`` with the wrong beam or draft.

What this module deliberately does not do is paper over a refusal.  A patch set that
arrives *trimmed* cannot be integrated as its basis surfaces (the untrimmed patch
carries sheet the hull does not have), so a hull whose census shows trimmed patches is
built, measured, and reported as ``TRIM_UNHANDLED`` with the over-cover it causes --
the number is the evidence for the next fix, not a hidden caveat.

Units are read from each file's global section, not assumed: the T2015 archives mix
metres (KCS, ONRT, DTC) and millimetres (KVLCC2, JBC), and a millimetre file read as
metres is a 1000x hull that still passes every self-consistency gate.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
import zipfile
import dataclasses
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

from tensorcad.fleet_baseline import ShipSpec, build_ship

#: Where the publisher archives and their extracted members live.  Nothing here is
#: redistributed (the publishers state no licence): the archive is fetched from
#: ``source_url`` and the member is extracted on demand, so a run in a fresh sandbox
#: needs the archives in this cache and nothing else.
HULL_CACHE = Path(os.environ.get("TENSORCAD_HULL_CACHE", "/root/.cache/hull_sources/extracted"))

SCHEMA = "tensorcad.fleet.official_geometry.v1"

#: Anchor acceptance: volume against the published displacement volume.
VOLUME_TOLERANCE = 0.01
#: Independent shape check: immersed block coefficient against the published one.
CB_TOLERANCE = 0.01


@dataclass(frozen=True)
class OfficialShip:
    """One hull whose geometry bytes came from the publisher, with its citation."""

    name: str
    iges_path: str
    iges_sha256: str
    source_url: str
    archive_member: str
    publisher: str
    license: str
    draft: float
    length_pp: float
    beam: float
    published_volume_m3: float
    published_cb: float | None = None
    published_lcb_percent: float | None = None
    #: Where the published numbers came from, and how they were read.
    volume_basis: str = ""
    particulars_page: str = ""
    particulars_page_sha256: str = ""
    particulars_locator: str = ""
    #: Published model scale of the bytes (1.0 full scale) and why it is that value.
    geometry_scale: float = 1.0
    scale_basis: str = "geometry file is full scale (publisher archive)"
    notes: str = ""

    def spec(self) -> ShipSpec:
        return ShipSpec(
            name=self.name,
            draft=self.draft,
            length_pp=self.length_pp,
            beam=self.beam,
            published_volume_m3=self.published_volume_m3,
            published_cb=self.published_cb,
            published_lcb_percent=self.published_lcb_percent,
            published_source=self.particulars_page or self.source_url,
            source_url=self.source_url,
            publisher=self.publisher,
            license=self.license,
            iges_sha256=self.iges_sha256,
            iges_path=self.iges_path,
            notes=self.notes,
            geometry_scale=self.geometry_scale,
            scale_basis=self.scale_basis,
        )


#: The fleet.  Every address below is the publisher's own archive; every digest is the
#: digest of the bytes this module reads, so the provenance is checkable rather than
#: asserted.  Particulars marked ``page-parsed`` were read out of the publisher's HTML
#: by :func:`tensorcad.fleet.official.parse_particulars` (re-runnable on the same page
#: bytes); the rest are recalled literature values and say so.
SPECS: dict[str, OfficialShip] = {
    "KCS": OfficialShip(
        name="KCS",
        iges_path="{cache}/KCS_hull_SVA.igs",
        iges_sha256="6f5556b2b5e58e14e721eaf438fc2f7874e051c1a858d77d7aa837af17db917d",
        source_url="https://www.t2015.nmri.go.jp/file/Geometry/IGES_files/kcs/KCS_hull_SVA.zip",
        archive_member="KCS_hull_SVA/KCS_hull_SVA.igs",
        publisher="NMRI Tokyo 2015 Workshop (t2015.nmri.go.jp)",
        license="not stated by publisher; benchmark archive, bytes not redistributed",
        draft=10.8,
        length_pp=230.0,
        beam=32.2,
        published_volume_m3=52030.0,
        published_cb=0.6505,
        published_lcb_percent=-1.48,
        volume_basis="page-parsed (NMRI T2015 KCS main particulars, full-scale column)",
        particulars_page="https://www.t2015.nmri.go.jp/kcs_gc.html",
        particulars_page_sha256="68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
        particulars_locator="Main particulars table, full-scale column",
        notes="bare hull below design waterline; 55 untrimmed type-128 patches",
    ),
    "JBC": OfficialShip(
        name="JBC",
        iges_path="{cache}/jbc.igs",
        iges_sha256="",
        source_url="https://www.t2015.nmri.go.jp/file/Geometry/IGES_files/jbc/JBC_IGES.zip",
        archive_member="jbc.igs",
        publisher="NMRI Tokyo 2015 Workshop (t2015.nmri.go.jp)",
        license="not stated by publisher; benchmark archive, bytes not redistributed",
        draft=16.5,
        length_pp=280.0,
        beam=45.0,
        published_volume_m3=178369.9,
        published_cb=0.858,
        volume_basis="page-parsed (NMRI T2015 JBC main particulars, DESIGNED FULL column)",
        particulars_page="https://www.t2015.nmri.go.jp/jbc_gc.html",
        particulars_page_sha256="fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
        particulars_locator="Main particulars table, DESIGNED FULL column",
        notes="millimetre file (unit flag 2); 45 untrimmed type-128 patches",
    ),
    "ONRT": OfficialShip(
        name="ONRT",
        iges_path="{cache}/onr5613-full_appendage_deckhous_mergedR.igs",
        iges_sha256="",
        source_url="https://www.t2015.nmri.go.jp/file/Geometry/IGES_files/onrt/onrt_wo_gap.zip",
        archive_member="onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
        publisher="NMRI Tokyo 2015 Workshop (t2015.nmri.go.jp)",
        license="not stated by publisher; benchmark archive, bytes not redistributed",
        draft=5.494,
        length_pp=154.0,
        beam=18.78,
        # The page publishes displacement *mass* 8507 t and Cb 0.535, not a volume.
        # Cb * LWL * BWL * T closes to 8500.1 m^3, and 8507 t / 1.0255 t/m^3 is 8295 m^3,
        # so the workshop's 8507 t is 1.000 t/m^3 (fresh water) rather than seawater.
        # The gate is therefore run against the page's own Cb closure, and the mass
        # reading is carried alongside as a cross-check with its density stated.
        published_volume_m3=8500.1,
        published_cb=0.535,
        volume_basis="page-parsed Cb*LWL*BWL*T (page gives displacement mass 8507 t, no volume)",
        particulars_page="https://www.t2015.nmri.go.jp/onrt_gc.html",
        particulars_page_sha256="5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
        particulars_locator="Main particulars table, Full scale column",
        notes="full hull with appendages and deckhouse merged; hull-only sampling keeps z in [0, T]",
    ),
    "KVLCC2": OfficialShip(
        name="KVLCC2",
        iges_path="{cache}/kvlcc2.igs",
        iges_sha256="",
        source_url="https://simman2014.dk/wp-content/uploads/2015/07/KVLCC2-Hull.zip",
        archive_member="kvlcc2.igs",
        publisher="SIMMAN 2014 archive (simman2014.dk), MOERI tanker",
        license="published workshop benchmark data, research use",
        draft=20.8,
        length_pp=320.0,
        beam=58.0,
        published_volume_m3=312622.0,
        published_cb=0.8098,
        published_lcb_percent=3.48,
        volume_basis="literature (SIMMAN KVLCC2 particulars; not page-parsed in this run)",
        notes="millimetre file; 19 type-128 patches, ALL bound to type-141 trim boundaries",
    ),
    "Series60": OfficialShip(
        name="Series60",
        iges_path="",
        iges_sha256="",
        source_url="",
        archive_member="",
        publisher="Todd (1963) Series 60 systematic series (DTMB report 1712)",
        license="methodical-series offset tables; no IGES is published",
        draft=5.33333,
        length_pp=100.0,
        beam=13.33333,
        published_volume_m3=0.6 * 100.0 * 13.33333 * 5.33333,
        published_cb=0.6,
        volume_basis="closure Cb*L*B*T of the series nominal point (Cb 0.60)",
        notes="a methodical series has offset tables, not a hull file: this entry has no official IGES by construction",
    ),
}


def materialise(ship: OfficialShip, cache: Path | None = None) -> Path | None:
    """Return the hull file, extracting it from a cached publisher archive if needed.

    The lookup order is: the extracted member under :data:`HULL_CACHE`, then any
    ``*.zip`` in the cache whose central directory contains a member with that
    basename (the archives are stored under their fetch-time names, which are not
    canonical).  Nothing is downloaded here; a missing archive is reported as absent
    so that a run without the cache says so instead of inventing geometry.
    """
    if not ship.iges_path:
        return None
    root = Path(cache or HULL_CACHE)
    target = Path(ship.iges_path.replace("{cache}", str(root)))
    if target.is_file():
        return target
    base = Path(ship.archive_member.replace("\\", "/")).name or target.name
    for archive in sorted(root.parent.glob("*.zip")) + sorted(root.glob("*.zip")):
        try:
            with zipfile.ZipFile(archive) as zf:
                member = next(
                    (n for n in zf.namelist() if Path(n.replace("\\", "/")).name == base), None
                )
                if member is None:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(zf.read(member))
                return target
        except (zipfile.BadZipFile, OSError):
            continue
    return None


def _euler_characteristic(mesh: Any) -> int:
    edges: set[tuple[int, int]] = set()
    for face in mesh.faces:
        a, b, c = face
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((u, v) if u < v else (v, u))
    return len(mesh.vertices) - len(edges) + len(mesh.faces)


def _rel_error(got: float, ref: float) -> float:
    return (got - ref) / ref if ref else float("nan")


def verify_ship(
    ship: OfficialShip,
    *,
    dx: float | None = None,
    dz: float | None = None,
) -> dict[str, Any]:
    """Build one hull from its official bytes and judge the mesh it produces."""
    record: dict[str, Any] = {"ship": ship.name, "schema": SCHEMA}
    record["provenance"] = {
        "source_url": ship.source_url,
        "archive_member": ship.archive_member,
        "publisher": ship.publisher,
        "license": ship.license,
        "geometry_sha256_declared": ship.iges_sha256,
        "particulars_page": ship.particulars_page,
        "particulars_page_sha256": ship.particulars_page_sha256,
        "particulars_locator": ship.particulars_locator,
        "volume_basis": ship.volume_basis,
        "scale_basis": ship.scale_basis,
        "notes": ship.notes,
    }
    path = materialise(ship)
    if path is None or not path.is_file():
        record["verdict"] = "NO_OFFICIAL_GEOMETRY"
        record["detail"] = (
            "no official hull file exists for this entry"
            if not ship.iges_path
            else f"geometry file absent: {ship.iges_path}"
        )
        return record

    blob = path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    record["provenance"]["geometry_sha256_measured"] = digest
    record["provenance"]["geometry_bytes"] = len(blob)
    if ship.iges_sha256 and digest != ship.iges_sha256:
        record["verdict"] = "DIGEST_MISMATCH"
        record["detail"] = f"declared {ship.iges_sha256} != measured {digest}"
        return record

    spec = ship.spec()
    # the citation keeps the placeholder-free archive address; the reader needs the
    # materialised path, which is the same file under the cache root.
    spec = dataclasses.replace(spec, iges_path=str(path))
    started = time.time()
    try:
        kwargs: dict[str, float] = {}
        if dx is not None:
            kwargs["dx"] = dx
        if dz is not None:
            kwargs["dz"] = dz
        mesh, info = build_ship(spec, **kwargs)  # type: ignore[arg-type]
    except Exception as exc:  # noqa: BLE001 - a refusal is a result, not a crash
        record["verdict"] = "WITHHELD"
        record["detail"] = f"{type(exc).__name__}: {exc}"
        record["seconds"] = round(time.time() - started, 1)
        return record
    record["seconds"] = round(time.time() - started, 1)

    stats = info["mesh"]
    derived = info.get("derived", {})
    census = info.get("census", {})
    volume = float(stats["volume_m3"])
    rel_volume = _rel_error(volume, ship.published_volume_m3)
    cb_measured = derived.get("cb")
    rel_cb = (
        _rel_error(float(cb_measured), ship.published_cb)
        if cb_measured is not None and ship.published_cb
        else None
    )
    trim_total = int(census.get("surfaces_type128", 0) or 0)
    trim_kept = int(census.get("basis_surfaces_untrimmed", 0) or 0)
    trim_dropped = trim_total - trim_kept

    record["census"] = census
    record["mesh"] = {
        "volume_m3": volume,
        "field_volume_m3": info.get("integrals", {}).get("volume_m3"),
        "triangle_count": stats["triangle_count"],
        "vertex_count": stats["vertex_count"],
        "is_watertight": bool(stats["is_watertight"]),
        "is_edge_manifold": bool(stats["is_edge_manifold"]),
        "boundary_edges": int(stats["boundary_edges"]),
        "zero_area_faces": int(stats["zero_area_faces"]),
        "euler_characteristic": _euler_characteristic(mesh),
        "extent_m": stats["extent_m"],
    }
    record["hydrostatics"] = {
        "cb_measured": cb_measured,
        "cb_published": ship.published_cb,
        "cb_rel_error_percent": None if rel_cb is None else round(100.0 * rel_cb, 3),
        "lcb_percent_measured": derived.get("lcb_percent_of_lpp"),
        "lcb_percent_published": ship.published_lcb_percent,
        "cm_measured": derived.get("cm"),
        "cwp_measured": derived.get("cwp"),
    }
    record["volume_vs_published"] = {
        "volume_m3": volume,
        "published_volume_m3": ship.published_volume_m3,
        "rel_error_percent": round(100.0 * rel_volume, 3),
        "tolerance_percent": 100.0 * VOLUME_TOLERANCE,
    }

    checks = {
        "watertight": bool(stats["is_watertight"]),
        "edge_manifold": bool(stats["is_edge_manifold"]),
        "no_boundary_edges": int(stats["boundary_edges"]) == 0,
        "no_zero_area_faces": int(stats["zero_area_faces"]) == 0,
        "genus_zero": _euler_characteristic(mesh) == 2,
        "volume_within_tolerance": abs(rel_volume) <= VOLUME_TOLERANCE,
        "cb_within_tolerance": (rel_cb is not None and abs(rel_cb) <= CB_TOLERANCE),
    }
    record["checks"] = checks

    if trim_dropped:
        record["verdict"] = "TRIM_UNHANDLED"
        record["detail"] = (
            f"{trim_dropped}/{trim_total} basis surfaces are bound to a trim boundary; "
            "they were integrated as whole basis surfaces, so the mesh over-covers the "
            "hull. The volume error below is the size of that over-cover, not mesh noise."
        )
    elif all(checks.values()):
        record["verdict"] = "VERIFIED"
        record["detail"] = "official bytes -> watertight genus-0 mesh within volume and Cb tolerance"
    else:
        failed = [name for name, ok in checks.items() if not ok]
        record["verdict"] = "FAILED"
        record["detail"] = "failed checks: " + ", ".join(failed)
    return record


def verify_fleet(
    out_dir: str | Path,
    *,
    only: Sequence[str] | None = None,
    specs: dict[str, OfficialShip] | None = None,
) -> dict[str, Any]:
    """Verify every hull, write ``official_geometry.json``, return the payload."""
    specs = specs or SPECS
    names = list(only) if only else list(specs)
    records = [verify_ship(specs[name]) for name in names]
    payload = {
        "schema": SCHEMA,
        "generated_by": "tensorcad.fleet.official_geometry",
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "tolerances": {"volume_percent": 100.0 * VOLUME_TOLERANCE, "cb_percent": 100.0 * CB_TOLERANCE},
        "verdict_counts": _counts(records),
        "ships": records,
    }
    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    (target / "official_geometry.json").write_text(json.dumps(payload, indent=1) + "\n")
    return payload


def _counts(records: Sequence[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        verdict = str(record.get("verdict", "UNKNOWN"))
        counts[verdict] = counts.get(verdict, 0) + 1
    return dict(sorted(counts.items()))


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.fleet.official_geometry",
        description="Build the baseline fleet from the publishers' own IGES bytes and gate it.",
    )
    parser.add_argument("--out", default="artifacts/baseline_fleet/official")
    parser.add_argument("--only", nargs="*", default=None, help="ship names (default: all)")
    parser.add_argument("--dx", type=float, default=None)
    parser.add_argument("--dz", type=float, default=None)
    args = parser.parse_args(argv)
    payload = verify_fleet(args.out, only=args.only, )
    compact = {record["ship"]: record.get("verdict") for record in payload["ships"]}
    print(json.dumps({"counts": payload["verdict_counts"], "verdicts": compact}, indent=1))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())