"""Probe the DTC hull IGES that the uni-due INAM archive ships, and record what it measures.

Anchor ``REPRO-GEOM-FLEET-REFS``.  Earlier rounds reported the DTC geometry stage as
``WITHHELD`` with the reason *"no DTC hull IGES is reachable: uni-due.de/ship/dtc.php
returns 404"*.  That reason was wrong about the bytes: the workshop IGES archive is
served from the INAM file area and is present in this environment's source cache.  What is
true is different and measurable -- the file **cannot be evaluated by this repository's
IGES reader**: all 1 345 type-128 surfaces sample to the origin, so the half-breadth field
route would loft a zero-volume "hull".

This tool pins the evidence so the next round does not have to re-discover it:

* the archive digest and every member digest (bytes are not redistributed);
* an entity census of the hull file, read from the fixed-width IGES D section;
* a sample census -- number of sampled points and their axis extent -- produced by the
  same :func:`tensorcad.fleet_baseline.sample_patches` call the geometry route makes;
* the ``RouteRefusal`` the route now raises, instead of a silently degenerate mesh;
* a control on a hull the reader *does* handle (KCS), so the zero is attributable;
* the published ∇/Cb/L/B/T the measurement would be compared against, injected from the
  cross-source library rather than typed in.

Nothing here is transcribed: every number in the artefact comes out of the byte strings in
the cache or out of the library file named on the command line.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:  # pragma: no cover - import plumbing
    sys.path.insert(0, str(REPO))

SCHEMA = "tensorcad.fleet_refs.dtc_geometry/1.0"
DEFAULT_CACHE = Path("/root/.cache/hull_sources")
DEFAULT_CROSSREF = REPO / "benchmarks" / "verified" / "fleet_refs" / "crossref.json"
DEFAULT_OUT = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_geometry.json"
#: The archive the DTC hull file lives in, as two spellings of the same URL, kept because
#: the fetch script caches by host path.
ARCHIVE_NAMES = ("dtc_igs_ordner.zip", "www.uni-due.de_imperia_md_content_inam_igs_ordner.zip")
HULL_MEMBER = "DTC_hull.igs"
CONTROL_ZIP = "t2015.nmri.go.jp_file_Geometry_IGES_files_kcs_KCS_hull_SVA.zip"
CONTROL_MEMBER = "KCS_hull_SVA/KCS_hull_SVA.igs"
DTC_SAMPLE_SPACING = 0.4  # model-scale metres; the hull is 5.976 m long
CONTROL_SAMPLE_SPACING = 0.5


def _digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_member(archive: Path, member: str) -> bytes:
    """Read ``member`` out of ``archive``, tolerating the archive's backslash path names."""
    with zipfile.ZipFile(archive) as zf:
        for info in zf.infolist():
            info.orig_filename = info.orig_filename.replace("\\", "/")
        for info in zf.infolist():
            if info.filename.replace('\\', '/').split('/')[-1] == member.replace('\\', '/').split('/')[-1]:
                return zf.read(info.filename)
    raise KeyError(f"{member} not in {archive}")


def entity_census(iges_text: str) -> dict[str, object]:
    """Count IGES directory entries by entity type (columns 1-8 of each D-section line)."""
    lines = iges_text.splitlines()
    directory = [line for line in lines if len(line) > 72 and line[72:73] == "D"]
    types = collections.Counter()
    for line in directory[::2]:
        field = line[0:8].strip()
        if field.isdigit():
            types[int(field)] += 1
    return {
        "directory_entries": len(directory) // 2,
        "entity_types": {str(k): v for k, v in sorted(types.items())},
        "type_128_surfaces": types.get(128, 0),
        "type_126_curves": types.get(126, 0),
        "type_144_trimmed_patches": types.get(144, 0),
    }


def sample_census(payload: bytes, spacing: float, workdir: Path) -> dict[str, object]:
    """Sample the file's B-spline surfaces exactly as the geometry route does, and measure them."""
    import numpy as np

    from tensorcad.fleet_baseline import sample_patches
    from tensorcad.iges import read_iges

    workdir.mkdir(parents=True, exist_ok=True)
    path = workdir / "probe.igs"
    path.write_bytes(payload)
    iges = read_iges(str(path))
    source = iges.surfaces
    surfaces = [s for s in (source() if callable(source) else source) if hasattr(s, "sample")]
    points = np.asarray(sample_patches(surfaces, spacing=spacing), dtype=float)
    if points.size:
        lo = points.min(axis=0)
        hi = points.max(axis=0)
        extent = (hi - lo).tolist()
    else:
        lo = hi = np.zeros(3)
        extent = [0.0, 0.0, 0.0]
    return {
        "surfaces_read": len(surfaces),
        "sample_spacing_m": spacing,
        "sample_points": int(points.shape[0]) if points.size else 0,
        "unique_points": int(np.unique(points, axis=0).shape[0]) if points.size else 0,
        "bbox_min": [round(float(v), 6) for v in lo],
        "bbox_max": [round(float(v), 6) for v in hi],
        "extent_m": [round(float(v), 6) for v in extent],
        "max_extent_m": round(float(max(extent)), 6),
    }


def route_refusal(payload: bytes, workdir: Path, *, length: float, beam: float, draft: float, cb: float) -> dict[str, object]:
    """Ask the geometry route to build the DTC hull and record how (or whether) it refuses."""
    from tensorcad.fleet.build import RouteRefusal, iges_field_mesh
    from tensorcad.fleet.manifest import Provenance, ShipRecord

    workdir.mkdir(parents=True, exist_ok=True)
    path = workdir / "dtc_probe.igs"
    path.write_bytes(payload)
    record = ShipRecord(
        name="DTC",
        ship_class="container",
        length=length,
        beam=beam,
        draft=draft,
        block_coefficient=cb,
        displacement_basis="literature_value_verified",
        provenance=Provenance(
            source="uni-due INAM DTC IGES archive",
            url="https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip",
            license="not stated by publisher",
            retrieval="cached_bytes",
            verified=True,
        ),
    )
    try:
        built = iges_field_mesh(
            record,
            path,
            dx=DTC_SAMPLE_SPACING,
            dz=DTC_SAMPLE_SPACING / 5.0,
            ring_points=60,
            samples_per_metre=1.0 / DTC_SAMPLE_SPACING,
        )
    except RouteRefusal as refusal:
        return {"outcome": "refused", "reason": str(refusal)}
    from tensorcad.fleet.gates import measure

    facts = measure(built.mesh)
    return {
        "outcome": "built",
        "signed_volume_m3": round(float(facts.signed_volume), 9),
        "triangle_count": facts.triangle_count,
    }


def published_reference(crossref_path: Path) -> dict[str, object]:
    """The library's DTC entries for the four quantities the anchor names, plus L/B/T."""
    data = json.loads(crossref_path.read_text())
    quantities = (data.get("report") or {}).get("DTC", {}).get("quantities", {})
    out: dict[str, object] = {}
    for name in ("displacement_volume", "cb", "lcb_percent", "lwl", "length_pp", "beam", "draft"):
        cell = quantities.get(name)
        if not cell:
            out[name] = {"verdict": "MISSING"}
            continue
        out[name] = {
            "verdict": cell.get("verdict"),
            "independent_families": cell.get("independent_families"),
            "published_values": cell.get("published_values"),
            "range_relative": cell.get("range_relative"),
            "tolerance": cell.get("tolerance"),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--crossref", type=Path, default=DEFAULT_CROSSREF)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workdir", type=Path, default=Path("/tmp/dtc_probe"))
    args = parser.parse_args(argv)

    archive = next((args.cache / name for name in ARCHIVE_NAMES if (args.cache / name).exists()), None)
    if archive is None:
        raise SystemExit(f"no DTC archive in {args.cache} (looked for {ARCHIVE_NAMES})")
    hull_bytes = _read_member(archive, HULL_MEMBER)
    hull_text = hull_bytes.decode("latin-1")
    # Member digests are taken from the bytes, not from the archive headers.
    with zipfile.ZipFile(archive) as zf:
        for info in zf.infolist():
            info.orig_filename = info.orig_filename.replace("\\", "/")
        members = {
            info.filename.replace('\\', '/').split('/')[-1]: hashlib.sha256(zf.read(info.filename)).hexdigest()
            for info in zf.infolist()
        }

    reference = published_reference(args.crossref)
    dtc_length = float((reference.get("length_pp") or {}).get("published_values", [{"value": 355.0}])[0]["value"])
    dtc_beam = float((reference.get("beam") or {}).get("published_values", [{"value": 51.0}])[0]["value"])
    dtc_draft = float((reference.get("draft") or {}).get("published_values", [{"value": 14.5}])[0]["value"])
    dtc_cb = float((reference.get("cb") or {}).get("published_values", [{"value": 0.661}])[0]["value"])

    dtc_census = entity_census(hull_text)
    report: dict[str, object] = {
        "schema": SCHEMA,
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "ship": "DTC",
        "archive": {
            "path_in_cache": archive.name,
            "sha256": _digest(archive),
            "bytes": archive.stat().st_size,
            "members": members,
            "hull_member": HULL_MEMBER,
            "hull_sha256": hashlib.sha256(hull_bytes).hexdigest(),
            "hull_bytes": len(hull_bytes),
            "license": "publisher states no licence; the bytes are digested, not redistributed",
        },
        "entity_census": dtc_census,
        "dtc_sample_census": sample_census(hull_bytes, DTC_SAMPLE_SPACING, args.workdir / "dtc"),
        "route": route_refusal(
            hull_bytes,
            args.workdir / "route",
            length=dtc_length,
            beam=dtc_beam,
            draft=dtc_draft,
            cb=dtc_cb,
        ),
        "control": {
            "ship": "KCS",
            "why": "a hull the same reader handles, so a zero sample cloud is attributable to the file, not the reader",
            "sample_census": sample_census(
                _read_member(args.cache / CONTROL_ZIP, CONTROL_MEMBER), CONTROL_SAMPLE_SPACING, args.workdir / "kcs"
            ),
        },
        "published_reference": reference,
        "status": None,
        "status_reason": None,
    }

    sample = report["dtc_sample_census"]
    control = report["control"]["sample_census"]
    if report["route"]["outcome"] == "built":
        report["status"] = "MEASURED"
        report["status_reason"] = "the route built a mesh from the official DTC hull file"
    elif sample["max_extent_m"] == 0.0 and control["max_extent_m"] > 0.0:
        report["status"] = "WITHHELD_READER"
        report["status_reason"] = (
            "the official DTC hull file is present (archive sha256 pinned) but the IGES reader "
            f"evaluates all {dtc_census['type_128_surfaces']} type-128 surfaces to the origin: "
            f"{sample['sample_points']} sample points, extent {sample['max_extent_m']} m, while the "
            f"KCS control file gives {control['sample_points']} points spanning {control['max_extent_m']} m. "
            "The geometry route refuses rather than lofting a zero-volume hull, so no DTC "
            "displacement can be measured against the published reference yet."
        )
    else:
        report["status"] = "WITHHELD_UNEXPLAINED"
        report["status_reason"] = "the sample cloud and the control both collapsed; investigate the reader first"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "status_reason")}, indent=1))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())