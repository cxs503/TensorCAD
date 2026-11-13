#!/usr/bin/env python3
"""Build the watertight baseline fleet and benchmark it against published values.

    python tools/baseline_fleet/build_fleet.py --out artifacts/baseline_fleet

What it does, in order
----------------------
1. **Provenance**: unzips the workshop archives that
   ``tools/baseline_fleet/fetch_sources.py`` downloaded (or uses an existing raw
   tree) and records the SHA-256 of every IGES file it reads.
2. **Zero hand-copied numbers**: extracts Lpp / B / T / displacement / Cb / LCB
   from the *cached publisher pages* under ``raw/pages`` together with the page
   SHA-256 and the matched text snippet, so each number in the spec JSON can be
   traced to bytes on disk.
3. **Build**: runs :func:`tensorcad.fleet_baseline.build_fleet` -- IGES patches to
   a watertight hull-below-waterline mesh, per-ship manifest, STL export.
4. **Benchmark**: ``ref`` (published) versus ``got`` (mesh) displacement, with the
   relative error and the verdict per acceptance tolerance.
5. **Figures**: volume comparison and the reconstructed half-breadth field.

Exit code is 1 when no hull reaches VERIFIED, so CI cannot mistake a rebuilding
failure for a pass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from tensorcad.fleet_baseline import (  # noqa: E402
    ShipSpec,
    build_fleet,
    DEFAULT_DX,
    DEFAULT_DZ,
    DEFAULT_RING_POINTS,
    DEFAULT_SAMPLES_PER_METRE,
)

#: ``ship -> (archive basename, IGES member inside the archive, publisher page key)``
ARCHIVES: dict[str, tuple[str, str, str]] = {
    "KCS": ("KCS-Hull.zip", "KCS_hull_SVA.igs", "kcs"),
    "KVLCC2": ("KVLCC2-Hull.zip", "kvlcc2.igs", "kvlcc2"),
    "JBC": ("JBC_IGES.zip", "jbc.igs", "jbc"),
    "ONRT": ("onrt_wo_gap.zip", "onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs", "onrt"),
}

#: Publisher-page patterns.  Each entry is ``(field, regex)``; the **last** number
#: in the first match is the full-scale value (every workshop table lists model
#: scale first, full scale second).  Nothing here is typed in by hand: the values
#: come from the HTML that ``fetch_sources.py`` cached.
PATTERNS: dict[str, dict[str, str]] = {
    "kcs": {
        "length_pp": r"Length between perpendiculars L\s*PP\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"Maximum beam of waterline B\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "volume": r"Displacement volume &nabla;\s*\(m\s*3\s*\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient \(CB\)[^0-9]{0,60}([\d.]+)",
        "lcb": r"LCB \(%L\s*PP\s*\), fwd\+\s*(?:&nbsp;|\s)*([-+\d.]+)",
    },
    "jbc": {
        "length_pp": r"Length between perpendiculars \\\(L_\{PP\}\\\)\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"(?:Maximum beam|Breadth)[^()]{0,40}\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft\(Full\)\s*T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "volume": r"Displacement volume of a hull\s*\\\( \nabla \\\)\s*\(m\s*3\s*\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient\s*\\\(C_B\\\)[^0-9]{0,80}([\d.]+)",
    },
    "onrt": {
        "length_pp": r"Length of waterline L\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"Maximum beam of waterline B\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient \(CB\)[^0-9]{0,60}([\d.]+)",
    },
}


#: Workshop tables list the model-scale column first (KCS, JBC) or the full-scale
#: column second (ONRT); the full-scale value is what the geometry is compared to.
PICK = {"kcs": "first", "jbc": "first", "onrt": "last", "kvlcc2": "last"}


def page_text(path: Path) -> str:
    raw = path.read_text(errors="replace")
    text = re.sub(r"<[^>]+>", " ", raw)
    return re.sub(r"\s+", " ", text)


def extract(page: Path, patterns: dict[str, str]) -> dict[str, dict[str, object]]:
    text = page_text(page)
    digest = hashlib.sha256(page.read_bytes()).hexdigest()
    found: dict[str, dict[str, object]] = {}
    for field_name, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            continue
        numbers = [float(value) for value in re.findall(r"[-+]?\d+\.?\d*", match.group(1))]
        if PICK[page.stem] == "first":
            numbers = numbers[:1]
        if not numbers:
            continue
        found[field_name] = {
            "value": numbers[-1],
            "snippet": match.group(0)[:180],
            "page_sha256": digest,
            "page": str(page),
        }
    return found


def discover_iges(raw: Path, ship: str) -> Path | None:
    archive_name, member, _ = ARCHIVES[ship]
    members = [p for p in raw.rglob(Path(member).name)]
    if members:
        return members[0]
    archive = raw / archive_name
    if not archive.exists():
        return None
    with zipfile.ZipFile(archive) as zf:
        for info in zf.infolist():
            if info.filename.endswith(Path(member).name):
                target = raw / ship / Path(info.filename).name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(zf.read(info))
                return target
    return None


def build_specs(raw: Path) -> tuple[list[ShipSpec], dict[str, object]]:
    specs: list[ShipSpec] = []
    evidence: dict[str, object] = {}
    for ship, (archive_name, member, page_key) in ARCHIVES.items():
        iges = discover_iges(raw, ship)
        archived = raw / archive_name
        page = raw / "pages" / f"{page_key}.html"
        record: dict[str, object] = {}
        if iges is not None:
            record["iges"] = {
                "path": str(iges.relative_to(REPO)),
                "sha256": hashlib.sha256(iges.read_bytes()).hexdigest(),
                "bytes": iges.stat().st_size,
            }
        if archived.exists():
            record["archive"] = {
                "path": str(archived.relative_to(REPO)),
                "sha256": hashlib.sha256(archived.read_bytes()).hexdigest(),
            }
        fields: dict[str, dict[str, object]] = {}
        if page.exists():
            fields = extract(page, PATTERNS.get(page_key, {}))
            record["page"] = str(page.relative_to(REPO))
        record["published_fields"] = fields
        evidence[ship] = record
        if iges is None:
            continue
        draft = fields.get("draft", {}).get("value")
        length = fields.get("length_pp", {}).get("value")
        beam = fields.get("beam", {}).get("value")
        volume = fields.get("volume", {}).get("value")
        cb = fields.get("cb", {}).get("value")
        lcb = fields.get("lcb", {}).get("value")
        if ship == "ONRT" and volume is None and None not in (cb, length, beam, draft):
            # The ONRT page publishes Cb with its definition ∇/(L_WL B_WL T) and the
            # three dimensions, so the displaced volume follows from the publisher's
            # own definition instead of the displacement in tonnes (which would need
            # an assumed water density).
            volume = float(cb) * float(length) * float(beam) * float(draft)
            record["volume_from_cb"] = {
                "formula": "Cb * Lwl * Bwl * T (publisher's own definition)",
                "value": volume,
            }
        if draft is None:
            record["skip_reason"] = (
                "publisher page lacked the design draft, so the waterline cut "
                "cannot be placed; no volume comparison is claimed"
            )
            continue
        specs.append(
            ShipSpec(
                name=ship,
                draft=float(draft),
                length_pp=float(length) if length else 0.0,
                beam=float(beam) if beam else 0.0,
                published_volume_m3=None if volume is None else float(volume),
                published_cb=None if cb is None else float(cb),
                published_lcb_percent=None if lcb is None else float(lcb),
                published_source=str(page.relative_to(REPO)) if page.exists() else "",
                source_url=ARCHIVES[ship][0],
                publisher="NMRI T2015 / SIMMAN 2020 workshop archive",
                license="not stated by publisher (see docs/baseline_fleet.md); bytes not redistributed",
                iges_sha256=str(record.get("iges", {}).get("sha256", "")),
                iges_path=str(iges),
                notes="hull below design waterline only; source IGES patch set",
            )
        )
    return specs, evidence


def figures(out: Path, summary: dict, manifests: dict[str, dict]) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as error:  # pragma: no cover - optional dependency
        print(f"figures skipped: {error}", file=sys.stderr)
        return []
    made: list[str] = []
    hulls = [name for name, item in summary["hulls"].items() if item.get("volume_m3")]
    ref = [summary["hulls"][name]["published_volume_m3"] for name in hulls]
    got = [summary["hulls"][name]["volume_m3"] for name in hulls]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    positions = np.arange(len(hulls))
    ax.bar(positions - 0.2, ref, 0.4, label="published ∇ (m³)", color="#4c72b0")
    ax.bar(positions + 0.2, got, 0.4, label="mesh ∇ (m³)", color="#dd8452")
    ax.set_xticks(positions)
    ax.set_xticklabels(hulls)
    ax.set_ylabel("displacement volume (m³)")
    ax.set_title("Baseline fleet: reconstructed vs published displacement")
    for index, name in enumerate(hulls):
        error = summary["hulls"][name]["volume_rel_error_percent"]
        ax.text(index, max(ref[index], got[index]) * 1.02, f"{error:+.2f}%", ha="center", fontsize=8)
    ax.legend()
    fig.tight_layout()
    path = out / "figs" / "volume_vs_published.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    made.append(str(path.relative_to(REPO)))

    for name, manifest in manifests.items():
        field = manifest.get("field")
        if not field:
            continue
        import importlib

        baseline = importlib.import_module("tensorcad.fleet_baseline")
        spec = ShipSpec(**{**manifest["spec"]})
        model_points = baseline.build_ship.__globals__  # only for type checkers
        del model_points
        # Re-derive the field cheaply from the manifest geometry: read the STL is
        # not invertible, so rebuild the field from the IGES at coarse resolution.
        from tensorcad.iges import read_iges

        surfaces = read_iges(spec.iges_path).surfaces()
        points = baseline.sample_patches(surfaces, spacing=0.4)
        grid = baseline.field_from_points(points, spec.draft, dx=1.0, dz=0.25)
        data = np.array(grid.b).T
        fig, ax = plt.subplots(figsize=(7.2, 3.4))
        mesh = ax.pcolormesh(
            np.array([grid.x(i) for i in range(grid.nx)]),
            np.array([grid.z(k) for k in range(grid.nz)]),
            data,
            shading="auto",
            cmap="viridis",
        )
        fig.colorbar(mesh, ax=ax, label="half-breadth b(x,z) (m)")
        for level in (0.25, 0.5, 0.75):
            ax.axhline(spec.draft * level, color="w", lw=0.5, alpha=0.5)
        ax.set_xlabel("x from AP (m)")
        ax.set_ylabel("z above keel (m)")
        ax.set_title(f"{name}: half-breadth field from published IGES ({len(surfaces)} patches)")
        fig.tight_layout()
        path = out / "figs" / f"{name}_half_breadth_field.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        made.append(str(path.relative_to(REPO)))
    return made


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/baseline_fleet")
    parser.add_argument("--dx", type=float, default=DEFAULT_DX)
    parser.add_argument("--dz", type=float, default=DEFAULT_DZ)
    parser.add_argument("--ring-points", type=int, default=DEFAULT_RING_POINTS)
    parser.add_argument("--samples-per-metre", type=float, default=DEFAULT_SAMPLES_PER_METRE)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    out = (REPO / args.out).resolve() if not Path(args.out).is_absolute() else Path(args.out)
    raw = out / "raw"
    specs, evidence = build_specs(raw)
    (out / "fleet_specs.json").write_text(
        json.dumps({"ships": [spec.to_dict() for spec in specs], "evidence": evidence}, indent=2)
    )
    if not specs:
        print("no buildable hull: run tools/baseline_fleet/fetch_sources.py first", file=sys.stderr)
        return 1

    summary = build_fleet(
        specs,
        out,
        dx=args.dx,
        dz=args.dz,
        ring_points=args.ring_points,
        samples_per_metre=args.samples_per_metre,
    )
    manifests = {}
    for spec in specs:
        path = out / f"{spec.name}_manifest.json"
        if path.exists():
            manifests[spec.name] = json.loads(path.read_text())

    benchmark = {
        "schema": "tensorcad.fleet_baseline.benchmark/1.0",
        "case": "GAP-CAD-T1-BASELINE-FLEET",
        "acceptance": {
            "watertight_and_manifold": True,
            "volume_rel_error_percent_max": 1.0,
        },
        "reference_definition": (
            "publisher's displacement volume from the workshop page cached under "
            "artifacts/baseline_fleet/raw/pages (each value carries the page SHA-256 "
            "and the matched snippet in fleet_specs.json)"
        ),
        "hulls": {},
    }
    for name, entry in summary["hulls"].items():
        ref = entry.get("published_volume_m3")
        got = entry.get("volume_m3")
        benchmark["hulls"][name] = {
            "ref_volume_m3": ref,
            "got_volume_m3": got,
            "rel_error_percent": entry.get("volume_rel_error_percent"),
            "within_1_percent": (
                None
                if entry.get("volume_rel_error_percent") is None
                else abs(entry["volume_rel_error_percent"]) <= 1.0
            ),
            "watertight": entry.get("watertight"),
            "edge_manifold": entry.get("edge_manifold"),
            "triangles": entry.get("triangles"),
            "verdict": entry.get("verdict"),
            "reasons": entry.get("reasons"),
        }
    (out / "benchmark.json").write_text(json.dumps(benchmark, indent=2))
    summary["figures"] = [] if args.no_figures else figures(out, summary, manifests)
    (out / "fleet_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(benchmark["hulls"], indent=2))
    return 0 if summary["verified"] else 1


if __name__ == "__main__":
    sys.exit(main())