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
import html
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
#: ``ship -> (archive basename, IGES member inside the archive, page keys)``.
#: Page keys are candidate publisher pages in preference order; the first one that is
#: present in the cache supplies the reference numbers.  ``simman_summary`` is the
#: SIMMAN 2020 "Ship Data" page whose Tables 1/2/3 carry the full-scale principal
#: particulars of KVLCC2 / KCS / ONRT, including the model scale and the displacement.
ARCHIVES: dict[str, tuple[str, str, tuple[str, ...]]] = {
    "KCS": ("KCS-Hull.zip", "KCS_hull_SVA.igs", ("t2015_kcs", "simman_summary")),
    "KVLCC2": ("KVLCC2-Hull.zip", "kvlcc2.igs", ("simman_summary",)),
    "JBC": ("JBC_IGES.zip", "jbc.igs", ("t2015_jbc",)),
    "ONRT": (
        "onrt_wo_gap.zip",
        "onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
        ("t2015_onrt", "simman_summary"),
    ),
}

#: Which model-space -> full-scale policy a hull's geometry file needs.
#: ``("full", None)``          the bytes are already full scale;
#: ``("inferred", None)``      no published model scale: normalise to the published L_PP
#:                             (recorded in the manifest; the volume gate becomes a
#:                             block-coefficient gate, see VolumeGateBasis);
#: ``("published", key)``      read the published model scale from that page's table.
GEOMETRY_SCALE_POLICY: dict[str, tuple[str, str | None]] = {
    "KCS": ("full", None),
    "KVLCC2": ("full", None),
    "JBC": ("inferred", None),
    "ONRT": ("published", "simman_summary"),
}

#: Publisher-page patterns.  Each entry is ``(field, regex)``; the **last** number in
#: the first match is the full-scale value on the T2015 pages (those tables list model
#: scale first, full scale second), and the **first** on the SIMMAN summary page (full
#: scale first).  Nothing here is typed in by hand: the values come from the HTML that
#: ``fetch_sources.py`` cached, and ``fleet_specs.json`` stores page sha256 + snippet.
PAGE_PATTERNS: dict[str, dict[str, str]] = {
    "t2015_kcs": {
        "length_pp": r"Length between perpendiculars L\s*PP\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"Maximum beam of waterline B\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "volume": r"Displacement volume\s*(?:∇|&nabla;)\s*\(m\s*3\s*\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient \(CB\)[^0-9]{0,60}([\d.]+)",
        "lcb": r"LCB \(%L\s*PP\s*\), fwd\+\s*(?:&nbsp;|\s)*([-+\d.]+)",
    },
    "t2015_jbc": {
        "length_pp": r"Length between perpendiculars[^0-9]{0,40}\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"Maximum beam[^0-9]{0,40}\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft\(Full\)\s*T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "volume": r"Displacement volume of a hull[^0-9]{0,40}\(m\s*3\s*\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient[^0-9]{0,80}([\d.]+)",
    },
    "t2015_onrt": {
        "length_pp": r"Length of waterline L\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "beam": r"Maximum beam of waterline B\s*WL\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "draft": r"Draft T\s*\(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
        "cb": r"Block coefficient \(CB\)[^0-9]{0,60}([\d.]+)",
    },
}

#: The SIMMAN 2020 summary page holds three stacked tables; each ship's row block is
#: sliced out before the patterns run so that "T (m)" cannot match a neighbour's table.
SUMMARY_SECTIONS: dict[str, tuple[str, str | None]] = {
    "KVLCC2": ("Table 1 Main Particulars of KVLCC2", "Table 2 Main Particulars of KCS"),
    "KCS": ("Table 2 Main Particulars of KCS", "Table 3 Main Particulars of ONRT"),
    "ONRT": ("Table 3 Main Particulars of ONRT", None),
}

SUMMARY_PATTERNS: dict[str, str] = {
    "length_pp": r"Lpp \(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
    "length_wl": r"Lwl \(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
    "beam": r"Bwl \(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
    "draft": r"T \(m\)\s*([\d.]+(?:\s*[\d.]+)*)",
    "volume": r"Displacement \(m3\)\s*([\d.]+(?:\s*[\d.]+)*)",
    "cb": r"\bCB\s*([\d.]+)",
    "lcb": r"LCB \(%\), fwd\+\s*([-+\d.]+)",
    "scale_row": r"Scale 1\s*([\d.]+(?:\s*[\d.]+)*)",
}

#: Workshop tables list the full-scale value in the first column on the summary page
#: and in the first (KCS, JBC) or second (ONRT) column on the T2015 pages.
PICK = {
    "t2015_kcs": "first",
    "t2015_jbc": "first",
    "t2015_onrt": "last",
    "simman_summary": "first",
}


def page_text(raw_html: str) -> str:
    text = re.sub(r"<script.*?</script>", " ", raw_html, flags=re.S)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text)


def locate_page(raw: Path, key: str) -> Path | None:
    """Publisher page cache lookup: ``raw/pages/<key>.html`` then the legacy flat name."""
    for candidate in (raw / "pages" / f"{key}.html", raw / f"{key}.html"):
        if candidate.is_file():
            return candidate
    return None


def slice_summary(text: str, ship: str) -> str:
    head, tail = SUMMARY_SECTIONS[ship]
    start = text.find(head)
    if start < 0:
        return ""
    end = text.find(tail, start) if tail else -1
    return text[start: end if end > 0 else len(text)]


def extract(text: str, patterns: dict[str, str], pick: str, digest: str,
            page: str) -> dict[str, dict[str, object]]:
    found: dict[str, dict[str, object]] = {}
    for field_name, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            continue
        numbers = [float(value) for value in re.findall(r"[-+]?\d+\.?\d*", match.group(1))]
        if field_name == "scale_row":
            # "Scale 1 48.9355 48.9355": the first number is the full-scale column
            # marker, the second is the model scale of the first model column.
            if len(numbers) < 2:
                continue
            found[field_name] = {
                "value": numbers[1],
                "all_numbers": numbers,
                "snippet": match.group(0)[:180],
                "page_sha256": digest,
                "page": page,
            }
            continue
        if pick == "first":
            numbers = numbers[:1]
        if not numbers:
            continue
        found[field_name] = {
            "value": numbers[-1],
            "snippet": match.group(0)[:180],
            "page_sha256": digest,
            "page": page,
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
    """Assemble one :class:`ShipSpec` per hull from cached bytes only.

    Every reference number travels with the page SHA-256 and the matched snippet, and
    every repair the assembler applies to a publisher table is recorded in
    ``evidence[ship]["corrections"]`` instead of being silently absorbed.
    """
    specs: list[ShipSpec] = []
    evidence: dict[str, object] = {}
    for ship, (archive_name, member, page_keys) in ARCHIVES.items():
        iges = discover_iges(raw, ship)
        archived = raw / archive_name
        record: dict[str, object] = {"corrections": []}
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
        page = None
        page_key = ""
        for candidate in page_keys:
            page = locate_page(raw, candidate)
            if page is not None:
                page_key = candidate
                break
        fields: dict[str, dict[str, object]] = {}
        if page is not None:
            digest = hashlib.sha256(page.read_bytes()).hexdigest()
            text = page_text(page.read_text(errors="replace"))
            if page_key == "simman_summary":
                text = slice_summary(text, ship)
                patterns = SUMMARY_PATTERNS
            else:
                patterns = PAGE_PATTERNS.get(page_key, {})
            fields = extract(text, patterns, PICK[page_key], digest, str(page.relative_to(REPO)))
            record["page"] = {
                "key": page_key,
                "path": str(page.relative_to(REPO)),
                "sha256": digest,
            }
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
            # The ONRT page publishes Cb with its definition nabla/(L_WL B_WL T) and the
            # three dimensions, so the displaced volume follows from the publisher's
            # own definition instead of the displacement in tonnes (which would need
            # an assumed water density).
            volume = float(cb) * float(length) * float(beam) * float(draft)
            record["volume_from_cb"] = {
                "formula": "Cb * Lwl * Bwl * T (publisher's own definition)",
                "value": volume,
            }
        # Decimal-shift repair for a publisher typo: the summary page prints KVLCC2's
        # B_WL as 580 m while its own Cb = nabla/(L_PP B_WL T) = 0.8098 implies 58.0 m.
        # Only a clean x10 discrepancy is corrected, and both the printed and the
        # derived value stay in the evidence file.
        if ship == "KVLCC2" and None not in (volume, cb, length, draft, beam):
            derived_beam = float(volume) / (float(cb) * float(length) * float(draft))
            for shift in (1.0, 10.0, 100.0):
                if abs(float(beam) / shift - derived_beam) <= 0.01 * derived_beam:
                    if shift != 1.0:
                        record["corrections"].append({
                            "field": "beam",
                            "as_published": float(beam),
                            "used": float(beam) / shift,
                            "rule": "publisher decimal typo: "
                                    f"nabla/(Cb*Lpp*T) = {derived_beam:.4f} m, printed {float(beam)} m",
                            "page_sha256": fields["beam"]["page_sha256"],
                        })
                        beam = float(beam) / shift
                    break
        if draft is None:
            record["skip_reason"] = (
                "publisher page lacked the design draft, so the waterline cut "
                "cannot be placed; no volume comparison is claimed"
            )
            continue
        policy, scale_page = GEOMETRY_SCALE_POLICY[ship]
        geometry_scale: float | None = 1.0
        scale_basis = "geometry file is full scale (policy: full)"
        if policy == "inferred":
            geometry_scale = None
            scale_basis = "no published model scale; normalised to the published L_PP"
        elif policy == "published":
            scale_fields = {}
            if scale_page and scale_page != page_key:
                other = locate_page(raw, scale_page)
                if other is not None:
                    other_digest = hashlib.sha256(other.read_bytes()).hexdigest()
                    scale_fields = extract(
                        slice_summary(page_text(other.read_text(errors="replace")), ship),
                        SUMMARY_PATTERNS, PICK[scale_page], other_digest,
                        str(other.relative_to(REPO)),
                    )
            else:
                scale_fields = fields
            if "scale_row" in scale_fields:
                geometry_scale = float(scale_fields["scale_row"]["value"])
                scale_basis = (
                    "published model scale from "
                    f"{scale_fields['scale_row']['page']}: {geometry_scale:.5f}"
                )
            else:
                geometry_scale = None
                scale_basis = "model scale not found on the cache; normalised to the published L_PP"
        specs.append(
            ShipSpec(
                name=ship,
                draft=float(draft),
                length_pp=float(length) if length else 0.0,
                beam=float(beam) if beam else 0.0,
                published_volume_m3=None if volume is None else float(volume),
                published_cb=None if cb is None else float(cb),
                published_lcb_percent=None if lcb is None else float(lcb),
                geometry_scale=geometry_scale,
                scale_basis=scale_basis,
                published_source=str(page.relative_to(REPO)) if page is not None else "",
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
    """Diagnostic figures: displacement comparison plus one half-breadth field per hull.

    Diagnostics only: a figure that cannot be drawn is reported and skipped, it never
    invalidates a hull's verdict (the verdict lives in the manifest and benchmark).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as error:  # pragma: no cover - optional dependency
        print(f"figures skipped: {error}", file=sys.stderr)
        return []
    from tensorcad import fleet_baseline as baseline

    made: list[str] = []
    hulls = [name for name, item in summary["hulls"].items() if item.get("volume_m3")]
    if hulls:
        ref = [summary["hulls"][name]["published_volume_m3"] for name in hulls]
        got = [summary["hulls"][name]["volume_m3"] for name in hulls]
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        positions = np.arange(len(hulls))
        ax.bar(positions - 0.2, ref, 0.4, label="published (m³)", color="#4c72b0")
        ax.bar(positions + 0.2, got, 0.4, label="reconstructed (m³)", color="#dd8452")
        ax.set_xticks(positions)
        ax.set_xticklabels(hulls)
        ax.set_ylabel("displacement volume below DWL (m³)")
        ax.set_title("Baseline fleet: reconstructed vs published displacement")
        for index, name in enumerate(hulls):
            error = summary["hulls"][name]["volume_rel_error_percent"]
            ax.text(index, max(ref[index], got[index]) * 1.02, f"{error:+.2f}%",
                    ha="center", fontsize=8)
        ax.legend()
        fig.tight_layout()
        path = out / "figs" / "volume_vs_published.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=140)
        plt.close(fig)
        made.append(str(path.relative_to(REPO)))

    for name, manifest in manifests.items():
        try:
            field = manifest.get("field")
            if not field:
                continue
            spec = baseline.ShipSpec(**{key: manifest["spec"][key] for key in
                                       ("name", "draft", "length_pp", "beam")})
            points = [(row[0], 0.0, row[1]) for row in field["breadth_samples"]]
            grid = baseline.field_from_points(points, spec.draft, dx=1.0, dz=0.25)
            fig, ax = plt.subplots(figsize=(7.2, 3.2))
            mesh = ax.pcolormesh(grid.x, grid.z, grid.breadth.T, shading="auto", cmap="viridis")
            ax.set_aspect("equal")
            ax.set_xlabel("x (m, aft -> fwd)")
            ax.set_ylabel("z (m above keel)")
            ax.set_title(f"{name}: reconstructed half-breadth field b(x, z)")
            fig.colorbar(mesh, ax=ax, label="b (m)")
            fig.tight_layout()
            path = out / "figs" / f"{name}_half_breadth_field.png"
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=140)
            plt.close(fig)
            made.append(str(path.relative_to(REPO)))
        except Exception as error:  # noqa: BLE001 - diagnostics only
            print(f"figure skipped for {name}: {type(error).__name__}: {error}", file=sys.stderr)
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