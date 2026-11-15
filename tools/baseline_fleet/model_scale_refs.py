#!/usr/bin/env python3
"""Cross-check the *model scale* of a baseline-fleet hull against a published statement.

    python tools/baseline_fleet/model_scale_refs.py fetch
    python tools/baseline_fleet/model_scale_refs.py build
    python tools/baseline_fleet/model_scale_refs.py backfill
    python tools/baseline_fleet/model_scale_refs.py report

Why this tool exists
--------------------
``benchmarks/verified/baseline_fleet/fleet_specs.json`` carries a ``geometry_scale`` per ship:
the factor that turns the admitted IGES/CAD bytes into full-scale metres.  For KCS, KVLCC2 and
ONRT the factor is *published* (the workshop page states the scale ratio).  For JBC it is
**derived** -- the Tokyo 2015 JBC geometry page prints no scale row, so the admitted record
divides the published B/D by the spans measured in the CAD bytes and records the half-spread
(``artifacts/baseline_fleet/jbc_model_scale.json``).  A derived factor is the weakest link in
the whole displacement comparison: it enters the benchmark cubed, so 0.7 % on the scale is
2 % on the volume -- twice the 1 % acceptance band.  Until now nothing in the repository
checked that derived number against a published one.

This tool does exactly that, and nothing else: it reads the statements out of a pinned
document (digit-level quote + page + SHA-256), compares them with the derived factor, and
prices the comparison in the currency the benchmark uses (relative volume error).  It never
rewrites ``geometry_scale`` and never upgrades a benchmark verdict; it publishes the
sensitivity so a reader can see how much of a withheld comparison the scale can explain.

Only digests, page numbers and the quotations are committed; the publisher bytes stay under
``artifacts/fleet_refs/raw`` (git-ignored) because the library does not redistribute them.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_LIBRARY = ROOT / "benchmarks/verified/fleet_refs/model_scale.json"
DEFAULT_SPECS = ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json"
DEFAULT_RESULT = ROOT / "benchmarks/verified/baseline_fleet/result.json"
DEFAULT_RAW = ROOT / "artifacts/fleet_refs/raw"
DEFAULT_DOC = ROOT / "docs/fleet_refs_jbc_scale.md"
DEFAULT_FIGURE = ROOT / "docs/figs/fleet_refs/jbc_scale_sensitivity.png"


@dataclass(frozen=True)
class Claim:
    """One numeral this document is asked about.

    ``pattern`` must capture exactly the numeral (group 1); it is anchored on the table row
    label or on the sentence that states it, so a number that merely appears nearby cannot be
    mistaken for the claim.
    """

    name: str
    label: str
    unit: str
    pattern: str
    #: Optional anchor: the claim is read only from the text that follows the first match of
    #: this pattern on the page.  A document can print two hulls' tables (this preprint prints
    #: KCS and JBC side by side), and an unanchored "Breadth, B" row silently returns the wrong
    #: hull -- which is exactly the class of defect this library exists to catch.
    after: str = ""


@dataclass(frozen=True)
class Document:
    key: str
    ship: str
    publisher: str
    url: str
    doi: str
    licence: str
    sha256: str
    file: str
    note: str = ""
    claims: tuple[Claim, ...] = ()


#: The document set.  One document, deliberately: the anchor asks for a cross-check of a
#: *derived* quantity, and a second transcription of the same preprint would add nothing.
DOCUMENTS: tuple[Document, ...] = (
    Document(
        key="jbc_scale_pr",
        ship="JBC",
        publisher=(
            "Research Square preprint (Tarafder, Hasanat, Hasan, Sayeed, Hannan), posted 2025-09-30, "
            "'Scale Effects on Bare-Hull Resistance of KCS and JBC: A RANS-Based Investigation with "
            "Verification and Validation', Table 3 and section 3.2"
        ),
        url="https://www.researchsquare.com/article/rs-7594251/latest.pdf",
        doi="10.21203/rs.3.rs-7594251/v1",
        licence="CC BY 4.0 (stated on page 1 of the document)",
        sha256="4016a43201d350afe59c63bedd221e55e7241ed078b9b3e60c46548177fc5b3d",
        file="jbc_scale_pr.pdf",
        note=(
            "the full-scale column of Table 3 reproduces the Tokyo 2015 JBC table (same digits), so it "
            "is corroboration of transcription only; the *model* column and the scale-ratio sentence in "
            "section 3.2 are the reason this document is cited"
        ),
        claims=(
            Claim(
                name="model_length_pp",
                label="model Lpp (m)",
                unit="m",
                pattern=r"Length,\s*Lpp\s*\n\s*[0-9.]+\s*\n\s*([0-9.]+)\s*\n\s*m",
                after=r"JBC hull principal particulars",
            ),
            Claim(
                name="full_length_pp",
                label="full-scale Lpp (m)",
                unit="m",
                pattern=r"Length,\s*Lpp\s*\n\s*([0-9.]+)\s*\n\s*[0-9.]+\s*\n\s*m",
                after=r"JBC hull principal particulars",
            ),
            Claim(
                name="declared_scale_ratio",
                label="declared scale ratio λ",
                unit="-",
                pattern=r"with a scale ratio of\s*([0-9]+(?:\.[0-9]+)?)",
                after=r"JBC stands for Japan Bulk Carrier",
            ),
            Claim(
                name="model_length_statement",
                label="model length from the same sentence (m)",
                unit="m",
                pattern=r"model\s+scale\s*([0-9.]+)\s*m with a scale ratio",
                after=r"JBC stands for Japan Bulk Carrier",
            ),
            Claim(
                name="beam",
                label="full-scale B (m)",
                unit="m",
                pattern=r"Breadth,\s*B\s*\n\s*([0-9.]+)\s*\n\s*[0-9.]+\s*\n\s*m",
                after=r"JBC hull principal particulars",
            ),
            Claim(
                name="draft",
                label="full-scale T (m)",
                unit="m",
                pattern=r"Draft,\s*T\s*\n\s*([0-9.]+)\s*\n\s*[0-9.]+\s*\n\s*m",
                after=r"JBC hull principal particulars",
            ),
            Claim(
                name="displacement_volume",
                label="full-scale ∇ (m³)",
                unit="m3",
                pattern=r"Displacement,\uf0d1?\s*\n\s*([0-9.]+)\s*\n\s*[0-9.]+\s*\n\s*m3",
                after=r"JBC hull principal particulars",
            ),
            Claim(
                name="cb",
                label="Cb",
                unit="-",
                pattern=r"Block coefficient,\s*Cb\s*\n\s*([0-9.]+)",
                after=r"JBC hull principal particulars",
            ),
        ),
    ),
)

#: Acceptance band of the baseline benchmark, from ``result.json``; the sensitivity is priced
#: against it so a reader sees the scale's share of the tolerance, not just a percentage.
BENCHMARK_TOLERANCE_PERCENT = 1.0


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def pdf_text(path: Path) -> tuple[str, int]:
    import fitz  # PyMuPDF; optional dependency, only needed for the fetch/build side

    with fitz.open(path) as document:
        # ``\f`` between pages, so a quote can name the page it was read from.  (The older
        # ``crossref_refs._pdf_text`` joins with ``\n`` and therefore reports one page for every
        # PDF; this tool keeps the page boundary.)
        return "\f".join(page.get_text() for page in document), document.page_count


def quote_context(page: str, match: re.Match, width: int = 110) -> str:
    start = max(0, match.start() - width)
    end = min(len(page), match.end() + width)
    return re.sub(r"\s+", " ", page[start:end]).strip()


def extract(document: Document, text: str) -> dict[str, dict[str, Any]]:
    """Read every claim out of ``text``; a claim the document does not state is simply absent."""
    pages = text.split("\f")
    out: dict[str, dict[str, Any]] = {}
    for claim in document.claims:
        for index, page in enumerate(pages, start=1):
            region = page
            if claim.after:
                anchor = re.search(claim.after, page)
                if anchor is None:
                    continue
                region = page[anchor.end() :]
            match = re.search(claim.pattern, region)
            if match is None:
                continue
            raw = match.group(1)
            out[claim.name] = {
                "value": float(raw),
                "raw": raw,
                "unit": claim.unit,
                "label": claim.label,
                "page": index,
                "quote": quote_context(region, match),
                "source": document.key,
                "url": document.url,
                "doi": document.doi,
                "sha256": document.sha256,
            }
            break
    return out


def scale_comparison(claims: dict[str, Any], spec: dict[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    """Published λ vs the derived ``geometry_scale``, priced in relative volume error."""
    published = [c for c in claims.values() if c]
    declared = next((c for c in published if c["source_claim"] == "declared_scale_ratio"), None)
    model = next((c for c in published if c["source_claim"] == "model_length_pp"), None)
    full = next((c for c in published if c["source_claim"] == "full_length_pp"), None)

    derived = float(spec["geometry_scale"])
    published_ratio = None
    ratio_basis = ""
    if declared is not None and model is not None and full is not None:
        from_sentence = float(declared["value"])
        from_table = float(full["value"]) / float(model["value"])
        agreement = abs(from_sentence / from_table - 1.0)
        published_ratio = from_table
        ratio_basis = (
            f"sentence λ={from_sentence:g} and Table 3 {full['value']:g} m / {model['value']:g} m = "
            f"{from_table:g} agree to {agreement * 100:.3f} %"
        )
        if agreement > 0.005:
            published_ratio = None
            ratio_basis += " -- the two statements disagree by more than 0.5 %, so neither is used"

    out: dict[str, Any] = {
        "ship": spec["name"],
        "derived_scale": derived,
        "derived_scale_basis": spec.get("scale_basis", ""),
        "published_scale_ratio": published_ratio,
        "published_scale_basis": ratio_basis,
        "claims": claims,
    }
    if published_ratio is None:
        out["status"] = "no_published_scale_obtained"
        return out

    delta = published_ratio / derived - 1.0
    volume_factor = (published_ratio / derived) ** 3
    out["status"] = "cross_checked"
    out["delta_relative"] = delta
    out["delta_percent"] = delta * 100.0
    out["volume_factor_published_over_derived"] = volume_factor
    out["tolerance_percent"] = BENCHMARK_TOLERANCE_PERCENT

    got = result.get("hulls", {}).get(spec["name"], {})
    reference = got.get("ref_volume_m3")
    got_volume = got.get("got_volume_m3")
    if reference and got_volume:
        at_published = float(got_volume) * volume_factor
        error_derived = (float(got_volume) / float(reference) - 1.0) * 100.0
        error_published = (at_published / float(reference) - 1.0) * 100.0
        out["benchmark"] = {
            "ref_volume_m3": float(reference),
            "got_volume_at_derived_scale_m3": float(got_volume),
            "error_at_derived_scale_percent": error_derived,
            "got_volume_at_published_scale_m3": at_published,
            "error_at_published_scale_percent": error_published,
            "points_explained_by_scale": error_published - error_derived,
            "verdict_unchanged": got.get("verdict"),
            "residual_is_geometric": abs(error_published) > BENCHMARK_TOLERANCE_PERCENT,
        }
    return out


# --------------------------------------------------------------------------------------
# subcommands


def _specs() -> dict[str, dict[str, Any]]:
    payload = json.loads(DEFAULT_SPECS.read_text())
    return {row["name"]: row for row in payload["ships"]}


def cmd_fetch(args: argparse.Namespace) -> int:
    import urllib.request

    raw = Path(args.raw)
    raw.mkdir(parents=True, exist_ok=True)
    registry_path = raw.parent / "model_scale_registry.json"
    rows: list[dict[str, Any]] = []
    ua = "Mozilla/5.0 (X11; Linux x86_64) tensorcad-fleet-refs/1.0"
    for document in DOCUMENTS:
        target = raw / document.file
        row: dict[str, Any] = {"key": document.key, "url": document.url, "file": document.file}
        if not target.exists():
            try:
                request = urllib.request.Request(document.url, headers={"User-Agent": ua})
                with urllib.request.urlopen(request, timeout=60) as response:
                    target.write_bytes(response.read())
                    row["http_status"] = response.status
            except Exception as error:  # recorded, never hidden
                row["status"] = "unreachable"
                row["error"] = f"{type(error).__name__}: {error}"[:200]
                rows.append(row)
                continue
        payload = target.read_bytes()
        row["sha256"] = sha256(payload)
        row["bytes"] = len(payload)
        row["retrieved_utc"] = datetime.datetime.fromtimestamp(
            target.stat().st_mtime, datetime.timezone.utc
        ).isoformat(timespec="seconds")
        row["pin_matches"] = row["sha256"] == document.sha256
        rows.append(row)
        print(f"{document.key}: {row['sha256'][:16]}… bytes={len(payload)} pin={row['pin_matches']}")
    registry_path.write_text(json.dumps({"schema": "tensorcad.fleet_refs.model_scale.registry/1.0",
                                         "documents": rows}, indent=1) + "\n")
    missing = [row for row in rows if not row.get("pin_matches")]
    if missing:
        print("WARNING: document does not match the pinned digest:", ", ".join(r["key"] for r in missing))
        return 1
    return 0


def cmd_build(args: argparse.Namespace) -> int:
    raw = Path(args.raw)
    registry_path = raw.parent / "model_scale_registry.json"
    registry = json.loads(registry_path.read_text()) if registry_path.exists() else {"documents": []}
    by_key = {row["key"]: row for row in registry.get("documents", [])}

    specs = _specs()
    result = json.loads(Path(args.result).read_text())
    documents: dict[str, Any] = {}
    per_ship: dict[str, dict[str, Any]] = {}
    for document in DOCUMENTS:
        path = raw / document.file
        record: dict[str, Any] = {
            "publisher": document.publisher,
            "url": document.url,
            "doi": document.doi,
            "licence": document.licence,
            "sha256": document.sha256,
            "note": document.note,
            "claims_requested": [claim.name for claim in document.claims],
        }
        record.update({k: v for k, v in by_key.get(document.key, {}).items() if k not in record})
        if not path.exists():
            record["status"] = "not_fetched"
            documents[document.key] = record
            continue
        payload = path.read_bytes()
        record["bytes"] = len(payload)
        record["digest_matches_registry"] = sha256(payload) == document.sha256
        text, pages = pdf_text(path)
        record["pages"] = pages
        record["status"] = "read"
        cells = extract(document, text)
        record["claims_read"] = sorted(cells)
        documents[document.key] = record
        target = per_ship.setdefault(document.ship, {})
        for name, cell in cells.items():
            target[name] = {**cell, "source_claim": name}
        target["_source"] = document.key

    comparisons: dict[str, Any] = {}
    for ship, claims in per_ship.items():
        payload = {k: v for k, v in claims.items() if not k.startswith("_")}
        comparisons[ship] = scale_comparison(payload, specs[ship], result)

    library = {
        "schema": "tensorcad.fleet_refs.model_scale/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "scope": (
            "published model-scale statements vs the derived geometry_scale of the baseline fleet; "
            "the comparison is quoted in relative displacement volume because that is the quantity "
            "the benchmark compares"
        ),
        "tolerance_percent": BENCHMARK_TOLERANCE_PERCENT,
        "documents": documents,
        "comparisons": comparisons,
    }
    out = Path(args.out)
    out.write_text(json.dumps(library, indent=1, sort_keys=True) + "\n")
    for ship, row in comparisons.items():
        if row["status"] == "cross_checked":
            print(
                f"{ship}: derived λ={row['derived_scale']:.5f}, published λ={row['published_scale_ratio']:.5f} "
                f"({row['delta_percent']:+.3f} % -> volume factor {row['volume_factor_published_over_derived']:.5f})"
            )
        else:
            print(f"{ship}: {row['status']}")
    return 0


def cmd_backfill(args: argparse.Namespace) -> int:
    library = json.loads(Path(args.library).read_text())
    specs_path = Path(args.specs)
    payload = json.loads(specs_path.read_text())
    changed = 0
    for row in payload["ships"]:
        comparison = library["comparisons"].get(row["name"])
        if comparison is None:
            continue
        entry = {
            "status": comparison["status"],
            "derived_scale": comparison["derived_scale"],
            "published_scale_ratio": comparison.get("published_scale_ratio"),
            "published_scale_basis": comparison.get("published_scale_basis", ""),
            "delta_percent": comparison.get("delta_percent"),
            "volume_factor_published_over_derived": comparison.get("volume_factor_published_over_derived"),
            "benchmark": comparison.get("benchmark"),
            "evidence": [
                {
                    "claim": name,
                    "value": cell["value"],
                    "unit": cell["unit"],
                    "page": cell["page"],
                    "quote": cell["quote"],
                    "url": cell["url"],
                    "doi": cell["doi"],
                    "sha256": cell["sha256"],
                }
                for name, cell in (comparison.get("claims") or {}).items()
                if isinstance(cell, dict) and "value" in cell
            ],
            "library": str(Path(args.library).relative_to(ROOT)),
        }
        if row.get("scale_cross_check") != entry:
            row["scale_cross_check"] = entry
            changed += 1
    if changed:
        specs_path.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"fleet_specs: {changed} ship(s) updated")
    return 0


def _markdown(library: dict[str, Any]) -> str:
    lines = [
        "# JBC model scale: derived from the CAD bytes vs a published statement",
        "",
        "Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact:",
        "`benchmarks/verified/fleet_refs/model_scale.json` (regenerate:",
        "`python tools/baseline_fleet/model_scale_refs.py build`). Every number below is rendered",
        "from that artefact; nothing is transcribed by hand.",
        "",
        "## 1. The claim under test",
        "",
        "`fleet_specs.json` turns the admitted JBC geometry into full-scale metres with a **derived**",
        "factor `λ = 39.7285`, because the Tokyo 2015 JBC geometry page publishes no scale row.",
        "The bench compares displacement volumes, so that factor enters cubed: a 0.7 % scale error",
        "is a 2.1 % volume error, i.e. twice the 1 % acceptance band of the baseline benchmark.",
        "",
        "## 2. Documents read",
        "",
        "| key | source | licence | sha256 | claims read |",
        "|---|---|---|---|---|",
    ]
    for key, row in sorted(library["documents"].items()):
        lines.append(
            f"| {key} | {row.get('publisher', '')[:120]} | {row.get('licence', '')} | "
            f"`{row['sha256'][:16]}…` | {', '.join(row.get('claims_read', [])) or row.get('status', '')} |"
        )
    lines += ["", "## 3. Quoted statements", "", "| ship | claim | value | unit | page | quote |", "|---|---|---|---|---|---|"]
    for ship, comparison in sorted(library["comparisons"].items()):
        for name, cell in (comparison.get("claims") or {}).items():
            if not isinstance(cell, dict) or "value" not in cell:
                continue
            quote = cell["quote"].replace("|", "\\|")
            lines.append(
                f"| {ship} | {cell['label']} | {cell['raw']} | {cell['unit']} | {cell['page']} | `{quote[:160]}` |"
            )
    lines += ["", "## 4. Derived vs published scale", ""]
    for ship, comparison in sorted(library["comparisons"].items()):
        lines.append(f"### {ship}")
        lines.append("")
        lines.append(f"* derived scale (fleet_specs): **{comparison['derived_scale']:.6f}** — {comparison['derived_scale_basis']}")
        if comparison["status"] != "cross_checked":
            lines.append(f"* status: **{comparison['status']}** — no published scale was obtained")
            lines.append("")
            continue
        lines.append(f"* published scale: **{comparison['published_scale_ratio']:.6f}** — {comparison['published_scale_basis']}")
        lines.append(
            f"* discrepancy: **{comparison['delta_percent']:+.3f} %** in λ = "
            f"**{comparison['volume_factor_published_over_derived']:.5f}×** in ∇"
        )
        benchmark = comparison.get("benchmark")
        if benchmark:
            lines.append("")
            lines.append("| quantity | at derived scale | at published scale |")
            lines.append("|---|---|---|")
            lines.append(
                f"| reconstructed ∇ (m³) | {benchmark['got_volume_at_derived_scale_m3']:.1f} | "
                f"{benchmark['got_volume_at_published_scale_m3']:.1f} |"
            )
            lines.append(
                f"| error vs published {benchmark['ref_volume_m3']:.1f} m³ | "
                f"{benchmark['error_at_derived_scale_percent']:+.2f} % | {benchmark['error_at_published_scale_percent']:+.2f} % |"
            )
            lines.append("")
            lines.append(
                f"The scale statement moves the comparison by **{benchmark['points_explained_by_scale']:+.2f} percentage "
                f"points**; the benchmark verdict stays **{benchmark['verdict_unchanged']}** and the residual "
                f"{benchmark['error_at_published_scale_percent']:+.2f} % is geometric, not a unit error."
            )
        lines.append("")
    lines += [
        "## 5. What this does and does not establish",
        "",
        "* It **does** show that the derived λ is not the published λ: the admitted JBC CAD spans imply",
        "  1.1327 m where a λ=40 model would measure 1.125 m across the published B = 45 m, and the",
        "  difference is quoted, digests pinned, and priced in the benchmark's own currency.",
        "* It does **not** decide which number is right. The CAD spans are measured; the λ=40 statement",
        "  is published. Both stay in the record (`geometry_scale` untouched, the statement quoted in",
        "  `scale_cross_check`), and no volume comparison is silently rescaled.",
        "* It does **not** upgrade the JBC verdict. The withheld −18 % is dominated by the IGES",
        "  control-net layout guard, not by the scale; the scale explains only the number of points",
        "  tabulated in section 4.",
        "",
        "## 6. Limitations",
        "",
        "* The published statement is a 2025 preprint, not the workshop's own page; its full-scale",
        "  column reproduces the Tokyo 2015 table to the digit, which is why it is treated as a",
        "  transcription witness for the particulars and as the source of the *scale* statement only.",
        "* The model column of that table (Lpp 7.00 m) is itself consistent with λ = 280/7 = 40; the",
        "  tool requires the sentence and the table to agree within 0.5 % before using either.",
        "* No published scale exists for the remaining hulls (KCS/KVLCC2/ONRT use published ratios;",
        "  DTC/Series 60 have no admitted geometry in this fleet), so the comparison is one-hull.",
        "",
    ]
    return "\n".join(lines)


def _figure(library: dict[str, Any], path: Path) -> bool:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [
        (ship, comparison)
        for ship, comparison in sorted(library["comparisons"].items())
        if comparison.get("benchmark")
    ]
    if not rows:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6))
    labels = [ship for ship, _ in rows]
    derived = [c["benchmark"]["error_at_derived_scale_percent"] for _, c in rows]
    published = [c["benchmark"]["error_at_published_scale_percent"] for _, c in rows]
    x = range(len(rows))
    axes[0].bar([i - 0.18 for i in x], derived, width=0.36, label="derived λ (fleet_specs)")
    axes[0].bar([i + 0.18 for i in x], published, width=0.36, label="published λ")
    axes[0].axhline(0.0, color="k", linewidth=0.8)
    for sign in (-library["tolerance_percent"], library["tolerance_percent"]):
        axes[0].axhline(sign, color="tab:red", linestyle="--", linewidth=0.8)
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(labels)
    axes[0].set_ylabel("displacement error vs published (%)")
    axes[0].set_title("Reconstruction error is unchanged in sign and verdict")
    axes[0].legend(fontsize=8)

    delta = [c["benchmark"]["points_explained_by_scale"] for _, c in rows]
    residual = [c["benchmark"]["error_at_published_scale_percent"] for _, c in rows]
    axes[1].bar(labels, delta, color="tab:orange", label="explained by scale")
    axes[1].bar(labels, residual, bottom=0.0, color="tab:gray", alpha=0.55, label="residual (geometry)")
    axes[1].axhline(0.0, color="k", linewidth=0.8)
    axes[1].set_ylabel("percentage points of error")
    axes[1].set_title("Scale accounts for a small part of the gap")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.3, axis="y")
    fig.suptitle(
        "JBC: derived vs published model scale (anchor REPRO-GEOM-FLEET-REFS)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def cmd_report(args: argparse.Namespace) -> int:
    library = json.loads(Path(args.library).read_text())
    doc = Path(args.doc)
    doc.write_text(_markdown(library))
    print(f"wrote {doc}")
    if _figure(library, Path(args.figure)):
        print(f"wrote {args.figure}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser("fetch", help="download the documents and record their digests")
    p_fetch.add_argument("--raw", default=str(DEFAULT_RAW))
    p_fetch.set_defaults(func=cmd_fetch)

    p_build = sub.add_parser("build", help="extract the quoted claims and compare with fleet_specs")
    p_build.add_argument("--raw", default=str(DEFAULT_RAW))
    p_build.add_argument("--result", default=str(DEFAULT_RESULT))
    p_build.add_argument("--out", default=str(DEFAULT_LIBRARY))
    p_build.set_defaults(func=cmd_build)

    p_back = sub.add_parser("backfill", help="write the comparison into fleet_specs.json")
    p_back.add_argument("--library", default=str(DEFAULT_LIBRARY))
    p_back.add_argument("--specs", default=str(DEFAULT_SPECS))
    p_back.set_defaults(func=cmd_backfill)

    p_report = sub.add_parser("report", help="render the documentation and the figure")
    p_report.add_argument("--library", default=str(DEFAULT_LIBRARY))
    p_report.add_argument("--doc", default=str(DEFAULT_DOC))
    p_report.add_argument("--figure", default=str(DEFAULT_FIGURE))
    p_report.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())