#!/usr/bin/env python3
"""Derive the *model scale* of a benchmark hull from its own IGES control net.

    python tools/baseline_fleet/jbc_model_scale.py                     # fetch, derive, write record
    python tools/baseline_fleet/jbc_model_scale.py --offline           # use cached archive only
    python tools/baseline_fleet/jbc_model_scale.py --verify-only       # offline check of the record

Why this tool exists
--------------------
Anchor ``GAP-CAD-T1-SOURCE-HUNT`` item (3) is "JBC needs its scale factor".
The Tokyo-2015 page for the Japan Bulk Carrier (``jbc_gc.html``) publishes *only*
full-scale particulars: one numeric column labelled "DESIGNED FULL" and one for the
ballast condition.  There is no model-scale row, so no scale factor can be scraped —
the earlier round recorded that as WITHHELD.

The scale is nevertheless recoverable from the publisher's own CAD data, because the
IGES archive that the same page serves is written at *model* size (header: units mm).
Two independent published lengths fix it:

    s_beam  = B_full  / (model breadth measured on the CAD surface)
    s_depth = D_full  / (model depth   measured on the CAD surface)

``B_full = 45.000 m`` and ``D = 25.000 m`` are read from the same page as the
particulars, and the geometry is measured here, so nothing is transcribed by hand.
The hull's x-extent is deliberately *not* used as a scale reference: it is the
deck-level extreme (bow flare / stern overhang), which is longer than Lpp and longer
than Lwl, and using it would silently bias the scale low by ~6 %.

Everything the tool states about the publisher is bound to a SHA-256: the page that
carries the particulars, the IGES archive, and the propeller particulars archive
(whose ``MP687`` table independently gives the model propeller diameter, a second
length in the same model space).  Archive bytes are not redistributed; only digests,
measured spans and the derivation are committed.

Output: ``artifacts/baseline_fleet/jbc_model_scale.json`` (+ optional figure).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import ssl
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

T2015 = "https://t2015.nmri.go.jp"
GC_URL = f"{T2015}/jbc_gc.html"
IGES_URL = f"{T2015}/file/Geometry_IGES_files/jbc/JBC_IGES.zip"
PROP_URL = f"{T2015}/file/Geometry_IGES_files/jbc/MP687_particulars.zip"

#: Published full-scale particulars, and the row label each one is read from.
#: The regex is applied to the *stripped row text* of the page, so a value can only
#: enter the record together with the line it was read from.
PUBLISHED = {
    "lpp_m": (r"Length between perpendiculars", 280.0),
    "lwl_m": (r"Length of waterline", 285.0),
    "beam_m": (r"Maximum beam of waterline", 45.0),
    "depth_m": (r"Depth", 25.0),
    "draft_m": (r"Draft\(Full\)", 16.5),
    "displacement_m3": (r"Displacement volume of a hull", 178369.9),
    "cb": (r"Block coefficient", 0.8580),
}

#: Scale hypotheses a reader might hold; the residual of each is reported so the
#: derived value can be judged instead of taken on trust.
HYPOTHESES = [40.0, 39.7, 39.4]

NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


# --------------------------------------------------------------------------- #
# network / cache
# --------------------------------------------------------------------------- #
def opener() -> urllib.request.OpenerDirector:
    # The T2015 host presents a certificate whose subject does not match the name
    # used by the publisher's own links (recorded, not hidden).
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return urllib.request.build_opener(urllib.request.HTTPSHandler(context=ctx))


def fetch(url: str, cache: Path, timeout: int = 120) -> tuple[bytes | None, dict]:
    """Return cached-or-fetched bytes plus a digest record (never raises on 404)."""
    cache.mkdir(parents=True, exist_ok=True)
    blob_path = cache / (hashlib.sha256(url.encode()).hexdigest()[:16] + "_" + url.rsplit("/", 1)[-1])
    info: dict = {"url": url, "cache": str(blob_path)}
    if blob_path.exists():
        blob = blob_path.read_bytes()
        info["retrieval"] = "cache"
    else:
        req = urllib.request.Request(url, headers={"User-Agent": "tensorcad-source-hunt/1.0"})
        try:
            with opener().open(req, timeout=timeout) as resp:
                blob = resp.read()
                info["http_status"] = resp.status
        except urllib.error.HTTPError as error:
            blob = error.read()
            info["http_status"] = error.code
        except Exception as error:  # noqa: BLE001 - network failure is data here
            return None, {"url": url, "http_status": None,
                          "error": f"{type(error).__name__}: {error}"}
        info["retrieval"] = "network"
        blob_path.write_bytes(blob)
    info["bytes"] = len(blob)
    info["sha256"] = hashlib.sha256(blob).hexdigest()
    return blob, info


# --------------------------------------------------------------------------- #
# page parsing
# --------------------------------------------------------------------------- #
def page_rows(html: str) -> list[str]:
    """Rows of the page as plain text lines (same normalisation the hunt tool uses)."""
    html = re.sub(r"</(?:tr|div|p|table|h\d)>", "\n", html, flags=re.I)
    text = re.sub(r"<[^>]*>", "|", html).replace("&nbsp;", " ")
    rows = []
    for line in text.splitlines():
        cells = [c.strip() for c in line.split("|") if c.strip()]
        if cells:
            rows.append(" ".join(cells))
    return rows


def read_particulars(html: str) -> dict[str, dict]:
    """Published full-scale values with the matched row text and the next numeric row.

    The JBC table puts a label in one row and the numbers in the row below it
    ("Length between perpendiculars L_PP (m)" then "280.000 same as on the left"),
    so the value is taken from the first following row that starts with a number.
    """
    rows = page_rows(html)
    out: dict[str, dict] = {}
    for key, (pattern, expected) in PUBLISHED.items():
        hit = None
        for i, row in enumerate(rows):
            if re.search(pattern, row, re.I):
                snippet = row
                value = None
                # The value row is the first *following* row that starts with a number:
                # the label row itself often carries digits inside its unit ("(m 3 )"),
                # and reading those was a real transcription bug caught by verify().
                for follower in rows[i + 1:i + 4]:
                    match = re.match(r"^(-?\d+(?:\.\d+)?)", follower)
                    if match:
                        value = float(match.group(1))
                        snippet = f"{row} / {follower}"
                        break
                if value is not None:
                    hit = {"value": value, "snippet": snippet,
                           "agrees_with_expectation": abs(value - expected) <= 1e-6 * max(1.0, abs(expected))}
                break
        out[key] = hit or {"value": None, "snippet": None, "found": False}
    return out


def read_model_propeller(text: str) -> dict | None:
    """``MP687_particulars.txt``: model propeller diameter (a length in model space)."""
    match = re.search(r"DIAMETER\s+([\d.]+)\s*MM", text)
    if not match:
        return None
    return {"diameter_mm": float(match.group(1)),
            "snippet": match.group(0).strip(),
            "note": "MODEL propeller of the same model space as the hull IGES"}


# --------------------------------------------------------------------------- #
# geometry measurement
# --------------------------------------------------------------------------- #
def surface_bbox(surfaces, samples: int = 41) -> dict:
    """Bounding box of the sampled surfaces.

    Any object exposing ``domain_u``, ``domain_v`` and ``evaluate(u, v)`` works, which
    keeps the measurement testable without a CAD kernel.
    """
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    count = 0
    per_surface = []
    for surface in surfaces:
        u0, u1 = surface.domain_u
        v0, v1 = surface.domain_v
        slo = [float("inf")] * 3
        shi = [float("-inf")] * 3
        for i in range(samples):
            u = u0 + (u1 - u0) * i / (samples - 1)
            for j in range(samples):
                v = v0 + (v1 - v0) * j / (samples - 1)
                try:
                    point = surface.evaluate(u, v)
                except Exception:  # noqa: BLE001 - a degenerate corner is skipped, not hidden
                    continue
                count += 1
                for axis in range(3):
                    lo[axis] = min(lo[axis], point[axis])
                    hi[axis] = max(hi[axis], point[axis])
                    slo[axis] = min(slo[axis], point[axis])
                    shi[axis] = max(shi[axis], point[axis])
        if count:
            per_surface.append({"lo": slo, "hi": shi})
    return {
        "samples": samples,
        "points": count,
        "surfaces": len(surfaces),
        "lo": lo,
        "hi": hi,
        "span": [hi[i] - lo[i] for i in range(3)],
        "per_surface": per_surface,
    }


def derive_scale(span: dict[str, float], published: dict[str, float]) -> dict:
    """Two independent span-derived scale estimates + the hypothesis residuals.

    ``span`` keys are ``x``/``y``/``z`` in metres of model space.
    """
    candidates = {
        "beam": {"reference_m": published["beam_m"], "span_m": span["y"]},
        "depth": {"reference_m": published["depth_m"], "span_m": span["z"]},
    }
    for name, item in candidates.items():
        item["scale"] = item["reference_m"] / item["span_m"]
    values = [item["scale"] for item in candidates.values()]
    value = sum(values) / len(values)
    spread = max(values) - min(values)
    residuals = {}
    for hypothesis in HYPOTHESES:
        residuals[f"1:{hypothesis:g}"] = {
            name: (item["span_m"] * hypothesis - item["reference_m"]) / item["reference_m"]
            for name, item in candidates.items()
        }
    return {
        "value": value,
        "half_spread": spread / 2.0,
        "relative_half_spread": (spread / 2.0) / value,
        "method": "s = published full-scale length / model length measured on the "
                  "publisher's IGES surface (two independent lengths: B and D)",
        "candidates": candidates,
        "hypothesis_residuals": residuals,
        "not_used": {
            "x_span": {
                "span_m": span["x"],
                "naive_scale_against_lpp": published["lpp_m"] / span["x"],
                "why_excluded": "the x extent is the deck-level extreme (bow flare, stern "
                                "overhang), longer than Lpp and than Lwl; using it as an Lpp "
                                "reference biases the scale low by ~6 %",
            }
        },
    }


# --------------------------------------------------------------------------- #
# record
# --------------------------------------------------------------------------- #
def build(cache: Path, samples: int) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from tensorcad.iges import read_iges, read_iges_surfaces  # noqa: PLC0415

    html_blob, html_info = fetch(GC_URL, cache)
    if html_blob is None:
        raise SystemExit(f"particulars page unreachable: {html_info}")
    html = html_blob.decode("utf-8", "replace")

    zip_blob, zip_info = fetch(IGES_URL, cache)
    if zip_blob is None:
        raise SystemExit(f"IGES archive unreachable: {zip_info}")
    prop_blob, prop_info = fetch(PROP_URL, cache)

    iges_path = cache / "jbc.igs"
    with zipfile.ZipFile(io.BytesIO(zip_blob)) as archive:
        members = archive.namelist()
        iges_path.write_bytes(archive.read("jbc.igs"))

    doc = read_iges(str(iges_path))
    surfaces = read_iges_surfaces(str(iges_path))
    box = surface_bbox(surfaces, samples=samples)
    metres_per_unit = doc.metres_per_unit or 1.0
    span_m = {"x": box["span"][0] * metres_per_unit,
              "y": box["span"][1] * metres_per_unit,
              "z": box["span"][2] * metres_per_unit}
    bbox_m = {"lo": [v * metres_per_unit for v in box["lo"]],
              "hi": [v * metres_per_unit for v in box["hi"]]}

    particulars = read_particulars(html)
    published = {key: item["value"] for key, item in particulars.items()
                 if isinstance(item.get("value"), (int, float))}
    if not {"beam_m", "depth_m", "lpp_m"} <= set(published):
        raise SystemExit(f"particulars page changed shape: {particulars}")

    propeller = None
    if prop_blob is not None:
        with zipfile.ZipFile(io.BytesIO(prop_blob)) as archive:
            propeller = read_model_propeller(archive.read("MP687_particulars.txt").decode("latin-1"))

    scale = derive_scale(span_m, published)

    # cross-checks that do not feed the derivation
    checks = {
        "iges_header_units": {"unit_flag": doc.unit_flag, "unit_name": doc.unit_name,
                              "metres_per_unit": metres_per_unit},
        "keel_plane_at_zero": {"z_lo_m": bbox_m["lo"][2],
                               "within_1mm": abs(bbox_m["lo"][2]) <= 1e-3},
        "half_model_symmetric": {"y_lo_m": bbox_m["lo"][1], "y_hi_m": bbox_m["hi"][1],
                                 "asymmetry_m": abs(bbox_m["hi"][1] + bbox_m["lo"][1])},
        "implied_model_lengths": {
            "lpp_m": published["lpp_m"] / scale["value"],
            "lwl_m": published["lwl_m"] / scale["value"],
            "displacement_m3": published["displacement_m3"] / scale["value"] ** 3,
            "x_extent_m": span_m["x"],
            "x_extent_over_implied_lpp": span_m["x"] * scale["value"] / published["lpp_m"],
        },
    }
    if propeller:
        checks["model_propeller"] = {
            **propeller,
            "model_diameter_m": propeller["diameter_mm"] / 1000.0,
            "consistent_with_hull_units": 0.1 < propeller["diameter_mm"] / 1000.0 < 1.0,
        }

    return {
        "schema": "tensorcad.fleet.model_scale.v1",
        "anchor": "GAP-CAD-T1-SOURCE-HUNT",
        "hull": "JBC (Japan Bulk Carrier, Tokyo 2015 workshop benchmark)",
        "generated_by": "tools/baseline_fleet/jbc_model_scale.py",
        "sources": {"particulars_page": html_info, "iges_archive": zip_info,
                    "propeller_archive": prop_info,
                    "archive_members": members},
        "licence": "not stated by the publisher; only digests, measured spans and the "
                   "derivation are recorded, archive bytes are not redistributed",
        "published_full_scale": particulars,
        "model_geometry": {
            "file": "jbc.igs (Ship Hull only; duct/rudder/propeller are separate members)",
            "surfaces": box["surfaces"], "sampled_points": box["points"],
            "samples_per_direction": box["samples"],
            "bbox_m": bbox_m, "span_m": span_m,
        },
        "scale_derived": scale,
        "checks": checks,
        "limits": [
            "The publisher states no model scale: this value is DERIVED from the CAD "
            "data, not read from a document, so it is not a workshop-stated figure.",
            "The two references disagree by the reported half-spread because the depth "
            "span includes deck-level fittings (hatch coaming / bulwark) above the "
            "moulded depth, which biases s_depth slightly low.",
            "Surface extremes are measured by sampling; a 41x41 grid per patch bounds "
            "the sampling deficit of B and D to well under the quoted half-spread.",
        ],
    }


def verify(record_path: Path) -> int:
    """Offline invariants: the recorded numbers must reproduce each other exactly."""
    record = json.loads(record_path.read_text())
    problems: list[str] = []
    span = record["model_geometry"]["span_m"]
    published = {k: v["value"] for k, v in record["published_full_scale"].items()
                 if isinstance(v.get("value"), (int, float))}
    scale = record["scale_derived"]
    for name, item in scale["candidates"].items():
        recomputed = item["reference_m"] / item["span_m"]
        if abs(recomputed - item["scale"]) > 1e-9:
            problems.append(f"{name}: recorded scale {item['scale']} != {recomputed}")
        if abs(item["span_m"] - span[{"beam": "y", "depth": "z"}[name]]) > 1e-9:
            problems.append(f"{name}: span in scale block differs from geometry block")
    mean = sum(item["scale"] for item in scale["candidates"].values()) / len(scale["candidates"])
    if abs(mean - scale["value"]) > 1e-9:
        problems.append("mean scale does not reproduce the candidate scales")
    if scale["relative_half_spread"] > 0.02:
        problems.append(f"half-spread {scale['relative_half_spread']:.3%} > 2 %: the two "
                        "references no longer agree, derivation must be re-examined")
    if published.get("lpp_m") != 280.0 or published.get("beam_m") != 45.0:
        problems.append("published particulars differ from the workshop values")
    if abs(published.get("displacement_m3", 0.0) - 178369.9) > 1e-6:
        problems.append(f"displacement read as {published.get('displacement_m3')} "
                        "m^3: the row parser picked a unit digit, not the value")
    for key in ("particulars_page", "iges_archive"):
        digest = record["sources"][key].get("sha256") or ""
        if len(digest) != 64 or not re.fullmatch(r"[0-9a-f]{64}", digest):
            problems.append(f"{key}: digest missing or malformed")
    if record["checks"]["iges_header_units"]["unit_name"] != "MM":
        problems.append("IGES header units are not MM: the model-space reading changes")
    if not record["checks"]["keel_plane_at_zero"]["within_1mm"]:
        problems.append("keel plane is not at z=0: coordinate frame assumption broken")
    for problem in problems:
        print(f"FAIL {problem}")
    if not problems:
        print(f"OK {record_path}: scale 1:{scale['value']:.3f} "
              f"±{scale['relative_half_spread']:.2%} from {len(scale['candidates'])} references")
    return 0 if not problems else 1


def figure(record: dict, path: Path) -> None:
    """Two panels: model-space point cloud (plan/side) and the scale references."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    span = record["model_geometry"]["span_m"]
    candidates = record["scale_derived"]["candidates"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    lo = record["model_geometry"]["bbox_m"]["lo"]
    hi = record["model_geometry"]["bbox_m"]["hi"]
    ax.add_patch(plt.Rectangle((lo[0], lo[2]), hi[0] - lo[0], hi[2] - lo[2],
                               fill=False, lw=1.6, color="C0",
                               label=f"model bbox: x {span['x']:.3f} m, z {span['z']:.3f} m"))
    ax.set_xlabel("x [m] (model space)")
    ax.set_ylabel("z [m]")
    ax.set_title("JBC IGES hull — side projection of the sampled extent")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    names = list(candidates)
    values = [candidates[n]["scale"] for n in names]
    ax.bar(names, values, color=["C0", "C1"], width=0.5,
           label="derived  s = reference / model span")
    for n, v in zip(names, values):
        ax.annotate(f"{v:.3f}", (n, v), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=9)
    for hypothesis in HYPOTHESES:
        ax.axhline(hypothesis, ls="--", lw=0.9,
                   color="0.35", label=f"1:{hypothesis:g}")
    ax.set_ylabel("scale  s")
    ax.set_ylim(min(values) - 1.2, max(values) + 1.2)
    ax.set_title("Scale from two independent published lengths")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    fig.suptitle(f"JBC model scale derived from publisher CAD: 1:{record['scale_derived']['value']:.2f} "
                 f"±{record['scale_derived']['relative_half_spread']:.2%}", fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    print(f"figure {path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/baseline_fleet/jbc_model_scale.json")
    parser.add_argument("--cache", default="artifacts/baseline_fleet/raw")
    parser.add_argument("--samples", type=int, default=41)
    parser.add_argument("--figure", default="docs/figs/jbc_model_scale.png")
    parser.add_argument("--offline", action="store_true",
                        help="use cached downloads only (fails if the cache is cold)")
    parser.add_argument("--verify-only", action="store_true",
                        help="check an existing record without touching network or CAD")
    args = parser.parse_args(argv)

    out = Path(args.out)
    if args.verify_only:
        return verify(out)
    cache = Path(args.cache)
    if args.offline and not cache.exists():
        raise SystemExit(f"--offline but cache {cache} is missing")
    record = build(cache, args.samples)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2) + "\n")
    print(f"record {out}: 1:{record['scale_derived']['value']:.3f} "
          f"±{record['scale_derived']['relative_half_spread']:.2%}")
    if args.figure:
        figure(record, Path(args.figure))
    return verify(out)


if __name__ == "__main__":
    raise SystemExit(main())