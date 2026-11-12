#!/usr/bin/env python3
"""Can the *mesh* route name the facet-scale defects the BRep route names?

The 2026-11-09 cross-route audit (``docs/audit_crosscheck.md``) measured that four of the
geomquality defect families -- ``SLIVER``, ``MICRO_EDGE``, ``CURVATURE_JUMP``,
``TANGENT_JUMP`` -- all collapsed into one mesh-side observation, "the shell is open", and
recorded the gap as *"口径缺失（如实记录，未粉饰）"*. This example answers the follow-up with
numbers: over the same 3-hull x 14-injection corpus, which of those families does
:mod:`tensorcad.geomquality.mesh_defects` actually name, and which does it refuse to name?

Writes ``results.json``, ``report.md`` and ``naming.png`` into ``--output`` and exits
non-zero if a pinned claim stops holding.

    python examples/mesh_route_naming/run_example.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geomquality import (  # noqa: E402
    baselines,
    defects,
    name_mesh_defects,
    repair_facet_defects,
)

#: The four families the cross-route audit could not name, plus the clean reference.
FAMILIES = ("CLEAN", "SLIVER", "MICRO_EDGE", "TANGENT_JUMP", "CURVATURE_JUMP", "ROUGHNESS")


def collect() -> dict:
    catalogue = baselines.specimens()
    cases = {case.code: case for case in defects.DEFECT_CASES}
    rows = []
    for name in sorted(catalogue):
        specimen = catalogue[name]()
        for family in FAMILIES:
            brep = specimen.brep if family == "CLEAN" else cases[family].inject(specimen.brep)
            mesh = brep.tessellate(16, 8)
            report = name_mesh_defects(mesh, declared_length=specimen.length)
            row = {
                "specimen": name,
                "family": family,
                "codes": list(report.codes),
                "faces": report.metrics["faces"],
                "components": report.metrics["components"],
                "micro_edge_share": report.metrics["micro_edge_share"],
                "sliver_share": report.metrics["sliver_share"],
                "open_edges": report.metrics["open_edges"],
            }
            if family in ("SLIVER", "MICRO_EDGE"):
                repair = repair_facet_defects(mesh, declared_length=specimen.length)
                row["repair"] = repair.to_dict()
            rows.append(row)
    return {"rows": rows, "specimens": sorted(catalogue), "families": list(FAMILIES)}


def claims(payload: dict) -> list[dict]:
    rows = payload["rows"]
    out = []

    def claim(name: str, ok: bool, detail: str) -> None:
        out.append({"claim": name, "passed": bool(ok), "detail": detail})

    def of(family: str) -> list[dict]:
        return [row for row in rows if row["family"] == family]

    clean = of("CLEAN")
    claim(
        "clean_references_are_named_clean",
        all(row["codes"] == [] for row in clean),
        "codes per clean hull: " + str({row["specimen"]: row["codes"] for row in clean}),
    )
    sliver = of("SLIVER")
    claim(
        "sliver_is_named",
        all("SLIVER" in row["codes"] for row in sliver),
        "codes: " + str({row["specimen"]: row["codes"] for row in sliver})
        + f"; sliver share {min(r['sliver_share'] for r in sliver):.3f}-"
        f"{max(r['sliver_share'] for r in sliver):.3f} against a 0.02 limit, clean "
        f"max {max(r['sliver_share'] for r in clean):.3f}",
    )
    micro = of("MICRO_EDGE")
    claim(
        "micro_edge_is_named",
        all("MICRO_EDGE" in row["codes"] for row in micro),
        "codes: " + str({row["specimen"]: row["codes"] for row in micro})
        + f"; micro-edge share {min(r['micro_edge_share'] for r in micro):.3f}-"
        f"{max(r['micro_edge_share'] for r in micro):.3f} against a 0.01 limit, clean "
        f"max {max(r['micro_edge_share'] for r in clean):.3f}",
    )
    unnamed = [row for row in rows if row["family"] in ("TANGENT_JUMP", "CURVATURE_JUMP", "ROUGHNESS")]
    claim(
        "seam_and_ripple_families_are_not_overclaimed",
        all(set(row["codes"]) <= {"ISOLATED_PATCH"} for row in unnamed),
        "codes: "
        + str({f"{row['specimen']}/{row['family']}": row["codes"] for row in unnamed}),
    )
    repairs = [row["repair"] for row in rows if "repair" in row and row["repair"]["accepted"]]
    claim(
        "an_accepted_repair_never_leaves_a_facet_code",
        all(repair["codes_after"] == [] for repair in repairs),
        f"{len(repairs)} of 6 repairs accepted; codes_after: "
        + str([repair["codes_after"] for repair in repairs]),
    )
    refusals = [row["repair"] for row in rows if "repair" in row and not row["repair"]["accepted"]]
    claim(
        "a_refused_repair_names_the_number_that_tripped_it",
        all(("moved" in r["reason"]) or ("open edges" in r["reason"]) or ("no facet-scale" in r["reason"])
            or ("every face" in r["reason"]) for r in refusals),
        f"{len(refusals)} refusals, e.g. " + (refusals[0]["reason"] if refusals else "none"),
    )
    return out


def plot(payload: dict, output: Path) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - environment dependent
        return [f"no figure: matplotlib unavailable ({exc.__class__.__name__})"]
    families = payload["families"]
    named = []
    for family in families:
        rows = [row for row in payload["rows"] if row["family"] == family]
        hits = sum(1 for row in rows if row["codes"])
        named.append(hits)
    figure, axis = plt.subplots(figsize=(7.0, 3.4))
    axis.bar(range(len(families)), named, color="#2b6cb0")
    axis.set_xticks(range(len(families)))
    axis.set_xticklabels(families, rotation=20, ha="right")
    axis.set_ylabel("hulls with >= 1 code (of 3)")
    axis.set_title("Mesh-route naming on the geomquality corpus")
    axis.set_ylim(0, 3.4)
    axis.grid(axis="y", alpha=0.3)
    figure.tight_layout()
    path = output / "naming.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return ["naming.png"]


def write_report(payload: dict, verdicts: list[dict], output: Path, figures: list[str]) -> str:
    rows = payload["rows"]
    lines = [
        "# Mesh-route naming on the geomquality corpus",
        "",
        "Generated by `examples/mesh_route_naming/run_example.py`; every number below is read",
        "from the same run's `results.json`.",
        "",
        "| checkpoint | " + " | ".join(payload["families"]) + " |",
        "|---|" + "---|" * len(payload["families"]),
    ]
    for specimen in payload["specimens"]:
        cells = []
        for family in payload["families"]:
            row = next(r for r in rows if r["specimen"] == specimen and r["family"] == family)
            cells.append(", ".join(row["codes"]) or "-")
        lines.append(f"| {specimen} | " + " | ".join(cells) + " |")
    lines += ["", "## Claims", ""]
    for item in verdicts:
        lines.append(f"- {'PASS' if item['passed'] else 'FAIL'} `{item['claim']}` — {item['detail']}")
    lines += ["", "## Figures", ""] + [f"- `{name}`" for name in figures]
    text = "\n".join(lines) + "\n"
    (output / "report.md").write_text(text, encoding="utf-8")
    return text


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="artifacts/mesh_route_naming")
    options = parser.parse_args(argv)
    output = Path(options.output)
    output.mkdir(parents=True, exist_ok=True)

    payload = collect()
    verdicts = claims(payload)
    figures = plot(payload, output)
    payload["claims"] = verdicts
    payload["figures"] = figures
    (output / "results.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    write_report(payload, verdicts, output, figures)

    failed = [item for item in verdicts if not item["passed"]]
    for item in verdicts:
        print(("PASS " if item["passed"] else "FAIL ") + item["claim"] + " :: " + item["detail"])
    print(f"wrote {output}/results.json, report.md, figures={figures}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
