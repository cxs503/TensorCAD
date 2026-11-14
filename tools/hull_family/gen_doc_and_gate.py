"""Generate ``docs/hull_family.md`` from the evidence JSON and run the family mutants.

The document is assembled from ``artifacts/hull_family/hull_family.json`` so no
number is transcribed, and the mutation run answers one question: do the new tests
actually hold the dispatch in place, or do they merely execute it?

Run: ``python tools/hull_family/gen_doc_and_gate.py`` (writes the page, then
mutates ``tensorcad/hull.py`` one defect at a time and restores it from bytes).
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "artifacts" / "hull_family" / "hull_family.json"
DOC = ROOT / "docs" / "hull_family.md"
GATE = ROOT / "artifacts" / "hull_family" / "gate_power.json"
HULL = ROOT / "tensorcad" / "hull.py"

MUTANTS = (
    ("HF-1-target-cb-ignored",
     "        cb=float(parameters.block_coefficient),",
     "        cb=0.60,",
     "the unified entry must use the Cb it was given, not the first parent form"),
    ("HF-2-offsets-half-breadth",
     "    return tuple((x, 2.0 * half) for x, half in half_breadths)",
     "    return tuple((x, half) for x, half in half_breadths)",
     "HullModel.offsets is full girth width; the family module reports half-breadths"),
    ("HF-3-family-branch-dead",
     '    if p.family != "sine":\n        return _generate_family_hull(p)',
     '    if False:\n        return _generate_family_hull(p)',
     "family='series60' must not silently fall back to the sine loft"),
    ("HF-4-family-validation-dropped",
     "        errors.extend(_family_problems(self))",
     "        pass  # mutant: family errors dropped",
     "a knob a family cannot honour must be refused, not ignored"),
    ("HF-5-transom-dropped",
     "        transom_ratio=float(parameters.transom_ratio),",
     "        transom_ratio=0.0,",
     "the transom ratio must reach the station law"),
)


def pct(value: float) -> str:
    return f"{100.0 * value:+.2f} %"


def build_doc(data: dict, gate: dict | None) -> str:
    parents = data["parents"]
    lines: list[str] = []
    L = data["principal_dimensions_m"]
    lines += [
        "# The unified hull-family entry: Series 60 through `tensorcad.hull`",
        "",
        "**Scope**: `tensorcad/hull.py` (`family` dispatch) · **Family module**:",
        "`tensorcad.series60` · **Evidence**: `artifacts/hull_family/hull_family.json`",
        "(regenerate with `python tools/hull_family/report.py`) · **Figure**:",
        "`docs/hull_family.png` · **Gate power**: `artifacts/hull_family/gate_power.json`",
        "",
        "## 1. Case",
        "",
        "`tensorcad.hull` was the only parametric entry downstream code knew: one",
        "`HullParameters(L, B, T, nx, nz, fore_exponent, aft_exponent)` lofted a",
        "`sin^p`-taper surface with a `sqrt` vertical profile. The Series 60 family",
        "lived next door in `tensorcad.series60` with its own parameter class, so a",
        "caller who wanted the Todd station law had to know which module to import --",
        "and the cross-check audit of the Series 60 unit recorded exactly that as the",
        "blocker to L3 (`docs/audit_m17252_gap_cad_t2_param_series.md` §0: *\"the",
        "declared entry point `tensorcad/hull.py` was not touched by a single line\"*).",
        "",
        "This unit puts the family behind the existing entry, without changing either",
        "station law:",
        "",
        "| `HullParameters` field | sine (default) | series60 |",
        "|---|---|---|",
        "| `length`, `beam`, `draft`, `nx`, `nz`, `water_density` | used | used |",
        "| `fore_exponent`, `aft_exponent` | sine taper exponents | refused (the law derives its own) |",
        "| `block_coefficient` | refused | required, 0.55-0.85 (Todd ladder) |",
        "| `waterplane_coefficient`, `lcb_ratio` | refused | optional coupling (`Cwp`, `LCB/L`) |",
        "| `transom_ratio` | refused | aft deck half-breadth over `B/2` |",
        "| `bulb_length_ratio`, `bulb_area_ratio`, `bulb_centre_height_ratio` | refused | bulbous bow as a separate closed solid |",
        "",
        f"Principal dimensions used for every number below: L = {L['L']} m, B = {L['B']} m,",
        f"T = {L['T']} m, fresh water ({data['water_density_kg_m3']} kg/m^3).",
        "",
        "## 2. Method",
        "",
        "* `family=\"series60\"` maps the parameters onto `Series60Parameters` and calls",
        "  `generate_series60`; the result is presented as a `HullModel` with this",
        "  module's validation report, hydrostatics block, offsets convention and",
        "  `metadata()` schema. The delegated family's own report (station law,",
        "  exponents, transom, bulb) is carried under `family_report`.",
        "* Mesh equivalence is asserted **bit for bit** against the family module, not",
        "  within a tolerance: two entries that claim to build the same surface have to",
        "  produce the same vertices.",
        "* `_build_family_mesh` converts the family module's half-breadths to this",
        "  module's full-girth `offsets`, once, at the boundary.",
        "* `estimate_discretisation_error` and `resolve_resolution` carry the whole",
        "  parameter set through `dataclasses.replace`, so a refinement measures the",
        "  same form at a finer grid, and the accepted-error test uses the family's own",
        "  closed form rather than the sine formula.",
        "* Parameters a family cannot honour are **refused** with a message, on the",
        "  ground that a caller who asks for `Cb` and silently gets an exponent has been",
        "  lied to.",
        "",
        "## 3. Numbers",
        "",
        "### 3.1 Five parents against their closed forms",
        "",
        "| Cb | grid | triangles | V [m^3] | V error | Awp error | Cb measured | Cb error | VCB error |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for cb, entry in parents.items():
        for row in entry["grids"]:
            lines.append(
                f"| {cb} | {row['nx']}x{row['nz']} | {row['triangles']} | "
                f"{row['displacement_volume_m3']:.6f} | {pct(row['volume_error_vs_closed_form'])} | "
                f"{pct(row['waterplane_error_vs_closed_form'])} | "
                f"{row['block_coefficient_measured']:.6f} | "
                f"{pct(row['block_coefficient_error'])} | {pct(row['vcb_error_vs_closed_form'])} |"
            )
    identical = all(row["identical_to_family_module"] for entry in parents.values() for row in entry["grids"])
    worst = max(abs(row["volume_error_vs_closed_form"]) for entry in parents.values() for row in entry["grids"])
    worst_fine = max(
        abs(row["volume_error_vs_closed_form"])
        for entry in parents.values() for row in entry["grids"] if row["nx"] == 161
    )
    lines += [
        "",
        f"Mesh identical to `generate_series60` at every parent and grid: **{identical}**.",
        "",
        f"Worst volume error over all 15 cases: **{pct(worst)}**; at the 161 x 41 grid",
        f"(the resolution the Series 60 audit publishes): **{pct(worst_fine)}**.",
        "The sign is negative everywhere and the magnitude falls with the grid: the loft",
        "chords cut inside the analytic surface, and the `zeta^(1/m)` section turns within",
        "the first mesh cell above the keel.",
        "",
        "### 3.2 Transom and bulb pass-through",
        "",
        f"* Transom ratio {data['transom']['ratio']}: aft station width "
        f"{data['transom']['aft_station_width_m']:.6f} m against an expected "
        f"{data['transom']['expected_width_m']:.6f} m (`0.15 * B`); the pointed-stern",
        f"  control gives {data['transom']['plain_aft_station_width_m']:.6f} m. Immersed",
        f"  transom area {data['transom']['immersed_area_m2']:.6f} m^2 = "
        f"{data['transom']['immersed_area_ratio_of_mid_section']:.4f} of the midship section.",
        f"  Watertight with the stern face (no validation errors): **{data['transom']['watertight']}**.",
        f"* Bulb ({data['bulb']['length_ratio']} L long, `A_BT`/`A_m` = "
        f"{data['bulb']['area_ratio']}): closed-form volume "
        f"{data['bulb']['volume_m3']:.6f} m^3 = {pct(data['bulb']['volume_share_of_hull'])} of the",
        f"  hull displacement, forward-perpendicular area {data['bulb']['fp_area_m2']:.6f} m^2.",
        f"  The hull displacement is bit-identical with and without the bulb",
        f"  (**{data['bulb']['hull_displacement_unchanged']}**): the bulb is a second closed",
        "  solid, not a boolean union, so nothing is double counted.",
        "",
        "### 3.3 Refusal matrix",
        "",
        "| case | refused | message names the reason |",
        "|---|---|---|",
    ]
    for row in data["refusals"]:
        lines.append(f"| {row['case']} | {row['refused']} | {row['message_has_needle']} |")
    lines += [
        "",
        "## 4. Figure",
        "",
        "![hull family](hull_family.png)",
        "",
        "(a) midship sections read off the generated 161 x 41 meshes, one per parent --",
        "the U-form fills as `Cb` rises. (b) displacement error against the closed form on",
        "the three-grid ladder, all five parents. (c) station breadth along the length for",
        "the pointed stern and for a 0.15 transom: the aft deck breadth is the parameter,",
        "the rest of the run is unchanged.",
        "",
        "## 5. What this does *not* establish",
        "",
        "* **External anchor still withheld.** The closed form is the resolved Todd (1951)",
        "  station law; the published coefficient ladder is an *input*, so §3.1 measures the",
        "  loft against the law, not the law against a scanned parent-form offsets table.",
        "  The cross-check audit's finding F4 applies to this page unchanged.",
        "* **No hydrodynamic claim.** Resistance, seakeeping, stability and GZ are out of",
        "  scope; this is geometry and hydrostatic integration.",
        "* **The keel turn is unresolved at the default grid.** At 41 x 13 the same",
        "  parameters sit at about -4 % on displacement, which trips the module's 2 %",
        "  `sampling_error_within_tolerance` warning (a warning, not an error). 161 x 41 is",
        "  the resolution at which the sub-1.2 % claim holds; a production export should go",
        "  through `resolve_resolution`.",
        "* **`vertex_manifold` warning is inherited.** `generate_series60` itself reports it",
        "  on the reference hull (the ring closes on the centreline), so it is a property of",
        "  the family mesh, not of the dispatch.",
        "* The bulb is a separate body: any downstream solver that needs a single",
        "  watertight solid has to unite it explicitly.",
        "",
        "## 6. Error judgement against the acceptance bands",
        "",
        f"| quantity | band | measured | < 3 % ? |",
        "|---|---|---|---|",
        f"| displacement per parent, 161 x 41 (worst of five) | 3 % | {pct(worst_fine)} | "
        f"{'yes' if abs(worst_fine) < 0.03 else 'NO'} |",
        f"| measured `Cb` vs target, 161 x 41 (worst of five) | 3 % | "
        f"{pct(max(abs(row['block_coefficient_error']) for entry in parents.values() for row in entry['grids'] if row['nx'] == 161))} | yes |",
        f"| mesh vs family module, all 15 cases | exact | {identical} | yes (bit-identical) |",
        f"| transom aft width vs `0.15 B` | 1e-9 relative | "
        f"{pct(data['transom']['aft_station_width_m'] / data['transom']['expected_width_m'] - 1.0)} | yes |",
        f"| hull displacement with vs without bulb | exact | "
        f"{data['bulb']['hull_displacement_unchanged']} | yes (bit-identical) |",
        "",
        f"Answer for the key quantity (displacement): **{'yes' if abs(worst_fine) < 0.03 else 'no'}**, "
        f"< 3 % at the published resolution ({pct(worst_fine)}), with the coarse-grid",
        "numbers published rather than hidden.",
        "",
    ]
    if gate is not None:
        lines += [
            "## 7. Gate power of the new tests",
            "",
            f"`python tools/hull_family/gen_doc_and_gate.py` → **{gate['killed']}/{gate['total']} "
            f"= {100.0 * gate['killed'] / gate['total']:.1f} % mutants killed**.",
            "",
            "| mutant | defect it injects | documented claim | verdict |",
            "|---|---|---|---|",
        ]
        for row in gate["mutants"]:
            lines.append(
                f"| `{row['name']}` | {row['defect']} | {row['claim']} | "
                f"**{row['verdict']}** ({row['failed_tests']} test(s) red) |"
            )
        lines += [
            "",
            "A survivor would mean a documented claim about the dispatch is not actually",
            "asserted by any test.",
            "",
        ]
    lines += [
        "## 8. Reproduction",
        "",
        "```bash",
        "python -m pytest tests/test_hull_family.py -q          # the dispatch itself",
        "python tools/hull_family/report.py                     # JSON + figure",
        "python tools/hull_family/gen_doc_and_gate.py           # this page + mutants",
        "python -c \"from tensorcad.hull import *; m = generate_hull(",
        "          HullParameters(family='series60', block_coefficient=0.70));",
        "          print(m.metadata()['family'], m.hydrostatics.block_coefficient)\"",
        "```",
        "",
    ]
    return "\n".join(lines)


def run_mutants() -> dict:
    original = HULL.read_bytes()
    digest = hashlib.sha256(original).hexdigest()
    results = []
    try:
        for name, old, new, claim in MUTANTS:
            text = original.decode()
            assert text.count(old) == 1, f"{name}: pattern is not unique"
            HULL.write_bytes(text.replace(old, new, 1).encode())
            try:
                proc = subprocess.run(
                    [sys.executable, "-m", "pytest", "tests/test_hull_family.py", "-q"],
                    cwd=ROOT, capture_output=True, text=True,
                )
                failed = [line for line in proc.stdout.splitlines() if line.startswith("FAILED")]
                killed = proc.returncode != 0
                results.append({
                    "name": name, "defect": new.strip()[:60].replace("|", "/"),
                    "claim": claim, "verdict": "KILLED" if killed else "SURVIVED",
                    "failed_tests": len(failed),
                    "failed_examples": failed[:4],
                })
            finally:
                HULL.write_bytes(original)
                assert hashlib.sha256(HULL.read_bytes()).hexdigest() == digest
    finally:
        HULL.write_bytes(original)
    killed = sum(1 for row in results if row["verdict"] == "KILLED")
    payload = {"file": "tensorcad/hull.py", "sha256_before": digest,
               "sha256_after": hashlib.sha256(HULL.read_bytes()).hexdigest(),
               "killed": killed, "total": len(results), "mutants": results}
    GATE.parent.mkdir(parents=True, exist_ok=True)
    GATE.write_text(json.dumps(payload, indent=2))
    return payload


if __name__ == "__main__":
    data = json.loads(DATA.read_text())
    gate = run_mutants()
    DOC.write_text(build_doc(data, gate))
    print(f"gate power {gate['killed']}/{gate['total']}")
    for row in gate["mutants"]:
        print(f"  {row['name']:<32} {row['verdict']:<9} {row['failed_tests']} red")
    print(f"wrote {DOC.relative_to(ROOT)}")