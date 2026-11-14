"""Product entry for the Series 60 parametric hull family.

``tensorcad.series60`` is the kernel; this module is the *entry* on it, so the
family can be generated, exported and judged without importing Python::

    tensorcad-hull series60 parent --cb 0.70 --cwp 0.753 --lcb 0.490 --gate --stl
    tensorcad-hull series60 family --published docs/data/series60_published_parents.json

Why an entry at all: the same audit that certified the kernel's closed forms
repeatedly found the capability unreachable from any product surface (not
exported, not in ``tensorcad.cli``, not in ``hull.py``), which is the difference
between "a module the tests import" and "a capability". The entry is thin on
purpose -- it computes nothing the library does not already expose, and it
carries no constants of its own (every reference is recomputed in-process from
:func:`tensorcad.series60.analytic_series60_hydrostatics`).

Exit codes follow the package convention: ``0`` when every claim the entry makes
holds, ``1`` when a model/gate check fails, ``2`` on an invalid invocation.

Artifact policy: ``--out DIR`` writes what a consumer can re-read
(``metadata.json``, ``offsets.csv``, ``hull.obj``, optionally ``hull.stl``); the
run artifacts under ``runs/`` are produced by this entry, so no number in a
report is typed by hand.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

from .export import write_obj, write_stl
from .gate import gate_mesh
from .series60 import (
    PARENT_CB,
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

#: Predeclared bar for "measured mesh vs the module's own closed form".  The hull
#: is integrated on the same grid the surface is lofted on, so this is a
#: discretisation error and 2 % is loose enough to be a claim with teeth.
DISCRETISATION_TOLERANCE = 0.02

#: Predeclared bar for "generated coefficient vs the published parent form".
PUBLISHED_TOLERANCE = 0.03

#: Measured-on-the-mesh quantity -> key of ``analytic_series60_hydrostatics``.
#: ``None`` marks a coefficient that has to be assembled from two measurements
#: (there is no single attribute for it on :class:`tensorcad.hydrostatics.Hydrostatics`).
MEASURED_KEYS: tuple[tuple[str, str], ...] = (
    ("displacement_volume", "displacement_volume"),
    ("waterplane_area", "waterplane_area"),
    ("wetted_surface_area", "wetted_surface_area"),
    ("longitudinal_centre_buoyancy", "longitudinal_centre_buoyancy"),
    ("longitudinal_centre_flotation", "longitudinal_centre_flotation"),
    # The mesh integrates in the hull frame (z = 0 at the waterline, keel at
    # z = -T); the analytic side carries both, so compare like with like.  Pairing
    # the z-frame measurement with the above-keel closed form reports -193 % on
    # every hull, which is a frame mismatch, not a discretisation error.
    ("vertical_centre_buoyancy", "vertical_centre_buoyancy_z"),
    ("block_coefficient", "block_coefficient"),
    ("prismatic_coefficient", "prismatic_coefficient"),
    ("maximum_section_area", "maximum_section_area"),
    # Cm and Cwp are quotients of two measured quantities; the analytic side has
    # them as first-class keys, so they are compared the same way.
    ("midship_coefficient", "midship_coefficient"),
    ("waterplane_coefficient", "waterplane_coefficient"),
)


def _parameters(args: argparse.Namespace) -> Series60Parameters:
    return Series60Parameters(
        length=args.length,
        beam=args.beam,
        draft=args.draft,
        nx=args.nx,
        nz=args.nz,
        water_density=args.density,
        cb=args.cb,
        cwp=args.cwp,
        lcb=args.lcb,
        transom_ratio=args.transom_ratio,
        bulb_length_ratio=args.bulb_length_ratio,
        bulb_area_ratio=args.bulb_area_ratio,
        bulb_centre_height_ratio=args.bulb_centre_height_ratio,
    )


def measured_hydrostatics(model: Any) -> dict[str, float]:
    """The measured record, with the two coefficients that are quotients added."""
    hydro = model.hydrostatics
    p = model.parameters
    values: dict[str, float] = {}
    derived = ("midship_coefficient", "waterplane_coefficient")
    for name, _ in MEASURED_KEYS:
        if name in derived:
            continue
        values[name] = float(getattr(hydro, name))
    values["midship_coefficient"] = (
        hydro.maximum_section_area / (p.beam * p.draft) if p.beam * p.draft else float("nan")
    )
    values["waterplane_coefficient"] = (
        hydro.waterplane_area / (p.length * p.beam) if p.length * p.beam else float("nan")
    )
    return values


def relative_errors(model: Any) -> dict[str, float | None]:
    """Measured mesh against the module's analytic reference, per quantity."""
    analytic = analytic_series60_hydrostatics(model.parameters)
    measured = measured_hydrostatics(model)
    errors: dict[str, float | None] = {}
    for name, key in MEASURED_KEYS:
        reference = analytic.get(key)
        value = measured.get(name)
        if reference is None or value is None or not math.isfinite(float(value)):
            errors[name] = None
            continue
        errors[name] = (float(value) - float(reference)) / float(reference) if reference else None
    return errors


def _write_offsets(model: Any, path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_m", "half_breadth_m"])
        for x, half_breadth in model.offsets:
            writer.writerow([f"{x:.12g}", f"{half_breadth:.12g}"])
    return path


def _export(model: Any, output: Path, *, stl: bool, document: dict[str, Any]) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    artifacts = {
        "mesh": str(write_obj(model.mesh, output / "hull.obj")),
        "offsets": str(_write_offsets(model, output / "offsets.csv")),
    }
    if stl:
        artifacts["stl"] = str(write_stl(model.mesh, output / "hull.stl"))
    metadata_path = output / "metadata.json"
    metadata_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    artifacts["metadata"] = str(metadata_path)
    return artifacts


def _gate(model: Any) -> dict[str, Any]:
    """The mesh gate, with the reference volume coming from outside the mesh."""
    analytic = analytic_series60_hydrostatics(model.parameters)
    report = gate_mesh(
        model.mesh,
        reference_volume=analytic["displacement_volume"],
        volume_tolerance=DISCRETISATION_TOLERANCE,
    )
    return report.to_dict()


def _parent_record(args: argparse.Namespace) -> tuple[int, dict[str, Any]]:
    p = _parameters(args)
    model = generate_series60(p)
    analytic = analytic_series60_hydrostatics(p)
    errors = relative_errors(model)
    gate = _gate(model) if args.gate else None
    document = model.metadata()
    document["analytic_reference"] = {key: value for key, value in analytic.items()}
    document["analytic_relative_errors"] = errors
    document["tolerance"] = {"discretisation": DISCRETISATION_TOLERANCE,
                             "published": PUBLISHED_TOLERANCE}
    if gate is not None:
        document["gate"] = gate

    finite = {key: value for key, value in errors.items() if value is not None}
    worst_key = max(finite, key=lambda key: abs(finite[key])) if finite else None
    document["discretisation"] = {
        "grid": [p.nx, p.nz],
        "worst_quantity": worst_key,
        "worst_relative_error": finite.get(worst_key) if worst_key else None,
        "within_tolerance": bool(finite) and all(
            abs(value) <= DISCRETISATION_TOLERANCE for value in finite.values()
        ),
    }

    artifacts: dict[str, str] = {}
    refused = gate is not None and not gate["passed"]
    if refused:
        # The gate's whole point is to refuse: a failed gate that still writes the
        # artifacts leaves a consumer with a file it cannot trust.
        document["artifacts"] = {}
        document["refused_by_gate"] = list(gate["failures"])
        return 1, document
    if not args.no_export:
        try:
            artifacts = _export(model, Path(args.output), stl=args.stl, document=document)
        except Exception as exc:  # noqa: BLE001 - the entry contract is an exit code
            print(f"tensorcad-hull series60: export failed: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            return 2, {"error": str(exc), "stage": "export"}
    document["artifacts"] = artifacts

    ok = model.validation.passed and bool(document["discretisation"]["within_tolerance"])
    if gate is not None:
        ok = ok and bool(gate["passed"])
    return (0 if ok else 1), document


def load_published(path: str | Path | None) -> dict[float, dict[str, float]]:
    """Published parent particulars keyed by Cb; ``{}`` when no table was given."""
    if path is None:
        return {}
    data = json.loads(Path(path).read_text())
    return {round(float(row["cb"]), 6): row for row in data["parents"]}


def _role(name: str) -> str:
    """What a comparison of *name* against the published table actually proves.

    The distinction matters and the earlier round got it wrong: ``cwp`` and
    ``lcb`` are handed to the generator as targets, so agreement is a **closure**
    (it proves the inverse solve converged, not that the hull is a Series 60);
    ``cp`` is the published ladder the exponents are read from, so its agreement
    is the ladder being read, not a prediction; only quantities that are *not*
    inputs can be called a comparison.
    """
    if name in ("waterplane_coefficient", "lcb_over_length"):
        return "input target (closure check)"
    if name in ("prismatic_coefficient", "midship_coefficient"):
        return "published ladder (input, not a prediction)"
    return "derived"


def _family_row(cb: float, published: dict[float, dict[str, float]],
                nx: int, nz: int, water_density: float) -> dict[str, Any]:
    published_row = published.get(round(cb, 6), {})
    # A published parent carries its own Cwp/LCB; a parent generated *from* the
    # table is the only honest way to ask "does the family reproduce it", so the
    # sweep drives the generator with the published targets when they exist.
    kwargs: dict[str, Any] = {"length": 10.0, "beam": 1.5, "draft": 0.6,
                              "nx": nx, "nz": nz, "water_density": water_density, "cb": cb}
    if published_row:
        kwargs["cwp"] = float(published_row["cwp"])
        kwargs["lcb"] = float(published_row["lcb_over_length"])
    p = Series60Parameters(**kwargs)
    model = generate_series60(p)
    analytic = analytic_series60_hydrostatics(p)
    measured = measured_hydrostatics(model)
    row: dict[str, Any] = {
        "cb_requested": cb,
        "grid": [nx, nz],
        "parameters": p.as_dict(),
        "closed_form": {
            "block_coefficient": analytic["block_coefficient"],
            "prismatic_coefficient": analytic["prismatic_coefficient"],
            "midship_coefficient": analytic["midship_coefficient"],
            "waterplane_coefficient": analytic["waterplane_coefficient"],
            "lcb_over_length": analytic["longitudinal_centre_buoyancy"] / p.length,
            "vcb_over_draft": analytic["vertical_centre_buoyancy"] / p.draft,
            "displacement_volume": analytic["displacement_volume"],
        },
        "measured": {
            "block_coefficient": measured["block_coefficient"],
            "prismatic_coefficient": measured["prismatic_coefficient"],
            "midship_coefficient": measured["midship_coefficient"],
            "waterplane_coefficient": measured["waterplane_coefficient"],
            "lcb_over_length": measured["longitudinal_centre_buoyancy"] / p.length,
            "vcb_over_draft": measured["vertical_centre_buoyancy"] / p.draft,
            "displacement_volume": measured["displacement_volume"],
        },
        "error_vs_published": {},
        "discretisation_relative_errors": relative_errors(model),
        "validation": model.validation.to_dict(),
    }
    if published_row:
        for name, key in (("block_coefficient", "cb"),
                          ("midship_coefficient", "cm"),
                          ("prismatic_coefficient", "cp"),
                          ("waterplane_coefficient", "cwp"),
                          ("lcb_over_length", "lcb_over_length")):
            reference = float(published_row[key])
            for source in ("closed_form", "measured"):
                got = float(row[source][name])
                row["error_vs_published"].setdefault(source, {})[name] = {
                    "published": reference,
                    "got": got,
                    "relative": (got - reference) / reference if reference else None,
                    "role": _role(name),
                }
    return row


def _family_record(args: argparse.Namespace) -> tuple[int, dict[str, Any]]:
    published = load_published(args.published)
    rows = [
        _family_row(cb, published, args.nx, args.nz, args.density)
        for cb in (args.cb_values or PARENT_CB)
    ]
    worst = 0.0
    worst_at = None
    for row in rows:
        for source, entries in row["error_vs_published"].items():
            for name, entry in entries.items():
                value = entry.get("relative")
                if value is None:
                    continue
                if abs(value) > abs(worst):
                    worst, worst_at = value, f"Cb={row['cb_requested']:.2f} {source} {name}"
    discretisation_worst = 0.0
    for row in rows:
        for name, value in row["discretisation_relative_errors"].items():
            if value is None:
                continue
            discretisation_worst = max(discretisation_worst, abs(value))
    document: dict[str, Any] = {
        "schema": "tensorcad.series60_family/2.0",
        "case": {
            "family": "Series 60 parent forms (Todd 1951)",
            "length_m": 10.0, "beam_m": 1.5, "draft_m": 0.6,
            "lb": 10.0 / 1.5, "bt": 1.5 / 0.6,
            "cb_values": [row["cb_requested"] for row in rows],
            "grid": [args.nx, args.nz],
            "water_density": args.density,
        },
        "published_source": str(args.published) if args.published else None,
        "tolerance": {"published": PUBLISHED_TOLERANCE,
                      "discretisation": DISCRETISATION_TOLERANCE},
        "role_legend": {
            "input target (closure check)": "handed to the generator as a target; agreement "
                                            "proves the inverse solve converged",
            "published ladder (input, not a prediction)": "the published Cp read to pin the "
                                                          "exponents; agreement is transcription",
            "derived": "not a target",
        },
        "rows": rows,
        "summary": {
            "worst_published_relative": worst,
            "worst_published_at": worst_at,
            "worst_discretisation_relative": discretisation_worst,
            "within_published_tolerance": abs(worst) <= PUBLISHED_TOLERANCE,
            "within_discretisation_tolerance": discretisation_worst <= DISCRETISATION_TOLERANCE,
        },
    }
    # The sweep's claim is the *published* comparison.  The discretisation bar is
    # resolution-dependent and reported per row (and by ``parent --gate``, which
    # is the entry that exports a mesh), so it does not decide this exit code.
    ok = bool(document["summary"]["within_published_tolerance"])
    return (0 if ok else 1), document


def _emit(document: dict[str, Any], *, json_out: bool, out: str | None, stem: str) -> None:
    if out:
        target = Path(out)
        target.mkdir(parents=True, exist_ok=True)
        (target / f"{stem}.json").write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        rows = document.get("rows")
        if rows:
            keys = ["cb_requested", "grid", "closed_form", "measured"]
            with (target / f"{stem}.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["cb", "nx", "nz", "Cb_model", "Cwp_model", "LCB_model_over_L",
                                 "Cm_model", "Cp_model", "VCB_model_over_T"])
                for row in rows:
                    measured = row["measured"]
                    writer.writerow([row["cb_requested"], row["grid"][0], row["grid"][1],
                                     f"{measured['block_coefficient']:.6f}",
                                     f"{measured['waterplane_coefficient']:.6f}",
                                     f"{measured['lcb_over_length']:.6f}",
                                     f"{measured['midship_coefficient']:.6f}",
                                     f"{measured['prismatic_coefficient']:.6f}",
                                     f"{measured['vcb_over_draft']:.6f}"])
            del keys
    if json_out:
        print(json.dumps(document, indent=2, sort_keys=True))
    else:
        summary = document.get("summary") or {}
        if summary:
            print(f"series60 family: {len(document['rows'])} parents, "
                  f"worst vs published {100.0 * summary['worst_published_relative']:+.3f} % "
                  f"({summary['worst_published_at']}), "
                  f"worst discretisation {100.0 * summary['worst_discretisation_relative']:+.3f} %",
                  file=sys.stderr)
        else:
            discretisation = document.get("discretisation", {})
            params = document.get("parameters", {})
            print(f"series60 parent: Cb={params.get('cb')} grid={discretisation.get('grid')} "
                  f"worst mesh-vs-closed-form "
                  f"{100.0 * float(discretisation.get('worst_relative_error') or 0.0):+.3f} % "
                  f"({discretisation.get('worst_quantity')})", file=sys.stderr)
            print(json.dumps(document.get("hydrostatics", {}), indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-hull series60",
        description="Generate and judge the Series 60 parent-form family.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(sub: argparse.ArgumentParser, *, defaults: bool) -> None:
        sub.add_argument("--length", type=float, default=10.0, help="length in metres")
        sub.add_argument("--beam", type=float, default=1.5, help="maximum moulded breadth in metres")
        sub.add_argument("--draft", type=float, default=0.6, help="keel-to-waterline draft in metres")
        sub.add_argument("--nx", type=int, default=161, help="longitudinal stations")
        sub.add_argument("--nz", type=int, default=41, help="vertical levels")
        sub.add_argument("--density", type=float, default=1025.0, help="water density in kg/m^3")
        sub.add_argument("--out", default=None, help="directory for the JSON/CSV record")
        sub.add_argument("--json", action="store_true", help="print the full record")
        if not defaults:
            return
        sub.add_argument("--transom-ratio", type=float, default=0.0,
                         help="aft-end deck-breadth ratio (0 = pointed stern)")
        sub.add_argument("--bulb-length-ratio", type=float, default=0.0,
                         help="bulb semi-axis ahead of the FP over L")
        sub.add_argument("--bulb-area-ratio", type=float, default=0.0,
                         help="bulb transverse area at the FP over the midship section area")
        sub.add_argument("--bulb-centre-height-ratio", type=float, default=0.35,
                         help="bulb centre height above the keel over the draft")

    parent = subparsers.add_parser("parent", help="one parent form, generated and exported")
    add_common(parent, defaults=True)
    parent.add_argument("--cb", type=float, default=0.70, help="block coefficient target")
    parent.add_argument("--cwp", type=float, default=None,
                        help="waterplane coefficient target (omit for the legacy hull)")
    parent.add_argument("--lcb", type=float, default=None,
                        help="LCB/L target measured from x=0 (needs --cwp)")
    parent.add_argument("--output", default="tensorcad-series60",
                        help="directory for hull.obj / offsets.csv / metadata.json")
    parent.add_argument("--no-export", action="store_true", help="report only, write no files")
    parent.add_argument("--stl", action="store_true", help="also write a binary STL")
    parent.add_argument("--gate", action="store_true",
                        help="run the mesh gate against the analytic volume and refuse a miss")

    family = subparsers.add_parser("family", help="the five Todd parents against the published table")
    add_common(family, defaults=False)
    family.add_argument("--published", default=None,
                        help="published parent particulars (JSON); omitted = closure only")
    family.add_argument("--cb", dest="cb_values", type=float, action="append", default=None,
                        help="restrict the sweep (repeatable); default = the five parents")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(list(sys.argv[1:] if argv is None else argv))
    try:
        if args.command == "parent":
            code, document = _parent_record(args)
            _emit(document, json_out=args.json, out=args.out, stem="series60_parent")
            return code
        code, document = _family_record(args)
        _emit(document, json_out=args.json, out=args.out, stem="five_cb_gradient")
        return code
    except ValueError as exc:
        print(f"tensorcad-hull series60: invalid model: {exc}", file=sys.stderr)
        if getattr(args, "json", False):
            print(json.dumps({"error": str(exc), "stage": "validate"}, sort_keys=True))
        return 2


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main(sys.argv[1:]))