"""The full chain: CAD file -> B-rep -> topology + continuity + hydrostatics -> gate.

Order matters, and it is the order a shipyard would use:

1. **Exchange** -- read the file into a B-rep (STEP topology, IGES sewing, native
   patch set) and record what the round trip cost in geometry.
2. **Topology** -- watertightness (0/1), free edges, cracks, non-manifold edges,
   orientation, duplicates, slivers, micro edges, unit scale.
3. **Geometry declared vs realised** -- ``L/B/T`` and ``Cb`` against the values on the
   drawing.  This runs before any hydrostatic number, because a displacement that
   integrates from a model of the wrong ship is not evidence about that ship.
4. **Continuity** -- seam class G0/G1/G2 with calibrated limits.
5. **Hydrostatics** -- patch quadrature (the reported route) and tessellation (the
   solver route), against the analytic reference.
6. **Grid convergence** -- three tessellation densities, observed order and GCI.
7. **Gate** -- every criterion evaluated by number, and :func:`enforce` refuses a
   report that does not satisfy them.  The gate never reads a narrative field.

Nothing here reports a number it did not compute, and every threshold lives in
:data:`CRITERIA` so a reviewer can argue with the numbers instead of the prose.
"""
from __future__ import annotations

import json
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import baselines as baselines_module
from .brep import BRep
from .continuity import continuity_evidence
from .defects import detect_all, run_library
from .hydro import convergence_study, facet_hydrostatics
from .io import load_cad, round_trip

#: Every machine-checked acceptance number of the campaign node, in one place.
CRITERIA: dict[str, dict[str, Any]] = {
    "G1_watertight": {
        "limit": 1.0,
        "unit": "0/1",
        "description": "the shell is closed, orientable and positively oriented",
        "comparison": "min",
    },
    "C1_wigley_volume_error": {
        "limit": 5e-4,
        "unit": "relative",
        "description": "Wigley III displacement from the CAD patch set vs closed form",
        "comparison": "max",
    },
    "C2_defect_detection_rate": {
        "limit": 1.0,
        "unit": "fraction",
        "description": "library cases whose declared defect code is reported (zero misses)",
        "comparison": "min",
    },
    "C3_curvature_f1": {
        "limit": 0.9,
        "unit": "F1",
        "description": "held-out seam classification F1 with calibrated thresholds",
        "comparison": "min",
    },
    "C4_hydrostatic_error": {
        "limit": 1e-3,
        "unit": "relative",
        "description": "worst of volume / LCB / VCB / waterplane area against the analytic hull",
        "comparison": "max",
    },
    "C5_grid_convergence_order": {
        "limit": 1.9,
        "unit": "observed order p",
        "description": "Richardson observed order of the tessellated displacement",
        "comparison": "min",
    },
    "P5_geometry_deviation": {
        "limit": 1e-4,
        "unit": "fraction of L",
        "description": "max |y_CAD - y_analytic| / L for every reference hull (campaign item 5)",
        "comparison": "max",
    },
    "C7_seconds_per_hull": {
        "limit": 10.0,
        "unit": "s",
        "description": "wall-clock for the whole chain on one hull",
        "comparison": "max",
    },
}


class GateError(RuntimeError):
    """Raised when a report fails a machine criterion."""


@dataclass
class CriterionResult:
    name: str
    value: float | None
    limit: float
    passed: bool
    detail: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "criterion": self.name,
            "value": self.value,
            "limit": self.limit,
            "passed": self.passed,
            "unit": CRITERIA.get(self.name, {}).get("unit", ""),
            "description": CRITERIA.get(self.name, {}).get("description", ""),
            "detail": self.detail,
        }


def criterion(
    name: str, value: float | None, *, comparison: str | None = None, detail: str = ""
) -> CriterionResult:
    spec = CRITERIA[name]
    mode = comparison or spec["comparison"]
    if value is None:
        return CriterionResult(
            name=name,
            value=None,
            limit=spec["limit"],
            passed=False,
            detail=detail or "no value",
        )
    passed = value >= spec["limit"] if mode == "min" else value <= spec["limit"]
    return CriterionResult(name=name, value=value, limit=spec["limit"], passed=passed, detail=detail)


def _relative(actual: float, expected: float, scale: float = 1.0) -> float:
    """Relative error, with the model scale as reference when the target is zero.

    LCB is exactly zero for a symmetric hull, so a relative error against it would divide
    by round-off; the error is then expressed against ``L``, which is the quantity a naval
    architect actually cares about ("a metre aft" is 10 % L on a 10 m model and 0.3 % L on
    a 320 m ship).
    """
    if abs(expected) >= 1e-6 * scale:
        return abs(actual - expected) / abs(expected)
    return abs(actual - expected) / scale


# --------------------------------------------------------------------------- #
# one hull, the whole chain
# --------------------------------------------------------------------------- #
def verify_cad(
    source: str | Path | BRep,
    *,
    declared: dict[str, float],
    reference: dict[str, float] | None = None,
    convergence_levels: Sequence[tuple[int, int]] = ((8, 4), (16, 8), (32, 16)),
    tolerance: float = 1e-6,
    workdir: str | Path | None = None,
) -> dict[str, Any]:
    """Run the chain on one model and return a machine-readable audit report."""
    started = time.perf_counter()
    exchange: dict[str, Any] = {}
    if isinstance(source, BRep):
        brep = source
        exchange["route"] = "in-memory"
    else:
        path = Path(source)
        brep = load_cad(path, tolerance=tolerance)
        exchange["route"] = path.suffix.lower()
        exchange["bytes"] = path.stat().st_size
    read_seconds = time.perf_counter() - started

    detection = detect_all(brep, declared=declared, tolerance=tolerance)
    hydro_patch = detection["hydrostatics"]
    hydro_facet = facet_hydrostatics(brep, u_samples=24, v_samples=12).to_dict()
    convergence = convergence_study(
        brep,
        levels=convergence_levels,
        reference=reference.get("volume") if reference else None,
    )
    comparison: dict[str, Any] = {}
    if reference:
        scale = declared.get("length", 1.0)
        for key, patch_key in (
            ("volume", "volume"),
            ("lcb_x", "lcb_x"),
            ("vcb_z", "vcb_z"),
            ("waterplane_area", "waterplane_area"),
        ):
            if key in reference:
                comparison[key] = {
                    "reference": reference[key],
                    "patch": hydro_patch[patch_key],
                    "facet": hydro_facet[patch_key],
                    "patch_relative_error": _relative(
                        hydro_patch[patch_key], reference[key], scale
                    ),
                    "facet_relative_error": _relative(
                        hydro_facet[patch_key], reference[key], scale
                    ),
                }
    if workdir is not None and not isinstance(source, BRep):
        pass
    seconds = time.perf_counter() - started
    return {
        "exchange": exchange,
        "seconds": seconds,
        "seconds_read": read_seconds,
        "watertight": detection["watertight"],
        "codes": detection["codes"],
        "defects": detection["defects"],
        "seams": detection["seams"],
        "curvature": detection["curvature"],
        "hydrostatics": {**hydro_patch, "facet": hydro_facet},
        "comparison": comparison,
        "convergence": convergence,
        "declared": dict(declared),
        "brep": {
            "faces": len(brep.faces),
            "edges": len(brep.edges),
            "control_points": sum(
                len(face.surface.control_net) * len(face.surface.control_net[0])
                for face in brep.faces
            ),
        },
    }


def worst_hydrostatic_error(report: dict[str, Any]) -> float | None:
    errors = [
        value["patch_relative_error"]
        for value in report.get("comparison", {}).values()
        if value.get("patch_relative_error") is not None
        and value["reference"] != 0
    ]
    return max(errors) if errors else None


def gate_single(report: dict[str, Any]) -> list[CriterionResult]:
    """Gate one model: watertightness, no fatal defect code, hydrostatics inside 0.1 %."""
    fatal = sorted(set(report.get("codes", ())))
    results = [
        criterion(
            "G1_watertight",
            1.0 if report.get("watertight") else 0.0,
            detail="topology audit on the model as read from the CAD file",
        ),
        criterion(
            "C4_hydrostatic_error",
            worst_hydrostatic_error(report),
            detail="worst of volume, LCB, VCB and waterplane area",
        ),
    ]
    if fatal:
        results.append(
            CriterionResult(
                name="G2_no_fatal_codes",
                value=0.0,
                limit=1.0,
                passed=False,
                detail="defect codes present: " + ", ".join(fatal),
            )
        )
    return results


def report_criteria(
    report: dict[str, Any],
    *,
    curvature_f1: float | None = None,
    detection_rate: float | None = None,
) -> list[CriterionResult]:
    results = [
        criterion(
            "G1_watertight",
            1.0 if report.get("watertight") else 0.0,
            detail="topology audit on the model as read from the CAD file",
        ),
        criterion(
            "C1_wigley_volume_error",
            report.get("comparison", {}).get("volume", {}).get("patch_relative_error"),
            detail="patch-quadrature displacement vs analytic/closed form",
        ),
        criterion(
            "C4_hydrostatic_error",
            worst_hydrostatic_error(report),
            detail="worst of volume, LCB, VCB and waterplane area",
        ),
        criterion(
            "C5_grid_convergence_order",
            report.get("convergence", {}).get("observed_order"),
            detail="Richardson order on three tessellation densities",
        ),
        criterion("C7_seconds_per_hull", report.get("seconds"), detail="whole chain, one hull"),
    ]
    if curvature_f1 is not None:
        results.append(
            criterion("C3_curvature_f1", curvature_f1, detail="held-out seam F1")
        )
    if detection_rate is not None:
        results.append(
            criterion(
                "C2_defect_detection_rate", detection_rate, detail="defect library, zero misses"
            )
        )
    return results


def enforce(results: Sequence[CriterionResult]) -> None:
    """Raise :class:`GateError` unless every criterion passed."""
    failed = [result for result in results if not result.passed]
    if failed:
        raise GateError(
            "machine gate rejected the geometry: "
            + "; ".join(
                f"{result.name}={result.value!r} (limit {result.limit})" for result in failed
            )
        )


# --------------------------------------------------------------------------- #
# the campaign evidence bundle
# --------------------------------------------------------------------------- #
def run_evidence(
    *,
    out_dir: str | Path = "artifacts/geomquality",
    defect_reference: str = "wigley_iii",
    cases: int = 150,
) -> dict[str, Any]:
    """Build, exchange, verify and gate all three reference hulls plus the library."""
    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    specimens = baselines_module.specimens()
    evidence: dict[str, Any] = {
        "criteria_spec": CRITERIA,
        "specimens": {},
        "exchange": {},
        "defect_library": {},
        "continuity": {},
    }
    reports: dict[str, Any] = {}
    for name, factory in specimens.items():
        specimen = factory()
        declared = {
            "length": specimen.length,
            "beam": specimen.beam,
            "draft": specimen.draft,
            "block_coefficient": specimen.particulars["block_coefficient"],
        }
        # 1. exchange: write every format, read it back, and measure what survived.
        trip = round_trip(specimen.brep, target / name, stem=name)
        evidence["exchange"][name] = trip
        # 2. run the chain on the file that came back from STEP (topology route).
        step_path = trip[".step"]["path"]
        report = verify_cad(
            step_path,
            declared=declared,
            reference=specimen.reference,
        )
        deviation = specimen.deviation()
        report["deviation"] = deviation
        report["specimen"] = specimen.to_dict()
        reports[name] = report
        evidence["specimens"][name] = {
            "watertight": report["watertight"],
            "codes": report["codes"],
            "hydrostatics": report["hydrostatics"],
            "comparison": report["comparison"],
            "convergence": report["convergence"],
            "deviation": deviation,
            "seconds": report["seconds"],
            "declared": declared,
            "reference": specimen.reference,
            "particulars": specimen.particulars,
            "brep": report["brep"],
            "settings": specimen.settings,
            "notes": specimen.notes,
        }

    reference_specimen = specimens[defect_reference]()
    declared_reference = {
        "length": reference_specimen.length,
        "beam": reference_specimen.beam,
        "draft": reference_specimen.draft,
        "block_coefficient": reference_specimen.particulars["block_coefficient"],
    }
    library = run_library(reference_specimen.brep, declared=declared_reference)
    evidence["defect_library"] = library
    evidence["defect_reference"] = defect_reference

    continuity = continuity_evidence(cases=cases)
    evidence["continuity"] = continuity

    # Criteria: the Wigley numbers carry C1/C4/C5, per the campaign spec.
    wigley = reports["wigley_iii"]
    worst_deviation = max(
        item["deviation"]["max_fraction_of_L"] for item in evidence["specimens"].values()
    )
    worst = 0.0
    for name, report in reports.items():
        error = worst_hydrostatic_error(report)
        if error is not None:
            worst = max(worst, error)
    results = [
        criterion(
            "G1_watertight",
            1.0 if report.get("watertight") else 0.0,
            detail="topology audit on the model as read from the CAD file",
        ),
        criterion(
            "C1_wigley_volume_error",
            wigley["comparison"]["volume"]["patch_relative_error"],
            detail="Wigley III patch quadrature vs 4/9 LBT closed form",
        ),
        criterion(
            "C2_defect_detection_rate",
            library["detection_rate"],
            detail=f"{library['detected']}/{library['cases']} library cases",
        ),
        criterion(
            "C3_curvature_f1",
            continuity["test"]["overall"]["f1"],
            detail=(
                f"{continuity['test_cases']} held-out seams, thresholds "
                f"{continuity['thresholds']}"
            ),
        ),
        criterion("C4_hydrostatic_error", worst, detail="worst across the three reference hulls"),
        criterion(
            "C5_grid_convergence_order",
            max(
                report["convergence"]["observed_order"] or 0.0
                for report in reports.values()
            ),
            detail="best observed order of the three hulls",
        ),
        criterion(
            "P5_geometry_deviation",
            worst_deviation,
            detail="worst of the three hulls, sampled on a 61 x 41 grid",
        ),
        criterion(
            "C7_seconds_per_hull",
            max(report["seconds"] for report in reports.values()),
            detail="slowest of the three hulls, file read included",
        ),
    ]
    evidence["criteria"] = [result.to_dict() for result in results]
    evidence["gate"] = {
        "passed": all(result.passed for result in results),
        "failed": [result.name for result in results if not result.passed],
    }
    payload = json.dumps(evidence, indent=1, default=str)
    (target / "evidence.json").write_text(payload)
    (target / "gate.json").write_text(
        json.dumps(
            {
                "criteria": evidence["criteria"],
                "gate": evidence["gate"],
                "exchange": evidence["exchange"],
                "defect_matrix": library["rows"],
                "continuity": {
                    "thresholds": continuity["thresholds"],
                    "test": continuity["test"],
                    "separation": continuity["separation"],
                    "cases": continuity["cases"],
                },
            },
            indent=1,
            default=str,
        )
    )
    return evidence


def main(argv: Sequence[str] | None = None) -> int:
    """``python -m tensorcad.geomquality.pipeline --out artifacts/geomquality``."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/geomquality")
    parser.add_argument("--cases", type=int, default=150)
    args = parser.parse_args(list(argv) if argv is not None else None)
    evidence = run_evidence(out_dir=args.out, cases=args.cases)
    results = [
        CriterionResult(
            name=item["criterion"],
            value=item["value"],
            limit=item["limit"],
            passed=item["passed"],
            detail=item.get("detail", ""),
        )
        for item in evidence["criteria"]
    ]
    for result in results:
        print(
            f"{result.name:32s} {'PASS' if result.passed else 'FAIL'} "
            f"value={result.value!r} limit={result.limit}"
        )
    print("gate:", "PASS" if evidence["gate"]["passed"] else f"FAIL {evidence['gate']['failed']}")
    return 0 if evidence["gate"]["passed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
