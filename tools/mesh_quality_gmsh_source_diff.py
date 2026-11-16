"""Gmsh differential, driven by the Gmsh *source* oracle.

`tools/mesh_quality_external_diff.py` is the general door: it compares
:mod:`tensorcad.meshquality` against any oracle that speaks the ``W <kind> <metric>
<index> <value>`` protocol. This script is the Gmsh arm of that door, run against
``gmsh_source_oracle`` -- Gmsh's own ``qualityMeasures.cpp`` /
``qualityMeasuresJacobian.cpp`` compiled from a Gmsh release source tree (see
``tools/mesh_quality_oracles/build_gmsh_source_oracle.sh``).

Two classes of comparison, kept apart so neither can borrow the other's credit:

``MAPPING``
    pairs that must agree to ``--tolerance``. A disagreement here is a failure and
    exits 1. These are the pairs whose *definitions* are the same in both trees.

``EXPLAINED``
    pairs that are measured but expected to differ, because the two trees answer
    slightly different questions. They are recorded with the measured difference and the
    reason, and they do **not** decide the exit code -- the anchor accepts "差异<1e-6
    或差异解释", and a pair that is quietly identical would mean the explanation is
    wrong. Keeping them separate is what stops an explained pair from being counted as
    a pass (and what stops a genuine regression in it from hiding behind the label).

Both tables were derived from executed numbers, not from reading the formulas; the
Gmsh unit for each kind is established by measuring Gmsh's own value on the *ideal*
element of that kind (see ``docs/m23119_gmsh_eigen_differential.md``).

Usage::

    python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin
    python tools/mesh_quality_gmsh_source_diff.py --oracle ... --require-oracle

Exit codes: 0 every ``MAPPING`` pair within 1e-6; 1 a ``MAPPING`` pair disagreed; 2 the
oracle did not run (and ``--require-oracle`` was given).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import mesh_quality_external_diff as drv  # noqa: E402

TOLERANCE = 1.0e-6

#: our metric -> Gmsh measure name, for pairs that assert agreement.
#:
#: ``tri``/``quad`` include ``minSICN``: Gmsh's ICN for a simplex/quadrilateral is the
#: corner quantity normalised to the ideal element of that kind, which is exactly the
#: definition of our ``min_sicn`` -- measured agreement is <= 1e-15 (see the doc).
#: ``hex`` carries only ``volume``: its ICN pair is in ``EXPLAINED`` (Gmsh minimises an
#: ideal-Jacobian bound over the whole cell, the library samples the physical corner
#: Jacobians), and Gmsh reports 0 for hex ``gamma``/``eta`` because ``MHexahedron``
#: overrides neither.
MAPPING: dict[str, dict[str, str]] = {
    "tri": {"gamma": "gamma", "eta": "eta", "min_sicn": "minSICN"},
    "quad": {"eta": "eta", "area": "area", "min_sicn": "minSICN"},
    "tet": {"gamma": "gamma", "eta": "eta", "volume": "volume"},
    "hex": {},
}

#: ours -> (Gmsh measure, why the two numbers are not expected to coincide).
#: Measured on the corpus; the recorded difference is the evidence for the sentence.
EXPLAINED: dict[str, dict[str, tuple[str, str]]] = {
    "tet": {
        "min_sicn": (
            "minSICN",
            "Gmsh's tetrahedral ICN is built from the element's ideal (reference) "
            "Jacobian and ideal gradients and is then minimised over the cell by Bezier "
            "subdivision (qualityMeasuresJacobian.cpp: minICNMeasure); our min_sicn "
            "samples the physical corner Jacobians normalised to the regular "
            "tetrahedron. Both equal 1.0 on the regular tetrahedron, so they are the "
            "same idea, but they are not the same number off it (measured spread on "
            "this corpus: see max_rel_diff).",
        ),
        "gmsh_icn": (
            "minSICN",
            "This pair is the one the library's docstring claims is an identity "
            "('Gmsh minSICN unit: 1 = right-corner tetrahedron, regular tetrahedron = "
            "0.7937'). Measured against Gmsh 4.15.2 that claim is inverted: Gmsh "
            "reports 1.0 for a regular tetrahedron and 0.8165 for Gmsh's own "
            "right-corner reference tetrahedron, while our gmsh_icn reports 0.7937 and "
            "1.0 for the same two cells. Recorded, not asserted, so the day the "
            "normalisation is fixed the number moves visibly instead of silently.",
        ),
    },
    "hex": {
        "volume": (
            "volume",
            "Two definitions of the same word. Gmsh integrates det J of the trilinear "
            "map over the reference cube with Gauss quadrature (MElement::getVolume, "
            "getIntegrationPoints(order+3)); the library applies the divergence theorem "
            "to the six faces split into 12 triangles (_hex_volume). Identical on the "
            "ideal cube (1.3e-15) and on any cell with planar faces, diverging once a "
            "face is warped -- the corpus's jittered 3x2x2 grid is exactly that case. "
            "Neither is wrong; the anchor's pair is only meaningful for planar faces.",
        ),
        "min_sicn": (
            "minSICN",
            "Gmsh minimises ICN over the cell (Bezier lower bound over ideal-Jacobian "
            "gradients); our min_sicn is the minimum over the 8 physical corner "
            "Jacobians. Identical on the ideal cube (both 1.0), up to a few percent "
            "apart on a jittered grid, in both directions.",
        ),
    },
}

#: Measures Gmsh itself cannot answer for a kind, with the observed value, so the gap is
#: evidence rather than silence.
NOT_COMPARABLE: dict[str, str] = {
    "gamma/eta (hex)": "MHexahedron overrides neither gammaShapeMeasure() nor "
                       "etaShapeMeasure(); Gmsh reports exactly 0.0 for every hex "
                       "(measured on hex_unit and hex_distorted). The hex arm's "
                       "external reference for those two measures is the Verdict/VTK "
                       "oracle, which tools/mesh_quality_external_diff.py drives.",
    "minSIGE (tri/quad/tet)": "Gmsh 4.15.2 returns nan for minSIGE on tri/quad/tet "
                              "corpora and a number only for hex; the library has no "
                              "SIGE metric, so nothing to compare.",
}

ARTIFACT = REPO / "artifacts" / "mesh_quality_external_diff" / "gmsh_source_diff.json"


def ours_for(kind: str, points, elems):
    """Per-element measures from the library, as ``{metric: [values]}``.

    ``METRIC_FN`` returns an :class:`tensorcad.meshquality.metrics.ElementQuality`, not a
    mapping, so ``.metrics`` has to be unwrapped before :func:`drv.compare` can index it.
    The scalar backend is requested explicitly: the driver compares against a scalar C++
    oracle, so a vectorised backend would be a second variable in the measurement.
    """
    fn = drv.METRIC_FN[kind]
    try:
        quality = fn(points, elems, backend="scalar")
    except TypeError:
        quality = fn(elems, points)
    return quality.metrics


def _summarise(cmp: dict, counted: bool) -> tuple[int, int, float]:
    """(pairs, elements_beyond_tol, worst_rel_diff) for one ``compare`` result."""
    pairs = 0
    beyond = 0
    worst = 0.0
    for info in cmp.values():
        if info.get("status") not in ("AGREE", "DIFFER"):
            continue
        pairs += 1
        worst = max(worst, info.get("max_rel_diff", 0.0) or 0.0)
        if counted:
            beyond += int(info.get("n_beyond_tol", 0) or 0)
    return pairs, beyond, worst


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--oracle", required=False, default=None)
    ap.add_argument("--require-oracle", action="store_true")
    ap.add_argument("--artifact", default=str(ARTIFACT))
    ap.add_argument("--tolerance", type=float, default=TOLERANCE)
    args = ap.parse_args(argv)

    tmpdir = Path(args.artifact).resolve().parent
    tmpdir.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "driver": "tools/mesh_quality_gmsh_source_diff.py",
        "oracle": "gmsh_source_oracle (Gmsh qualityMeasures.cpp + "
                  "qualityMeasuresJacobian.cpp, source build)",
        "tolerance": args.tolerance,
        "meshes": {},
    }

    disagreed = 0
    compared_pairs = 0
    explained_pairs = 0
    ran_any = False

    for name, (kind, points, elems) in sorted(drv.meshes().items()):
        if kind not in MAPPING and kind not in EXPLAINED:
            report["meshes"][name] = {"kind": kind, "status": "SKIPPED",
                                      "why": "no Gmsh measure for this kind; see "
                                             "NOT_COMPARABLE in this driver"}
            continue
        mesh_file = tmpdir / f"_{name}.mkm"
        drv.write_mesh_file(mesh_file, kind, points, elems)
        status, external = drv.run_oracle(args.oracle, mesh_file, "W")
        mesh_file.unlink(missing_ok=True)
        entry: dict = {"kind": kind, "elements": len(elems), "status": status,
                       "not_comparable": dict(NOT_COMPARABLE)}
        if status != "RAN" or kind not in external:
            entry["why"] = "oracle did not run"
            report["meshes"][name] = entry
            continue
        ran_any = True
        ours = ours_for(kind, points, elems)
        cmp = drv.compare(ours, external[kind], MAPPING.get(kind, {}))
        entry["metrics"] = cmp
        pairs, beyond, worst = _summarise(cmp, counted=True)
        compared_pairs += pairs
        disagreed += beyond
        entry["worst_rel_diff"] = worst

        explained = {}
        for ours_metric, (gmsh_metric, why) in EXPLAINED.get(kind, {}).items():
            if ours_metric not in ours or gmsh_metric not in external[kind]:
                continue
            one = drv.compare(ours, external[kind], {ours_metric: gmsh_metric})
            info = one[ours_metric]
            info["why"] = why
            explained[ours_metric] = info
            explained_pairs += 1
        if explained:
            entry["explained"] = explained
        report["meshes"][name] = entry

    report["pairs_compared"] = compared_pairs
    report["pairs_explained"] = explained_pairs
    report["elements_beyond_tolerance"] = disagreed
    report["oracle_ran"] = ran_any
    Path(args.artifact).parent.mkdir(parents=True, exist_ok=True)
    Path(args.artifact).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    for name, entry in report["meshes"].items():
        if entry.get("status") != "RAN":
            print(f"{name:16s} {entry['status']}")
            continue
        for metric, info in entry.get("metrics", {}).items():
            print(f"{name:16s} {metric:10s} {info.get('status'):9s} "
                  f"worst_rel={info.get('max_rel_diff', float('nan')):.3e} "
                  f"beyond_tol={info.get('n_beyond_tol', 0)}/{entry['elements']}")
        for metric, info in entry.get("explained", {}).items():
            print(f"{name:16s} {metric:10s} {'EXPLAINED':9s} "
                  f"worst_rel={info.get('max_rel_diff', float('nan')):.3e} "
                  f"(definitional, not a gate)")
    print(f"pairs compared {compared_pairs}, beyond {args.tolerance:g}: {disagreed}; "
          f"explained {explained_pairs}")

    if not ran_any:
        return 2 if args.require_oracle else 0
    return 1 if disagreed else 0


if __name__ == "__main__":
    raise SystemExit(main())