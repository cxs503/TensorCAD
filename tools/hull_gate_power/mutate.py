#!/usr/bin/env python3
"""Gate-power audit of the parametric hull generator's closed forms.

Continues the 2026-11-04 claim-gate-power audit
(``docs/analysis/2026-11-04-claim-gate-power.md``), which explicitly left
``hull.py`` uncovered: "the generator's own closed forms are shared constants --
a mutant of the exponent would have to be paired with the analytic side to be
meaningful". This harness performs that pairing: one single-defect mutant at a
time, injected by exact-string replacement into ``tensorcad/hull.py`` (the
generator) or ``tensorcad/hydrostatics_analytic.py`` (the analytic reference),
with the test files that guard the documented claim run against each mutant.

A mutant the suite does not notice is a documented number with no teeth.

Method (same discipline as the 2026-11-04 audit):

* one defect at a time, patched by exact-string replacement (a pattern that is
  not unique is refused, not patched);
* a clean-tree guard aborts if ``git diff -- tensorcad`` is non-empty, because
  an interrupted run would leave one mutant in place;
* the file is restored from its original bytes with a SHA-256 check before the
  next mutant starts;
* each mutant runs only the test files that guard its claim, so a kill is
  attributable; baselines (unmutated) are recorded for the same file sets.

Usage::

    python3 tools/hull_gate_power/mutate.py --out results.json [--only HUL-1,ANA-2]

The run is resumable: mutants already present in ``--out`` are skipped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HULL = ROOT / "tensorcad" / "hull.py"
ANA = ROOT / "tensorcad" / "hydrostatics_analytic.py"
SERIES = ROOT / "tensorcad" / "series60.py"
SERIES_CLI = ROOT / "tensorcad" / "series60_cli.py"

GEN_TESTS = ["tests/test_hull.py", "tests/test_analytic_hydrostatics.py"]
ANA_TESTS = ["tests/test_hydrostatics_analytic.py", "tests/test_analytic_hydrostatics.py", "tests/test_hull.py"]
#: The transom/bulb claims are guarded by the two Series 60 test modules: the first
#: pins the new ends, the second pins that they did not move the legacy hull.
SERIES_TESTS = ["tests/test_series60_bulb_transom.py", "tests/test_series60.py"]
#: The product entry (``tensorcad-hull series60``) has its own claims: that the
#: command reaches the family, that the exported document carries the analytic
#: reference, that the gate refuses instead of writing, and that the published
#: comparison labels a closure as a closure.
SERIES_CLI_TESTS = ["tests/test_series60_cli.py"]
#: The family-wide comparison table (frame pairing, derived coefficients) is
#: guarded by the parameter-box sweep, not by any single-hull test: the defects
#: it pins (a `-193 %` VCB frame mismatch, two `nan` coefficients) pass every
#: test that only looks at one hull.
SERIES_CONFORMANCE_TESTS = ["tests/test_series60_family_conformance.py"]

MUTANTS: list[dict] = [
    {
        "id": "HUL-1-taper-exponent",
        "kind": "model",
        "file": "hull",
        "old": "        _breadth_exponent(parameters, index),\n    )",
        "new": "        0.8,  # mutant HUL-1\n    )",
        "claim": "mesh volume 7.446554 m^3 at 41x13 vs analytic 7.512698 (-0.88 %): the sampled surface must be the surface the reference integrates",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-2-vertical-taper",
        "kind": "model",
        "file": "hull",
        "old": "    return math.sqrt(level / (parameters.nz - 1))",
        "new": "    return level / (parameters.nz - 1)  # mutant HUL-2",
        "claim": "the sqrt vertical profile (the 2/3 volume factor of B*L*T*I*2/3)",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-3-end-zeroing",
        "kind": "gate",
        "file": "hull",
        "old": "    if index == 0 or index == parameters.nx - 1:\n        return 0.0",
        "new": "    if False:  # mutant HUL-3\n        return 0.0",
        "claim": "the end stations close on the centreline (math.sin(pi) is not exactly zero in floating point)",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-4-midship-split",
        "kind": "gate",
        "file": "hull",
        "old": "    if 2 * index >= parameters.nx - 1:",
        "new": "    if 2 * index > parameters.nx - 1:  # mutant HUL-4",
        "claim": "midship belongs to both halves; the tie must not bend the surface",
        "tests": GEN_TESTS,
    },
    {
        "id": "ANA-1-volume-factor",
        "kind": "value",
        "file": "ana",
        "old": "    volume = beam * length * draft * taper * (2.0 / 3.0)",
        "new": "    volume = beam * length * draft * taper * (1.0 / 2.0)  # mutant ANA-1",
        "claim": "displaced volume 7.512698 m^3 = B*L*T*I*2/3",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-2-vcb",
        "kind": "value",
        "file": "ana",
        "old": "        vertical_centre_buoyancy=-0.4 * draft,",
        "new": "        vertical_centre_buoyancy=-0.5 * draft,  # mutant ANA-2",
        "claim": "VCB = -0.4*T (0.6*T above the keel)",
        "claim_note": "the z-weighted mean of zeta**(1/2) on [-T, 0]",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-3-waterplane-half-percent",
        "kind": "tolerance",
        "file": "ana",
        "old": "    waterplane = beam * length * taper",
        "new": "    waterplane = beam * length * taper * 0.995  # mutant ANA-3",
        "claim": "waterplane area 14.061669 m^2 at 41x13; the fine grid is within half a percent on every quantity",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-4-taper-arg",
        "kind": "value",
        "file": "ana",
        "old": "math.gamma(0.5 * (exponent + 1.0))",
        "new": "math.gamma(0.5 * (exponent + 1.5))  # mutant ANA-4",
        "claim": "taper_integral is a Beta integral: e=0 -> 1, e=1 -> 2/pi, e=2 -> 1/2",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-5-volume-half-percent",
        "kind": "tolerance",
        "file": "ana",
        "old": "    volume = beam * length * draft * taper * (2.0 / 3.0)",
        "new": "    volume = beam * length * draft * taper * (2.0 / 3.0) * 1.005  # mutant ANA-5",
        "claim": "the documented -0.880 % volume grid error, pinned to 5e-5",
        "tests": ANA_TESTS,
    },
    {
        "id": "HUL-5-taper-integral-const",
        "kind": "value",
        "file": "hull",
        "old": "TAPER_INTEGRAL = (2.0 / math.pi) * math.sqrt(math.pi) * math.gamma(0.85) / (2.0 * math.gamma(1.35))",
        "new": "TAPER_INTEGRAL = (1.0 / math.pi) * math.sqrt(math.pi) * math.gamma(0.85) / (2.0 * math.gamma(1.35))  # mutant HUL-5",
        "claim": "TAPER_INTEGRAL = int_0^1 sin(pi*x)^0.7 dx, cross-derived from the volume identity",
        "tests": ANA_TESTS,
    },
    # ------------------------------------------------------------------ Series 60 ends
    # A mutant nobody kills is a documented number with no teeth. These target the
    # transom blend, the bulb closed forms and the exported bulb solid.
    {
        "id": "S60-1-transom-blend-area",
        "kind": "model",
        "file": "series",
        "old": "        keep * h + transom * 0.5,",
        "new": "        keep * h + transom * 0.55,  # mutant S60-1 (transom share of the aft area)",
        "claim": "transom_ratio blends the aft taper as (1-tau) f + tau, so the stern "
                 "share of the waterplane is tau/2 and Cwp/LCB/Cb still close",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-2-transom-blend-keep",
        "kind": "model",
        "file": "series",
        "old": "    keep = 1.0 - transom",
        "new": "    keep = 1.0 + transom  # mutant S60-2",
        "claim": "the base taper is *thinned* as the stern fills, not amplified",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-3-transom-first-moment",
        "kind": "model",
        "file": "series",
        "old": "        keep * m + transom * 0.125,",
        "new": "        keep * m + transom * 0.25,  # mutant S60-3 (int_0^1/2 xi dxi = 1/8)",
        "claim": "the transom's own first moment is int_0^1/2 xi dxi = 1/8",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-4-transom-square-moment",
        "kind": "model",
        "file": "series",
        "old": "        keep * keep * mq + 2.0 * keep * transom * m + transom * transom * 0.125,",
        "new": "        keep * keep * mq + 2.0 * keep * transom * m + transom * transom * 0.25,  # mutant S60-4",
        "claim": "the squared blend's transom term weights int xi dxi with 1/8, so the "
                 "section taper gain and Cb stay on target",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-5-deck-taper-half",
        "kind": "model",
        "file": "series",
        "old": "        if self.transom_ratio > 0.0 and xi < 0.5:",
        "new": "        if self.transom_ratio > 0.0 and xi < 0.6:  # mutant S60-5",
        "claim": "the transom applies to the aft half only, xi < 1/2",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-6-transom-area-factor",
        "kind": "model",
        "file": "series",
        "old": "    return p.beam * law.deck_taper(0.0) * p.draft * law.fullness(0.0)",
        "new": "    return 0.5 * p.beam * law.deck_taper(0.0) * p.draft * law.fullness(0.0)  # mutant S60-6",
        "claim": "the immersed stern section area is B f T Cm (port and starboard), and "
                 "the mesh's own stern faces must reproduce it",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-7-single-exponent-bracket",
        "kind": "model",
        "file": "series",
        "old": "        if cwp(mid) < cwp_target:",
        "new": "        if cwp(mid) > cwp_target:  # mutant S60-7",
        "claim": "Cwp rises with the longitudinal exponent; the bisection direction",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-8-bulb-volume-cubic",
        "kind": "model",
        "file": "series",
        "old": "    return 0.5 * math.pi * a * b * ((h + b) - (h ** 3 + b ** 3) / (3.0 * b * b))",
        "new": "    return 0.5 * math.pi * a * b * ((h + b) - (h ** 3 + b ** 3) / (2.0 * b * b))  # mutant S60-8",
        "claim": "the submerged half-spheroid volume (pi a b / 2)[(h+b) - (h^3+b^3)/(3b^2)], "
                 "which collapses to 2/3 pi a b^2 and to pi a b^2 / 3",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-9-bulb-submerged-clip",
        "kind": "model",
        "file": "series",
        "old": "    if waterline_offset >= b:",
        "new": "    if waterline_offset >= 2.0 * b:  # mutant S60-9",
        "claim": "a bulb whose axis is more than one radius below the waterline is fully "
                 "submerged: its volume stops depending on the offset",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-10-bulb-centroid-k",
        "kind": "model",
        "file": "series",
        "old": "    k = 0.125 * (s * (5.0 - 2.0 * s * s) * root + 3.0 * math.asin(s)) + 3.0 * math.pi / 16.0",
        "new": "    k = 0.125 * (s * (5.0 - 2.0 * s * s) * root + 3.0 * math.asin(s)) + 3.0 * math.pi / 8.0  # mutant S60-10",
        "claim": "the submerged bulb's x centroid, which must collapse to 3a/8 when the "
                 "bulb is fully submerged",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-11-bulb-fp-segment",
        "kind": "model",
        "file": "series",
        "old": "        fp_submerged = b * b * (theta - math.sin(theta) * math.cos(theta))",
        "new": "        fp_submerged = b * b * (theta - math.sin(theta))  # mutant S60-11",
        "claim": "the circular-segment area of the FP ellipse below the waterline",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-12-bulb-waterplane",
        "kind": "model",
        "file": "series",
        "old": "        waterplane = 0.5 * math.pi * a * b * (1.0 - s * s)",
        "new": "        waterplane = 0.5 * math.pi * a * b * (1.0 - s * s * 0.9)  # mutant S60-12",
        "claim": "the waterplane patch a piercing bulb adds is a *half* ellipse",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-13-bulb-wetted-arc",
        "kind": "model",
        "file": "series",
        "old": "        arc = 2.0 * (math.pi - math.acos(h / b))",
        "new": "        arc = 2.0 * (math.pi - math.acos(h / b)) * 0.9  # mutant S60-13",
        "claim": "the submerged arc of a bulb ring is 2(pi - acos(h/r))",
        "tests": SERIES_TESTS,
    },
    {
        "id": "S60-14-bulb-wetted-stretch",
        "kind": "model",
        "file": "series",
        "old": "        stretch = math.sqrt(a ** 2 + b ** 2 * s * s)",
        "new": "        stretch = math.sqrt(a ** 2 + b ** 2 * s * s) * 0.9  # mutant S60-14",
        "claim": "the axial stretch enters the bulb area element as r sqrt(r'^2 + a^2)",
        "tests": SERIES_TESTS,
    },
    # ------------------------------------------------------------------ Series 60 CLI
    {
        "id": "S60-15-entry-ignores-cb",
        "kind": "gate",
        "file": "series_cli",
        "old": "        cb = params.cb",
        "new": "        cb = 0.7  # mutant S60-15",
        "claim": "the CLI honours the block coefficient from the request",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-16-entry-vcb-frame",
        "kind": "value",
        "file": "series_cli",
        "old": "        vcb = hydro.vertical_centre_buoyancy",
        "new": "        vcb = -0.5 * params.draft  # mutant S60-16",
        "claim": "the VCB frame is correct (not -0.5*T)",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-17-entry-self-referential-gate",
        "kind": "gate",
        "file": "series_cli",
        "old": "        if hydro.block_coefficient < 0.5 or hydro.block_coefficient > 0.85:",
        "new": "        if hydro.block_coefficient < 0.4 or hydro.block_coefficient > 0.9:  # mutant S60-17",
        "claim": "the CLI gate refuses out-of-range Cb",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-18-entry-closure-as-prediction",
        "kind": "gate",
        "file": "series_cli",
        "old": "        role = 'closure'",
        "new": "        role = 'external'  # mutant S60-18",
        "claim": "the CLI labels the comparison as a closure, not external",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-19-derive-midship-coefficient",
        "kind": "value",
        "file": "series_cli",
        "old": "        cm = hydro.maximum_section_area / (params.beam * params.draft)",
        "new": "        cm = hydro.block_coefficient / hydro.prismatic_coefficient * 1.1  # mutant S60-19",
        "claim": "the CLI derives Cm from the mesh, not from the Cp ladder",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-20-derive-waterplane-coefficient",
        "kind": "value",
        "file": "series_cli",
        "old": "        cwp = hydro.waterplane_area / (params.length * params.beam)",
        "new": "        cwp = hydro.prismatic_coefficient * 0.95  # mutant S60-20",
        "claim": "the CLI derives Cwp from the mesh, not from the Cp ladder",
        "tests": SERIES_CLI_TESTS,
    },
    {
        "id": "S60-21-comparison-is-self-referential",
        "kind": "gate",
        "file": "series_cli",
        "old": "        if abs(err) > 0.03:",
        "new": "        if abs(err) > 0.1:  # mutant S60-21",
        "claim": "the CLI comparison bar is 3%, not 10%",
        "tests": SERIES_CLI_TESTS,
    },
]

#: Mutable modules that the harness can patch.
MUTABLE = {
    "hull": HULL,
    "ana": ANA,
    "series": SERIES,
    "series_cli": SERIES_CLI,
}


def stray_markers(paths: list[Path] | None = None) -> list[str]:
    """Scan files for ``# mutant`` markers.

    Returns a list of ``"filename:line: content"`` strings for every line that
    contains ``# mutant``.  If ``paths`` is ``None``, scans all files in
    ``MUTABLE``.
    """
    if paths is None:
        paths = list(MUTABLE.values())
    offenders: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        for i, line in enumerate(path.read_text().splitlines(), start=1):
            if "# mutant" in line:
                offenders.append(f"{path.name}:{i}: {line.strip()}")
    return offenders
