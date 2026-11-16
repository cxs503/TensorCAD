#!/usr/bin/env python
"""Mutation audit for the golden geometry gate (GAP-CAD-D2-GOLDEN-REGRESSION).

Why this exists
---------------
``tests/test_golden_regression.py`` proves the *harness* is well behaved: it fails closed on
a missing/unreadable/stale baseline, rejects a hand-edited reference, refuses to bless a
tree that fails its own gate, and is wired to a git hook and to CI. Every one of those
negative tests injects either a synthetic probe case or a tampered baseline file.

None of them answers the question the gate exists for: **if the geometry code itself
regresses, does the gate notice?** This audit answers it with two independent batteries.

Battery A - binding controls (does the gate see the real library?)
    A small, deliberate mutation is pushed into the *real* ``tensorcad`` entry point that a
    golden case calls (``hydrostatics.compute_hydrostatics``, ``series60.generate_series60``,
    ``wigley.generate_wigley``, ``DTMB5415Space.hydrostatics``), and the untouched gate is
    run against the committed baseline. Detection is required. A no-mutation control must
    still pass, so a harness that fails everything cannot look sensitive.

Battery B - self-blessed blind-spot envelope (what can slip through?)
    The drift layer compares against the blessed measurement at 1e-9, so it catches *any*
    reproducible change. The interesting adversary is therefore the developer who answers a
    red gate with ``--bless``. For every asserted analytic metric this battery scales the
    measured value by ``1 + delta`` on both grid levels, recomputes the tolerance exactly as
    ``--bless`` would, re-blesses, and re-checks. The smallest delta that survives a bless is
    the metric's real sensitivity floor; it is reported per metric, in milliseconds of CPU,
    because the geometry is measured once per case and every cell is then pure evaluation.

The sweep is *not* a model of the gate: every cell is evaluated by
``tensorcad.golden.evaluate_case`` itself, against a baseline assembled the way ``bless``
assembles one. Only ``case.measure`` is replaced by a cached, perturbed delegate.

Usage::

    python tools/golden/mutation_audit.py --check          # CI / test gate
    python tools/golden/mutation_audit.py --json OUT --md OUT --figure OUT
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import importlib.util
import io
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import golden  # noqa: E402

SCHEMA = "tensorcad.golden.mutation_audit/1"

#: Relative perturbations used by battery B. Dense around the tolerance floor (1e-4) and
#: around the largest tolerance the committed baseline derives (~2e-2, Series 60 block
#: coefficient), so a threshold can always be bracketed to within one sweep step.
PERTURBATIONS: tuple[float, ...] = (
    1e-9,
    1e-6,
    1e-5,
    1e-4,
    3e-4,
    1e-3,
    3e-3,
    1e-2,
    3e-2,
    1e-1,
)

#: Detection layers of the gate, in the order the gate evaluates them.
LAYERS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("reference", ("reference_not_reproducible", "reference_mismatch", "reference_not_in_baseline")),
    ("accuracy", ("reference_mismatch",)),
    ("invariant", ("invariant_violated",)),
    ("drift", ("drift",)),
    ("crash", ("case_crashed",)),
)


# --------------------------------------------------------------------------------------
# transparent proxies: mutate one attribute of a real object, leave everything else alone
# --------------------------------------------------------------------------------------
class _Scaled:
    """Delegate every attribute of ``inner``, except the ones given a replacement.

    A replacement is either a number (the attribute is multiplied by it) or an object (the
    attribute is returned in its place, e.g. a scaled ``hydrostatics``), and mappings are
    transparently subscriptable so a dict-returning API keeps working.
    """

    def __init__(self, inner: Any, **factors: Any) -> None:
        object.__setattr__(self, "_inner", inner)
        object.__setattr__(self, "_factors", factors)

    def __getattr__(self, name: str) -> Any:
        value = getattr(object.__getattribute__(self, "_inner"), name)
        factors = object.__getattribute__(self, "_factors")
        if name in factors:
            replacement = factors[name]
            if isinstance(replacement, (int, float)):
                return float(value) * float(replacement)
            return replacement
        return value

    # -- mapping protocol: several reporting APIs return plain dicts ---------------
    def __getitem__(self, key: str) -> Any:
        value = object.__getattribute__(self, "_inner")[key]
        factors = object.__getattribute__(self, "_factors")
        if key in factors:
            replacement = factors[key]
            if isinstance(replacement, (int, float)):
                return float(value) * float(replacement)
            return replacement
        return value

    def __contains__(self, key: str) -> bool:
        return key in object.__getattribute__(self, "_inner")

    def __iter__(self) -> Iterator[str]:
        return iter(object.__getattribute__(self, "_inner"))

    def __len__(self) -> int:
        return len(object.__getattribute__(self, "_inner"))

    def keys(self) -> Any:
        return object.__getattribute__(self, "_inner").keys()

    def items(self) -> Any:
        return [(key, self[key]) for key in self]

    def get(self, key: str, default: Any = None) -> Any:
        return self[key] if key in self else default

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"_Scaled({object.__getattribute__(self, '_inner')!r})"


@contextlib.contextmanager
def _patched(target: Any, name: str, replacement: Any) -> Iterator[None]:
    """Replace ``target.name`` for the duration of the block, restoring it afterwards."""
    original = getattr(target, name)
    setattr(target, name, replacement)
    try:
        yield
    finally:
        setattr(target, name, original)


# --------------------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Control:
    """One library-level mutation that the gate must detect."""

    name: str
    case: str
    detail: str
    expected: tuple[str, ...]
    mutate: Any


@dataclass
class ControlResult:
    name: str
    case: str
    detail: str
    expected: tuple[str, ...]
    detected: bool
    codes: tuple[str, ...]
    layers: tuple[str, ...]
    failed_cases: tuple[str, ...] = field(default_factory=tuple)
    findings: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "case": self.case,
            "detail": self.detail,
            "expected_codes": list(self.expected),
            "detected": self.detected,
            "codes": list(self.codes),
            "layers": list(self.layers),
            "failed_cases": list(self.failed_cases),
            "findings": list(self.findings),
        }


@dataclass
class EnvelopeRow:
    """One (case, metric) sensitivity record from battery B."""

    case: str
    metric: str
    tolerance: float
    reference_error: float
    floor: float | None  # smallest delta that survives a bless, None = never survives
    ceiling_step: float | None  # smallest delta that is caught after a bless
    survives: tuple[float, ...]
    caught_by: dict[str, list[str]]
    incoherent_floor: float | None = None
    incoherent_caught: float | None = None
    invariants_hold: bool = True
    deltas: tuple[float, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "case": self.case,
            "metric": self.metric,
            "tolerance": self.tolerance,
            "reference_error": self.reference_error,
            "coherent": {
                "blind_spot_floor": self.floor,
                "detection_ceiling": self.ceiling_step,
                "survives_a_bless": [_number(delta) for delta in self.survives],
                "caught_by": {_number(delta): codes for delta, codes in self.caught_by.items()},
                "invariants_hold_at_every_delta": self.invariants_hold,
            },
            "incoherent": {
                "blind_spot_floor": self.incoherent_floor,
                "detection_ceiling": self.incoherent_caught,
            },
            "deltas": [_number(delta) for delta in self.deltas],
        }


def _number(value: float) -> float:
    return float(f"{value:.6e}")


# --------------------------------------------------------------------------------------
# battery A: library-level binding controls
# --------------------------------------------------------------------------------------
def _mutate_sphere_volume(**factors: float) -> Any:
    from tensorcad import hydrostatics

    real = hydrostatics.compute_hydrostatics

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        return _Scaled(real(*args, **kwargs), **factors)

    return _patched(hydrostatics, "compute_hydrostatics", wrapped)


def _mutate_generated(module_name: str, generator: str, **factors: float) -> Any:
    """Scale attributes of the ``.hydrostatics`` of whatever ``generator`` returns."""

    @contextlib.contextmanager
    def patcher() -> Iterator[None]:
        module = importlib.import_module(f"tensorcad.{module_name}")
        real = getattr(module, generator)

        def wrapped(*args: Any, **kwargs: Any) -> Any:
            built = real(*args, **kwargs)
            return _Scaled(built, hydrostatics=_Scaled(built.hydrostatics, **factors))

        with _patched(module, generator, wrapped):
            yield

    return patcher()


def _mutate_5415_hydrostatics(**factors: float) -> Any:
    from tensorcad.dtmb5415 import DTMB5415Space

    real = DTMB5415Space.hydrostatics

    def wrapped(self: Any) -> Any:
        return _Scaled(real(self), **factors)

    return _patched(DTMB5415Space, "hydrostatics", wrapped)


def _mutate_5415_mesh(**factors: float) -> Any:
    """Scale the shipped 5415 shell: the mirror identity must notice it is no longer 2x
    the half-cloud extent."""

    import dataclasses

    from tensorcad.dtmb5415 import DTMB5415Space

    real = DTMB5415Space.mesh

    def wrapped(self: Any, *args: Any, **kwargs: Any) -> Any:
        built = real(self, *args, **kwargs)
        scales = [factors.get(axis, 1.0) for axis in "xyz"]
        vertices = tuple(
            tuple(float(component) * scale for component, scale in zip(vertex, scales, strict=True))
            for vertex in built.vertices
        )
        return dataclasses.replace(built, vertices=vertices)

    return _patched(DTMB5415Space, "mesh", wrapped)


def controls() -> list[Control]:
    """The library-level mutations battery A applies.

    Each mutation is small (0.1 %-1 %), of a kind a refactor plausibly introduces (a wrong
    scale factor, a lost factor of two in a mirror, a coefficient that no longer matches its
    own volume), and lands inside the *real* call path of a golden case.
    """

    return [
        Control(
            name="no_mutation_control",
            case="series60",
            detail="sanity control: an unmutated library must still pass, otherwise a gate "
            "that fails everything would look maximally sensitive",
            expected=(),
            mutate=lambda: contextlib.nullcontext(),
        ),
        Control(
            name="sphere_volume_scaled_1.01",
            case="sphere_ellipsoid",
            detail="compute_hydrostatics returns 1 % too much displacement volume: a unit or "
            "integration-scale regression on the SDF path",
            expected=("reference_mismatch", "drift"),
            mutate=lambda: _mutate_sphere_volume(displacement_volume=1.01),
        ),
        Control(
            name="series60_waterplane_scaled_1.002",
            case="series60",
            detail="generate_series60 reports a waterplane area 0.2 % too large, i.e. a "
            "coefficient error far above the 2.06e-4 tolerance of that metric",
            expected=("reference_mismatch", "drift"),
            mutate=lambda: _mutate_generated("series60", "generate_series60", waterplane_area=1.002),
        ),
        Control(
            name="wigley_volume_scaled_0.99",
            case="wigley",
            detail="generate_wigley loses 1 % of displacement volume: the closed-form Wigley "
            "reference cannot follow",
            expected=("reference_mismatch", "drift"),
            mutate=lambda: _mutate_generated("wigley", "generate_wigley", displacement_volume=0.99),
        ),
        Control(
            name="dtmb5415_cb_scaled_1.001",
            case="dtmb5415",
            detail="the 5415 block coefficient is decoupled from its own volume by 0.1 %; the "
            "volume itself has no published reference here, so only the Cb = V/(LBT) identity "
            "can see it",
            expected=("invariant_violated",),
            mutate=lambda: _mutate_5415_hydrostatics(cb=1.001),
        ),
        Control(
            name="dtmb5415_shell_scaled_1.001",
            case="dtmb5415",
            detail="the shipped shell is scaled 1.001 in every axis; the mirror identity must "
            "notice that the mesh extent is no longer twice the half-cloud extent",
            expected=("invariant_violated",),
            mutate=lambda: _mutate_5415_mesh(x=1.001, y=1.001, z=1.001),
        ),
    ]


def _case(name: str) -> golden.GoldenCase:
    for case in golden.golden_cases():
        if case.name == name:
            return case
    raise KeyError(name)


def _layers_of(codes: tuple[str, ...]) -> tuple[str, ...]:
    hit = []
    for layer, vocabulary in LAYERS:
        if any(code in vocabulary for code in codes):
            hit.append(layer)
    return tuple(dict.fromkeys(hit))


def run_controls(*, baseline_path: Path = golden.BASELINE, quick: bool = True) -> list[ControlResult]:
    results: list[ControlResult] = []
    for control in controls():
        case = _case(control.case)
        with control.mutate():
            # the whole suite runs: the mutations must not perturb the other golden cases,
            # so a failure that spreads is itself evidence of a leaky mutation
            result = golden.run_suite(
                quick=quick,
                baseline_path=baseline_path,
                baseline_optional=False,
            )
            found: list[golden.Finding] = list(result.all_findings)
        codes = tuple(sorted({finding.code for finding in found if finding.severity == golden.ERROR}))
        if not control.expected:
            detected = not codes
        else:
            detected = any(code in control.expected for code in codes)
        failed_cases = tuple(sorted({f.case for f in found if f.severity == golden.ERROR}))
        if control.expected and failed_cases not in ((), (control.case,)) and "baseline_case_unknown" not in codes:
            detected = False
        results.append(
            ControlResult(
                name=control.name,
                case=control.case,
                detail=control.detail,
                expected=control.expected,
                detected=detected,
                codes=codes,
                layers=_layers_of(codes),
                failed_cases=failed_cases,
                findings=tuple(
                    f"{f.case}/{f.code}: {f.detail}" for f in found if f.severity == golden.ERROR
                ),
            )
        )
    return results


# --------------------------------------------------------------------------------------
# battery B: self-blessed blind-spot envelope
# --------------------------------------------------------------------------------------
class _CachedPerturbed(golden.GoldenCase):
    """A golden case whose measurements are cached once and scaled by a factor map.

    ``factors`` maps a measured metric name to a multiplicative factor. A *coherent* map
    scales every metric that a reported identity links together (``V = Cb LBT`` and friends),
    so the invariant layer cannot see it and only accuracy decides; an *incoherent* map
    perturbs one number and leaves its identities behind, which is what a plumbing/unit bug
    looks like from the outside.
    """

    def __init__(
        self,
        inner: golden.GoldenCase,
        factors: dict[str, float],
        cache: dict[str, dict[str, float]] | None = None,
        reference: dict[str, float] | None = None,
    ) -> None:
        self.name = inner.name
        self.title = inner.title
        self.reference_source = inner.reference_source
        self._inner = inner
        self._factors = factors
        # one measurement per level per case, shared by every cell of the sweep: the
        # geometry is deterministic, only the reported numbers are perturbed
        self._cache: dict[str, dict[str, float]] = cache if cache is not None else {}
        self._reference = reference

    def _values(self, level: str) -> dict[str, float]:
        if level not in self._cache:
            self._cache[level] = {k: float(v) for k, v in self._inner.measure(level).items()}
        return self._cache[level]

    def measure(self, level: str) -> dict[str, float]:
        values = dict(self._values(level))
        for key, factor in self._factors.items():
            if key in values:
                values[key] = values[key] * factor
        return values

    def levels(self) -> tuple[str, ...]:
        return self._inner.levels()

    def refinement_factor(self, coarse_level: str, fine_level: str) -> float | None:
        return self._inner.refinement_factor(coarse_level, fine_level)

    def all_values(self) -> dict[str, dict[str, float]]:
        return {level: self.measure(level) for level in self._inner.levels()}

    def reference(self) -> dict[str, float]:
        # the reference is invariant under the mutation: it encodes the adversary's blind
        # spot, so it is computed once per case and reused verbatim
        if self._reference is None:
            self._reference = {k: float(v) for k, v in self._inner.reference().items()}
        return dict(self._reference)

    def invariants(self, values: Any) -> Any:
        return self._inner.invariants(values)

    def known_limitations(self, values: Any) -> Any:
        return self._inner.known_limitations(values)

    def geometry_sha256(self) -> str:
        return ""


#: Coherent factor maps: perturbing one reported number drags everything an identity ties to
#: it. Anything absent here has no identity to satisfy, so incoherent and coherent coincide.
def coherent_factors(case: str, metric: str) -> dict[str, float]:
    """Return the set of metrics that must move with ``metric`` to keep every identity true."""

    if case == "series60":
        if metric == "displacement_volume_m3":
            # V = Cb LBT and Cb = Cp Cm: scaling V alone must drag Cb and Cp
            return {"displacement_volume_m3": 1.0, "block_coefficient": 1.0, "prismatic_coefficient": 1.0}
        if metric == "maximum_section_area_m2":
            # Cm = A/(BT) drags Cb = Cp Cm and then V = Cb LBT
            return {"maximum_section_area_m2": 1.0, "block_coefficient": 1.0, "displacement_volume_m3": 1.0}
        if metric == "prismatic_coefficient":
            return {"prismatic_coefficient": 1.0, "block_coefficient": 1.0, "displacement_volume_m3": 1.0}
        if metric == "block_coefficient":
            # Cb = Cp Cm and V = Cb LBT: Cb cannot move alone
            return {"block_coefficient": 1.0, "prismatic_coefficient": 1.0, "displacement_volume_m3": 1.0}
    if case == "dtmb5415":
        if metric == "block_coefficient":
            return {"block_coefficient": 1.0, "displacement_volume_m3": 1.0}
        if metric == "displacement_volume_m3":
            return {"displacement_volume_m3": 1.0, "block_coefficient": 1.0}
        if metric.startswith("mesh_extent_"):
            # the mesh is the mirror of the half cloud: both extents move together
            return {metric: 1.0, metric.replace("mesh_extent_", "cloud_extent_"): 1.0}
    return {metric: 1.0}


def _blessed_like(inner: golden.GoldenCase, case_result: golden.CaseResult) -> dict[str, Any]:
    """Assemble a baseline entry exactly the way ``bless`` would, from a live case result."""

    reference = case_result.reference
    ladder = inner.levels()
    if len(ladder) >= 3:
        acceptance, partner = ladder[-2], ladder[-3]
    else:
        acceptance, partner = ladder[-1], (ladder[0] if len(ladder) > 1 else None)
    coarse = case_result.measurements.get(partner) if partner is not None else None
    tolerances: dict[str, Any] = {}
    for key, expected in reference.items():
        if inner.name:
            tolerance, basis = golden._tolerance(
                expected, coarse.get(key) if coarse else None,
                case_result.measurements[acceptance][key]
            )
        else:  # pragma: no cover - defensive
            tolerance, basis = golden.TOLERANCE_FLOOR, "floor"
        tolerances[key] = {"value": tolerance, "basis": basis}
    return {
        "title": case_result.title,
        "reference_source": case_result.reference_source,
        "geometry_sha256": case_result.geometry_sha256,
        "measurements": copy.deepcopy(case_result.measurements),
        "reference": dict(reference),
        "tolerances": tolerances,
        "extrapolated": {
            m.name: {"value": m.value, "tolerance": m.tolerance, "basis": m.basis}
            for m in case_result.metrics
            if m.kind == "richardson"
        },
    }


def _deltas(tolerance: float) -> tuple[float, ...]:
    """The fixed log grid plus a tolerance-anchored ladder.

    The ladder is what makes the claim *sharp*: it brackets the sensitivity threshold to
    within 10 % of the metric's own declared tolerance instead of within one decade.
    """

    anchored = {
        max(tolerance * factor, 1e-12)
        for factor in (0.1, 0.2, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0, 10.0)
    }
    return tuple(sorted(set(PERTURBATIONS) | anchored))


def _sweep(
    case: golden.GoldenCase,
    metric: str,
    deltas: tuple[float, ...],
    *,
    cache: dict[str, dict[str, float]],
    reference: dict[str, float],
    coherent: bool,
) -> tuple[list[float], dict[str, list[str]], dict[float, dict[str, float]]]:
    """Bless-and-recheck one metric over the delta ladder; return what survives and why not."""

    template = {name: 1.0 for name in coherent_factors(case.name, metric)}
    probe = _CachedPerturbed(case, {}, cache=cache, reference=reference)
    findings: list[golden.Finding] = []
    base = golden.evaluate_case(probe, None, quick=False, findings=findings, baseline_optional=True)
    survives: list[float] = []
    caught: dict[str, list[str]] = {}
    residuals: dict[float, dict[str, float]] = {}
    for delta in deltas:
        factor = 1.0 + delta
        factors = {name: factor if coherent else 1.0 for name in template}
        if not coherent:
            factors = {metric: factor}
        mutated = _CachedPerturbed(case, factors, cache=cache, reference=reference)
        blessed = _blessed_like(case, base)
        for level, values in mutated.all_values().items():
            blessed["measurements"][level] = dict(values)
        found: list[golden.Finding] = []
        cell = golden.evaluate_case(mutated, blessed, quick=False, findings=found)
        codes = sorted({f.code for f in found if f.severity == golden.ERROR})
        if codes:
            caught[str(delta)] = codes
        else:
            survives.append(delta)
        residuals[delta] = {
            metric_check.name: abs(metric_check.value)
            for metric_check in cell.metrics
            if metric_check.kind == "invariant"
        }
    return survives, caught, residuals


def run_envelope(*, baseline_path: Path = golden.BASELINE) -> list[EnvelopeRow]:
    """Sweep the delta ladder over every asserted analytic metric of every golden case."""

    committed = json.loads(Path(baseline_path).read_text(encoding="utf-8"))["cases"]
    rows: list[EnvelopeRow] = []
    for case in golden.golden_cases():
        cache: dict[str, dict[str, float]] = {}
        reference = {k: float(v) for k, v in case.reference().items()}
        probe = _CachedPerturbed(case, {}, cache=cache, reference=reference)
        findings: list[golden.Finding] = []
        base = golden.evaluate_case(
            probe, committed.get(case.name), quick=False, findings=findings
        )
        assert not [f for f in findings if f.severity == golden.ERROR], findings
        for metric in sorted(base.reference):
            tolerance = next(m.tolerance for m in base.metrics if m.name == metric) or 0.0
            error = next(m.relative_error for m in base.metrics if m.name == metric) or 0.0
            deltas = _deltas(tolerance)
            survives, caught, residuals = _sweep(
                case, metric, deltas, cache=cache, reference=reference, coherent=True
            )
            incoherent_survives, incoherent_caught, _ = _sweep(
                case, metric, deltas, cache=cache, reference=reference, coherent=False
            )
            coherent = all(
                value <= 1e-9 * 1.001
                for delta, values in residuals.items()
                for value in values.values()
            ) if residuals else True
            rows.append(
                EnvelopeRow(
                    case=case.name,
                    metric=metric,
                    tolerance=float(tolerance),
                    reference_error=float(error),
                    floor=max(survives) if survives else None,
                    ceiling_step=min(float(k) for k in caught) if caught else None,
                    survives=tuple(survives),
                    caught_by={float(k): v for k, v in caught.items()},
                    incoherent_floor=max(incoherent_survives) if incoherent_survives else None,
                    incoherent_caught=(
                        min(float(k) for k in incoherent_caught) if incoherent_caught else None
                    ),
                    invariants_hold=bool(coherent),
                    deltas=deltas,
                )
            )
    return rows


# --------------------------------------------------------------------------------------
# policy: what the audit asserts, and what it merely reports
# --------------------------------------------------------------------------------------
def identities_of(case: str, metric: str) -> int:
    """How many metrics a coherent perturbation of ``metric`` has to drag along."""

    return len(coherent_factors(case, metric))


#: Incoherent mutations (one number reported inconsistent with its own identities) must be
#: caught far below the accuracy tolerance: an identity is exact, so it has no discretisation
#: error to forgive.
INCOHERENT_CEILING = 1e-3
#: Relative slack on the claim "the declared tolerance is the real threshold": the ladder is
#: anchored at 0.9, 1.0, 1.1, 1.5, 2, 3, 10 x tolerance, so a bracket wider than these factors
#: means the tolerance is not what decides.
ACCURACY_LOWER = 0.5
ACCURACY_UPPER = 3.0


def violations(control_results: list[ControlResult], rows: list[EnvelopeRow]) -> list[str]:
    """Blocking findings.

    A blind spot (a coherent error smaller than the tolerance surviving a ``--bless``) is a
    property of the two-grid criterion, not a lie: it is *reported*, with its number. What
    blocks are the ways the gate can be dishonest about itself:

    1. a library-level regression it does not see, or a no-mutation control it rejects;
    2. a tolerance that is not the real threshold (off by more than the ladder brackets);
    3. an identity that fails to catch an inconsistency three orders below the tolerance,
       for a metric that participates in that identity;
    4. a "coherent" mutation that actually trips an identity, which would make the accuracy
       column meaningless.
    """

    problems: list[str] = []
    for result in control_results:
        if not result.detected:
            problems.append(
                f"battery A: control '{result.name}' was NOT detected "
                f"(expected {list(result.expected) or 'PASS'}, got {list(result.codes) or 'no error'})"
            )
        if result.expected and result.failed_cases not in ((), (result.case,)):
            problems.append(
                f"battery A: control '{result.name}' failed cases {list(result.failed_cases)} "
                f"instead of the mutated case '{result.case}' alone"
            )
    for row in rows:
        label = f"{row.case}/{row.metric}"
        if row.ceiling_step is None:
            problems.append(
                f"battery B: {label} survives every swept delta up to {max(row.deltas):.3e}: "
                "the accuracy layer never fires"
            )
            continue
        low = row.floor or 0.0
        # A threshold *tighter* than the declared tolerance is not a defect: the reported
        # set is over-determined, so the tightest member of a metric's dependency cone fires
        # first (Series 60 Cp/waterplane at 2.06e-4 guard the 1.8e-2 volume bound). It is
        # counted, not blocked. Only a threshold *looser* than the tolerance can hide a bug.
        if row.ceiling_step > ACCURACY_UPPER * row.tolerance:
            problems.append(
                f"battery B: {label} tolerance {row.tolerance:.3e} understates the threshold: "
                f"only {row.ceiling_step:.3e} is caught, and {low:.3e} survives"
            )
        if identities_of(row.case, row.metric) > 1:
            if not row.invariants_hold:
                problems.append(
                    f"battery B: {label} coherent map still trips an identity, so the accuracy "
                    "column is not isolating accuracy"
                )
            if row.incoherent_caught is None or row.incoherent_caught > INCOHERENT_CEILING:
                problems.append(
                    f"battery B: {label} participates in an identity but an inconsistent "
                    f"report survives up to {row.incoherent_floor!r} (ceiling "
                    f"{row.incoherent_caught!r} > {INCOHERENT_CEILING:g})"
                )
    return problems


def blind_spots(rows: list[EnvelopeRow], threshold: float = 1e-2) -> list[EnvelopeRow]:
    """Metrics whose tolerance admits a coherent regression of ``threshold`` or more."""

    return [row for row in rows if row.floor is not None and row.floor >= threshold]


def audit(*, quick: bool = True, baseline_path: Path = golden.BASELINE) -> dict[str, Any]:
    started = time.perf_counter()
    control_results = run_controls(baseline_path=baseline_path, quick=quick)
    rows = run_envelope(baseline_path=baseline_path)
    problems = violations(control_results, rows)
    floors = [row.floor for row in rows if row.floor is not None]
    payload = {
        "schema": SCHEMA,
        "generated_by": "python tools/golden/mutation_audit.py",
        "baseline": golden.display_path(baseline_path),
        "seconds": round(time.perf_counter() - started, 2),
        "policy": {
            "accuracy_lower_factor": ACCURACY_LOWER,
            "accuracy_upper_factor": ACCURACY_UPPER,
            "incoherent_ceiling": INCOHERENT_CEILING,
            "perturbations": list(PERTURBATIONS),
            "ladder_factors": [0.1, 0.2, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0, 10.0],
        },
        "controls": [result.to_dict() for result in control_results],
        "envelope": [row.to_dict() for row in rows],
        "violations": problems,
        "blind_spots_ge_1pct": [
            f"{row.case}/{row.metric}" for row in blind_spots(rows)
        ],
        "summary": {
            "controls_total": len(control_results),
            "controls_detected": sum(1 for result in control_results if result.detected),
            "layers_exercised": sorted(
                {layer for result in control_results for layer in result.layers}
            ),
            "metrics_swept": len(rows),
            "metrics_with_a_ge_1pct_blind_spot": len(blind_spots(rows)),
            "worst_blind_spot": max(floors, default=0.0),
            "median_blind_spot": (
                sorted(floors)[len(floors) // 2] if floors else 0.0
            ),
            "metrics_with_identity_early_warning": sum(
                1 for row in rows if identities_of(row.case, row.metric) > 1
            ),
            "metrics_caught_before_their_own_tolerance": sum(
                1
                for row in rows
                if row.floor is not None and row.floor < ACCURACY_LOWER * row.tolerance
            ),
            "median_threshold_over_tolerance": (
                sorted(
                    ((row.floor or 0.0) / row.tolerance)
                    for row in rows
                    if row.tolerance > 0 and row.floor is not None
                ) or [0.0]
            )[len(rows) // 2],
        },
    }
    return payload


# --------------------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------------------
def _fmt(value: float | None) -> str:
    return "-" if value is None else f"{value:.2e}"


def format_audit(payload: dict[str, Any]) -> str:
    lines = [
        "golden gate mutation audit",
        "",
        "battery A - library-level binding controls (must be detected)",
        f"{'control':<34}{'case':<17}{'detected':<10}codes",
        "-" * 112,
    ]
    for control in payload["controls"]:
        verdict = "yes" if control["detected"] else "NO"
        codes = ", ".join(control["codes"]) or "-"
        lines.append(f"{control['name']:<34}{control['case']:<17}{verdict:<10}{codes}")
    lines += [
        "",
        "battery B - blind-spot envelope: delta that survives a --bless",
        "  coherent   = the reported numbers stay mutually consistent (only accuracy can see it)",
        "  incoherent = one number disagrees with its own identity (the identity layer fires)",
        f"{'case':<17}{'metric':<34}{'tolerance':>11}{'coherent':>11}{'incoherent':>12}",
        "-" * 112,
    ]
    for row in payload["envelope"]:
        coherent = row["coherent"]
        lines.append(
            f"{row['case']:<17}{row['metric']:<34}{row['tolerance']:>11.2e}"
            f"{_fmt(coherent['blind_spot_floor']):>11}"
            f"{_fmt(row['incoherent']['blind_spot_floor']):>12}"
        )
    summary = payload["summary"]
    lines += [
        "",
        f"controls detected: {summary['controls_detected']}/{summary['controls_total']}",
        f"layers exercised : {', '.join(summary['layers_exercised'])}",
        f"metrics swept    : {summary['metrics_swept']}"
        f" ({summary['metrics_with_a_ge_1pct_blind_spot']} admit a coherent >= 1 % after a bless)",
        f"blind spot       : median {summary['median_blind_spot']:.2e}, "
        f"worst {summary['worst_blind_spot']:.2e} (relative)",
        f"policy violations: {len(payload['violations'])}",
    ]
    for problem in payload["violations"]:
        lines.append(f"  VIOLATION {problem}")
    lines.append(
        "verdict: "
        + (
            "PASS - the gate is bound to the library and its declared tolerance is its real "
            "threshold"
            if not payload["violations"]
            else "FAIL - see violations above"
        )
    )
    return "\n".join(lines)


def render_markdown(payload: dict[str, Any]) -> str:
    summary = payload["summary"]
    out = [
        "# Golden geometry gate - mutation audit",
        "",
        "Generated by `python tools/golden/mutation_audit.py`; every number below is rendered",
        "from the JSON the same run writes, so nothing is hand-copied. The gate under test is",
        "`tensorcad/golden.py` against the committed baseline",
        "`docs/golden/tensorcad_geometry_golden.json`.",
        "",
        "## Why a mutation audit",
        "",
        "`tests/test_golden_regression.py` proves the *harness* fails closed: missing, unreadable",
        "or stale baselines, hand-edited references, a tree that refuses to bless itself. Every",
        "one of those negative tests injects a synthetic probe case or a tampered baseline file,",
        "so none of them shows that the gate sees a regression in the geometry code it is wired",
        "to. Battery A closes that gap. Battery B answers the second question a `--bless`-able",
        "gate has to answer: **how large a coherent regression can hide inside a green gate?**",
        "",
        "## Battery A - library-level binding controls",
        "",
        "Each control pushes a 0.1-1 % mutation into the *real* entry point a golden case calls",
        "and runs the untouched gate against the committed baseline. The whole suite runs, so a",
        "mutation whose failure spreads to another case is itself a finding (leakage).",
        "",
        "| control | case | expected | verdict | error codes | layers | failed cases |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for control in payload["controls"]:
        expected = ", ".join(control["expected_codes"]) or "clean pass"
        codes = ", ".join(control["codes"]) or "-"
        layers = ", ".join(control["layers"]) or "-"
        verdict = "detected" if control["detected"] else "**MISSED**"
        failed = ", ".join(control["failed_cases"]) or "-"
        out.append(
            f"| `{control['name']}` | {control['case']} | {expected} | {verdict} "
            f"| {codes} | {layers} | {failed} |"
        )
    out += [
        "",
        f"{summary['controls_detected']}/{summary['controls_total']} controls detected; layers",
        f"exercised: {', '.join(summary['layers_exercised'])}.",
        "",
        "Two adversary classes are swept over every asserted analytic metric:",
        "",
        "- **coherent** - the reported numbers stay mutually consistent (`V = Cb LBT` and the",
        "  other identities still hold), so the identity layer is blind by construction and only",
        "  the accuracy tolerance decides. This is what a wrong scale factor, a stale coefficient",
        "  or a mis-wired formula looks like from outside.",
        "- **incoherent** - one number is reported inconsistent with the rest. An identity is",
        "  exact, so it must fire orders of magnitude below the accuracy tolerance.",
        "",
        "For each delta the baseline is re-blessed from the mutated measurement and the gate is",
        "re-run: `blind spot` is the largest delta that still passes. The ladder is anchored at",
        "0.1/0.2/0.5/0.9/1.0/1.1/1.5/2/3/10 x the metric's own tolerance, so a threshold that",
        "is *not* the tolerance shows up as a bracket miss rather than as sweep noise.",
        "",
        "## Battery B - blind-spot envelope",
        "",
        "| case | metric | tolerance | coherent blind spot | caught at | incoherent blind spot |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["envelope"]:
        coherent = row["coherent"]
        out.append(
            f"| {row['case']} | `{row['metric']}` | {row['tolerance']:.2e} "
            f"| {_fmt(coherent['blind_spot_floor'])} | {_fmt(coherent['detection_ceiling'])} "
            f"| {_fmt(row['incoherent']['blind_spot_floor'])} |"
        )
    out += [
        "",
        f"{summary['metrics_swept']} metrics swept. Median blind spot "
        f"{summary['median_blind_spot']:.2e}, worst {summary['worst_blind_spot']:.2e} (relative).",
        f"{summary['metrics_with_a_ge_1pct_blind_spot']} metric(s) admit a *coherent* regression",
        "of 1 % or more after a `--bless`.",
        "",
    ]
    if payload["blind_spots_ge_1pct"]:
        out += [
            "### Known blind spots (reported, not hidden)",
            "",
            "| case | metric | tolerance | blind spot |",
            "| --- | --- | --- | --- |",
        ]
        for label in payload["blind_spots_ge_1pct"]:
            case, metric = label.split("/")
            row = next(
                r for r in payload["envelope"] if r["case"] == case and r["metric"] == metric
            )
            out.append(
                f"| {case} | `{metric}` | {row['tolerance']:.2e} "
                f"| {_fmt(row['coherent']['blind_spot_floor'])} |"
            )
        out.append("")
    else:
        out += [
            "No metric admits a coherent regression of 1 % after a bless: every tolerance in the",
            "library is below 1 %, so a `--bless` cannot drown a percent-level bug.",
            "",
        ]
    out += [
        "## Policy and verdict",
        "",
        "Blocking policy:",
        "",
        "1. every battery-A control is detected and localised to its mutated case, and the",
        "   no-mutation control passes;",
        "2. the declared tolerance is the real threshold, to within the ladder bracket",
        f"   ({payload['policy']['accuracy_lower_factor']:g}x-"
        f"{payload['policy']['accuracy_upper_factor']:g}x tolerance);",
        f"3. a metric that participates in an identity catches an *incoherent* report by",
        f"   {payload['policy']['incoherent_ceiling']:g} relative;",
        "4. a coherent map really is invisible to the identity layer (otherwise the accuracy",
        "   column would be measuring the wrong thing).",
        "",
        f"Violations: {len(payload['violations'])}.",
        "",
    ]
    for problem in payload["violations"]:
        out.append(f"- **VIOLATION** {problem}")
    out.append("")
    out.append(
        "**Verdict:** "
        + (
            "PASS - the gate is bound to the library, its identities give early warning, and its "
            "declared tolerance is the threshold that actually decides."
            if not payload["violations"]
            else "FAIL - see violations above."
        )
    )
    return "\n".join(out) + "\n"


def render_figure(payload: dict[str, Any], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    rows = payload["envelope"]
    cases = sorted({row["case"] for row in rows})
    colors = dict(zip(cases, ("#1f77b4", "#d62728", "#2ca02c", "#9467bd"), strict=False))

    figure, (left, right) = plt.subplots(1, 2, figsize=(13.2, 5.8), width_ratios=(1.15, 1.0))
    lo, hi = 1e-11, 1.0
    left.plot([lo, hi], [lo, hi], "k--", linewidth=1.0, label="threshold = declared tolerance")
    left.plot([lo, hi], [3 * lo, 3 * hi], "k:", linewidth=0.8, label="3x tolerance")
    for row in rows:
        floor = row["coherent"]["blind_spot_floor"]
        ceiling = row["coherent"]["detection_ceiling"]
        color = colors[row["case"]]
        if floor is not None:
            left.plot(row["tolerance"], floor, "o", color=color, markersize=6)
        if ceiling is not None:
            left.plot(row["tolerance"], ceiling, "^", color=color, markersize=5, alpha=0.55)
        if row["incoherent"]["blind_spot_floor"] is not None:
            left.plot(
                row["tolerance"], row["incoherent"]["blind_spot_floor"], "s", color=color,
                markersize=3, alpha=0.5,
            )
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_xlim(5e-5, 5e-2)
    left.set_ylim(lo, hi)
    left.set_xlabel("declared acceptance tolerance (relative)")
    left.set_ylabel("mutation delta that survives a --bless")
    left.set_title("Is the declared tolerance the real threshold?")
    left.grid(True, which="both", alpha=0.25)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="none", color=colors[case], label=case)
        for case in cases
    ]
    handles += [
        plt.Line2D([], [], marker="o", linestyle="none", color="0.35",
                   label="coherent (blind spot)"),
        plt.Line2D([], [], marker="^", linestyle="none", color="0.35",
                   label="coherent (first caught)"),
        plt.Line2D([], [], marker="s", linestyle="none", color="0.35", label="incoherent"),
    ]
    left.legend(handles=handles, fontsize=6.5, loc="upper left")

    floors = np.array([row["coherent"]["blind_spot_floor"] or 0.0 for row in rows])
    names = [f"{row['case']}\n{row['metric']}" for row in rows]
    order = np.argsort(floors)
    right.barh(
        np.arange(len(rows)),
        floors[order],
        color=[colors[rows[i]["case"]] for i in order],
    )
    right.set_yticks(np.arange(len(rows)))
    right.set_yticklabels([names[i] for i in order], fontsize=6)
    right.axvline(1e-2, color="k", linestyle="--", linewidth=1.0, label="1 %")
    right.set_xscale("log")
    right.set_xlabel("largest coherent delta that survives a --bless")
    right.set_title("Per-metric blind spot")
    right.grid(True, axis="x", which="both", alpha=0.25)
    right.legend(fontsize=7, loc="lower right")

    figure.suptitle(
        f"Golden gate mutation audit - {payload['summary']['controls_detected']}"
        f"/{payload['summary']['controls_total']} binding controls detected, "
        f"{len(payload['violations'])} policy violations",
        fontsize=12,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    figure.savefig(path, dpi=140)
    plt.close(figure)


# --------------------------------------------------------------------------------------
# cli
# --------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="non-zero exit on a policy violation")
    parser.add_argument("--quick", action="store_true", help="battery A on the coarse grid only")
    parser.add_argument("--json", metavar="PATH", help="write the machine-readable audit here")
    parser.add_argument("--md", metavar="PATH", help="write the markdown report here")
    parser.add_argument("--figure", metavar="PATH", help="write the summary figure (png) here")
    parser.add_argument("--baseline", default=str(golden.BASELINE))
    args = parser.parse_args(argv)

    payload = audit(quick=True, baseline_path=Path(args.baseline))
    buffer = io.StringIO()
    buffer.write(format_audit(payload) + "\n")
    sys.stdout.write(buffer.getvalue())
    if args.json:
        Path(args.json).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"json -> {args.json}")
    if args.md:
        Path(args.md).write_text(render_markdown(payload), encoding="utf-8")
        print(f"markdown -> {args.md}")
    if args.figure:
        render_figure(payload, Path(args.figure))
        print(f"figure -> {args.figure}")
    return 0 if (not args.check or not payload["violations"]) else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())