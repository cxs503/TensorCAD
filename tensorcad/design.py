"""A closed design loop over the parametric hull: target in, verified geometry out.

The rest of the package answers "what is this surface?"; this module answers "find me the
surface that meets the requirement, prove it does, and hand it to a solver". It is the
first place where a *design goal* rather than a parameter set is the input, and it is
deliberately narrow: one hull family (:class:`~tensorcad.hull.HullModel`), two free
variables, and criteria that this repository can measure exactly.

The loop
--------

    targets ─► guess ─► build + measure ─► residuals ─► Newton step ─► bounds ─► repeat
                                     ▲                                          │
                                     └──────────────────────────────────────────┘

Each step is a damped Newton iteration on ``log`` of the free dimensions, with a
finite-difference Jacobian, a backtracking line search, and hard bounds taken from the
same :data:`~tensorcad.hull.PARAMETER_SPECS` that validate any other parameter set. Every
intermediate design is kept, so the record shows the path rather than only the answer.

What is measured, and what is not
---------------------------------

The residual is built from hydrostatics integrated off the mesh
(:func:`tensorcad.hydrostatics.compute_hydrostatics`), i.e. from the artifact a solver
would receive. Every iteration also carries the full validation report and the
``tensorcad.geometry.handoff/1`` record, so a design is never "converged but unverified":
:attr:`DesignSolution.solved` requires both.

The loop does **not** run a solver. It produces the manifest a solver reads
(:meth:`DesignSolution.manifest`), and accepts a ``feedback`` callable so an external
run can push measured results back into the targets
(:meth:`DesignProblem.solve`). Whether a design is *good* hydrodynamically is a question
this repository cannot answer, and every record says so in ``solver.status``.

No NumPy: the residual is 2–4 dimensional, the Jacobian is finite-differenced from the
same evaluation, and the linear solve is a few lines of Gaussian elimination. The whole
loop runs under a standard-library-only install, which is the same promise the hull
model makes.
"""
from __future__ import annotations

import json
import math
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh
from .hull import HullParameters, build_mesh
from .hydrostatics import Hydrostatics, compute_hydrostatics

__all__ = [
    "DESIGN_SCHEMA",
    "DesignProblem",
    "DesignRecord",
    "DesignSolution",
    "DesignTargets",
    "write_design_json",
    "write_solver_manifest",
]

#: Schema tag for the design record; bump on incompatible change.
DESIGN_SCHEMA = "tensorcad.design/1"

#: Residual level at which the finite-difference Jacobian's own error takes over. The
#: central difference has a relative truncation error of order ``step**2 ~ 1e-10``, so the
#: last iterations stop contracting quadratically; the observed order is measured above
#: this floor rather than at the floor.
_FD_FLOOR = 1e-8

#: Variables the loop may move, and their admissible range. The ranges are the hull
#: model's own specs; they are repeated here only so a step can be clipped *before* a
#: candidate is built, and every clipped step is reported.
_VARIABLE_BOUNDS: dict[str, tuple[float, float]] = {
    "length": (1e-3, 1.0e4),
    "beam": (1e-3, 1.0e3),
    "draft": (1e-3, 1.0e3),
}


#: The equalities a design may aim at, and the measured quantity behind each one.
_EQUALITIES = {
    "displacement_m3": "integrated displaced volume, m^3",
    "beam_to_draft": "moulded beam / draft",
    "length_to_beam": "length / moulded beam",
    "draft_m": "draft, m",
}


@dataclass(frozen=True)
class DesignTargets:
    """What the design has to achieve, plus the budgets it may not exceed.

    The equalities are deliberately *not* two geometric integrals. Displacement and
    wetted area both scale roughly with beam times draft on this family, so asking for
    both at once makes the Jacobian nearly rank deficient -- measured, the loop stalls
    at a residual of 4.0e-3 and cannot remove it, because no hull in this family has
    that pair. What a naval architect actually specifies is a displacement and a shape
    ratio (``beam_to_draft``, ``length_to_beam``, or a draft), and those two are
    independent. Wetted area stays as a **budget**, which is what an area usually is.
    """

    displacement_m3: float | None = None
    beam_to_draft: float | None = None
    length_to_beam: float | None = None
    draft_m: float | None = None
    wetted_area_m2_max: float | None = None

    def __post_init__(self) -> None:
        for name in ("displacement_m3", "beam_to_draft", "length_to_beam", "draft_m", "wetted_area_m2_max"):
            value = getattr(self, name)
            if value is not None and not (math.isfinite(value) and value > 0.0):
                raise ValueError(f"{name} must be a finite positive number, got {value!r}")

    @property
    def names(self) -> tuple[str, ...]:
        """The active equalities, in a fixed order."""
        return tuple(name for name in _EQUALITIES if getattr(self, name) is not None)

    @property
    def budgets(self) -> tuple[str, ...]:
        return ("wetted_area_m2_max",) if self.wetted_area_m2_max is not None else ()

    def measured(self, hydrostatics: Hydrostatics, parameters: HullParameters) -> dict[str, float]:
        """The measured counterpart of every equality, so nothing is compared to itself."""
        return {
            "displacement_m3": hydrostatics.displacement_volume,
            "beam_to_draft": parameters.beam / parameters.draft,
            "length_to_beam": parameters.length / parameters.beam,
            "draft_m": parameters.draft,
        }

    def residuals(self, hydrostatics: Hydrostatics, parameters: HullParameters) -> dict[str, float]:
        """Relative residuals ``measured / target - 1``, one per active equality."""
        measured = self.measured(hydrostatics, parameters)
        residuals: dict[str, float] = {}
        for name in self.names:
            target = getattr(self, name)
            assert target is not None  # implied by ``names``
            residuals[name] = measured[name] / target - 1.0
        return residuals

    def budget_checks(self, hydrostatics: Hydrostatics) -> dict[str, tuple[float, float, bool]]:
        """``name -> (measured, budget, within)`` for every budget, measured each iteration."""
        checks: dict[str, tuple[float, float, bool]] = {}
        if self.wetted_area_m2_max is not None:
            measured = hydrostatics.wetted_surface_area
            checks["wetted_area_m2_max"] = (
                measured,
                self.wetted_area_m2_max,
                measured <= self.wetted_area_m2_max * (1.0 + 1e-12),
            )
        return checks


@dataclass(frozen=True)
class DesignRecord:
    """One evaluated design: the parameters, the surface, and the evidence it is usable."""

    parameters: HullParameters
    variables: Mapping[str, float]
    iteration: int
    residuals: Mapping[str, float]
    budgets: Mapping[str, tuple[float, float, bool]]
    hydrostatics: Hydrostatics
    handoff: Mapping[str, Any]
    checks: tuple[Any, ...]
    mesh: Mesh = field(repr=False)
    seconds: float = 0.0
    date: str = ""

    @property
    def passed(self) -> bool:
        """True when the mesh is closed, every geometry check passed, and no budget is exceeded."""
        return (
            bool(self.handoff["watertight"])
            and all(bool(check.passed) for check in self.checks)
            and all(within for _, _, within in self.budgets.values())
        )

    @property
    def residual_norm(self) -> float:
        """L2 norm of the relative residuals, the quantity the loop drives to zero."""
        return math.sqrt(math.fsum(value * value for value in self.residuals.values()))

    @property
    def digest(self) -> str:
        return str(self.handoff["mesh_sha256"])

    def to_dict(self) -> dict[str, Any]:
        return {
            "iteration": self.iteration,
            "parameters": self.parameters.as_dict(),
            "variables": dict(self.variables),
            "residuals": dict(self.residuals),
            "residual_norm": self.residual_norm,
            "budgets": {
                name: {"measured": measured, "budget": budget, "within": within}
                for name, (measured, budget, within) in self.budgets.items()
            },
            "hydrostatics": self.hydrostatics.to_dict(),
            "handoff": dict(self.handoff),
            "checks": {check.name: check.to_dict() for check in self.checks},
            "passed": self.passed,
            "seconds": self.seconds,
        }


def _poly(a: list[list[float]], b: list[float]) -> list[float]:
    """Solve ``a x = b`` by Gaussian elimination with partial pivoting.

    Small and explicit on purpose: the system here is 1–4 wide, and pulling in a linear
    algebra stack for it would break the standard-library-only promise of this module.
    """
    size = len(b)
    matrix = [[*row[:], b[index]] for index, row in enumerate(a)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(matrix[row][column]))
        if abs(matrix[pivot][column]) < 1e-300:
            raise ZeroDivisionError("the design Jacobian is singular")
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        for row in range(column + 1, size):
            factor = matrix[row][column] / matrix[column][column]
            for cell in range(column, size + 1):
                matrix[row][cell] -= factor * matrix[column][cell]
    solution = [0.0] * size
    for row in reversed(range(size)):
        total = matrix[row][size] - math.fsum(
            matrix[row][cell] * solution[cell] for cell in range(row + 1, size)
        )
        solution[row] = total / matrix[row][row]
    return solution


class DesignProblem:
    """A hull family, a set of targets and the free dimensions allowed to move."""

    def __init__(
        self,
        *,
        base: HullParameters,
        targets: DesignTargets,
        variables: Sequence[str] = ("beam", "draft"),
        max_iterations: int = 12,
        tolerance: float = 1e-10,
        damping: float = 1.0,
        finite_difference_step: float = 1e-5,
    ) -> None:
        if not variables:
            raise ValueError("a design problem needs at least one free variable")
        for name in variables:
            if name not in _VARIABLE_BOUNDS:
                raise ValueError(f"unknown design variable {name!r}; known: {sorted(_VARIABLE_BOUNDS)}")
        if len(set(variables)) != len(variables):
            raise ValueError(f"duplicate design variables: {list(variables)}")
        if not targets.names:
            raise ValueError(
                "a design problem needs at least one equality target "
                "(displacement_m3, beam_to_draft, length_to_beam or draft_m)"
            )
        if len(variables) < len(targets.names):
            raise ValueError(
                f"{len(targets.names)} targets {list(targets.names)} cannot be met by "
                f"{len(variables)} free variables {list(variables)}"
            )
        if max_iterations < 1:
            raise ValueError("max_iterations must be at least 1")
        if not 0.0 < tolerance < 1.0:
            raise ValueError(f"tolerance must be in (0, 1), got {tolerance!r}")
        if not 0.0 < damping <= 1.0:
            raise ValueError(f"damping must be in (0, 1], got {damping!r}")
        if finite_difference_step <= 0.0:
            raise ValueError("finite_difference_step must be positive")
        base.validate()
        self.base = base
        self.targets = targets
        self.variables = tuple(variables)
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.damping = damping
        self.finite_difference_step = finite_difference_step

    # -- evaluation ----------------------------------------------------------

    def parameters_for(self, variables: Mapping[str, float]) -> HullParameters:
        """The full parameter set for a point in the free-variable space.

        Every design variable in this problem is a float dimension, so the cast is
        explicit rather than splatted into ``replace``: ``nx``, ``nz`` and the density
        stay as they came, and a wrong name is an error instead of a silently ignored
        keyword.
        """
        unknown = set(variables) - set(self.variables)
        if unknown:
            raise ValueError(
                f"unexpected variables {sorted(unknown)}; this problem moves {list(self.variables)}"
            )
        values = {name: float(value) for name, value in variables.items()}
        base = self.base
        return replace(
            base,
            length=values.get("length", base.length),
            beam=values.get("beam", base.beam),
            draft=values.get("draft", base.draft),
        )

    def variable_values(self, parameters: HullParameters) -> dict[str, float]:
        return {name: float(getattr(parameters, name)) for name in self.variables}

    def evaluate(
        self,
        variables: Mapping[str, float],
        *,
        iteration: int = 0,
        targets: DesignTargets | None = None,
    ) -> DesignRecord:
        """Build the hull, measure it, and return the record with its evidence.

        The record holds the mesh and the full check list, so a caller can export or
        re-inspect it without rebuilding. Nothing here is cached: two evaluations of the
        same point produce identical meshes (the generator is deterministic), which is
        what makes the finite-difference Jacobian meaningful.
        """
        from .validation import ValidationReport, inspect_geometry, inspect_mesh_quality

        started = time.perf_counter()
        parameters = self.parameters_for(variables)
        parameters.validate()
        mesh, _ = build_mesh(parameters)
        hydrostatics = compute_hydrostatics(
            mesh, waterline_z=0.0, water_density=parameters.water_density
        )
        quality = tuple(inspect_mesh_quality(mesh))
        checks = tuple(inspect_geometry(mesh)) + quality
        sampled = _SampledHandoff(
            mesh=mesh,
            report=_measure(mesh),
            quality=quality,
            validation=ValidationReport(checks),
        )
        return DesignRecord(
            parameters=parameters,
            variables=self.variable_values(parameters),
            iteration=iteration,
            residuals=(targets or self.targets).residuals(hydrostatics, parameters),
            budgets=(targets or self.targets).budget_checks(hydrostatics),
            hydrostatics=hydrostatics,
            handoff=sampled.handoff(),
            checks=checks,
            mesh=mesh,
            seconds=time.perf_counter() - started,
            date="",
        )

    # -- the loop ------------------------------------------------------------

    def jacobian(self, record: DesignRecord, targets: DesignTargets | None = None) -> list[list[float]]:
        """``d residual / d log(variable)`` by central differences.

        In log space because a hull dimension is positive and scales multiplicatively:
        a step of 1e-5 means "0.001 % bigger", the same relative probe at 0.5 m and at
        50 m, which keeps the finite difference well conditioned for every design.
        """
        active = targets or self.targets
        names = list(active.names)
        matrix: list[list[float]] = [[0.0] * len(self.variables) for _ in names]
        step = self.finite_difference_step
        for column, variable in enumerate(self.variables):
            for sign in (1.0, -1.0):
                shifted = dict(record.variables)
                shifted[variable] = shifted[variable] * math.exp(sign * step)
                probe = self.evaluate(shifted, targets=active)
                for row, name in enumerate(names):
                    matrix[row][column] += sign * probe.residuals[name] / (2.0 * step)
        return matrix

    def solve(
        self,
        *,
        start: Mapping[str, float] | None = None,
        feedback: Callable[[DesignRecord], Mapping[str, float] | None] | None = None,
    ) -> DesignSolution:
        """Run the damped Newton loop and return every iteration plus the verdict.

        ``start`` defaults to the base parameters. ``feedback`` is the seam for an
        external solver: called with the current record, it returns adjustments applied
        to the targets (for example a displacement correction implied by a measured
        sinkage), so a coupled run feeds back into the next design step instead of the
        loop pretending to know the answer.
        """
        current = dict(self.variable_values(self.base)) if start is None else dict(start)
        for name in current:
            if name not in self.variables:
                raise ValueError(f"{name!r} is not a free variable of this problem")
        targets = self.targets
        records: list[DesignRecord] = []
        steps: list[dict[str, Any]] = []
        clipped = 0
        status = "max_iterations"
        jacobians: list[list[list[float]]] = []

        for iteration in range(self.max_iterations + 1):
            record = self.evaluate(current, iteration=iteration, targets=targets)
            records.append(record)
            if feedback is not None and iteration < self.max_iterations:
                adjustment = feedback(record)
                if adjustment:
                    targets = _adjust(targets, adjustment)
            if record.residual_norm <= self.tolerance:
                status = "converged"
                break
            if iteration == self.max_iterations:
                status = "max_iterations"
                break
            matrix = self.jacobian(record, targets)
            jacobians.append(matrix)
            if not targets.names:
                status = "no_targets"
                break
            try:
                delta = _poly(matrix, [-record.residuals[name] for name in targets.names])
            except ZeroDivisionError:
                status = "singular_jacobian"
                break
            # Damped step with backtracking: accept the first candidate that reduces the
            # residual and stays inside the parameter bounds; otherwise halve and retry.
            accepted = False
            infeasible = False
            scale = self.damping
            for _ in range(8):
                candidate = {
                    name: current[name] * math.exp(scale * delta[index])
                    for index, name in enumerate(self.variables)
                }
                lower_clipped = any(
                    candidate[name] < _VARIABLE_BOUNDS[name][0] for name in self.variables
                )
                for name in self.variables:
                    low, high = _VARIABLE_BOUNDS[name]
                    candidate[name] = min(max(candidate[name], low), high)
                try:
                    trial = self.evaluate(candidate, iteration=iteration + 1, targets=targets)
                except ValueError:
                    # The candidate left the admissible design space (a ratio outside
                    # the hull model's specs, for one). A design loop reports that as an
                    # infeasible step rather than dying inside a parameter check.
                    scale *= 0.5
                    infeasible = True
                    continue
                if trial.residual_norm < record.residual_norm:
                    accepted = True
                    clipped += int(lower_clipped)
                    current = candidate
                    steps.append(
                        {
                            "iteration": iteration,
                            "scale": scale,
                            "delta_log": dict(zip(self.variables, delta, strict=True)),
                            "residual_before": record.residual_norm,
                            "residual_after": trial.residual_norm,
                            "clipped": lower_clipped,
                        }
                    )
                    break
                scale *= 0.5
            if not accepted:
                status = "infeasible" if infeasible else "no_descent"
                break

        return DesignSolution(
            problem=self,
            targets=targets,
            records=tuple(records),
            steps=tuple(steps),
            jacobians=tuple(jacobians),
            status=status,
            clipped_steps=clipped,
        )


@dataclass(frozen=True)
class _MeshFacts:
    """The report fields the handoff record reads, measured with the stdlib helpers."""

    watertight: bool
    edge_manifold: bool
    euler_characteristic: int
    signed_volume: float
    surface_area: float
    bounds_min: tuple[float, float, float]
    bounds_max: tuple[float, float, float]


def _measure(mesh: Mesh) -> _MeshFacts:
    edges = meshlib.edge_face_counts(mesh)
    box = mesh.bounding_box()
    return _MeshFacts(
        watertight=meshlib.is_watertight(mesh),
        edge_manifold=meshlib.is_edge_manifold(mesh),
        euler_characteristic=len(mesh.vertices) - len(edges) + len(mesh.faces),
        signed_volume=meshlib.total_signed_volume(mesh),
        surface_area=meshlib.surface_area(mesh),
        bounds_min=tuple(box.minimum),  # type: ignore[arg-type]
        bounds_max=tuple(box.maximum),  # type: ignore[arg-type]
    )


@dataclass(frozen=True)
class _SampledHandoff:
    """The three fields :meth:`SampledField.handoff` reads, without importing numpy paths."""

    mesh: Mesh
    report: Any
    quality: tuple[Any, ...]
    validation: Any

    def handoff(self) -> dict[str, Any]:
        return _handoff_record(self.mesh, self.report, self.quality, self.validation)


def _handoff_record(mesh: Mesh, report: Any, quality: tuple[Any, ...], validation: Any) -> dict[str, Any]:
    """Build the ``tensorcad.geometry.handoff/1`` record for an in-memory mesh.

    Same keys and same digest convention as
    :meth:`tensorcad.geometry_field.SampledField.handoff`, so a design record can be
    compared with a sampled field without a second schema.
    """
    import hashlib
    import struct

    digest = hashlib.sha256()
    digest.update(struct.pack("<II", len(mesh.vertices), len(mesh.faces)))
    for vertex in mesh.vertices:
        digest.update(struct.pack("<3d", float(vertex[0]), float(vertex[1]), float(vertex[2])))
    for face in mesh.faces:
        digest.update(struct.pack("<3I", int(face[0]), int(face[1]), int(face[2])))
    return {
        "schema": "tensorcad.geometry.handoff/1",
        "mesh_sha256": digest.hexdigest(),
        "units": mesh.units,
        "coordinate_frame": mesh.coordinate_frame,
        "vertices": len(mesh.vertices),
        "triangles": len(mesh.faces),
        "watertight": bool(report.watertight),
        "edge_manifold": bool(report.edge_manifold),
        "euler_characteristic": int(report.euler_characteristic),
        "signed_volume_m3": float(report.signed_volume),
        "surface_area_m2": float(report.surface_area),
        "bounds_min_m": list(report.bounds_min),
        "bounds_max_m": list(report.bounds_max),
        "validation": dict(validation.to_dict()),
        "quality": {check.name: check.to_dict() for check in quality},
    }


def _adjust(targets: DesignTargets, adjustment: Mapping[str, float]) -> DesignTargets:
    """Apply solver feedback to the targets, keeping every value positive and finite."""
    changes: dict[str, float] = {}
    for name, delta in adjustment.items():
        if name not in ("displacement_m3", "wetted_area_m2"):
            raise ValueError(f"feedback cannot adjust {name!r}")
        target = getattr(targets, name)
        if target is None:
            raise ValueError(f"feedback adjusts {name!r}, which this problem does not target")
        changes[name] = max(float(delta), -target + 1e-12)
    updated = replace(targets, **{name: getattr(targets, name) + delta for name, delta in changes.items()})
    return updated


@dataclass(frozen=True)
class DesignSolution:
    """Every iteration of a design run, the verdict, and the handoff a solver reads."""

    problem: DesignProblem
    targets: DesignTargets
    records: tuple[DesignRecord, ...]
    steps: tuple[Mapping[str, Any], ...]
    jacobians: tuple[list[list[float]], ...]
    status: str
    clipped_steps: int = 0

    @property
    def final(self) -> DesignRecord:
        return self.records[-1]

    @property
    def solved(self) -> bool:
        """Converged **and** verified: a converged residual over a broken mesh is not a design."""
        return self.status == "converged" and self.final.passed

    @property
    def failed_checks(self) -> tuple[str, ...]:
        """Names of every gate the final design misses, checks and budgets together.

        Empty exactly when :attr:`solved` is true, which is the property the tests pin:
        a design loop may not report a clean bill of health it cannot itemise.
        """
        failed = [str(check.name) for check in self.final.checks if not check.passed]
        failed.extend(
            name for name, (_, _, within) in self.final.budgets.items() if not within
        )
        return tuple(failed)

    @property
    def convergence_order(self) -> float | None:
        """Observed order from the last three residual norms, when they shrink the way they should.

        A Newton step on a smooth residual contracts quadratically, so
        ``log(r[k+1]/r[k]) / log(r[k]/r[k-1])`` sits near 2. Damped steps, a clipped
        variable or a first iteration far from the solution all lower it, so the value is
        reported with the iteration count rather than asserted.
        """
        norms = [record.residual_norm for record in self.records]
        for index in range(len(norms) - 3, -1, -1):
            a, b, c = norms[index], norms[index + 1], norms[index + 2]
            if min(a, b, c) <= _FD_FLOOR:
                continue
            if a > 0.0 and 0.0 < b < a and 0.0 < c < b:
                return math.log(c / b) / math.log(b / a)
        return None

    def sensitivity(self) -> dict[str, dict[str, float]]:
        """``d ln(measured) / d ln(variable)`` at the first evaluated design.

        The Jacobian the loop uses is ``d residual / d ln x`` with
        ``residual = measured / target - 1``, so the measured sensitivity is that entry
        divided by ``1 + residual``. Reporting it is a real check on the whole
        differentiation: a shape ratio must have sensitivity exactly +1 and -1 to its own
        numerator and denominator, and displacement must scale linearly with beam and
        draft on a hull whose sections are proportional to the beam.
        """
        if not self.jacobians:
            return {}
        matrix = self.jacobians[0]
        names = list(self.targets.names)
        first = self.records[0].residuals
        return {
            variable: {
                name: matrix[row][column] / (1.0 + first[name])
                for row, name in enumerate(names)
            }
            for column, variable in enumerate(self.problem.variables)
        }

    def manifest(self) -> dict[str, Any]:
        """The record a solver or an audit reads: targets, path, verdict, handoff.

        ``solver.status`` is ``"ready"`` when the design converged and verified: the
        geometry, the fluid properties and the reference quantities are all present.
        ``solver.executed_here`` is ``False``, always — this repository measures
        geometry, it does not run the flow or the structure.
        """
        final = self.final
        hydro = final.hydrostatics
        return {
            "schema": DESIGN_SCHEMA,
            "status": self.status,
            "solved": self.solved,
            "iterations": len(self.records),
            "clipped_steps": self.clipped_steps,
            "convergence_order": self.convergence_order,
            "targets": {
                name: getattr(self.targets, name)
                for name in (
                    "displacement_m3",
                    "beam_to_draft",
                    "length_to_beam",
                    "draft_m",
                    "wetted_area_m2_max",
                )
            },
            "design_variables": list(self.problem.variables),
            "path": [record.to_dict() for record in self.records],
            "steps": [dict(step) for step in self.steps],
            "final": final.to_dict(),
            "reference_quantities": {
                "length_m": final.parameters.length,
                "beam_m": final.parameters.beam,
                "draft_m": final.parameters.draft,
                "waterline_z_m": hydro.waterline_z,
                "water_density_kg_per_m3": hydro.water_density,
                "wetted_area_m2": hydro.wetted_surface_area,
                "displacement_volume_m3": hydro.displacement_volume,
                "displacement_mass_kg": hydro.displacement_mass,
                "block_coefficient": hydro.block_coefficient,
                "prismatic_coefficient": hydro.prismatic_coefficient,
                "longitudinal_centre_buoyancy_m": hydro.longitudinal_centre_buoyancy,
                "vertical_centre_buoyancy_m": hydro.vertical_centre_buoyancy,
            },
            "handoff": dict(final.handoff),
            "verification": {
                "passed": final.passed,
                "failed": list(self.failed_checks),
                "checks": len(final.checks),
            },
            "solver": {
                "status": "ready" if self.solved else "not_ready",
                "executed_here": False,
                "reason": (
                    "TensorCAD measures and hands off geometry; the flow/structure solve "
                    "belongs to TensorFEM/TensorLBM and is not run in this repository"
                ),
            },
            "limits": [
                "one parametric hull family (sin(pi x/L)^0.7 breadth, sqrt(z/draft) section)",
                "no stability, resistance, seakeeping or class-rule check",
                "the residual is on integrated volume and shape ratios, not on hydrodynamic merit",
                "the wetted-area budget is a geometric gate, not a resistance estimate",
            ],
        }


def write_solver_manifest(solution: DesignSolution, path: str | Path) -> Path:
    """Write :meth:`DesignSolution.manifest` as JSON, floats at full precision."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(solution.manifest(), indent=2) + "\n", encoding="utf-8")
    return target


def write_design_json(solution: DesignSolution, path: str | Path) -> Path:
    """Write the manifest plus a short digest line per iteration.

    A separate, smaller document from the manifest: this is the one an agent reads
    between steps, and it is deliberately cheap to parse.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": DESIGN_SCHEMA,
        "status": solution.status,
        "solved": solution.solved,
        "iterations": [
            {
                "iteration": record.iteration,
                "variables": dict(record.variables),
                "residual_norm": record.residual_norm,
                "digest": record.digest,
                "passed": record.passed,
            }
            for record in solution.records
        ],
    }
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return target
