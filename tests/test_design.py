"""The design loop: does it find the design, and does it refuse to bless a bad one?

The tests are behavioural rather than numerical-where-they-can-be: the strongest check
available here is that the loop **recovers a hull whose measurements it was given**. Take
a reference hull, measure it, hand its displacement and shape ratios to a problem that
starts somewhere else, and the answer must be the reference hull again. That fails if the
Jacobian, the step, the hydrostatics or the mesh are wrong, and it needs no tolerance that
was fitted to the code.

The rest pins the verdicts, because a design loop that only reports success is not a
design loop: converged-and-verified, converged-but-rejected-by-the-validation-layer,
budget exceeded, iteration budget exhausted, and the argument checks.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import struct

import pytest

from tensorcad.design import (
    DESIGN_SCHEMA,
    DesignProblem,
    DesignTargets,
    write_design_json,
    write_solver_manifest,
)
from tensorcad.hull import HullParameters

REFERENCE = HullParameters(length=12.0, beam=2.4, draft=0.9)


def _measure(parameters: HullParameters):
    problem = DesignProblem(base=parameters, targets=DesignTargets(displacement_m3=1.0))
    return problem.evaluate(problem.variable_values(parameters))


@pytest.fixture(scope="module")
def reference_record():
    return _measure(REFERENCE)


def _targets_from(record) -> DesignTargets:
    """The measured displacement and shape ratios of a hull, as design targets."""
    return DesignTargets(
        displacement_m3=record.hydrostatics.displacement_volume,
        beam_to_draft=record.parameters.beam / record.parameters.draft,
        length_to_beam=record.parameters.length / record.parameters.beam,
    )


def test_a_target_built_from_a_hull_is_satisfied_by_that_hull(reference_record) -> None:
    """Self-consistency: evaluating the reference against its own measurements is exact."""
    problem = DesignProblem(
        base=REFERENCE,
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
    )
    record = problem.evaluate(problem.variable_values(REFERENCE))
    assert record.residuals == {"displacement_m3": 0.0, "beam_to_draft": 0.0, "length_to_beam": 0.0}
    assert record.residual_norm == 0.0


def test_the_loop_recovers_a_known_hull_from_its_measurements(reference_record) -> None:
    """The core acceptance test: three variables in, the reference hull out."""
    problem = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
    )
    solution = problem.solve()
    assert solution.status == "converged"
    assert solution.solved
    final = solution.final.parameters
    for name in ("length", "beam", "draft"):
        recovered = getattr(final, name)
        expected = getattr(REFERENCE, name)
        assert recovered == pytest.approx(expected, rel=1e-7), name
    assert solution.final.hydrostatics.displacement_volume == pytest.approx(
        reference_record.hydrostatics.displacement_volume, rel=1e-9
    )


def test_the_convergence_rate_is_quadratic(reference_record) -> None:
    """A Newton step on a smooth residual doubles its correct digits; measure it, don't claim it."""
    problem = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
    )
    solution = problem.solve()
    order = solution.convergence_order
    assert order is not None
    assert 1.8 < order < 2.3, order
    norms = [record.residual_norm for record in solution.records]
    assert all(later < earlier for earlier, later in itertools.pairwise(norms))


def test_the_finite_difference_jacobian_measures_the_geometry_not_the_grid(reference_record) -> None:
    """A shape ratio has sensitivity exactly +1 and -1; displacement is linear in beam and draft."""
    problem = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=reference_record.hydrostatics.displacement_volume,
            beam_to_draft=REFERENCE.beam / REFERENCE.draft,
        ),
        variables=("beam", "draft"),
    )
    solution = problem.solve()
    sensitivity = solution.sensitivity()
    assert sensitivity["beam"]["beam_to_draft"] == pytest.approx(1.0, abs=1e-9)
    assert sensitivity["draft"]["beam_to_draft"] == pytest.approx(-1.0, abs=1e-9)
    assert sensitivity["beam"]["displacement_m3"] == pytest.approx(1.0, abs=1e-6)
    assert sensitivity["draft"]["displacement_m3"] == pytest.approx(1.0, abs=1e-6)


def test_two_solves_of_the_same_problem_agree_bit_for_bit(reference_record) -> None:
    """The generator is deterministic, so a design run is replayable from its targets."""
    targets = _targets_from(reference_record)
    base = HullParameters(length=10.0, beam=2.0, draft=0.8)
    first = DesignProblem(base=base, targets=targets, variables=("length", "beam", "draft")).solve()
    second = DesignProblem(base=base, targets=targets, variables=("length", "beam", "draft")).solve()
    assert first.final.digest == second.final.digest
    assert first.final.variables == second.final.variables
    def without_timing(solution):
        return [
            {key: value for key, value in record.to_dict().items() if key != "seconds"}
            for record in solution.records
        ]

    assert without_timing(first) == without_timing(second)


def test_a_solved_design_can_be_re_used_as_the_start_of_the_next_one(reference_record) -> None:
    """State is carried, not rebuilt: reopening a design and changing a target converges."""
    problem = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=reference_record.hydrostatics.displacement_volume,
            beam_to_draft=REFERENCE.beam / REFERENCE.draft,
        ),
    )
    first = problem.solve()
    assert first.solved
    # Same hull family, 20 % more displacement, starting from the first answer.
    heavier = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=1.2 * reference_record.hydrostatics.displacement_volume,
            beam_to_draft=REFERENCE.beam / REFERENCE.draft,
        ),
    ).solve(start=first.final.variables)
    assert heavier.solved
    assert heavier.final.hydrostatics.displacement_volume == pytest.approx(
        1.2 * reference_record.hydrostatics.displacement_volume, rel=1e-9
    )
    assert heavier.final.parameters.beam > first.final.parameters.beam


def test_the_wetted_area_is_a_budget_and_a_tight_one_is_enforced(reference_record) -> None:
    """A budget is measured every iteration and blocks the verdict, not the convergence."""
    measured = reference_record.hydrostatics.wetted_surface_area
    tight = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=reference_record.hydrostatics.displacement_volume,
            beam_to_draft=REFERENCE.beam / REFERENCE.draft,
            wetted_area_m2_max=0.5 * measured,
        ),
    ).solve()
    assert tight.status == "converged"  # the residual really did reach zero
    assert not tight.solved  # but the design is not deliverable
    assert tight.final.budgets["wetted_area_m2_max"][2] is False
    assert tight.failed_checks == ("wetted_area_m2_max",)
    assert tight.manifest()["verification"]["failed"] == ["wetted_area_m2_max"]

    loose = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=reference_record.hydrostatics.displacement_volume,
            beam_to_draft=REFERENCE.beam / REFERENCE.draft,
            wetted_area_m2_max=2.0 * measured,
        ),
    ).solve()
    assert loose.solved
    assert loose.final.budgets["wetted_area_m2_max"][2] is True


def test_a_converged_design_that_fails_its_own_mesh_checks_is_not_delivered() -> None:
    """The verdict is the conjunction: converging is not the same as being usable.

    A 10 x 10 m hull with a 4.1 m draft hits the target displacement and the target
    beam/draft ratio exactly, and the validation layer then rejects its triangles
    (1.38 deg minimum angle at 41.5 aspect ratio). The loop reports
    ``unverified_geometry`` and refuses to call it solved.
    """
    problem = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=DesignTargets(displacement_m3=200.0, beam_to_draft=2.5),
    )
    solution = problem.solve()
    assert solution.status == "converged"
    assert not solution.solved
    assert solution.failed_checks  # named, not implied
    assert solution.final.residual_norm < 1e-9  # it *did* hit the target
    assert solution.final.handoff["watertight"] is True
    failed = {check.name for check in solution.final.checks if not check.passed}
    assert failed, "the record must name what failed"
    assert any("sliver" in name or "angle" in name for name in failed), failed
    assert solution.manifest()["solver"]["status"] == "not_ready"


def test_an_exhausted_iteration_budget_is_reported_as_such(reference_record) -> None:
    problem = DesignProblem(
        base=HullParameters(length=9.0, beam=1.6, draft=0.6),
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
        max_iterations=2,
    )
    solution = problem.solve()
    assert solution.status == "max_iterations"
    assert not solution.solved
    assert len(solution.records) == 3
    assert solution.final.residual_norm > 1e-6


def test_solver_feedback_moves_the_targets_and_is_recorded(reference_record) -> None:
    """The seam for a coupled run: feedback adjusts the target, the manifest shows both."""
    original = reference_record.hydrostatics.displacement_volume
    starts = {"calls": 0}

    def feedback(record):
        starts["calls"] += 1
        if starts["calls"] > 1:
            return None
        return {"displacement_m3": 0.05 * original}

    problem = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(displacement_m3=original, beam_to_draft=REFERENCE.beam / REFERENCE.draft),
    )
    solution = problem.solve(feedback=feedback)
    assert solution.solved
    assert starts["calls"] >= 1
    assert solution.targets.displacement_m3 == pytest.approx(1.05 * original)
    assert solution.final.hydrostatics.displacement_volume == pytest.approx(1.05 * original, rel=1e-9)
    assert solution.manifest()["targets"]["displacement_m3"] == pytest.approx(1.05 * original)


def test_feedback_cannot_adjust_something_that_is_not_a_target(reference_record) -> None:
    problem = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=DesignTargets(
            displacement_m3=reference_record.hydrostatics.displacement_volume,
            beam_to_draft=3.0,
        ),
    )
    with pytest.raises(ValueError, match="cannot adjust"):
        problem.solve(feedback=lambda record: {"drag_n": 1.0})


def test_the_handoff_record_uses_the_shared_digest_convention(reference_record) -> None:
    """A design record must be comparable with any other handoff, so it hashes the same way."""
    mesh = reference_record.mesh
    digest = hashlib.sha256()
    digest.update(struct.pack("<II", len(mesh.vertices), len(mesh.faces)))
    for vertex in mesh.vertices:
        digest.update(struct.pack("<3d", *vertex))
    for face in mesh.faces:
        digest.update(struct.pack("<3I", *face))
    assert reference_record.digest == digest.hexdigest()
    assert reference_record.handoff["schema"] == "tensorcad.geometry.handoff/1"
    assert reference_record.handoff["mesh_sha256"] == digest.hexdigest()
    assert reference_record.handoff["watertight"] is True
    assert reference_record.handoff["euler_characteristic"] == 2


def test_the_manifest_is_json_ready_and_states_the_solver_was_not_run(reference_record) -> None:
    solution = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
    ).solve()
    manifest = solution.manifest()
    assert manifest["schema"] == DESIGN_SCHEMA
    assert manifest["solved"] is True
    assert manifest["solver"]["status"] == "ready"
    assert manifest["solver"]["executed_here"] is False
    assert manifest["solver"]["reason"]
    assert manifest["reference_quantities"]["water_density_kg_per_m3"] == 1025.0
    assert len(manifest["path"]) == len(solution.records)
    assert json.loads(json.dumps(manifest)) == manifest


def test_both_writers_produce_the_documents_they_promise(reference_record, tmp_path) -> None:
    solution = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=_targets_from(reference_record),
        variables=("length", "beam", "draft"),
    ).solve()
    manifest_path = write_solver_manifest(solution, tmp_path / "nested" / "manifest.json")
    design_path = write_design_json(solution, tmp_path / "design.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    digest = json.loads(design_path.read_text(encoding="utf-8"))
    assert manifest["schema"] == DESIGN_SCHEMA
    assert len(digest["iterations"]) == len(solution.records)
    assert digest["iterations"][-1]["digest"] == solution.final.digest
    assert digest["iterations"][-1]["passed"] is True


def test_a_problem_without_targets_is_refused() -> None:
    with pytest.raises(ValueError, match="at least one equality target"):
        DesignProblem(base=HullParameters(), targets=DesignTargets())


def test_an_underdetermined_problem_is_refused() -> None:
    """Two equalities cannot be met by one knob, and saying so early beats stalling."""
    with pytest.raises(ValueError, match="cannot be met by"):
        DesignProblem(
            base=HullParameters(),
            targets=DesignTargets(displacement_m3=7.0, beam_to_draft=2.5),
            variables=("beam",),
        )


def test_the_argument_checks_name_the_problem() -> None:
    base = HullParameters()
    targets = DesignTargets(displacement_m3=7.0)
    with pytest.raises(ValueError, match="unknown design variable"):
        DesignProblem(base=base, targets=targets, variables=("skeg",))
    with pytest.raises(ValueError, match="duplicate"):
        DesignProblem(base=base, targets=targets, variables=("beam", "beam"))
    with pytest.raises(ValueError, match="at least one free variable"):
        DesignProblem(base=base, targets=targets, variables=())
    with pytest.raises(ValueError, match="max_iterations"):
        DesignProblem(base=base, targets=targets, max_iterations=0)
    with pytest.raises(ValueError, match="tolerance"):
        DesignProblem(base=base, targets=targets, tolerance=2.0)
    with pytest.raises(ValueError, match="damping"):
        DesignProblem(base=base, targets=targets, damping=0.0)
    with pytest.raises(ValueError, match="finite_difference_step"):
        DesignProblem(base=base, targets=targets, finite_difference_step=0.0)
    with pytest.raises(ValueError, match="displacement_m3"):
        DesignTargets(displacement_m3=-1.0)


def test_evaluating_an_unknown_variable_is_refused() -> None:
    problem = DesignProblem(base=HullParameters(), targets=DesignTargets(displacement_m3=7.0))
    with pytest.raises(ValueError, match="unexpected variables"):
        problem.evaluate({"beam": 2.0, "draft": 0.8, "length": 10.0})


def test_a_start_outside_the_free_variables_is_refused() -> None:
    problem = DesignProblem(base=HullParameters(), targets=DesignTargets(displacement_m3=7.0))
    with pytest.raises(ValueError, match="not a free variable"):
        problem.solve(start={"beam": 2.0, "draft": 0.8, "length": 10.0})
