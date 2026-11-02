"""Tests for resolution control driven by the analytic displacement.

The acceptance quantity is :func:`tensorcad.hull.analytic_volume` — a closed-form
integral of the analytic surface, which the sampled mesh never sees. Refining a
mesh against its own refinement would show the numbers moving without showing
that they are right; these tests pin the result to the external integral instead.

The expected facet/vertex counts are derived by hand from the documented
topology (one closed ring of ``2*nz-1`` points per station, ``nz`` distinct points
on the two pointed end sections, two degenerate end-cap triangles dropped) and
are therefore independent of the generator's bookkeeping.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import (
    HullParameters,
    analytic_volume,
    generate_hull,
    resolve_resolution,
)


def _triangles(nx: int, nz: int) -> int:
    ring = 2 * nz - 1
    return 2 * ring * (nx - 1) - 2


def _vertices(nx: int, nz: int) -> int:
    return (nx - 2) * (2 * nz - 1) + 2 * nz


def test_plan_refines_until_the_analytic_error_is_within_target() -> None:
    start = HullParameters(nx=11, nz=5)
    plan = resolve_resolution(start, target_volume_error=0.01)
    assert plan.accepted
    assert abs(plan.relative_error) <= 0.01
    assert plan.parameters.nx > start.nx and plan.parameters.nz > start.nz
    assert plan.refinements >= 1
    # Bookkeeping is checked against the hand-derived topology, not against the
    # plan's own fields alone.
    assert plan.triangles == _triangles(plan.parameters.nx, plan.parameters.nz)
    assert plan.vertices == _vertices(plan.parameters.nx, plan.parameters.nz)


def test_plan_agrees_with_a_fresh_generation_and_the_analytic_integral() -> None:
    plan = resolve_resolution(HullParameters(nx=11, nz=5), target_volume_error=0.01)
    rebuilt = generate_hull(plan.parameters)
    exact = analytic_volume(
        plan.parameters.length, plan.parameters.beam, plan.parameters.draft
    )
    assert rebuilt.hydrostatics is not None
    measured = rebuilt.hydrostatics.displacement_volume
    assert plan.relative_error == pytest.approx(measured / exact - 1.0, abs=1e-12)
    assert abs(measured / exact - 1.0) <= plan.target_volume_error


def test_plan_leaves_an_already_adequate_grid_alone() -> None:
    """Refinement must not be a habit: an accurate grid is returned unchanged."""
    start = HullParameters(nx=81, nz=21)
    plan = resolve_resolution(start, target_volume_error=0.005)
    assert plan.accepted
    assert plan.refinements == 0
    assert plan.parameters == start
    assert plan.triangles == 6558  # hand-derived from the topology, matches earlier audits


def test_plan_reports_an_honest_failure_when_the_budget_is_exhausted() -> None:
    plan = resolve_resolution(
        HullParameters(nx=11, nz=5), target_volume_error=1e-4, max_triangles=1000
    )
    assert not plan.accepted
    assert abs(plan.relative_error) > plan.target_volume_error
    # It stopped at the last grid the budget allowed, not at a grid it never built.
    assert plan.triangles == _triangles(plan.parameters.nx, plan.parameters.nz)
    assert plan.triangles <= 1000
    assert plan.relative_error < -0.01


def test_plan_rejects_arguments_that_could_not_terminate() -> None:
    with pytest.raises(ValueError, match="target_volume_error"):
        resolve_resolution(HullParameters(), target_volume_error=0.0)
    with pytest.raises(ValueError, match="target_volume_error"):
        resolve_resolution(HullParameters(), target_volume_error=math.nan)
    with pytest.raises(ValueError, match="growth"):
        resolve_resolution(HullParameters(), growth=1.0)
    with pytest.raises(ValueError, match="growth"):
        resolve_resolution(HullParameters(), growth=math.inf)
    with pytest.raises(ValueError, match="max_triangles"):
        resolve_resolution(HullParameters(), max_triangles=0)


def test_plan_serialises_its_parameters_and_evidence() -> None:
    plan = resolve_resolution(HullParameters(nx=11, nz=5), target_volume_error=0.01)
    payload = plan.to_dict()
    assert payload["schema"] == "tensorcad.hull.resolution/1"
    assert payload["accepted"] is True
    assert payload["parameters"]["nx"] == plan.parameters.nx
    assert payload["target_volume_error"] == 0.01
    assert payload["triangles"] == plan.triangles
