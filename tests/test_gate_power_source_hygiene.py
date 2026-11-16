"""The source tree must never carry an injected mutation.

Why this file exists: on 2026-11-16 ``tensorcad/hydrostatics_analytic.py`` on main
contained

    vertical_centre_buoyancy=-0.5 * draft,  # mutant ANA-2

-- a gate-power mutant injected by ``tools/hull_gate_power/mutate.py`` into a run
that was killed before it could restore the file.  The leftover was then committed
("m22596: matured mission deliverables") and six tests in
``tests/test_hydrostatics_analytic.py`` went red on main.  The harness restored the
file only on a clean exit, and nothing ever looked at the tree for leftovers.

So the invariant is checked here, cheaply, by reading the sources: a file that says
"# mutant" is broken on purpose, and shipping it silently swaps a correct closed
form for a wrong one.  The companion behaviours (refuse to start on a dirty tree,
restore on SIGTERM/SIGINT/SIGHUP) are in ``stray_markers`` and the signal handlers
of the harness.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tools" / "hull_gate_power" / "mutate.py"
ANA = ROOT / "tensorcad" / "hydrostatics_analytic.py"


@pytest.fixture(scope="module")
def harness():
    spec = importlib.util.spec_from_file_location("mutate", HARNESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_no_shipped_source_carries_an_injected_mutation(harness):
    offenders = harness.stray_markers()
    assert offenders == [], (
        "a mutation marker is still in the source tree -- a killed run left it "
        "behind and it would ship as a wrong closed form: " + "; ".join(offenders)
    )


def test_the_scan_finds_a_marker_wherever_it_is(tmp_path):
    spec = importlib.util.spec_from_file_location("mutate", HARNESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    dirty = tmp_path / "dirty_module.py"
    dirty.write_text("x = 1  # fine\ny = -0.5  # mutant ANA-2\n")
    clean = tmp_path / "clean_module.py"
    clean.write_text("x = 1  # fine\n")
    assert [line for line in module.stray_markers([dirty, clean])] == [
        "dirty_module.py:2: y = -0.5  # mutant ANA-2"
    ]


def test_the_scan_covers_every_mutable_module(harness):
    covered = {path.name for path in harness.MUTABLE.values()}
    assert covered == {"hull.py", "hydrostatics_analytic.py", "series60.py",
                       "series60_cli.py"}


def test_the_analytic_vertical_centre_of_buoyancy_is_the_closed_form():
    """The value the ANA-2 mutant broke: -0.4*T, not -0.5*T.

    z-weighted mean of zeta**(1/2) on [-T, 0] is (-4T/15)/(2/3) = -0.4 T; the
    mutant shipped -0.5 T.  Pinned analytically here so the constant has to be
    derived, not copied.
    """
    from tensorcad.hull import HullParameters
    from tensorcad.hydrostatics_analytic import analytic_hydrostatics

    draft = 0.6
    result = analytic_hydrostatics(
        HullParameters(length=10.0, beam=1.5, draft=draft))
    vcb = result.vertical_centre_buoyancy
    assert math.isclose(vcb, -0.4 * draft, rel_tol=1e-12), (
        f"analytic VCB is {vcb}, expected -0.4*T = {-0.4 * draft} "
        "(the ANA-2 mutant shipped -0.5*T)"
    )
    # And the closed form must be the one the comment derives: int z*sqrt(1+z/T).
    integral = sum(
        (z * math.sqrt(1.0 + z / draft)) for z in
        [(-draft) * (i + 0.5) / 20000 for i in range(20000)]
    ) * (draft / 20000)
    area = sum(
        math.sqrt(1.0 + z / draft) for z in
        [(-draft) * (i + 0.5) / 20000 for i in range(20000)]
    ) * (draft / 20000)
    assert math.isclose(vcb, integral / area, rel_tol=2e-3)