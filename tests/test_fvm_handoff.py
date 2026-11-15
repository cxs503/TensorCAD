"""Tests for the watertight-solid -> FVM consumption route.

The claims under test, in the order the module makes them:

* a CAD solid written as gmsh ``C3D4`` survives the contract verifier and reaches
  TensorFVM's own mesh object **from the file**;
* the finite-volume operator reproduces a linear (harmonic) field to round-off on
  a box *and* on a curved hull solid -- that is a patch test, not a residual;
* a non-trivial closed-form harmonic field is reproduced within a stated
  tolerance, and the error falls under refinement;
* the z-up frame is worth something: the same table read as y-up is wrong by
  orders of magnitude more;
* refusals are refusals: a moved node and a non-``z-up`` manifest raise, and an
  absent solver library is a ``status``, never a pass.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

pytest.importorskip("numpy")
pytest.importorskip("tensorfvm.native_mesh3d", reason="TensorFVM checkout is required for the FVM leg")

from tensorcad import fvm_handoff as fh  # noqa: E402
from tensorcad.solver_interface import HandoffError  # noqa: E402

BARGE = dict(length=2.0, beam=0.4, draft=0.3, nx=9, nz=5, layers=2)
HULL = dict(length=10.0, beam=1.0, draft=0.6, nx=13, nz=7, layers=2)


@pytest.fixture(scope="module")
def barge_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("barge") / "bundle"
    fh.barge_bundle(out, **BARGE)
    return out


@pytest.fixture(scope="module")
def hull_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("hull") / "bundle"
    fh.hull_volume_bundle(out, **HULL)
    return out


def test_patch_test_is_machine_precision_on_a_box(barge_bundle: Path) -> None:
    """A consistent FVM operator must return a linear field exactly."""
    result = fh.fvm_patch_from_bundle(barge_bundle)
    assert result.status == "ran", result.error
    assert result.metrics["relative_linf"] < 1e-9, result.metrics
    assert result.metrics["relative_l2"] < 1e-12
    assert not result.gates_failed, result.gates_failed


def test_patch_test_on_a_curved_hull_solid(hull_bundle: Path) -> None:
    """The same exactness claim on a ship-shaped body, not on a box."""
    result = fh.fvm_patch_from_bundle(hull_bundle)
    assert result.status == "ran", result.error
    metrics = result.metrics
    assert metrics["cells"] > 1000
    assert metrics["relative_linf"] < 1e-8, metrics
    # skew cells limit the round-off level, so the residual ledger is the sharper witness
    assert metrics["linear_relative_residual"] < 1e-12
    assert metrics["global_relative_balance"] < 1e-10


def test_the_wall_surface_is_closed_and_genus_zero(barge_bundle: Path) -> None:
    report = fh._watertight(barge_bundle)
    assert report["bad_edges"] == 0
    assert report["duplicate_directed_edges"] == 0
    assert report["euler_characteristic"] == 2
    assert report["outward"] is True


def test_harmonic_reference_is_met_and_is_a_real_discretisation_error(barge_bundle: Path) -> None:
    result = fh.fvm_harmonic_from_bundle(barge_bundle, length=BARGE["length"], beam=BARGE["beam"], draft=BARGE["draft"])
    assert result.status == "ran", result.error
    metrics = result.metrics
    # a non-trivial field: the error must be visibly above round-off and below the band
    assert 1e-3 < metrics["relative_l2"] < 0.15, metrics
    assert not result.gates_failed, result.gates_failed
    # the same solve must also satisfy the solver's own conservation ledger
    assert metrics["linear_relative_residual"] < 1e-12
    assert metrics["global_relative_balance"] < 1e-10


def test_refinement_reduces_the_harmonic_error(tmp_path: Path) -> None:
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    fh.barge_bundle(coarse, **{**BARGE, "nx": 9, "nz": 5, "layers": 2})
    fh.barge_bundle(fine, **{**BARGE, "nx": 17, "nz": 9, "layers": 3})
    a = fh.fvm_harmonic_from_bundle(coarse, length=BARGE["length"], beam=BARGE["beam"], draft=BARGE["draft"])
    b = fh.fvm_harmonic_from_bundle(fine, length=BARGE["length"], beam=BARGE["beam"], draft=BARGE["draft"])
    assert a.status == b.status == "ran"
    assert b.metrics["cells"] > a.metrics["cells"]
    assert b.metrics["relative_l2"] < a.metrics["relative_l2"], (a.metrics["relative_l2"], b.metrics["relative_l2"])


def test_the_frame_contract_has_a_measurable_price(barge_bundle: Path) -> None:
    """Reading the z-up table as y-up is the mistake the 5415 rounds paid for."""
    result = fh.frame_counterexample(barge_bundle)
    assert result.status == "ran", result.error
    metrics = result.metrics
    assert metrics["frame"] == "z-up"
    assert metrics["zup_relative_linf"] < 1e-9
    assert metrics["yup_relative_linf"] > 1e-2
    assert metrics["frame_price_ratio"] > 1e6


def test_a_moved_node_is_refused_before_any_solve(barge_bundle: Path, tmp_path: Path) -> None:
    tampered = tmp_path / "tampered"
    tampered.mkdir()
    for name in ("handoff.json", "volume.msh"):
        (tampered / name).write_text((barge_bundle / name).read_text())
    lines = (tampered / "volume.msh").read_text().splitlines()
    for index, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit() and "." in parts[1]:
            parts[1] = f"{float(parts[1]) + 1e-6:.17g}"
            lines[index] = " ".join(parts)
            break
    (tampered / "volume.msh").write_text("\n".join(lines) + "\n")
    with pytest.raises(HandoffError):
        fh._open_bundle(tampered)


def test_a_y_up_manifest_is_refused(barge_bundle: Path, tmp_path: Path) -> None:
    tampered = tmp_path / "yup"
    tampered.mkdir()
    for name in ("handoff.json", "volume.msh"):
        (tampered / name).write_text((barge_bundle / name).read_text())
    manifest = json.loads((tampered / "handoff.json").read_text())
    manifest["frame"] = {**manifest["frame"], "name": "y-up", "up_axis": "y"}
    (tampered / "handoff.json").write_text(json.dumps(manifest))
    with pytest.raises(HandoffError):
        fh._open_bundle(tampered)


def test_an_absent_solver_is_a_status_not_a_pass(barge_bundle: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    real_load = fh._load

    def fake_load(name: str):
        if name.startswith("tensorfvm"):
            raise fh.SolverUnavailableError(f"{name}: simulated absence")
        return real_load(name)

    monkeypatch.setattr(fh, "_load", fake_load)
    patch = fh.fvm_patch_from_bundle(barge_bundle)
    harmonic = fh.fvm_harmonic_from_bundle(barge_bundle, length=2.0, beam=0.4, draft=0.3)
    frame = fh.frame_counterexample(barge_bundle)
    for result in (patch, harmonic, frame):
        assert result.status == "unavailable"
        assert "simulated absence" in result.error
        assert result.gates == ()


def test_run_all_cases_writes_a_results_table(tmp_path: Path) -> None:
    payload = fh.run_all_cases(
        tmp_path / "artifacts",
        barge_resolutions=((9, 5, 2),),
        hull_resolutions=((13, 7, 2),),
    )
    out = tmp_path / "artifacts"
    written = json.loads((out / "results.json").read_text())
    assert written["schema"] == "tensorcad.fvm_handoff/1"
    assert payload["all_gates_passed"] is True, [c["gates_failed"] for c in payload["cases"]]
    assert written["statuses"] == ["ran"]
    assert {c["route"] for c in written["cases"]} == {"volume_msh"}
    summary = (out / "summary.md").read_text()
    assert "PASS" in summary and "FAIL" not in summary


def test_harmonic_field_really_satisfies_laplace() -> None:
    """The reference must be harmonic, else the comparison means nothing."""
    import numpy as np

    length, beam, draft = 2.0, 0.4, 0.3
    h = 1e-5
    point = np.array([[0.37 * length, 0.41 * beam, -0.5 * draft]])
    laplacian = 0.0
    for axis in range(3):
        for sign in (+1, -1):
            shifted = point.copy()
            shifted[0, axis] += sign * h
            laplacian += float(fh.exact_harmonic(shifted, length=length, beam=beam, draft=draft)[0])
    laplacian -= 6.0 * float(fh.exact_harmonic(point, length=length, beam=beam, draft=draft)[0])
    laplacian /= h * h
    scale = abs(float(fh.exact_harmonic(point, length=length, beam=beam, draft=draft)[0]))
    assert abs(laplacian) / (scale * math.pi**2) < 1e-3


def test_a_collapsed_keel_column_is_refused_rather_than_silently_dropped(tmp_path: Path) -> None:
    """The documented limit: a zero-breadth keel line yields degenerate faces.

    Recorded as a test so the finding cannot quietly disappear: the CAD mesher
    accepts a collapsed column, the tets it splits into are degenerate, and
    TensorFVM's auditor refuses the mesh.  The handoff must therefore emit a
    solid with a flat of the bottom (``hull_section_grid`` does).
    """
    from tensorcad.volume_mesh import SectionGrid, build_volume_mesh
    from tensorcad.solver_interface import write_bundle_volume

    length, beam, draft = 10.0, 1.0, 0.6
    nx, nz = 9, 5
    xs = tuple(length * i / (nx - 1) for i in range(nx))
    zs = tuple(-draft + draft * j / (nz - 1) for j in range(nz))
    collapsed = SectionGrid(
        xs=xs,
        zs=zs,
        breadth=lambda i, j: 0.5 * beam * (1.0 - (2.0 * xs[i] / length - 1.0) ** 2) * (1.0 - (zs[j] / draft) ** 2),
    )
    mesh = build_volume_mesh(collapsed, first_layer=0.006, growth=1.2, layers=2, split_tets=True)
    out = tmp_path / "collapsed"
    write_bundle_volume(out, vertices=mesh.vertices, tets=mesh.tets, case="collapsed_keel")
    with pytest.raises(Exception) as excinfo:
        fh._open_bundle(out)
    assert "degenerate" in str(excinfo.value).lower() or "nonpositive" in str(excinfo.value).lower()