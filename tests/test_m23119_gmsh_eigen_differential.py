"""Gmsh differential with Gmsh's Jacobian measures switched on (`HAVE_EIGEN`).

Background (see ``docs/m23119_gmsh_eigen_differential.md``): the Gmsh arm of the anchor
``GAP-CAD-D1-MESH-QUALITY-LIB`` used to answer ``nan``/``0``/``DBL_MAX`` for every
``minSICN``/``minSIGE`` because Gmsh's ``CondNumBasis`` inverts a matrix through
``fullMatrix.h``, which is ``#if defined(HAVE_EIGEN)``. Those numbers are not an oracle
gap that can be explained away -- they are a *silent* one: the binary exits 0 and the
differential reads it as data.

These tests pin the three properties that keep that from coming back, without needing a
Gmsh source tree or Eigen at test time:

* the build switches on both config switches it depends on and states which way each
  went (``HAVE_MESH`` for the shape measures, ``HAVE_EIGEN`` for the Jacobian ones);
* the driver keeps "must agree" and "expected to differ, here is why" in two disjoint
  tables -- a pair may not sit in both, and every explained pair carries a reason;
* the recorded run is a real one: the committed artifact must still show the four element
  kinds compared, zero disagreements, and the two *measured* deviations that falsify
  claims made earlier in this repository (``gmsh_icn`` is not Gmsh's unit; hex ``volume``
  is a triangulated volume, not ``∫det J``).
"""
from __future__ import annotations

import json
import importlib.util
import subprocess
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ORACLES = REPO / "tools" / "mesh_quality_oracles"
BUILD = ORACLES / "build_gmsh_source_oracle.sh"
DRIVER = REPO / "tools" / "mesh_quality_gmsh_source_diff.py"
ARTIFACT = REPO / "artifacts" / "mesh_quality_external_diff" / "gmsh_source_diff.json"


def _driver():
    spec = importlib.util.spec_from_file_location("gmsh_source_diff_m23119", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def test_build_switches_on_have_eigen_as_well_as_have_mesh():
    """`HAVE_EIGEN` decides whether minSICN/minSIGE are numbers or nan/0/DBL_MAX."""
    script = BUILD.read_text()
    assert 'DEFINE_ON="HAVE_MESH"' in script, "HAVE_MESH is the base switch"
    assert 'DEFINE_ON="$DEFINE_ON,HAVE_EIGEN"' in script, (
        "Eigen must be added to the config defines when it is found; without it Gmsh "
        "answers nan for every ICN and the oracle still exits 0"
    )
    for candidate in ("TP_EIGEN_INCLUDE", "contrib/eigen", "/usr/include/eigen3"):
        assert candidate in script, candidate
    assert "HAVE_EIGEN=" in script, "the build must report which way the switch went"


def test_build_compiles_in_parallel_so_the_eigen_build_fits_a_command_window():
    script = BUILD.read_text()
    assert "xargs -P" in script and "TP_JOBS" in script, (
        "a serial build is ~70 s without Eigen and >2 min with it; the object loop must "
        "be parallel or the build gets abandoned half-way"
    )
    assert "qualityMeasuresJacobian.cpp" in script, (
        "the Jacobian measures (minSICN/minSIGE) live in that translation unit"
    )


def test_build_rejects_a_tree_without_the_jacobian_measures(tmp_path):
    """A tree missing qualityMeasuresJacobian.cpp must be refused, not silently thinned."""
    sh = shutil.which("sh")
    if sh is None:
        pytest.skip("no /bin/sh on this host")
    fake = tmp_path / "tree"
    (fake / "src" / "geo").mkdir(parents=True)
    (fake / "src" / "mesh").mkdir(parents=True)
    (fake / "src" / "mesh" / "qualityMeasures.cpp").write_text("// stub\n")
    completed = subprocess.run([sh, str(BUILD), str(fake), str(tmp_path / "out")],
                               capture_output=True, text=True, timeout=120)
    assert completed.returncode == 2, (completed.returncode, completed.stdout, completed.stderr)
    assert "qualityMeasuresJacobian.cpp missing" in completed.stderr


def test_agree_and_explained_tables_are_disjoint_and_reasons_are_given():
    drv = _driver()
    for kind, mapping in drv.MAPPING.items():
        explained = drv.EXPLAINED.get(kind, {})
        assert not set(mapping) & set(explained), (
            f"{kind}: {set(mapping) & set(explained)} appear both as 'must agree' and as "
            "'explained difference'; one pair, one status"
        )
    for kind, explained in drv.EXPLAINED.items():
        for metric, (gmsh_name, why) in explained.items():
            assert gmsh_name and isinstance(gmsh_name, str), (kind, metric)
            assert len(why.split()) > 12, (
                f"{kind}.{metric} has no explanatory sentence; an unexplained difference "
                "is a finding, not an explanation"
            )


def test_every_mapped_metric_exists_in_the_library_under_that_kind():
    """The tables must name real metrics -- a typo currently reads as 'no such pair'."""
    import mesh_quality_external_diff as ext  # noqa: PLC0415  (tools/ path added below)

    drv = _driver()
    corpus = ext.meshes()
    by_kind = {}
    for name, (kind, points, elems) in corpus.items():
        by_kind.setdefault(kind, (points, elems))
    for kind, mapping in list(drv.MAPPING.items()) + list(drv.EXPLAINED.items()):
        points, elems = by_kind[kind]
        ours = drv.ours_for(kind, points, elems)
        for metric in mapping:
            assert metric in ours, f"{kind}: no such library metric {metric!r}"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="recorded differential absent")
def test_recorded_run_compares_four_kinds_without_a_disagreement():
    report = json.loads(ARTIFACT.read_text())
    assert report["oracle_ran"] is True
    assert report["pairs_compared"] >= 12, report["pairs_compared"]
    assert report["elements_beyond_tolerance"] == 0, (
        "a recorded disagreement in a 'must agree' pair is a failed run, not a baseline"
    )
    kinds = set()
    for entry in report["meshes"].values():
        if entry.get("status") != "RAN":
            continue
        if entry.get("metrics") or entry.get("explained"):
            kinds.add(entry["kind"])
    assert kinds == {"tri", "quad", "tet", "hex"}, kinds


@pytest.mark.skipif(not ARTIFACT.exists(), reason="recorded differential absent")
def test_recorded_run_keeps_the_two_falsified_claims_visible():
    """Deviations that contradict earlier documentation must stay in the record.

    * tet ``gmsh_icn`` vs Gmsh ``minSICN``: ``metrics.py`` claimed the two were the same
      quantity in the same unit. Measured, Gmsh's unit is the *regular* tetrahedron
      (Gmsh regular = 1.0, right-corner = 0.8165) while ours is the right-corner one, and
      the ratio is not even constant (>= 35 % on the Kuhn corpus, >= 92 % on slivers).
    * hex ``volume``: the library triangulates the faces (divergence theorem), Gmsh
      integrates ``det J`` of the trilinear map. Agreement on planar faces only.
    """
    report = json.loads(ARTIFACT.read_text())
    tet = report["meshes"]["tet_cube6"]["explained"]
    assert tet["gmsh_icn"]["max_rel_diff"] > 1e-3, (
        "if gmsh_icn ever matches Gmsh's minSICN, the docstring claim becomes true and "
        "this test (and the explanation) must be updated deliberately"
    )
    hexd = report["meshes"]["hex_distorted"]["explained"]
    assert hexd["volume"]["max_rel_diff"] > 1e-2, hexd["volume"]
    hexu = report["meshes"]["hex_unit"]["explained"]
    assert hexu["volume"]["max_rel_diff"] < 1e-12, (
        "on planar faces the two volume definitions coincide; a change here means the "
        "library's hex volume moved"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))