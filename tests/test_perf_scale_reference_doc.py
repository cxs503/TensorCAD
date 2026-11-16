"""Guards for ``docs/perf_scale_reference.md`` (GAP-CAD-D4-PERF-SCALE).

The refererence document is the anchor's synthesis of every committed perf
ladder.  Its tables carried ``perf-evidence`` markers for a year while nothing
rendered them: no renderer was registered under the ``ref-*`` names, and no
script or test read the file, so the markers promised programmatic injection
while the tables were typed by hand.  Three of the four had drifted from their
artifacts -- the single-body table stopped at 3.3e6 triangles while the artifact
reaches 1.3e7, the mesh table had silently dropped its ``not_run`` rows, and the
slope table stated ``load = 1.17`` where the ladder measures 1.08.

These tests are the guard that was missing.  They assert *invariants*, not that
the renderer runs: every block regenerates byte-identically from the artifacts,
every committed ladder is cited, the acceptance sub-tasks are actually covered in
the 10^6 band, and closed debts are not re-stated as open.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import perf_evidence as pe

ROOT = Path(__file__).resolve().parents[1]
PERF = ROOT / "benchmarks" / "perf_scale"
DOC = ROOT / "docs" / "perf_scale_reference.md"
SCRIPT = ROOT / "scripts" / "perf_scale_reference_doc.py"

#: Artifact paths the document is rendered from, keyed exactly as the script.
ARTIFACT_PATHS = {
    "mesh": PERF / "artifacts" / "ladder_mesh.json",
    "voxel": PERF / "artifacts" / "ladder_voxel.json",
    "fleet3pt": PERF / "ladder_fleet_10m_3pt.json",
    "single10m": PERF / "artifacts" / "ladder_single_10m.json",
    "import_single10m": PERF / "artifacts" / "ladder_import_single_10m.json",
    "import_fleet": PERF / "artifacts" / "ladder_import_fleet_10m.json",
}


def _artifacts() -> dict[str, dict]:
    return {key: pe.load_artifact(path) for key, path in ARTIFACT_PATHS.items()}


def _doc() -> str:
    return DOC.read_text(encoding="utf-8")


def _table(name: str) -> list[list[str]]:
    """Cells of a rendered table, header excluded."""
    body = pe.extract_block(_doc(), name)
    assert body is not None, f"{name} markers missing from the document"
    rows = [line for line in body.splitlines() if line.startswith("|")]
    assert len(rows) >= 2, f"{name} has no data rows"
    return [[c.strip() for c in row.strip("|").split("|")] for row in rows[2:]]


# ---------------------------------------------------------------------------
# Drift: the tables are what the artifacts render to
# ---------------------------------------------------------------------------

def test_every_reference_table_regenerates_from_artifacts() -> None:
    """The committed document *is* the render of the committed artifacts."""
    problems = pe.check_document(_doc(), _artifacts(), names=pe.REFERENCE_TABLE_NAMES)
    assert problems == [], problems


def test_all_reference_markers_are_registered_renderers() -> None:
    """A marker with no registered renderer is a hand-written table in disguise."""
    text = _doc()
    for name in pe.REFERENCE_TABLE_NAMES:
        assert pe.extract_block(text, name) is not None, f"{name}: markers missing"
        assert name in pe.ALL_TABLE_SPECS, f"{name}: not registered in ALL_TABLE_SPECS"


def test_registration_refresh_covers_late_defined_tables() -> None:
    """Tables added after the first ``_register`` call must still resolve.

    ``ALL_TABLE_SPECS`` is filled by ``update`` at call time, so a spec added to
    ``REFERENCE_SPECS`` afterwards is invisible unless the registry is refreshed
    -- the failure mode that made ``--write`` raise ``KeyError`` on the first run.
    """
    assert "ref-exponents-table" in pe.ALL_TABLE_SPECS
    assert "ref-error-table" in pe.ALL_TABLE_SPECS
    assert set(pe.REFERENCE_TABLE_NAMES) <= set(pe.ALL_TABLE_SPECS)


def test_positional_artifact_order_matches_specs() -> None:
    """``render_block`` splats artifacts positionally; order is load-bearing."""
    for name in ("ref-coverage-table", "ref-exponents-table", "ref-error-table"):
        assert pe.REFERENCE_SPECS[name][0] == pe.REFERENCE_ARTIFACT_ORDER, name


def test_script_check_passes_on_committed_doc() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "ok" in result.stdout


def test_script_check_detects_a_hand_edited_number(tmp_path: Path) -> None:
    """A single edited digit must be reported, not tolerated."""
    text = _doc()
    body = pe.extract_block(text, "ref-coverage-table")
    assert body is not None
    first_data_row = [ln for ln in body.splitlines() if ln.startswith("| `")][0]
    edited = first_data_row.replace("8,199,998", "9,199,998")
    assert edited != first_data_row
    doc = tmp_path / "reference.md"
    doc.write_text(text.replace(first_data_row, edited), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check", "--doc", str(doc)],
        cwd=str(ROOT), capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 1
    assert "DRIFT" in result.stdout


# ---------------------------------------------------------------------------
# Coverage: the acceptance sub-tasks really are measured in the 10^6 band
# ---------------------------------------------------------------------------

def _coverage() -> dict[str, list[str]]:
    return {row[0].strip("`"): row for row in _table("ref-coverage-table")}


@pytest.mark.parametrize(
    "subtask",
    # The four sub-tasks named by the ledger's acceptance text, plus the
    # accelerated point-location path that is the only one reaching 10^6.
    ["load", "point_locate_accel", "volume", "voxelize"],
)
def test_acceptance_subtask_reaches_the_1e6_band(subtask: str) -> None:
    row = _coverage()[subtask]
    assert row[1] == "三角面" or subtask == "voxelize", row
    assert row[4] == "✅", f"{subtask} does not reach 10^6: {row}"
    assert int(row[2].replace(",", "")) >= 10**6


@pytest.mark.parametrize("subtask", ["obj_read", "volume", "point_locate_accel", "subdivide"])
def test_paths_that_cross_1e7_are_recorded(subtask: str) -> None:
    row = _coverage()[subtask]
    assert row[5] == "✅", f"{subtask} does not reach 10^7: {row}"
    assert int(row[2].replace(",", "")) >= 10**7


def test_generator_load_ceiling_is_stated_not_hidden() -> None:
    """``load`` (generator) is structurally short of 10^7: the table says so.

    If a future ladder pushes the generator past 10^7 this test starts failing,
    which is the intended signal: the document's prose would then be wrong.
    """
    row = _coverage()["load"]
    assert row[5] == "—", row
    assert 8.19e6 <= int(row[2].replace(",", "")) < 10**7


def test_brute_force_point_location_is_not_claimed_in_the_1e6_band() -> None:
    row = _coverage()["point_locate"]
    assert row[4] == "—", row


def test_every_committed_ladder_is_cited_by_the_document() -> None:
    """A new measurement cannot land without the synthesis document noticing."""
    text = _doc()
    measured = [
        p for p in sorted(PERF.rglob("*.json"))
        if p.name not in {"baseline.json"} and ".provenance" not in p.name
        and "mechanism_smoke" not in p.name
    ]
    assert measured, "no perf ladders found"
    missing = []
    for path in measured:
        rel = path.relative_to(PERF)
        if path.stem not in text and str(rel) not in text:
            missing.append(str(rel))
    assert missing == [], f"ladders not cited in {DOC.name}: {missing}"


# ---------------------------------------------------------------------------
# Accuracy: every ladder carries a verdict against the 3% standard
# ---------------------------------------------------------------------------

def test_every_volume_ladder_passes_the_3_percent_standard() -> None:
    rows = _table("ref-error-table")
    judged = 0
    for row in rows:
        ladder, _count, _worst, verdict, _basis = row[:5]
        if verdict == "not_run":
            continue
        judged += 1
        assert verdict == "✅", f"{ladder} is RED: {row}"
        assert float(row[2]) < pe.REFERENCE_ERROR_LIMIT, row
    assert judged >= 5, f"only {judged} ladders judged"


def test_error_table_reads_percent_notes_as_percent() -> None:
    """``rel_err=-1.842%`` is 1.8%, not 184% (the unit bug this table hit)."""
    assert pe._volume_rel_err({"note": "rel_err=-1.842%"}) == pytest.approx(0.01842)
    assert pe._volume_rel_err({"note": "rel_err=0.003%"}) == pytest.approx(3e-5)
    # A structured field is already a fraction and wins over the note.
    assert pe._volume_rel_err({"volume_rel_err": -0.0003, "note": "rel_err=-9%"}) == pytest.approx(0.0003)
    assert pe._volume_rel_err({"note": "no error recorded"}) is None


# ---------------------------------------------------------------------------
# Honesty: closed debts are not re-stated as open
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "stale",
    [
        "10⁷ 单体档未到件",   # closed: ladder_single_10m.json reaches 1.32e7
        "文件导入路径未测",   # closed: ladder_import_{single,fleet}_10m.json
        "单体细分到 3.3M tri",
    ],
)
def test_closed_debts_are_not_listed_as_open(stale: str) -> None:
    assert stale not in _doc(), f"stale debt claim re-appears: {stale}"


def test_debt_section_still_names_the_real_gaps() -> None:
    """Removing stale items must not empty the section: honesty both ways."""
    text = _doc()
    section = text.split("## 7. 局限与欠账", 1)[1]
    assert "CI 基线" in section
    assert "float64" in section or "内存" in section
    assert "抖动" in section


def test_exponents_table_refuses_a_slope_from_two_points() -> None:
    """Two points fit any line perfectly; the table must print ``not_run``."""
    rows = _table("ref-exponents-table")
    two_point = [r for r in rows if r[2] == "2"]
    assert two_point, "expected at least one two-rung ladder"
    for row in two_point:
        assert row[4] == "not_run", row


def test_exponents_table_fits_the_measured_load_slope() -> None:
    """The stale 1.17 belonged to the fleet ladder; the single ladder is 1.08."""
    rows = {r[0] + "|" + r[1]: r for r in _table("ref-exponents-table")}
    assert rows["ladder_mesh|load"][4] == "1.08"
    assert rows["ladder_fleet_10m_3pt|load"][4] == "1.17"


def test_single_body_table_reaches_the_1e7_band() -> None:
    """The synthesis document must show the 10^7 rung it claims in prose."""
    triangles = [
        int(row[2].replace(",", "")) for row in _table("ref-single-table")
    ]
    assert max(triangles) >= 10**7, max(triangles)


def test_import_table_joins_both_scenarios() -> None:
    """The file path is covered as a join of single + fleet, both present."""
    scenarios = {row[0] for row in _table("ref-import-table")}
    assert scenarios == {"单体", "多体"}, scenarios
    subtasks = {row[1] for row in _table("ref-import-table")}
    assert {"obj_write", "obj_read"} <= subtasks