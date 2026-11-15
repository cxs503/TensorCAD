"""``docs/hull_family.md`` is generated: hold it to the evidence it is generated from.

The page publishes 15 hydrostatic rows, five derived ``< 3 % ?`` verdicts, a
transom/bulb pass-through, a refusal matrix and a gate-power table.  Not one of
those numbers is typed into the markdown -- all of them are injected by
``tools/hull_family/gen_doc_and_gate.py:build_doc`` from
``docs/hull_family/evidence.json`` (``tools/hull_family/report.py``) and
``docs/hull_family/gate.json`` (the mutation run).  That construction is worth
nothing unless the two JSON files are *in the tree* and the page still equals what
they produce, so this module asserts exactly that:

1. both files sit at the paths the page cites, and those paths are in ``docs/`` --
   ``artifacts/`` is scratch and a merge drops it, which is how an earlier revision
   of this page came to cite JSON that was not in the repository at all;
2. ``build_doc(evidence, gate)`` reproduces the committed page byte for byte, so a
   stale published verdict (or a number re-typed by hand) is red;
3. the figure on the page is the figure the evidence describes;
4. the gate evidence was produced against *this* ``tensorcad/hull.py``, and every
   mutant died -- a survivor would mean a documented claim is not asserted;
5. the five verdicts the page publishes are what the evidence actually says.

Regenerating the page (``python tools/hull_family/report.py`` and
``python tools/hull_family/gen_doc_and_gate.py``) is what makes these tests pass
after a deliberate change; the tests are the reason the page cannot drift silently.

The mutation patterns in ``gen_doc_and_gate.py`` are matched against
``tensorcad/hull.py``; editing that file's family dispatch therefore requires a
re-run of the gate, which is exactly the coupling this module pins down (test 4
fails with the digest of the edited file against the recorded one).
"""
from __future__ import annotations

import hashlib
import json
import struct
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "hull_family"))

import gen_doc_and_gate as tool  # noqa: E402

PAGE = ROOT / "docs" / "hull_family.md"
EVIDENCE = ROOT / "docs" / "hull_family" / "evidence.json"
GATE = ROOT / "docs" / "hull_family" / "gate.json"


@pytest.fixture(scope="module")
def evidence() -> dict:
    return json.loads(EVIDENCE.read_text())


@pytest.fixture(scope="module")
def gate() -> dict:
    return json.loads(GATE.read_text())


def test_evidence_lives_in_the_tree_beside_the_page() -> None:
    page = PAGE.read_text()
    for path in (EVIDENCE, GATE):
        assert path.is_file(), f"{path.relative_to(ROOT)} is cited by the page but absent"
        cited = str(path.relative_to(ROOT))
        assert f"`{cited}`" in page, f"the page does not point at {cited}"
        assert path.relative_to(ROOT).parts[0] == "docs", (
            "evidence must live where a merge keeps it, not under artifacts/")


def test_the_tools_read_the_same_paths_the_page_cites() -> None:
    assert tool.DATA == EVIDENCE
    assert tool.GATE == GATE
    assert EVIDENCE.relative_to(ROOT) == Path(tool.DATA.relative_to(ROOT))


def test_page_is_reproduced_from_the_evidence(evidence: dict, gate: dict) -> None:
    committed = PAGE.read_text()
    regenerated = tool.build_doc(evidence, gate)
    if committed != regenerated:
        diff = [
            (old, new) for old, new in zip(committed.splitlines(), regenerated.splitlines())
            if old != new
        ]
        pytest.fail(
            "docs/hull_family.md is stale: it differs from what the evidence produces "
            f"({len(diff)} line(s), first: {diff[:2] if diff else 'length differs'}). "
            "Re-run tools/hull_family/gen_doc_and_gate.py."
        )


def test_figure_is_the_one_the_evidence_describes(evidence: dict) -> None:
    record = evidence["figure"]
    png = ROOT / record["png"]
    assert png.is_file() and record["png"] == "docs/hull_family.png"
    data = png.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert record["dpi"] == 150
    width, height = record["pixels"]
    assert data[16:24] == struct.pack(">II", width, height), (
        "the recorded pixel size is not the PNG's own IHDR size: the evidence describes a "
        "different figure than the one on disk")
    assert width >= 1200 and height >= 400, (
        f"the figure is {width} x {height} px, too small to read the three panels")

    matplotlib = pytest.importorskip("matplotlib")
    if matplotlib.__version__ == record["matplotlib"]:
        assert hashlib.sha256(data).hexdigest() == record["sha256"], (
            "docs/hull_family.png is not the figure recorded in the evidence")
    else:
        pytest.skip(
            f"sha256 is renderer-specific (evidence: matplotlib {record['matplotlib']}, "
            f"here {matplotlib.__version__}); pixel geometry is checked above")


def test_gate_evidence_is_about_this_hull_module(gate: dict) -> None:
    digest = hashlib.sha256((ROOT / "tensorcad" / "hull.py").read_bytes()).hexdigest()
    assert gate["file"] == "tensorcad/hull.py"
    assert gate["sha256_before"] == digest, (
        "tensorcad/hull.py changed since the gate ran: re-run "
        "tools/hull_family/gen_doc_and_gate.py (its mutant patterns may no longer match)")
    assert gate["sha256_after"] == digest, "the mutation run did not restore hull.py"
    assert gate["total"] == len(tool.MUTANTS) == 5
    assert gate["killed"] == gate["total"], (
        "a mutant survived: a claim the page documents is not asserted by any test -- "
        + ", ".join(row["name"] for row in gate["mutants"] if row["verdict"] != "KILLED"))


def test_published_verdicts_are_what_the_evidence_says(evidence: dict) -> None:
    rows = [row for entry in evidence["parents"].values() for row in entry["grids"]]
    assert len(rows) == 15  # five parents x three grids, as the table shows

    # The page's band table is the published resolution (161 x 41); the coarse grids
    # are printed beside it as the documented limitation (~4 % at 41 x 13, which trips
    # the module's 2 % sampling warning).  So the bands are asserted where they are
    # claimed, and the coarse grid must still be in the evidence -- not silently dropped.
    fine = [row for row in rows if row["nx"] == 161]
    assert len(fine) == 5
    worst_fine = max(abs(row["volume_error_vs_closed_form"]) for row in fine)
    assert worst_fine < 0.03, f"displacement error {100 * worst_fine:.2f} % breaches the band"
    assert max(abs(row["block_coefficient_error"]) for row in fine) < 0.03
    assert max(abs(row["waterplane_error_vs_closed_form"]) for row in fine) < 0.03
    assert max(abs(row["vcb_error_vs_closed_form"]) for row in fine) < 0.03
    coarse = [row for row in rows if row["nx"] == 41]
    assert len(coarse) == 5 and max(
        abs(row["volume_error_vs_closed_form"]) for row in coarse) > 0.02, (
        "the page documents a ~4 % coarse-grid error; the evidence no longer shows one")

    # bit-for-bit delegation, transom width and the bulb's non-interference
    assert all(row["identical_to_family_module"] for row in rows)
    transom = evidence["transom"]
    assert transom["watertight"] is True
    assert abs(transom["aft_station_width_m"] / transom["expected_width_m"] - 1.0) < 1e-9
    bulb = evidence["bulb"]
    assert bulb["hull_displacement_unchanged"] is True

    # the refusal matrix: every listed case is refused and says why
    assert len(evidence["refusals"]) >= 6
    for row in evidence["refusals"]:
        assert row["refused"] is True, row["case"]
        assert row["message_has_needle"] is True, row["case"]