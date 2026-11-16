"""End-to-end lock for the baseline-fleet geometry diagnostics (m20163 / m20487).

``docs/fleet_geometry_kernel_diagnostics.md`` and
``benchmarks/verified/baseline_fleet/geometry_diagnostics.json`` are the evidence
record for the type-128 control-net reader: which layout every patch was read as,
why that reading is admissible, and what the file's units are.  Two things can go
wrong with a record like that, and neither is caught by a test that needs the
publisher bytes (they are not redistributable, so everything that touches them is
skipped in CI):

* **the record drifts from the code.**  It did: the committed JSON said JBC "43
  weight_block + 2 interleaved, 4 ambiguous" and ONRT "106 + 26, 26 ambiguous",
  while the reader in ``main`` reports 45 and 132 ``weight_block`` and 0 ambiguous
  on the same bytes, because m20487 added the non-positive-weight rejection.
* **the record cannot be tied to the bytes it describes.**  It could not: every
  ``sha256`` in the JSON was ``null``, because :func:`tensorcad.iges.read_iges`
  never computed one, so a reader had no way to check that the counts were
  measured on the file the publishers ship.

So the tests here run the tool on *synthetic* files whose bytes this repository
owns (a format-faithful 80-column writer below), which makes the reader's layout
and units behaviour verifiable offline, and they gate the committed JSON against
itself (counts, histograms, SHA-256) and against the document that cites it.

Synthetic files, no publisher bytes, so the suite runs offline.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from tensorcad.fleet_baseline import partition_surfaces
from tensorcad.iges import read_iges

REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "benchmarks" / "verified" / "baseline_fleet" / "geometry_diagnostics.json"
DOC = REPO / "docs" / "fleet_geometry_kernel_diagnostics.md"


def _tool():
    """Import ``tools/baseline_fleet/diagnose_sources.py`` (no package there)."""
    spec = importlib.util.spec_from_file_location(
        "baseline_fleet_diagnose_sources", REPO / "tools" / "baseline_fleet" / "diagnose_sources.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# a format-faithful IGES writer (fixed 80-column sections, real DE/PD pointers)
# --------------------------------------------------------------------------- #
def _record(text: str, letter: str, sequence: int) -> str:
    return text[:72].ljust(72) + letter + f"{sequence:7d}"


def _parameter_lines(tokens: list[str], limit: int = 63) -> list[str]:
    """Pack comma-terminated tokens into parameter lines of at most ``limit`` chars.

    Each line must end on a delimiter: the reader takes columns 1-64 as the data
    and treats an unterminated line as continuing into columns 65-72.
    """
    lines: list[str] = []
    current = ""
    for index, token in enumerate(tokens):
        piece = token + (";" if index == len(tokens) - 1 else ",")
        if current and len(current) + len(piece) > limit:
            lines.append(current)
            current = ""
        current += piece
    if current:
        lines.append(current)
    return lines


def _global_lines(payload: str, limit: int = 70) -> list[str]:
    """Split the global section into delimiter-terminated lines of at most ``limit``."""
    pieces = payload.split(",")
    lines: list[str] = []
    current = ""
    for index, piece in enumerate(pieces):
        token = piece + ("," if index < len(pieces) - 1 else "")
        if current and len(current) + len(token) > limit:
            lines.append(current)
            current = ""
        current += token
    if current:
        lines.append(current)
    return lines


def write_iges(path: Path, entities: list[tuple[int, list[str]]], global_payload: str) -> Path:
    """Write a minimal IGES 5.3 file carrying ``entities`` = ``[(type, tokens)]``.

    The directory entry puts the parameter-line count in columns 25-32 (field 4
    of the second record) and the parameter data at the pointer in columns 9-16.
    """
    directory: list[str] = []
    parameters: list[str] = []
    parameter_sequence = 1
    for index, (entity_type, tokens) in enumerate(entities):
        chunks = _parameter_lines(tokens)
        pointer = parameter_sequence
        for chunk in chunks:
            parameters.append(_record(chunk, "P", parameter_sequence))
            parameter_sequence += 1
        directory_sequence = 2 * index + 1
        first = (
            f"{entity_type:8d}{pointer:8d}{0:8d}{0:8d}{0:8d}"
            f"{0:8d}{0:8d}{0:8d}{'00000000':>8s}"
        )
        second = (
            f"{entity_type:8d}{0:8d}{0:8d}{len(chunks):8d}{0:8d}"
            f"{0:8d}{0:8d}{0:8d}{'00000000':>8s}"
        )
        directory.append(_record(first, "D", directory_sequence))
        directory.append(_record(second, "D", directory_sequence + 1))
    global_lines = [
        _record(chunk, "G", index + 1)
        for index, chunk in enumerate(_global_lines(global_payload))
    ]
    lines = [
        _record("S      1", "S", 1),
        *global_lines,
        *directory,
        *parameters,
        _record("S      1", "T", 1),
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


#: The Rhino convention: explicit ``1H,`` / ``1H;`` delimiter fields, unit flag 2
#: (millimetres) and units name ``MM``.
GLOBAL_MM = (
    "1H,,1H;,4HTEST,4HT.igs,8Htensorca,8Htensorca,32,38,6,38,15,"
    "4HTEST,1.0,2,2HMM,1000,1.0,13H20261114.1200,1.0E-07,"
)
#: The KCS convention: delimiter fields left empty, unit flag 6 (metres).
GLOBAL_M = (
    "1H,,1H;,4HTEST,4HT.igs,8Htensorca,8Htensorca,32,38,6,38,15,"
    "4HTEST,1.0,6,1HM,1.0,1.0,13H20261114.1200,1.0E-07,"
)


def surface_tokens(kind: str, coords, weights=None, *, degree: int = 1, trailing=()) -> list[str]:
    """A type-128 record: header, knot vectors, then the control-point tail."""
    count = len(coords)
    side = int(round(count ** 0.5))
    count_u = count_v = side if side * side == count else count // 2 or 2
    if count_u * count_v != count:  # 2 x n nets
        count_u, count_v = 2, count // 2
    rational = "0" if kind != "polynomial" else "1"
    tokens = [
        "128", str(count_u - 1), str(count_v - 1), str(degree), str(degree),
        "0", "0", rational, "0", "0",
    ]
    knots = ["0", "0", "1", "1"]
    if degree == 2:
        knots = ["0", "0", "0", "1", "1", "1"]
    tokens += knots * 2
    if kind == "interleaved":
        for point, weight in zip(coords, weights or [1.0] * count):
            tokens += [repr(point[0]), repr(point[1]), repr(point[2]), repr(weight)]
    elif kind == "weight_block":
        tokens += [repr(weight) for weight in (weights or [1.0] * count)]
        for point in coords:
            tokens += [repr(point[0]), repr(point[1]), repr(point[2])]
    else:
        for point in coords:
            tokens += [repr(point[0]), repr(point[1]), repr(point[2])]
    tokens += [repr(value) for value in trailing]
    return tokens


#: A 2x2 patch of the ONRT kind: no interior control point, so the smoothness score
#: cannot decide the layout -- only the weight evidence can.
DEGENERATE = [(-9.0, 0.0, -5.0), (-9.0, 18.0, -5.0), (-6.0, 0.0, -5.0), (-6.0, 18.0, -5.0)]


@pytest.fixture()
def synthetic(tmp_path):
    path = write_iges(
        tmp_path / "synth.igs",
        [
            # 2x2 NMRI-style weight block: the tail opens with 1, 1, 1, 1 weights.
            (128, surface_tokens("weight_block", DEGENERATE, trailing=(0.0, 1.0, 0.0, 1.0))),
            # 2x2 IGES-specification interleaved quadruples.
            (128, surface_tokens("interleaved", [(0.0, 0.0, 0.0), (0.0, 9.0, 0.0),
                                                 (3.0, 0.0, 0.0), (3.0, 9.0, 0.0)])),
            # genuinely rational 2x2: varying weights, so the weight prior is silent
            (128, surface_tokens("weight_block", [(0.0, 0.0, 0.0), (0.0, 4.0, 0.0),
                                                  (2.0, 0.0, 0.0), (2.0, 4.0, 0.0)],
                                 weights=[1.0, 0.9, 1.1, 1.0])),
        ],
        GLOBAL_MM,
    )
    return path


# --------------------------------------------------------------------------- #
# 1. the reader, end to end from bytes: layout and units
# --------------------------------------------------------------------------- #
def test_sha256_binds_the_parsed_model_to_the_bytes_it_came_from(synthetic):
    model = read_iges(str(synthetic))
    assert model.sha256 == hashlib.sha256(synthetic.read_bytes()).hexdigest()
    assert len(model.sha256) == 64


def test_millimetre_global_section_is_read_as_millimetres(synthetic):
    model = read_iges(str(synthetic))
    assert (model.unit_flag, model.unit_name, model.units) == (2, "MM", "mm")
    assert model.metres_per_unit == pytest.approx(1e-3)
    assert model.unit_source == "positional"


def test_metre_global_section_stays_metres(tmp_path):
    path = write_iges(tmp_path / "m.igs", [(128, surface_tokens("weight_block", DEGENERATE))], GLOBAL_M)
    model = read_iges(str(path))
    assert (model.unit_flag, model.unit_name, model.metres_per_unit) == (6, "M", 1.0)


def test_degenerate_weight_block_patch_is_not_misread_as_interleaved(synthetic):
    """The m20163 defect: ``interleaved`` used to win by tuple sort order.

    Reading the weight block ``1, 1, 1, 1`` as ``x, y, z, w`` invents a control
    point at (1, 1, 1), i.e. a unit-cube net on an 18 m-beam hull.
    """
    surfaces = read_iges(str(synthetic)).surfaces()
    block = surfaces[0]
    assert block.layout == "weight_block"
    assert block.layout_ambiguous is False
    assert [tuple(round(value, 12) for value in point) for point in block.control] == [
        tuple(round(value, 12) for value in point) for point in DEGENERATE
    ]
    assert block.weights == (1.0, 1.0, 1.0, 1.0)


def test_interleaved_patch_is_still_read_as_interleaved(synthetic):
    interleaved = read_iges(str(synthetic)).surfaces()[1]
    assert interleaved.layout == "interleaved"
    assert interleaved.layout_ambiguous is False


def test_genuinely_rational_patch_keeps_its_varying_weights(synthetic):
    rational = read_iges(str(synthetic)).surfaces()[2]
    assert rational.layout == "weight_block"
    assert rational.weights == (1.0, 0.9, 1.1, 1.0)
    assert rational.layout_ambiguous is False, "the competing reading needs a non-positive weight"


def test_a_true_tie_is_flagged_and_quarantined_instead_of_silently_chosen(tmp_path):
    """Fail closed: if both readings are admissible and perfectly smooth, say so.

    A positive, *non-constant* weight block makes the ``interleaved`` reading
    admissible too (its weight slots then hold coordinates, all positive here) and
    neither reading has constant weights, so the weight prior is silent; on a 2x2 net
    both are perfectly smooth, so the record genuinely cannot decide it.
    The patch must come back ambiguous and be dropped by the partition rather than
    meshed from a guess.
    """
    coords = [(1.0, 1.0, 1.0), (1.0, 2.0, 1.0), (2.0, 3.0, 1.0), (2.0, 4.0, 1.0)]
    path = write_iges(
        tmp_path / "tie.igs",
        [(128, surface_tokens("weight_block", coords, weights=[1.0, 0.9, 1.1, 1.0]))],
        GLOBAL_M,
    )
    surface = read_iges(str(path)).surfaces()[0]
    assert surface.layout_ambiguous is True
    usable, quarantined = partition_surfaces([surface])
    assert len(usable) == 0
    assert "ambiguous control-point layout" in quarantined[0]["reason"]


# --------------------------------------------------------------------------- #
# 2. the diagnostics tool: what a reviewer reads
# --------------------------------------------------------------------------- #
def test_diagnose_reports_units_layout_evidence_and_provenance(synthetic, tmp_path):
    tool = _tool()
    report = tool.diagnose("SYNTH", synthetic, root=tmp_path)
    assert report["file"] == "synth.igs", "the record must be portable, not a worktree path"
    assert report["sha256"] == hashlib.sha256(synthetic.read_bytes()).hexdigest()
    assert (report["unit_flag"], report["units"]) == (2, "MM")
    assert report["metres_per_unit"] == pytest.approx(1e-3)
    assert report["patch_count"] == 3
    assert report["layout_histogram"] == {"weight_block": 2, "interleaved": 1}
    assert report["ambiguous_patches"] == 0
    assert report["usable_patches"] == 3 and report["quarantined_patches"] == 0
    assert report["transformation_matrices"] == 0
    assert report["trailing_values_max"] == 4, "writer noise after the control data is counted"
    # Every axis is reported in metres, so the 0.001 factor is visible in the numbers.
    spans = [box["span"] for box in report["control_cloud_axes_m"]]
    assert spans[0] == pytest.approx(0.012, abs=1e-9)
    assert spans[1] == pytest.approx(0.018, abs=1e-9)


def test_diagnose_explains_why_each_patch_was_not_ambiguous(synthetic, tmp_path):
    tool = _tool()
    patches = tool.diagnose("SYNTH", synthetic, root=tmp_path)["patches"]
    block, interleaved, rational = patches
    assert block["layout_evidence"] == "unique admissible reading"
    assert block["admissible_layouts"] == ["weight_block"]
    assert block["weights_uniform"] is True
    assert interleaved["layout_evidence"] == "unique admissible reading"
    assert interleaved["admissible_layouts"] == ["interleaved"]
    assert rational["weights_uniform"] is False
    assert rational["weight_count"] == 3


def test_main_writes_the_json_and_names_a_missing_hull(tmp_path):
    tool = _tool()
    raw = tmp_path / "raw"
    (raw / "HAVE").mkdir(parents=True)
    write_iges(raw / "HAVE" / "hull.igs", [(128, surface_tokens("weight_block", DEGENERATE))], GLOBAL_MM)
    specs = tmp_path / "fleet_specs.json"
    specs.write_text(json.dumps([
        {"name": "HAVE", "draft": 5.0, "length_pp": 100.0, "beam": 20.0,
         "published_volume_m3": None, "geometry_scale": 1.0,
         "iges_path": "HAVE/hull.igs"},
        {"name": "GONE", "draft": 5.0, "length_pp": 100.0, "beam": 20.0,
         "published_volume_m3": None, "geometry_scale": 1.0,
         "iges_path": "GONE/hull.igs"},
    ]))
    out = tmp_path / "diag.json"
    assert tool.main(["--raw", str(raw), "--specs", str(specs), "--json", str(out)]) == 0
    report = json.loads(out.read_text())
    assert report["schema"].startswith("tensorcad.baseline_fleet.diagnostics/")
    assert report["hulls"]["HAVE"]["patch_count"] == 1
    assert report["hulls"]["GONE"]["status"] == "NOT_FETCHED"
    assert "fetch_sources.py" in report["hulls"]["GONE"]["remedy"]


# --------------------------------------------------------------------------- #
# 3. the committed record: self-consistent, bound to bytes, and cited by the doc
# --------------------------------------------------------------------------- #
def test_committed_evidence_is_present_and_internally_consistent():
    report = json.loads(EVIDENCE.read_text())
    assert report["schema"] == "tensorcad.baseline_fleet.diagnostics/1.1"
    assert set(report["hulls"]) == {"KCS", "KVLCC2", "JBC", "ONRT"}
    for name, hull in report["hulls"].items():
        assert hull.get("status") != "NOT_FETCHED", name
        assert sum(hull["layout_histogram"].values()) == hull["patch_count"], name
        assert len(hull["patches"]) == hull["patch_count"], name
        assert hull["usable_patches"] + hull["quarantined_patches"] == hull["patch_count"], name
        assert hull["ambiguous_patches"] == hull["quarantined_patches"], name
        assert sum(hull["layout_evidence_histogram"].values()) == hull["patch_count"], name
        assert not Path(hull["file"]).is_absolute(), f"{name}: record must be portable"


def test_committed_evidence_is_bound_to_the_source_bytes():
    """Every hull carries the SHA-256 of the file it was measured on.

    This is the gate that the original record failed: ``sha256`` was ``null`` for
    all four hulls because the reader never computed one, so nothing tied the
    counts to the publisher's bytes.
    """
    report = json.loads(EVIDENCE.read_text())
    for name, hull in report["hulls"].items():
        digest = hull["sha256"]
        assert isinstance(digest, str) and len(digest) == 64, name
        int(digest, 16)
        assert digest.rstrip("0") != "", name


def test_uniform_weights_field_does_not_claim_the_weight_test_decided():
    """The field is named for what it measures.

    It used to be ``decided_by_weight_test`` and was computed as "the weights are
    constant", which counted the 55 KCS polynomial patches as decided by a test
    that never ran on them (a polynomial record has no weight block to test).
    """
    report = json.loads(EVIDENCE.read_text())
    kcs = report["hulls"]["KCS"]
    assert "patches_decided_by_weight_test" not in kcs
    assert kcs["patches_with_uniform_weights"] == kcs["patch_count"]
    assert kcs["layout_evidence_histogram"].get("weight prior (constant positive weights)", 0) == 0


def test_document_cites_the_digests_of_the_recorded_files():
    """The document and the JSON must describe the same bytes."""
    report = json.loads(EVIDENCE.read_text())
    text = DOC.read_text()
    for name, hull in report["hulls"].items():
        assert hull["sha256"][:16] in text, f"{name}: digest {hull['sha256'][:16]} not cited in the doc"