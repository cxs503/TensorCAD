"""Entity 142 is read by field position: the ``IT`` flag is not a pointer.

IGES 5.3 section 4.25 defines the record as ``142, CRV, SPTR, BPTR, CPTR, IT``.  The
reader used to resolve *every* integer token that named an existing curve, so a preference
flag of ``1`` became a reference to directory entry 1.  On the KVLCC2 archive that is a
type-126 curve of another face, 28.4 m from the surface being trimmed: every trimmed patch
gained a phantom third loop, its real parameter-space image was vetoed against that
phantom (measured disagreement 1.33e1 model units) and the face was refused.

These tests pin the field map on the two layouts that actually occur in the archives, and
on the real KVLCC2 file when the asset is on the machine.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from tensorcad.geomquality.iges_trim import DirectoryEntry, _fields_of_142, read_iges_trimmed

KVLCC2 = Path("/root/tp-agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs")

SURFACE = 128
COMPOSITE = 102
CURVE = 126
PLANE = 108


def _entry(sequence: int, entity_type: int, tokens: list[str]) -> DirectoryEntry:
    return DirectoryEntry(
        sequence=sequence,
        entity_type=entity_type,
        form=0,
        transformation=0,
        level=0,
        status="00000000",
        parameter_lines=1,
        tokens=tuple(tokens),
        raw="",
    )


def _document(entries: dict[int, DirectoryEntry]) -> SimpleNamespace:
    """``_fields_of_142`` only needs ``document.entries.get``."""
    return SimpleNamespace(entries=entries)


def _builder(document) -> SimpleNamespace:
    return SimpleNamespace(document=document)


def test_kvlcc2_layout_the_it_flag_is_not_a_pointer() -> None:
    """``142,0,5,33,19,1,0,0``: CRV=0, SPTR=5, BPTR=33, CPTR=19, IT=1 (DE 35)."""
    entries = {
        5: _entry(5, SURFACE, ["11", "5", "5"]),
        33: _entry(33, COMPOSITE, ["6", "21", "23"]),
        19: _entry(19, COMPOSITE, ["6", "7", "9"]),
        1: _entry(1, CURVE, ["1", "1", "0"]),  # the phantom: another face's curve
    }
    base, fields = _fields_of_142(
        _builder(_document(entries)), _entry(35, 142, ["0", "5", "33", "19", "1", "0", "0"])
    )
    assert base == 0
    assert fields == {"CRV": 0, "SPTR": 5, "BPTR": 33, "CPTR": 19, "IT": 1}
    assert 1 not in (fields["CRV"], fields["BPTR"], fields["CPTR"])


def test_kp505_layout_has_no_parameter_curve_and_a_composite() -> None:
    """``142,0,1,0,11,2,0,0``: CRV=0, SPTR=1 (plane), BPTR=0, CPTR=11, IT=2 (DE 13)."""
    entries = {
        1: _entry(1, PLANE, ["1.0", "0.0", "0.0", "-28.5"]),
        11: _entry(11, COMPOSITE, ["2", "5", "9"]),
    }
    _, fields = _fields_of_142(
        _builder(_document(entries)), _entry(13, 142, ["0", "1", "0", "11", "2", "0", "0"])
    )
    assert fields == {"CRV": 0, "SPTR": 1, "BPTR": 0, "CPTR": 11, "IT": 2}


def test_record_that_keeps_the_entity_type_token_aligns_on_the_surface() -> None:
    """A writer that does not strip the leading ``142`` is one field longer."""
    entries = {
        5: _entry(5, SURFACE, ["11"]),
        33: _entry(33, COMPOSITE, ["6"]),
        19: _entry(19, COMPOSITE, ["6"]),
    }
    base, fields = _fields_of_142(
        _builder(_document(entries)),
        _entry(35, 142, ["142", "0", "5", "33", "19", "1", "0", "0"]),
    )
    assert base == 1
    assert fields["SPTR"] == 5
    assert fields["BPTR"] == 33
    assert fields["CPTR"] == 19


def test_trailing_padding_that_names_a_curve_is_still_not_a_pointer() -> None:
    """Exporter padding (``,0,0`` on these files) must not be resolved either way."""
    entries = {7: _entry(7, CURVE, ["1"]), 8: _entry(8, CURVE, ["1"])}
    _, fields = _fields_of_142(
        _builder(_document(entries)), _entry(9, 142, ["0", "0", "7", "8", "0", "7", "8"])
    )
    # fields past IT are outside the record's definition: they are not returned at all
    assert fields["BPTR"] == 7 and fields["CPTR"] == 8
    assert fields["IT"] == 0


@pytest.mark.skipif(not KVLCC2.exists(), reason="KVLCC2 IGES asset is not redistributable")
def test_kvlcc2_every_trimmed_face_has_exactly_one_loop() -> None:
    """15 of 19 faces used to carry three loops (the phantom third one included)."""
    model = read_iges_trimmed(str(KVLCC2))
    assert len(model.faces) == 19
    assert [len(face.loops) for face in model.faces] == [1] * 19
    assert all(face.trim_routes() == ["parametric"] for face in model.faces)
    # the cross-check is reported, never used as a veto: the record's 3-D composite is a
    # different curve on most patches, and the parametric trim is kept anyway
    assert all(value is None or value >= 0.0 for face in model.faces for value in face.cross_checks)