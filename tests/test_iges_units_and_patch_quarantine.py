"""Global-section units and control-point-layout regressions (m17800).

Both defects surfaced while grading the baseline fleet against the publishers' own
tables:

* the units flag was located by a content scan that returns the first real number
  followed by a small integer.  Global fields 7-11 (precision magnitudes) look exactly
  like that, so the JBC file -- millimetres, flag 2 -- was read as metres (flag 6), a
  silent factor of 1000 in every length, area and volume;
* one ambiguous 3x3 patch (JBC DE 386, a 3.4 m x 0.5 m flat piece) raised, which made
  all 45 patches of that hull unreadable.

The fixtures are synthetic global-section payloads (no publisher bytes needed); the last
class of tests runs against the cached workshop files when they are present and is
skipped otherwise.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.fleet_baseline import partition_surfaces
from tensorcad.iges import (
    BSplineSurface,
    _global_units,
    _split_control,
    hollerith_fields,
    read_iges,
)

RAW = Path("artifacts/baseline_fleet/raw")

RAW_GLOBALS = {
    "JBC": RAW / "JBC" / "jbc.igs",
    "KCS": RAW / "KCS" / "KCS_hull_SVA.igs",
    "KVLCC2": RAW / "KVLCC2" / "kvlcc2.igs",
    "ONRT": RAW / "ONRT" / "onr5613-full_appendage_deckhous_mergedR.igs",
}
have = {name: path.exists() for name, path in RAW_GLOBALS.items()}
need = pytest.mark.skipif(not all(have.values()), reason="cached workshop IGES absent")


def _global_fields(name: str) -> tuple[str, ...]:
    """Tokenise a cached file's global section (ground truth: the publisher's bytes)."""
    from tensorcad.iges import _split_sections

    lines = RAW_GLOBALS[name].read_text(errors="replace").splitlines()
    sections = _split_sections(lines)
    return hollerith_fields("".join(line[0:72] for line in sections["G"]))


@need
def test_units_flag_lands_on_the_same_token_for_both_writer_conventions():
    """KCS leaves the two delimiter fields empty, the Rhino files write ``1H,``/``1H;``.

    A naive ``str.split(",")`` puts those two conventions one token apart, which is what
    the old content scan tripped over; Hollerith tokenising must align them.
    """
    kcs = _global_fields("KCS")
    jbc = _global_fields("JBC")
    assert (kcs[13], kcs[14]) == ("6", "M")
    assert (jbc[13], jbc[14]) == ("2", "MM")


@need
def test_the_old_content_scan_would_have_read_jbc_as_metres():
    """Regression: JBC is millimetres, and the historical fallback still says metres."""
    from tensorcad.iges import _unit_flag_content_scan

    fields = _global_fields("JBC")
    assert _unit_flag_content_scan(fields) == 6  # what shipped before m17800
    flag, name, scale, source, warnings = _global_units(fields)
    assert (flag, name, scale, source, warnings) == (2, "MM", 1.0, "positional", [])


@need
def test_catomatic_file_is_millimetres_with_unit_name_cross_check():
    flag, name, _scale, source, _warn = _global_units(_global_fields("KVLCC2"))
    assert (flag, name, source) == (2, "MM", "positional")


@need
def test_units_name_overrides_a_contradicting_flag():
    fields = list(_global_fields("JBC"))
    fields[13] = "6"  # a writer that declares metres but names MM
    flag, name, _scale, source, warnings = _global_units(tuple(fields))
    assert (flag, name, source) == (2, "MM", "name_overrides_flag")
    assert warnings and "disagrees" in warnings[0]


@need
def test_model_space_scale_is_read_and_applied():
    fields = list(_global_fields("KCS"))
    fields[12] = "2.0D0"
    _flag, _name, scale, _source, _warn = _global_units(tuple(fields))
    assert scale == 2.0


@need
def test_ambiguous_layout_is_reported_not_raised():
    """The JBC DE-386 tail against the real bytes: 46 values for 9 control points.

    Before m17800 this raised IgesError and made all 45 patches of the hull unreadable.
    """
    from tensorcad.iges import _split_control

    model = read_iges(str(RAW_GLOBALS["JBC"]))
    values = model.raw_parameters[386]
    count_u = int(values[1]) + 1
    count_v = int(values[2]) + 1
    cursor = 10 + (count_u + int(values[3]) + 1) + (count_v + int(values[4]) + 1)
    body = values[cursor:]
    control, weights, layout, ambiguous = _split_control(
        body, points=count_u * count_v, rational=int(values[7]) == 0,
        count_u=count_u, count_v=count_v,
    )
    assert len(control) == len(weights) == count_u * count_v
    assert layout in {"interleaved", "weight_block", "polynomial"}
    assert ambiguous is True


def _surface(ambiguous: bool) -> BSplineSurface:
    return BSplineSurface(
        degree_u=1, degree_v=1, knots_u=(0.0, 0.0, 1.0, 1.0), knots_v=(0.0, 0.0, 1.0, 1.0),
        control=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
        weights=(1.0, 1.0, 1.0, 1.0), count_u=2, count_v=2,
        layout="weight_block", layout_ambiguous=ambiguous,
    )


def test_partition_surfaces_quarantines_only_the_ambiguous_patch():
    usable, quarantined = partition_surfaces([_surface(False), _surface(True), _surface(False)])
    assert len(usable) == 2 and len(quarantined) == 1
    assert quarantined[0]["patch_index"] == 1
    assert "ambiguous" in quarantined[0]["reason"]


@pytest.mark.skipif(not (RAW / "JBC" / "jbc.igs").exists(), reason="cached workshop IGES absent")
def test_cached_jbc_is_millimetres_and_parses_completely():
    model = read_iges(str(RAW / "JBC" / "jbc.igs"))
    assert (model.unit_flag, model.unit_name) == (2, "MM")
    assert model.metres_per_unit == pytest.approx(1e-3)
    surfaces = model.surfaces()  # must not raise on the ambiguous patch
    assert len(surfaces) == 45
    usable, quarantined = partition_surfaces(surfaces)
    assert (len(usable), len(quarantined)) == (44, 1)


@pytest.mark.skipif(not (RAW / "KCS" / "KCS_hull_SVA.igs").exists(), reason="cached IGES absent")
def test_cached_kcs_stays_metres():
    model = read_iges(str(RAW / "KCS" / "KCS_hull_SVA.igs"))
    assert (model.unit_flag, model.unit_name, model.metres_per_unit) == (6, "M", 1.0)
    assert len(partition_surfaces(model.surfaces())[1]) == 0
