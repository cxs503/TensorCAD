"""``hull_characteristics`` must not halve a model that carries both sides.

Before this fix, ``hull_characteristics`` integrated ``envelope_breadth = max |y|`` at
every station whatever the model was, and then multiplied by ``symmetry_factor``.  For a
half-breadth sheet that is right (one side, doubled).  For a **both-sides** model
``symmetry_factor == 1`` and ``max |y|`` is exactly *half* the breadth, so A(x), the
displacement and the beam all came out 2x low, while ``Cb = volume / (Lwl * beam * T)``
stayed correct because both its numerator and its denominator shrank by the same two.
That is why the defect survived the existing suite: the only full-model test asserted the
*wetted area* (which never used ``envelope_breadth``), and Cb looked perfect.

The reference here is closed form: the box below has V = L*B*T, B = 4 m, Cb = Cm = 1 by
construction, so the numbers are not another implementation's opinion.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.hull_surface import (  # noqa: E402
    classify_model,
    envelope_breadth,
    hull_characteristics,
    validate_against_particulars,
    _Slabs,
)

L, B, T = 10.0, 4.0, 2.0


def _full_box():
    """Both sides: y in [-B/2, +B/2]."""
    vertices = [(0.0, -B / 2, 0.0), (L, -B / 2, 0.0), (L, B / 2, 0.0), (0.0, B / 2, 0.0),
                (0.0, -B / 2, T), (L, -B / 2, T), (L, B / 2, T), (0.0, B / 2, T)]
    faces = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
             (2, 3, 7), (2, 7, 6), (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5)]
    return vertices, faces


def _half_box():
    """One side: y in [0, B/2], the publisher half-breadth convention."""
    vertices = [(0.0, 0.0, 0.0), (L, 0.0, 0.0), (L, B / 2, 0.0), (0.0, B / 2, 0.0),
                (0.0, 0.0, T), (L, 0.0, T), (L, B / 2, T), (0.0, B / 2, T)]
    faces = [(0, 1, 2), (0, 2, 3), (4, 6, 5), (4, 7, 6), (1, 5, 4), (1, 4, 0),
             (2, 6, 7), (2, 7, 3), (1, 2, 6), (1, 6, 5)]
    return vertices, faces


def test_envelope_breadth_returns_the_two_sided_span_for_full_models():
    vertices, faces = _full_box()
    slabs = _Slabs.build(vertices, faces)
    station = L / 2
    assert envelope_breadth(slabs, station, T / 2) == pytest.approx(B / 2)  # one side
    assert envelope_breadth(slabs, station, T / 2, full=True) == pytest.approx(B)


def test_full_box_displacement_beam_and_coefficients_are_exact():
    vertices, faces = _full_box()
    characteristic = hull_characteristics(vertices, faces, waterline=T)
    assert characteristic.symmetry_factor == 1
    assert characteristic.volume == pytest.approx(L * B * T, rel=1e-12)
    assert characteristic.beam == pytest.approx(B, rel=1e-9)
    assert characteristic.lwl == pytest.approx(L, rel=1e-9)
    assert characteristic.cb == pytest.approx(1.0, rel=1e-9)
    assert characteristic.cm == pytest.approx(1.0, rel=1e-9)


def test_forcing_half_breadth_readings_reproduces_the_halving():
    """One keyword, one measurement: the pre-fix value is a flag away, not a claim."""
    vertices, faces = _full_box()
    before = hull_characteristics(vertices, faces, waterline=T, full_breadth=False)
    after = hull_characteristics(vertices, faces, waterline=T)
    assert before.volume == pytest.approx(after.volume / 2.0, rel=1e-12)
    assert before.beam == pytest.approx(after.beam / 2.0, rel=1e-12)
    # ... and Cb cannot see any of it, which is how the defect hid.
    assert before.cb == pytest.approx(after.cb, rel=1e-12)
    assert before.cm == pytest.approx(after.cm, rel=1e-12)


def test_half_breadth_path_is_unchanged():
    vertices, faces = _half_box()
    assert classify_model(vertices)["half_breadth"] is True
    characteristic = hull_characteristics(vertices, faces, waterline=T)
    assert characteristic.symmetry_factor == 2
    assert characteristic.volume == pytest.approx(L * B * T, rel=1e-12)
    assert characteristic.beam == pytest.approx(B, rel=1e-9)
    # forcing the flag on a one-sided model is a no-op (y_min is 0 up to weld noise)
    forced = hull_characteristics(vertices, faces, waterline=T, full_breadth=True)
    assert forced.volume == pytest.approx(characteristic.volume, rel=1e-9)
    assert forced.beam == pytest.approx(characteristic.beam, rel=1e-9)


def test_a_model_living_entirely_on_positive_y_is_declared_half_breadth():
    """The classifier's convention, and the escape hatch when it is wrong for you.

    ``classify_model`` decides "half-breadth" from ``y_min >= 0`` alone.  Vertices cannot
    tell a one-sided sheet from a both-sides model that happens to lie at ``y > 0``, and
    no ship model is of the second kind, so the convention stays.  What a caller with an
    unusual frame needs is the override: ``full_breadth=True`` reads the two-sided span
    whatever the convention says.
    """
    vertices, faces = _full_box()
    shifted = [(x, y + 3.0, z) for (x, y, z) in vertices]
    assert classify_model(shifted)["half_breadth"] is True  # convention, not truth
    forced = hull_characteristics(shifted, faces, waterline=T, full_breadth=True)
    # forced full span: A(x) = B*T, symmetry factor stays the classifier's 2 -> 2x too big,
    # which is why the switch is a diagnostic and not a default for offset frames.
    assert forced.volume == pytest.approx(2.0 * L * B * T, rel=1e-12)


def test_validate_against_particulars_passes_a_full_model_now():
    vertices, faces = _full_box()
    characteristic = hull_characteristics(vertices, faces, waterline=T)
    verdict = validate_against_particulars(
        characteristic,
        {"volume": L * B * T, "cb": 1.0, "lpp": L, "beam": B, "source": "closed form"})
    assert verdict["quantities"]["volume"]["error_percent"] == pytest.approx(0.0, abs=1e-9)
    assert verdict["passed"] is True
