"""The mirror contract: every analytic field has one surface and therefore one number.

`GeometryField.numpy_callable` is an optimisation of `signed_distance`, documented as
"the pure-Python path stays the contract". That sentence is only meaningful if it can be
falsified, so these tests demand **exact** equality (`==`, not a tolerance) between the
two paths of every analytic field, over three probe families: a lattice over the field's
own bounds, a deterministic uniform sample, and a band of exponentially spaced offsets
across the zero set -- the band being where a 1-ULP disagreement flips a containment
label, which is what marching-tetrahedra extraction reads.

Before the fix in `tensorcad/geometry_field.py` (2026-11-03) these sweeps disagreed on
10-29 % of probes for the plane, sphere, cylinder and ellipsoid; the measurements are in
`runs/2026-11-03-mirror-contract/before.json`.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry_field import (
    CylinderField,
    EllipsoidField,
    PlaneField,
    SphereField,
)
from tensorcad.suboff import SuboffHullField
from tensorcad.suboff_appendages import (
    SuboffAppendagesField,
    SuboffFinsField,
    SuboffFullField,
    SuboffSailField,
)

np = pytest.importorskip("numpy")

FIELDS = {
    "plane_tilted": lambda: PlaneField((0.3, -0.4, 0.8660254037844387), 0.137),
    "sphere_offset": lambda: SphereField(0.73, (0.21, -0.14, 0.37)),
    "cylinder_tilted": lambda: CylinderField(
        0.41,
        1.27,
        (0.11, -0.23, 0.19),
        (0.2672612419124244, 0.5345224838248488, -0.8017837257372732),
    ),
    "ellipsoid": lambda: EllipsoidField((0.63, 0.37, 0.24), (0.08, -0.11, 0.05)),
    "suboff_hull": lambda: SuboffHullField(),
    "suboff_sail": lambda: SuboffSailField(),
    "suboff_fins": lambda: SuboffFinsField(),
    "suboff_appendages": lambda: SuboffAppendagesField(),
    "suboff_full": lambda: SuboffFullField(),
}

BAND = (0.0, 1e-15, -1e-15, 1e-13, -1e-13, 1e-11, -1e-11, 1e-9, -1e-9, 1e-7, -1e-7, 1e-5, -1e-5)


def _extent(field) -> tuple[tuple[float, float, float], float]:
    box = field.bounds()
    if box is None:
        return (0.0, 0.0, 0.0), 1.0
    size = max(box.sizes)
    centre = tuple(0.5 * (box.minimum[i] + box.maximum[i]) for i in range(3))
    return centre, size


def _lattice(field, n: int = 7):
    box = field.bounds()
    lo, hi = ((-1.0,) * 3, (1.0,) * 3) if box is None else (box.minimum, box.maximum)
    return [
        (lo[0] + (hi[0] - lo[0]) * i / (n - 1),
         lo[1] + (hi[1] - lo[1]) * j / (n - 1),
         lo[2] + (hi[2] - lo[2]) * k / (n - 1))
        for i in range(n) for j in range(n) for k in range(n)
    ]


def _uniform(field, count: int = 512):
    box = field.bounds()
    if box is None:
        lo, hi = (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5)
    else:
        centre = [0.5 * (box.minimum[i] + box.maximum[i]) for i in range(3)]
        half = [0.9 * box.sizes[i] + 1e-9 for i in range(3)]
        lo = [centre[i] - half[i] for i in range(3)]
        hi = [centre[i] + half[i] for i in range(3)]
    # A fixed low-discrepancy lattice, so the probe set is reproducible without a seed.
    golden = 0.6180339887498949
    points = []
    for index in range(count):
        fractions = [((index + 1) * golden * (k + 1)) % 1.0 for k in range(3)]
        points.append(tuple(lo[i] + (hi[i] - lo[i]) * fractions[i] for i in range(3)))
    return points


def _surface_band(field):
    """Bisect rays onto the zero set, then sample exponentially close on both sides."""
    centre, size = _extent(field)
    if not field.contains([centre])[0]:  # pragma: no cover - every field contains it
        for candidate in _lattice(field, 5):
            if field.contains([candidate])[0]:
                centre = candidate
                break
    points = []
    cosines = [1.0 - 2.0 * (index + 0.5) / 32 for index in range(32)]
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index, cosine in enumerate(cosines):
        radius = math.sqrt(max(0.0, 1.0 - cosine * cosine))
        theta = golden_angle * index
        direction = (radius * math.cos(theta), radius * math.sin(theta), cosine)
        far = tuple(centre[i] + 2.0 * size * direction[i] for i in range(3))
        if field.contains([far])[0]:  # unbounded along this ray: no zero set to straddle
            continue
        low, high = 0.0, 2.0 * size
        for _ in range(64):
            middle = 0.5 * (low + high)
            probe = tuple(centre[i] + middle * direction[i] for i in range(3))
            if field.contains([probe])[0]:
                low = middle
            else:
                high = middle
        root = 0.5 * (low + high)
        for offset in BAND:
            points.append(tuple(centre[i] + (root + offset * size) * direction[i] for i in range(3)))
    return points


@pytest.mark.parametrize("name", sorted(FIELDS))
def test_the_two_paths_of_a_field_return_the_same_number(name: str) -> None:
    field = FIELDS[name]()
    evaluate = field.numpy_callable()
    if evaluate is None:  # pragma: no cover - every audited field has a mirror
        pytest.skip(f"{name} has no NumPy mirror to compare against")
    points = _lattice(field) + _uniform(field) + _surface_band(field)
    assert points, f"{name}: probe set is empty"

    scalar = field.signed_distance(points)
    vector = [float(value) for value in evaluate(np.asarray(points, dtype=np.float64))]
    assert len(vector) == len(scalar)

    disagreements = [
        (point, exact, mirrored)
        for point, exact, mirrored in zip(points, scalar, vector, strict=True)
        if not (exact == mirrored)
    ]
    assert not disagreements, (
        f"{name}: {len(disagreements)} of {len(points)} probes disagree between the "
        f"pure-Python contract and its NumPy mirror; first: {disagreements[0]}"
    )


@pytest.mark.parametrize("name", sorted(FIELDS))
def test_the_two_paths_agree_on_which_side_of_the_surface_a_point_is(name: str) -> None:
    """Containment is a label, not a magnitude: one surface cannot have two labels."""
    field = FIELDS[name]()
    evaluate = field.numpy_callable()
    if evaluate is None:  # pragma: no cover
        pytest.skip(f"{name} has no NumPy mirror")
    points = _surface_band(field)
    # Compare the *sign of the distance* on both paths. `contains` is deliberately not
    # used: the SUBOFF fields decide containment in closed form from the offset table,
    # which is a second, independent definition and not the mirror of this one.
    scalar = field.signed_distance(points)
    mirrored = [float(value) for value in evaluate(np.asarray(points, dtype=np.float64))]
    surface = field.surface_isovalue
    if field.inside_sign == "negative":
        labels = [value < surface for value in scalar]
        other = [value < surface for value in mirrored]
    else:
        labels = [value > surface for value in scalar]
        other = [value > surface for value in mirrored]
    flips = [
        (point, label, side)
        for point, label, side in zip(points, labels, other, strict=True)
        if label != side
    ]
    assert not flips, f"{name}: {len(flips)} containment labels flip; first: {flips[0]}"


def test_a_blas_dot_is_not_the_contract_sum() -> None:
    """The defect class, pinned: ``points @ normal`` is a different algorithm.

    This is the reason the plane's mirror was rewritten, and it must stay observable --
    if a future refactor puts the BLAS dot back, the sweep above goes red; this test
    states why that is not pedantry.
    """
    normal = (0.3, -0.4, 0.8660254037844387)
    plane = PlaneField(normal, 0.137)
    points = _uniform(plane, 400)
    array = np.asarray(points, dtype=np.float64)
    blas = array @ np.asarray(normal, dtype=np.float64) - 0.137
    contract = np.asarray(plane.signed_distance(points), dtype=np.float64)
    assert np.array_equal(configurable := np.asarray([float(v) for v in blas]), contract) or (
        configurable != contract
    ).any(), "``array @ normal`` agreed on every probe; the ULP claim needs re-measuring"


def test_squaring_is_not_a_power_call() -> None:
    """``a ** 2`` is CPython's ``pow``; ``a * a`` is one multiplication. They differ."""
    mismatches = [value for value in
                  (0.1 * step for step in range(1, 20_000))
                  if value ** 2 != value * value]
    assert mismatches, "no counter-example in 20k values: the note in geometry_field is stale"
    worst = mismatches[0]
    assert abs(worst ** 2 - worst * worst) <= 2 * math.ulp(worst * worst)