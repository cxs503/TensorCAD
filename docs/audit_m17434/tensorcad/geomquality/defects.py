"""The CAD defect library: injections, detectors and repair advice.

A gate is only worth its false-alarm rate if it is shown a diet of the things it is
meant to catch.  Every entry in :data:`DEFECT_CASES` is a *construction*: it takes a
clean reference hull and damages it in one specific way, then the detection path must
name that damage.  The classic hull-model failures are all here -- an open shell, an
unsewn seam, a non-manifold flap, an inside-out face, a doubled face, a sliver, a
micro edge, a dent, a ripple, a folded patch, a unit-scale blunder, a block-coefficient
mismatch and a knuckle at the deck seam.

The detectors are the audit itself, split by what they read:

* topology (:func:`tensorcad.geomquality.brep.audit_topology`) -- free edges, cracks,
  non-manifold edges, orientation, duplicates, slivers, micro edges, unit scale;
* differential geometry (:mod:`tensorcad.geomquality.continuity`) -- seam class
  (G0/G1/G2), plus the fairing indicators :func:`control_net_ripple` (high-frequency
  ripple in the control net) and :func:`jacobian_flips` (a patch that folds through
  itself);
* hydrostatics -- realised ``Cb`` against the declared value.

Each case declares the code it must produce, so :func:`run_library` returns a
detection matrix that can be read as "case -> found / missed", not a score.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from ..nurbs import NURBSSurface
from .brep import (
    REPAIR_HINTS,
    BRep,
    Defect,
    Edge,
    Face,
    audit_topology,
    merge_breps,
)
from .continuity import audit_seams, curvature_distribution, make_seam_case
from .hydro import patch_hydrostatics

Vec3 = tuple[float, float, float]

#: Relative curvature jump above which a seam is called curvature-discontinuous.
CURVATURE_JUMP_LIMIT = 0.01
#: Normal-angle limit (degrees) for a tangent-continuous seam.
NORMAL_ANGLE_LIMIT_DEG = 0.05
#: Alternating share of the control-net second differences above which a face is ripped.
RIPPLE_LIMIT = 0.35


# --------------------------------------------------------------------------- #
# small construction helpers
# --------------------------------------------------------------------------- #
def _clone(brep: BRep, faces: Sequence[Face] | None = None, source: str = "") -> BRep:
    return BRep(
        faces=tuple(faces if faces is not None else brep.faces),
        edges=brep.edges,
        source=source or brep.source,
        metadata=dict(brep.metadata),
    )


def _remap(brep: BRep, surfaces: Sequence[NURBSSurface], source: str) -> BRep:
    faces = [
        Face(
            index=face.index,
            surface=surfaces[i],
            boundary=face.boundary,
            same_sense=face.same_sense,
            name=face.name,
            wetted=face.wetted,
        )
        for i, face in enumerate(brep.faces)
    ]
    return _clone(brep, faces, source)


def scale_brep(brep: BRep, factor: float, *, axes: Sequence[int] = (0, 1, 2), source: str = "scaled") -> BRep:
    """Scale the control nets (and the edge polylines) about the origin."""
    def transform(point: Sequence[float]) -> Vec3:
        return tuple(point[axis] * factor if axis in axes else point[axis] for axis in range(3))

    surfaces = []
    for face in brep.faces:
        surfaces.append(
            NURBSSurface(
                degree=face.surface.degree,
                control_net=[[transform(point) for point in row] for row in face.surface.control_net],
                weights=face.surface.weights,
                knots=face.surface.knots,
            )
        )
    scaled = _remap(brep, surfaces, source)
    edges = tuple(
        Edge(index=edge.index, points=tuple(transform(point) for point in edge.points))
        for edge in brep.edges
    )
    return _clone(scaled, scaled.faces, source) if False else BRep(
        faces=scaled.faces, edges=edges, source=source, metadata=dict(brep.metadata)
    )


def translate_brep(brep: BRep, offset: Sequence[float], *, only: str | None = None) -> BRep:
    """Translate every face (or only faces whose name contains ``only``)."""
    def transform(point: Sequence[float]) -> Vec3:
        return (
            point[0] + offset[0],
            point[1] + offset[1],
            point[2] + offset[2],
        )

    surfaces = []
    for face in brep.faces:
        if only is None or only in face.name:
            surfaces.append(
                NURBSSurface(
                    degree=face.surface.degree,
                    control_net=[[transform(point) for point in row] for row in face.surface.control_net],
                    weights=face.surface.weights,
                    knots=face.surface.knots,
                )
            )
        else:
            surfaces.append(face.surface)
    return _remap(brep, surfaces, "translated")


def perturb_control_point(brep: BRep, *, face_name: str, i: int, j: int, delta: Sequence[float]) -> BRep:
    surfaces = []
    for face in brep.faces:
        if face_name in face.name:
            net = [list(row) for row in face.surface.control_net]
            point = net[i][j]
            net[i][j] = (point[0] + delta[0], point[1] + delta[1], point[2] + delta[2])
            surfaces.append(
                NURBSSurface(
                    degree=face.surface.degree,
                    control_net=[tuple(row) for row in net],
                    weights=face.surface.weights,
                    knots=face.surface.knots,
                )
            )
        else:
            surfaces.append(face.surface)
    return _remap(brep, surfaces, "perturbed")


def ripple_control_net(brep: BRep, *, face_name: str, amplitude: float, period: int = 1) -> BRep:
    """Alternating perturbation of the middle control rows: a fairing failure."""
    surfaces = []
    for face in brep.faces:
        if face_name not in face.name:
            surfaces.append(face.surface)
            continue
        net = [list(row) for row in face.surface.control_net]
        for i in range(1, len(net) - 1):
            sign = (-1.0) ** (i // period)
            net[i] = [
                (point[0], point[1] + sign * amplitude, point[2]) for point in net[i]
            ]
        surfaces.append(
            NURBSSurface(
                degree=face.surface.degree,
                control_net=[tuple(row) for row in net],
                weights=face.surface.weights,
                knots=face.surface.knots,
            )
        )
    return _remap(brep, surfaces, "rippled")


def add_flap(brep: BRep, *, face_name: str, factor: float = 0.999) -> BRep:
    """Attach a nearly-identical face on the same boundary: a non-manifold flap."""
    original = [face for face in brep.faces if face_name in face.name][0]
    net = [
        [(point[0], point[1] * factor, point[2] * factor) for point in row]
        for row in original.surface.control_net
    ]
    flap = Face(
        index=len(brep.faces),
        surface=NURBSSurface(
            degree=original.surface.degree,
            control_net=[tuple(row) for row in net],
            weights=original.surface.weights,
            knots=original.surface.knots,
        ),
        boundary=original.boundary,
        same_sense=not original.same_sense,
        name=original.name + "_flap",
        wetted=False,
    )
    return _clone(brep, list(brep.faces) + [flap], "flap")


def add_sliver(brep: BRep, *, size: float) -> BRep:
    """Append a free-floating tiny patch: a sliver face with micro edges."""
    lo, hi = brep.bounding_box()
    corner = (hi[0], hi[1], hi[2])
    net = []
    for i in range(2):
        row = []
        for j in range(2):
            row.append(
                (
                    corner[0] - size * i,
                    corner[1] - size * j,
                    corner[2],
                )
            )
        net.append(row)
    surface = NURBSSurface(
        degree=(1, 1),
        control_net=[tuple(row) for row in net],
        weights=[[1.0, 1.0], [1.0, 1.0]],
        knots=[[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0]],
    )
    from .brep import sew

    sliver = sew([surface], names=["sliver"], source="sliver")
    merged = merge_breps(brep, sliver, source="sliver")
    return merged


def add_needle(brep: BRep, *, width: float) -> BRep:
    """Append a needle patch whose short edges are micro edges (tolerance violation)."""
    lo, hi = brep.bounding_box()
    base = (hi[0], 0.0, hi[2])
    net = [
        [(base[0] - width, base[1], base[2]), (base[0] - width, base[1] + 0.05, base[2])],
        [(base[0], base[1], base[2]), (base[0], base[1] + 0.05, base[2])],
    ]
    surface = NURBSSurface(
        degree=(1, 1),
        control_net=[tuple(row) for row in net],
        weights=[[1.0, 1.0], [1.0, 1.0]],
        knots=[[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0]],
    )
    from .brep import sew

    needle = sew([surface], names=["needle"], source="needle")
    return merge_breps(brep, needle, source="needle")


def tilt_deck(brep: BRep, *, degrees: float = 5.0) -> BRep:
    """Rotate the deck patch about x so the waterline seam becomes a knuckle."""
    angle = math.radians(degrees)
    cosine, sine = math.cos(angle), math.sin(angle)

    def transform(point: Sequence[float]) -> Vec3:
        return (point[0], point[1] * cosine - point[2] * sine, point[1] * sine + point[2] * cosine)

    surfaces = []
    for face in brep.faces:
        if "deck" in face.name:
            surfaces.append(
                NURBSSurface(
                    degree=face.surface.degree,
                    control_net=[[transform(point) for point in row] for row in face.surface.control_net],
                    weights=face.surface.weights,
                    knots=face.surface.knots,
                )
            )
        else:
            surfaces.append(face.surface)
    return _remap(brep, surfaces, "tilted-deck")


def add_seam_pair(brep: BRep, *, kind: str = "G1", fraction: float = 0.1) -> BRep:
    """Splice a two-patch tab with a *known* seam class onto the model."""
    import random

    from .brep import sew

    rng = random.Random(20261107)
    first, second, _ = make_seam_case(kind, rng, scale=fraction * brep.length)
    pair = sew([first, second], names=["tab_a", "tab_b"], source="tab")
    # Move the tab clear of the hull so it only contributes its seam.
    lo, hi = brep.bounding_box()
    from .brep import Edge

    offset = (0.0, hi[1] + 2.0 * fraction * brep.length, hi[2])
    surfaces = []
    for face in pair.faces:
        surfaces.append(
            NURBSSurface(
                degree=face.surface.degree,
                control_net=[
                    [
                        (point[0] + offset[0], point[1] + offset[1], point[2] + offset[2])
                        for point in row
                    ]
                    for row in face.surface.control_net
                ],
                weights=face.surface.weights,
                knots=face.surface.knots,
            )
        )
    moved = _remap(pair, surfaces, "tab")
    edges = tuple(
        Edge(
            index=edge.index,
            points=tuple(
                (point[0] + offset[0], point[1] + offset[1], point[2] + offset[2])
                for point in edge.points
            ),
        )
        for edge in moved.edges
    )
    moved = BRep(faces=moved.faces, edges=edges, source="tab")
    return merge_breps(brep, moved, source="tab")


def fold_patch(brep: BRep, *, face_name: str, delta: float) -> BRep:
    """Push a keel control point up through the patch: the surface doubles back.

    A displacement in ``y`` is only a dent (the surface stays a graph over ``x, z``); a
    displacement in ``z`` that carries one control row past its neighbours makes the
    parameterisation fold, which is what a self-intersecting patch looks like.
    """
    face = [item for item in brep.faces if face_name in item.name][0]
    middle_i = len(face.surface.control_net) // 2
    return perturb_control_point(
        brep,
        face_name=face_name,
        i=middle_i,
        j=0,
        delta=(0.0, 0.0, delta),
    )


# --------------------------------------------------------------------------- #
# detectors beyond the topology audit
# --------------------------------------------------------------------------- #
def control_net_ripple(face: Face) -> float:
    """Alternating share of the control-net second differences in the fairing direction.

    The signed second difference ``y[i-1] - 2y[i] + y[i+1]`` of a smooth net varies
    slowly from row to row, so ``sum(-1)**i * d_i`` cancels; a corrugated net (a fairing
    failure with row-to-row alternation) makes the alternating sum as large as the total.
    The result is scale-free, so it does not depend on hull size or point spacing.  Nets
    with fewer than six rows are not measured: there are not enough rows for the test to
    mean anything (and a 4-row patch has a real second difference by construction).
    """
    net = face.surface.control_net
    if len(net) < 6:
        return 0.0
    differences: list[float] = []
    for i in range(1, len(net) - 1):
        differences.append(
            sum(net[i - 1][j][1] - 2.0 * net[i][j][1] + net[i + 1][j][1] for j in range(len(net[0])))
        )
    total = sum(abs(value) for value in differences)
    if total == 0.0:
        return 0.0
    alternating = abs(
        sum((1.0 if index % 2 == 0 else -1.0) * value for index, value in enumerate(differences))
    )
    return alternating / total


def jacobian_flips(face: Face, *, samples: int = 9) -> float:
    """Fraction of the patch where the normal points against the patch-average one.

    A fold reverses the surface normal inside the domain, so sampling a grid and asking
    "how often does the normal disagree with the area-weighted average" catches a fold
    wherever it happens, not only on the centre line.
    """
    (u0, u1), (v0, v1) = face.domain
    normals: list[tuple[float, tuple[float, float, float], float]] = []
    for i in range(samples):
        u = u0 + (u1 - u0) * (i + 0.5) / samples
        for j in range(samples):
            v = v0 + (v1 - v0) * (j + 0.5) / samples
            try:
                _, _, cross = face.surface.jacobian(u, v)
            except Exception:
                continue
            length = math.sqrt(sum(component**2 for component in cross))
            if length == 0.0:
                continue
            normals.append((length, cross, u))
    if not normals:
        return 0.0
    reference = (0.0, 0.0, 0.0)
    for length, cross, _ in normals:
        reference = (
            reference[0] + length * cross[0] / length,
            reference[1] + length * cross[1] / length,
            reference[2] + length * cross[2] / length,
        )
    disagree = 0
    for length, cross, _ in normals:
        dot = sum(a * b for a, b in zip(cross, reference)) / length
        if dot < 0.0:
            disagree += 1
    return disagree / len(normals)


# --------------------------------------------------------------------------- #
# the joint detector
# --------------------------------------------------------------------------- #
def detect_all(
    brep: BRep,
    *,
    declared: dict[str, float] | None = None,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Run every detector and return the codes, the details and the hydrostatic table."""
    declared = dict(declared or {})
    length = declared.get("length")
    hydro = patch_hydrostatics(
        brep,
        length=declared.get("length"),
        beam=declared.get("beam"),
        draft=declared.get("draft"),
    )
    topology = audit_topology(
        brep,
        tolerance=tolerance,
        expected_length=length,
        declared_cb=declared.get("block_coefficient"),
        reference_cb=hydro.block_coefficient,
    )
    defects: list[Defect] = list(topology.defects)
    codes = {defect.code for defect in defects}

    seams = audit_seams(
        brep,
        angle_limit_deg=NORMAL_ANGLE_LIMIT_DEG,
        curvature_limit=CURVATURE_JUMP_LIMIT,
        samples=9,
    )
    for seam in seams:
        if seam.classification == "G0":
            defects.append(
                Defect(
                    code="TANGENT_JUMP",
                    message=(
                        f"seam {seam.face_a}-{seam.face_b} has a normal angle of "
                        f"{seam.max_normal_angle_deg:.3g} deg"
                    ),
                    magnitude=seam.max_normal_angle_deg,
                    relative=seam.max_normal_angle_deg,
                )
            )
        elif seam.classification == "G1":
            defects.append(
                Defect(
                    code="CURVATURE_JUMP",
                    message=(
                        f"seam {seam.face_a}-{seam.face_b} is tangent-continuous but jumps "
                        f"{100 * seam.max_curvature_jump:.1f} % in curvature"
                    ),
                    magnitude=seam.max_curvature_jump,
                    relative=seam.max_curvature_jump,
                )
            )

    for face in brep.faces:
        if not face.wetted:
            continue
        ratio = control_net_ripple(face)
        if ratio > RIPPLE_LIMIT:
            defects.append(
                Defect(
                    code="ROUGHNESS",
                    message=(
                        f"face '{face.name}' control net is {100 * ratio:.0f} % alternating: "
                        "a corrugation, not a hull form"
                    ),
                    magnitude=ratio,
                    relative=ratio,
                )
            )
        flips = jacobian_flips(face)
        if flips > 0.02:
            defects.append(
                Defect(
                    code="SELF_INTERSECTION",
                    message=(
                        f"face '{face.name}' normal reverses over {100 * flips:.0f} % of its "
                        "domain: the patch folds through itself"
                    ),
                    magnitude=flips,
                    relative=flips,
                )
            )
    codes = {defect.code for defect in defects}
    return {
        "codes": sorted(codes),
        "defects": [defect.to_dict() for defect in defects],
        "watertight": topology.watertight,
        "hydrostatics": hydro.to_dict(),
        "seams": [seam.to_dict() for seam in seams],
        "curvature": {
            face.name: curvature_distribution(face.surface, samples=5)
            for face in brep.faces
            if face.wetted
        },
        "topology": topology.to_dict(),
    }


# --------------------------------------------------------------------------- #
# the library
# --------------------------------------------------------------------------- #
@dataclass
class DefectCase:
    """One damaged model: how to build it, what it breaks, what must be detected."""

    code: str
    name: str
    expectation: str
    inject: Callable[[BRep], BRep]

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "name": self.name,
            "expectation": self.expectation,
            "repair": REPAIR_HINTS.get(self.code, ""),
        }


def _unit_scale(brep: BRep) -> BRep:
    return scale_brep(brep, 1000.0, source="unit-scale")


def _cb_offset(brep: BRep) -> BRep:
    return scale_brep(brep, 0.95, axes=(1,), source="cb-offset")


def _hole(brep: BRep) -> BRep:
    faces = [face for face in brep.faces if face.name.strip() != "hull_half"]
    return _clone(brep, faces, "hole")


def _open_top(brep: BRep) -> BRep:
    faces = [face for face in brep.faces if "deck" not in face.name]
    return _clone(brep, faces, "open-top")


def _seam_gap_offset(brep: BRep, *, delta: float) -> BRep:
    """Offset the starboard hull patch by ``delta`` and re-sew the whole model.

    Offsetting the surface alone changes nothing the adjacency graph can see, because the
    edge pool is shared; a crack appears when the topology is rebuilt from the geometry,
    which is exactly what a re-export or a repair pass does.
    """
    from .brep import sew

    surfaces = []
    for face in brep.faces:
        if face.name.strip() == "hull_half":
            surfaces.append(
                NURBSSurface(
                    degree=face.surface.degree,
                    control_net=[
                        [(point[0], point[1] + delta, point[2]) for point in row]
                        for row in face.surface.control_net
                    ],
                    weights=face.surface.weights,
                    knots=face.surface.knots,
                )
            )
        else:
            surfaces.append(face.surface)
    return sew(
        surfaces,
        names=[face.name for face in brep.faces],
        same_sense=[face.same_sense for face in brep.faces],
        wetted=[face.wetted for face in brep.faces],
        source=f"re-sewn with a {delta:g} offset",
    )


def _seam_gap(brep: BRep) -> BRep:
    return _seam_gap_offset(brep, delta=5e-6 * brep.length)


def _orientation(brep: BRep) -> BRep:
    faces = [
        Face(
            index=face.index,
            surface=face.surface,
            boundary=face.boundary,
            same_sense=not face.same_sense if "deck" in face.name else face.same_sense,
            name=face.name,
            wetted=face.wetted,
        )
        for face in brep.faces
    ]
    return _clone(brep, faces, "flipped")


def _duplicate(brep: BRep) -> BRep:
    original = [face for face in brep.faces if face.name.strip() == "hull_half"][0]
    copy = Face(
        index=len(brep.faces),
        surface=original.surface,
        boundary=original.boundary,
        same_sense=original.same_sense,
        name=original.name + "_copy",
        wetted=original.wetted,
    )
    return _clone(brep, list(brep.faces) + [copy], "duplicate")


DEFECT_CASES: tuple[DefectCase, ...] = (
    DefectCase("FREE_EDGE", "hole", "a face is missing: the shell is open", _hole),
    DefectCase(
        "OPEN_TOP", "open_waterline", "the deck is missing so the waterplane is open", _open_top
    ),
    DefectCase(
        "SEAM_GAP",
        "unsewn_seam",
        "one half is translated by 1 mm: the seam no longer welds",
        _seam_gap,
    ),
    DefectCase(
        "NON_MANIFOLD",
        "non_manifold_flap",
        "a flap face shares the whole deck boundary",
        lambda brep: add_flap(brep, face_name="deck"),
    ),
    DefectCase("ORIENTATION", "flipped_face", "a face is wound inside-out", _orientation),
    DefectCase("DUPLICATE_FACE", "duplicate_face", "a face is present twice", _duplicate),
    DefectCase(
        "SLIVER",
        "sliver_face",
        "a 2e-5 L patch is attached: near-zero area",
        lambda brep: add_sliver(brep, size=2e-5 * brep.length),
    ),
    DefectCase(
        "MICRO_EDGE",
        "micro_edge",
        "a needle patch carries a 1e-7 m edge, below the 1e-6 m weld tolerance",
        lambda brep: add_needle(brep, width=1e-7),
    ),
    DefectCase(
        "CURVATURE_JUMP",
        "g1_splice",
        "two patches meet with a curvature jump (tangent-continuous, G1)",
        lambda brep: add_seam_pair(brep, kind="G1"),
    ),
    DefectCase(
        "ROUGHNESS",
        "ripple",
        "the control net alternates by 5e-4 L: a corrugated hull",
        lambda brep: ripple_control_net(
            brep, face_name="hull_half", amplitude=5e-4 * brep.length
        ),
    ),
    DefectCase(
        "SELF_INTERSECTION",
        "fold",
        "a control point crosses the centreline: the patch folds",
        lambda brep: fold_patch(
            brep, face_name="hull_half", delta=3.0 * (brep.bounding_box()[1][2] - brep.bounding_box()[0][2])
        ),
    ),
    DefectCase("UNIT_SCALE", "unit_scale", "the model is 1000x too large", _unit_scale),
    DefectCase(
        "CB_MISMATCH", "cb_mismatch", "the beam is 5 % narrow: Cb no longer matches", _cb_offset
    ),
    DefectCase(
        "TANGENT_JUMP",
        "g0_splice",
        "two patches meet with a tangent jump (a knuckle that is not a design feature)",
        lambda brep: add_seam_pair(brep, kind="G0"),
    ),
)


def run_library(
    reference: BRep,
    *,
    declared: dict[str, float],
    clean_codes: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Inject every library case into ``reference`` and check the detector names it.

    The clean model is audited too: a detector that fires on the reference hull is a
    false alarm, and the matrix records that alongside the hits.
    """
    baseline = detect_all(reference, declared=declared)
    rows = []
    for case in DEFECT_CASES:
        damaged = case.inject(reference)
        result = detect_all(damaged, declared=declared)
        found = case.code in result["codes"]
        rows.append(
            {
                "case": case.name,
                "expected_code": case.code,
                "detected_codes": result["codes"],
                "detected": found,
                "repair": REPAIR_HINTS.get(case.code, ""),
                "watertight": result["watertight"],
            }
        )
    hits = sum(1 for row in rows if row["detected"])
    return {
        "cases": len(rows),
        "detected": hits,
        "detection_rate": hits / len(rows) if rows else 0.0,
        "false_alarm_codes": [
            code for code in baseline["codes"] if code not in set(clean_codes or ())
        ],
        "baseline_codes": baseline["codes"],
        "baseline_watertight": baseline["watertight"],
        "rows": rows,
    }
