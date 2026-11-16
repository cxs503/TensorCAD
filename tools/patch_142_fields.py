"""One-shot patcher: replace ``_curves_of_142`` with a field-positional reader.

Kept in the tree so the change can be re-derived from the archived function text, and so
the before/after can be reproduced without hand editing (see
``docs/gap_cad_t1_iges_pipe_round5.md``).  Run from the repository root:

    python tools/patch_142_fields.py

It is idempotent: if the new reader is already present it does nothing.
"""
from __future__ import annotations

import sys
from pathlib import Path

TARGET = Path(__file__).resolve().parents[1] / "tensorcad/geomquality/iges_trim.py"

HEAD = '''def _curves_of_142(builder: _Builder, entry: DirectoryEntry,
                   surface: Surface | None = None) -> list[BoundaryLoop]:
'''

NEW = '''#: Type-142 parameter fields, IGES 5.3 section 4.25, in file order.
_142_FIELDS = ("CRV", "SPTR", "BPTR", "CPTR", "IT")
#: Entity types ``_Builder.surface`` can build from.
_SURFACE_ENTITY_TYPES = (108, 120, 122, 128)
#: Entity types ``_Builder.curve`` can build from (composites included).
_CURVE_ENTITY_TYPES = (100, 102, 106, 110, 126)


def _fields_of_142(builder: _Builder, entry: DirectoryEntry) -> tuple[int, dict[str, int]]:
    """A type-142 record read **by field position**, as IGES 5.3 section 4.25 defines it.

    ``CRV`` (model-space curve), ``SPTR`` (surface), ``BPTR`` (the same curve in the
    parameter space of that surface), ``CPTR`` (a composite carrying the preferred
    representation) and finally ``IT`` -- a 0/1 flag selecting which of the two is
    preferred.  ``IT`` is a **flag, not a pointer**.

    An earlier revision scanned every integer token of the record and resolved each one
    that happened to name an existing curve, so a flag of ``1`` became a reference to
    directory entry 1 -- on the KVLCC2 archive a type-126 curve belonging to another
    face, 28.4 m away from the surface being trimmed.  Every trimmed patch then acquired
    a phantom loop, its genuine parameter-space image was vetoed against that phantom,
    and the face was refused after minutes of futile inversion.  Two independent writers
    agree on the layout and both pad the record with trailing zeros that are not
    pointers either::

        KVLCC2  DE 35: 142, 0, 5, 33, 19, 1, 0, 0
        KP505   DE 13: 142, 0, 1,  0, 11, 2, 0, 0

    ``parse_iges`` strips a leading entity-type token when the writer repeats it, so
    ``CRV`` is normally the first remaining parameter; a record that keeps the type is
    one field longer.  ``SPTR`` decides the alignment, because it is the only field whose
    role is unmistakable -- it has to name a surface.
    """
    document = builder.document

    def is_surface(position: int) -> bool:
        try:
            value = entry.integer(position)
        except IgesError:
            return False
        target = document.entries.get(value) if value > 0 else None
        return target is not None and target.entity_type in _SURFACE_ENTITY_TYPES

    base = 1 if (is_surface(2) and not is_surface(1)) else 0
    fields: dict[str, int] = {}
    for offset, name in enumerate(_142_FIELDS):
        position = base + offset
        value = 0
        if position < len(entry.tokens):
            try:
                value = entry.integer(position)
            except IgesError:
                value = 0
        fields[name] = value if value > 0 else 0
    return base, fields


def _curves_of_142(builder: _Builder, entry: DirectoryEntry,
                   surface: Surface | None = None) -> list[BoundaryLoop]:
    """Resolve a curve-on-surface record (142) into the loops that trim ``surface``.

    The record is read by field position (:func:`_fields_of_142`) and the fields carry
    their spec meaning: ``BPTR`` **is** the curve in the parameter space of ``SPTR``, so
    when it exists and lies inside the surface's knot domain it is consumed as an exact
    parametric loop -- no inversion, no branch ambiguity, no minutes-long walk.

    The other representation (``CRV`` in model space, or the 3-D children of the
    ``CPTR`` composite) is used as a *cross-check*, and a disagreement is **reported,
    not used as a veto**: on the KVLCC2 archive that 3-D composite is not the same curve
    (measured disagreement 1.33e1 model units on DE 35, 1.41e1 on DE 107), so vetoing the
    parameter-space image against it discarded the only usable trim and left the face
    with two nonsense loops that the mesh stage then refused.  When the two
    representations *do* agree, the 3-D copy is dropped exactly as before, with the same
    refusal text and the measured agreement.

    Records without a usable ``BPTR`` -- ``KP505`` writes ``BPTR = 0`` and puts both
    representations in the ``CPTR`` composite -- are classified by measurement, as
    before: a reference whose coordinates fall inside the domain and that reproduces
    another reference's 3-D curve to within the chord tolerance is the parameter-space
    image and is used; a reference inside the domain that disagrees is refused; and a
    record whose only reference sits inside the domain keeps the 3-D route, because that
    case cannot be classified without an inversion probe.
    """
    document = builder.document
    base, fields = _fields_of_142(builder, entry)

    def curve_like(value: int) -> bool:
        target = document.entries.get(value) if value > 0 else None
        return target is not None and target.entity_type in _CURVE_ENTITY_TYPES

    references: list[int] = []
    for name in ("CRV", "BPTR", "CPTR"):
        value = fields[name]
        if curve_like(value) and value not in references:
            references.append(value)
    if not references:
        raise IgesError(
            f"curve-on-surface entity {entry.sequence} does not reference any curve this "
            f"reader can evaluate (fields {dict(fields)}, tokens: {list(entry.tokens)})"
        )

    # ---- the spec-declared parameter-space curve wins when it lies in the domain ----
    bptr = fields["BPTR"] if curve_like(fields["BPTR"]) else 0
    if surface is not None and bptr:
        loop = None
        try:
            loop = _loop_from_curve(builder, bptr)
        except IgesError as exc:
            builder.refused[bptr] = str(exc)
        if loop is not None and len(loop) >= 2 and _uv_in_domain(surface, loop):
            uv = [(point[0], point[1]) for point in loop]
            limit = max(2.0 * builder.tolerance, 1e-12)
            cross_check: float | None = None
            for other in [value for value in references if value != bptr]:
                try:
                    candidate = _loop_from_curve(builder, other)
                except IgesError as exc:
                    builder.refused[other] = str(exc)
                    continue
                if len(candidate) < 2:
                    continue
                deviation = _uv_deviation(surface, uv, candidate)
                cross_check = deviation if cross_check is None else min(cross_check, deviation)
                if deviation <= limit:
                    builder.refused[other] = (
                        f"3-D curve {other} of 142 entity {entry.sequence} is the same "
                        f"trim as its parameter-space image {bptr} "
                        f"(point-to-polyline agreement {deviation:.3e} <= {limit:.3e} "
                        f"model units); the image is used because it is exact"
                    )
                else:
                    builder.refused[other] = (
                        f"3-D curve {other} of 142 entity {entry.sequence} is not the same "
                        f"curve as the parameter-space image {bptr} declared by the "
                        f"record's BPTR field (point-to-polyline disagreement "
                        f"{deviation:.3e} model units > {limit:.3e}); the parameter-space "
                        f"curve is the trim, the 3-D curve is kept as a cross-check only"
                    )
            closed, gap = _close_uv(surface, uv)
            if len(closed) >= 3:
                image = tuple(surface.evaluate(u, v) for u, v in closed)
                builder.consumed.add(entry.sequence)
                return [BoundaryLoop(points=image, uv=closed, route="parametric",
                                     residual=0.0, cross_check=cross_check, source=bptr)]
            builder.refused[bptr] = (
                f"parameter-space curve {bptr} of 142 entity {entry.sequence} collapses to "
                f"{len(closed)} parameters (closing gap {gap:.3e} model units)"
            )

    # ---- otherwise: classify the record's references by measurement ----
    sampled: list[tuple[int, tuple[Vec3, ...]]] = []
    for reference in references:
        try:
            loop = _loop_from_curve(builder, reference)
        except IgesError as exc:
            builder.refused[reference] = str(exc)
            continue
        if len(loop) >= 2:
            sampled.append((reference, loop))
    if not sampled:
        raise IgesError(
            f"curve-on-surface entity {entry.sequence} references only curves that were refused"
        )
    if surface is None:
        return [BoundaryLoop(points=points, source=reference) for reference, points in sampled]
'''


def main() -> int:
    text = TARGET.read_text()
    if "_fields_of_142" in text:
        print("already patched")
        return 0
    start = text.index(HEAD)
    tail_marker = '    gate = 1e-6 * builder.extent\n'
    anchor = text.index(tail_marker, start)
    patched = text[:start] + NEW + "\n" + text[anchor:]
    TARGET.write_text(patched)
    print("patched", TARGET)
    return 0


if __name__ == "__main__":
    sys.exit(main())