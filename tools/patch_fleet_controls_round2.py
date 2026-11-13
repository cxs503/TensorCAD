"""One-shot patch: add the mis-transcribed-particulars negative control (the DTC defect)."""

from pathlib import Path

CONTROLS = Path("tensorcad/fleet/controls.py")
text = CONTROLS.read_text()

if "mis_transcribed_beam" not in text:
    text = text.replace(
        '"""\n\nfrom __future__ import annotations',
        "* :func:`mis_transcribed_beam` — the defect the round-2 cross-check actually found: a\n"
        "  published displacement quoted next to a beam that was transcribed wrongly (DTC carried\n"
        "  B = 45.8 m against a cited 173467 m^3, 11.3 % off its own Cb*L*B*T), which must fail\n"
        "  ``G6b_reference_consistency`` before any geometry is built.\n"
        '"""\n\nfrom __future__ import annotations',
        1,
    )
    text += '''

def mis_transcribed_beam(base: ShipRecord, *, beam: float = 45.8) -> ShipRecord:
    """The same record with one particular transcribed wrongly, volume left as cited.

    This is not a synthetic trap: it is the DTC entry as it stood before the round-2
    cross-check (B = 45.8 m next to a cited 173467 m^3), and it is the reason the gate
    exists. A record like this passed every pre-existing gate because it had no cited
    displacement and its reference was the closure, i.e. it was checked against itself.
    """
    return _clone(
        base,
        name=f"{base.name}_mis_transcribed",
        beam=beam,
        notes=base.notes + f" [negative control: beam mis-transcribed to {beam} m, volume kept as cited]",
        provenance=Provenance(
            source="negative control derived from " + base.provenance.source,
            url=base.provenance.url,
            license=base.provenance.license,
            retrieval="derived_control",
            verified=False,
        ),
    )
'''
CONTROLS.write_text(text)

RUN = Path("tensorcad/fleet/run.py")
text = RUN.read_text()
if "mis_transcribed_particulars" not in text:
    text = text.replace(
        "    placeholder_box,\n    wrong_volume_record,",
        "    mis_transcribed_beam,\n    placeholder_box,\n    wrong_volume_record,",
        1,
    )
    control = '''
    # 5. a cited displacement quoted next to a mis-transcribed particular: the DTC defect
    dtc = next((ship for ship in manifest.ships if ship.name == "DTC"), base)
    mistyped = mis_transcribed_beam(dtc, beam=dtc.beam * 45.8 / 51.0)
    gate = gate_reference_consistency(mistyped)
    controls.append(
        {
            "control": "mis_transcribed_particulars",
            "expected": "rejected",
            "rejected": not gate.passed,
            "failed_gates": [] if gate.passed else [gate.name],
            "detail": gate.detail,
            "cited_volume": gate.values["cited_volume"],
            "closure_volume": gate.values["closure_volume"],
            "relative_error": gate.values["relative_error"],
        }
    )
    return controls
'''
    text = text.replace("    return controls\n", control, 1)
RUN.write_text(text)
print("controls patched:", "mis_transcribed_beam" in CONTROLS.read_text(), "mis_transcribed_particulars" in RUN.read_text())