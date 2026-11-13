"""One-shot patch: add the G6b reference-triple consistency gate to the fleet gate set."""

from pathlib import Path

TARGET = Path("tensorcad/fleet/gates.py")
text = TARGET.read_text()

if "G6b_reference_consistency" not in text:
    text = text.replace(
        "#: Internal consistency of a declared (Cb, L, B, T, volume) set.\nPARTICULARS_TOLERANCE = 0.03",
        "#: Internal consistency of a declared (Cb, L, B, T, volume) set.\nPARTICULARS_TOLERANCE = 0.03\n"
        "#: Agreement required between a *cited* displacement volume and Cb*L*B*T. Seven\n"
        "#: significant digits are carried by the published sets, so a transcription error in any\n"
        "#: of L, B, T, Cb or the volume shows up far above this band, while round-off of the\n"
        "#: published particulars stays below it.\nREFERENCE_CONSISTENCY_TOLERANCE = 0.005",
        1,
    )

    GATE = '''
def gate_reference_consistency(
    record: ShipRecord, tolerance: float = REFERENCE_CONSISTENCY_TOLERANCE
) -> GateResult:
    """A cited displacement volume must agree with the record's own Cb*L*B*T.

    This is the gate that catches a *transcription* error rather than a bad mesh: if a
    published displacement is quoted next to published L, B, T and Cb, then the four
    numbers over-determine each other, and any one of them entered wrongly makes the
    set inconsistent. A closure reference is not a check at all (it is defined as
    Cb*L*B*T), and the gate says so instead of reporting a vacuous pass.
    """
    closure = record.closed_form_volume()
    if record.displacement_volume is None:
        return GateResult(
            "G6b_reference_consistency",
            True,
            f"reference is the closure Cb*L*B*T = {closure:.6g} m^3: no independent displacement is "
            f"cited, so this gate cannot check anything (G7 then only proves the mesh matches the "
            f"record's own particulars)",
            {
                "cited_volume": None,
                "closure_volume": closure,
                "relative_error": None,
                "tolerance": tolerance,
                "vacuous": True,
            },
        )
    cited = float(record.displacement_volume)
    rel = (cited - closure) / cited
    return GateResult(
        "G6b_reference_consistency",
        abs(rel) <= tolerance,
        f"cited V = {cited:.6g} m^3 vs Cb*L*B*T = {closure:.6g} m^3: {rel * 100:+.4f}% "
        f"(tolerance {tolerance * 100:g}%); a disagreement beyond this band means L, B, T, Cb or the "
        f"cited volume was mis-transcribed",
        {
            "cited_volume": cited,
            "closure_volume": closure,
            "relative_error": rel,
            "tolerance": tolerance,
            "vacuous": False,
            "basis": record.displacement_basis,
        },
    )


'''
    text = text.replace("def gate_volume(", GATE.lstrip("\n") + "def gate_volume(", 1)
    text = text.replace(
        "        gate_particulars_consistent(record),\n        gate_volume(",
        "        gate_particulars_consistent(record),\n"
        "        gate_reference_consistency(record),\n"
        "        gate_volume(",
        1,
    )

TARGET.write_text(text)
print("gates patched:", "G6b_reference_consistency" in TARGET.read_text())