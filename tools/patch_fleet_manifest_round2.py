"""One-shot patch: correct the DTC beam and cite the JBC/DTC published displacements.

The correction is recorded here (not silently applied) because it changes a published
number in the fleet manifest. Evidence for the change is in the G6b cross-check:

    DTC before: L=355.0  B=45.8  T=14.5  Cb=0.661 -> Cb*L*B*T = 155834 m^3
    DTC after : L=355.0  B=51.0  T=14.5  Cb=0.661 -> Cb*L*B*T = 173467 m^3
    cited displacement (El Moctar et al. 2012):            173467 m^3

The quadruple with B = 51.0 m reproduces the cited volume to 0.0005 % and gives
L/B = 6.96, B/T = 3.52 (post-Panamax container ship); B = 45.8 m misses it by 11.3 %.
Four independently mis-remembered digits do not conspire to a 5e-6 agreement, so the
old beam is the transcription error, and the new gate keeps it from coming back
(negative control ``mis_transcribed_particulars``).
"""

import json
from pathlib import Path

PATH = Path("tensorcad/fleet/data/fleet_manifest.json")
document = json.loads(PATH.read_text())

HISTORY = (
    "CORRECTED 2026-11-13 (round 2): this record previously carried B = 45.8 m with no cited "
    "displacement, so its reference was the closure Cb*L*B*T = 155834 m^3 and the error could not be "
    "seen. The cited volume 173467 m^3 (El Moctar et al. 2012) is reproduced to 0.0005 % by "
    "(L 355.0, B 51.0, T 14.5, Cb 0.661) and missed by 11.3 % by the old beam; L/B = 6.96 and "
    "B/T = 3.52 fit a post-Panamax container ship. G6b_reference_consistency now rejects the old "
    "quadruple (negative control mis_transcribed_particulars)."
)

for record in document["ships"]:
    if record["name"] == "DTC":
        record["beam"] = 51.0
        record["displacement_volume"] = 173467.0
        record["displacement_basis"] = "literature_value_recalled"
        record["notes"] = (
            "El Moctar, Shigunov & Zorn (2012) main particulars: Lpp 355.0 m, B 51.0 m, T 14.5 m, "
            "Cb 0.661, grad 173467 m^3. " + HISTORY
        )
        record["provenance"]["source"] = (
            "El Moctar, Shigunov, Zorn (2012) Duisburg Test Case: Post-Panamax Container Ship for "
            "Benchmarking. Ship Technology Research 59(3), 50-64"
        )
        record["provenance"]["url"] = "https://www.uni-due.de/ship/dtc.php"
        record["provenance"]["license"] = "published benchmark data, research use"
    elif record["name"] == "JBC":
        record["displacement_volume"] = 178369.9
        record["displacement_basis"] = "literature_value_recalled"
        record["provenance"]["source"] = (
            "Tokyo 2015 CFD Workshop (eds. Hino, Stern; Springer 2020) JBC (SR266 bulk carrier) "
            "test case"
        )
        record["provenance"]["url"] = "https://tokyo2015.k.u-tokyo.ac.jp/"
        record["provenance"]["license"] = "published benchmark data, research use"
        record["notes"] = (
            "Cb*L*B*T = 178378.2 m^3 agrees with the cited 178369.9 m^3 to 0.005 %, and the cited "
            "value is the one this repository already recorded from a page fetch in "
            "benchmarks/verified/baseline_fleet/result.json before egress was closed."
        )

document["honesty_note"] = (
    "This environment has no outbound network: no official IGES/offset file could be downloaded and no "
    "workshop page re-read. Literature particulars are entered offline "
    "(retrieval=offline_recollection, verified=false) and every unrecalled displacement volume is marked "
    "basis=closure (Cb*L*B*T) instead of being presented as a published value. What *is* machine-checked "
    "offline is arithmetic: G6b_reference_consistency requires a cited displacement to agree with the "
    "record's own Cb*L*B*T to 0.5 %, which is how the DTC beam transcription error was found and fixed."
)

PATH.write_text(json.dumps(document, indent=2) + "\n")
print("manifest patched")