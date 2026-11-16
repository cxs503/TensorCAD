"""Regenerate (or verify) every number that ``docs/lines_fairing/lines_fairing.md`` reports.

The document claims "from ``results.json``, not hand-copied"; this script is what makes the
claim checkable.  Three tables are produced -- the plain-energy sweep, the per-family breakdown
at the 0.8 % B/2 budget, and the design-feature comparison -- plus the ``feature_rows.png``
figure.  Run from the repository root::

    python docs/lines_fairing/feature_evidence.py            # rewrite the tables + figure
    python docs/lines_fairing/feature_evidence.py --check     # exit 1 if the document has drifted

``--check`` fails when a generated row is missing from the document or the freshly rendered
figure differs from the committed PNG, so a stale table (or a figure drawn from an older
``results.json``) cannot be merged by accident.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCUMENT = ROOT / "docs/lines_fairing/lines_fairing.md"
RESULTS = ROOT / "docs/lines_fairing/results.json"
FIGURE = ROOT / "docs/lines_fairing/feature_rows.png"

LABELS = {
    "wigley_iii": "Wigley III",
    "series60_cb060": "Series 60 $C_B$=0.60",
    "kvlcc2_class": "KVLCC2-class",
}
COLOURS = {"wigley_iii": "tab:blue", "series60_cb060": "tab:orange", "kvlcc2_class": "tab:green"}


def _sum(flips: dict[str, int]) -> int:
    return sum(flips.values())


def plain_table(cases: list[dict]) -> list[str]:
    lines = [
        "| hull form | flips of the analytic table | flips with jitter | tolerance "
        "| deviation (m) | flips after | accepted | deviation argmax (waterline) "
        "| rms to reference (m) | dVolume |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for case in cases:
        analytic = _sum(case["analytic_flips"])
        jittered = _sum(case["unfair_flips"])
        for entry in case["sweeps"]["off"]:
            reference = entry["reference"]
            volume = entry["hydrostatics"]["relative_change"]["volume"]
            lines.append(
                f"| `{case['form']}` | {analytic} | {jittered} "
                f"| {entry['tolerance_fraction_of_half_beam'] * 100:.1f} % B/2 "
                f"| {entry['deviation_max_m']:.4f} | {entry['flips_after']} "
                f"| {entry['accepted']} | j={entry['deviation_argmax']['waterline_index']} "
                f"| {reference['rms_before_m']:.5f} -> {reference['rms_after_m']:.5f} "
                f"| {volume * 100:+.3f} % |"
            )
    return lines


def family_table(cases: list[dict]) -> list[str]:
    lines = [
        "| hull form | waterlines | stations | buttocks |",
        "|---|---|---|---|",
    ]
    for case in cases:
        entry = next(
            entry for entry in case["sweeps"]["off"]
            if abs(entry["tolerance_fraction_of_half_beam"] - 0.008) < 1e-12
        )
        cells = [
            f"{entry['flips_by_direction'][key]['before']} -> "
            f"{entry['flips_by_direction'][key]['after']}"
            for key in ("waterlines", "stations", "buttocks")
        ]
        lines.append(f"| `{case['form']}` | " + " | ".join(cells) + " |")
    return lines


def feature_table(cases: list[dict]) -> list[str]:
    lines = [
        "| hull form | budget [% B/2] | flips before | flips, plain | flips, features "
        "| rms before [m] | rms, plain [m] | rms, features [m] | feature rows |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for case in cases:
        for plain, featured in zip(case["sweeps"]["off"], case["sweeps"]["auto"]):
            rows = ", ".join(str(index) for index in featured["feature_waterlines"]) or "none"
            lines.append(
                f"| `{case['form']}` "
                f"| {plain['tolerance_fraction_of_half_beam'] * 100:.1f} "
                f"| {plain['flips_before']} | {plain['flips_after']} | {featured['flips_after']} "
                f"| {plain['reference']['rms_before_m']:.5f} "
                f"| {plain['reference']['rms_after_m']:.5f} "
                f"| {featured['reference']['rms_after_m']:.5f} | {rows} |"
            )
    return lines


def render_figure(cases: list[dict], target: Path) -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (ax_flips, ax_rms) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    for case in cases:
        form = case["form"]
        plain = case["sweeps"]["off"]
        featured = case["sweeps"]["auto"]
        budgets = [entry["tolerance_fraction_of_half_beam"] * 100.0 for entry in plain]
        ax_flips.plot(budgets, [entry["flips_after"] for entry in plain], "o--",
                      color=COLOURS[form], lw=1.2, alpha=0.85,
                      label=f"{LABELS[form]}: plain energy")
        ax_flips.plot(budgets, [entry["flips_after"] for entry in featured], "o-",
                      color=COLOURS[form], lw=1.8,
                      label=f"{LABELS[form]}: feature rows")
        ax_rms.plot(budgets, [entry["reference"]["rms_after_m"] for entry in plain], "o--",
                    color=COLOURS[form], lw=1.2, alpha=0.85,
                    label=f"{LABELS[form]}: plain energy")
        ax_rms.plot(budgets, [entry["reference"]["rms_after_m"] for entry in featured], "o-",
                    color=COLOURS[form], lw=1.8,
                    label=f"{LABELS[form]}: feature rows")
        ax_rms.axhline(plain[0]["reference"]["rms_before_m"], color=COLOURS[form], lw=0.8,
                       ls=":", alpha=0.9)

    ax_flips.set_xlabel("deviation budget [% of B/2]")
    ax_flips.set_ylabel("curvature sign flips after fairing")
    ax_flips.set_title("flips (lower is better), same jittered tables")
    ax_flips.legend(fontsize=7)
    ax_flips.grid(alpha=0.3)

    ax_rms.set_xlabel("deviation budget [% of B/2]")
    ax_rms.set_ylabel("rms distance to the analytic table [m]")
    ax_rms.set_title("distance to the manual-fairing stand-in (dotted = before fairing)")
    ax_rms.legend(fontsize=7)
    ax_rms.grid(alpha=0.3)

    figure.suptitle("design-feature rows (near-keel band of a power-law section) in the fairing")
    figure.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=150)
    plt.close(figure)
    return target.read_bytes()


def _replace_block(text: str, signature: str, block: list[str]) -> str:
    """Replace the table whose header starts with ``signature`` (up to the next blank line)."""
    start = text.index(signature)
    end = text.index("\n\n", start)
    return text[:start] + "\n".join(block) + text[end:]


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    check = "--check" in argv
    cases = json.loads(RESULTS.read_text(encoding="utf-8"))["cases"]
    text = DOCUMENT.read_text(encoding="utf-8")
    blocks = (
        ("| hull form | flips of the analytic table", plain_table(cases)),
        ("| hull form | waterlines | stations", family_table(cases)),
        ("| hull form | budget [% B/2]", feature_table(cases)),
    )

    if check:
        problems: list[str] = []
        for entry in cases:
            for mode in ("off", "auto"):
                for row in entry["sweeps"][mode]:
                    if mode == "auto" and row["feature_waterlines"] and not row["within_tolerance"]:
                        problems.append(
                            f"{entry['form']} @ {row['tolerance_fraction_of_half_beam']}: "
                            "a feature-row run left the tolerance"
                        )
        for signature, block in blocks:
            missing = [line for line in block[2:] if line not in text]
            if signature not in text:
                problems.append(f"table with signature {signature!r} is missing")
            elif missing:
                problems.append(f"stale table {signature!r}: {len(missing)} row(s) differ")
        original = FIGURE.read_bytes()
        rendered = render_figure(cases, FIGURE)
        if original != rendered:
            problems.append("feature_rows.png does not match a fresh render of results.json")
        for problem in problems:
            print("FAIL:", problem)
        print("feature rows / tables:", "OK" if not problems else f"{len(problems)} problem(s)")
        return 1 if problems else 0

    for signature, block in blocks:
        text = _replace_block(text, signature, block)
    DOCUMENT.write_text(text, encoding="utf-8")
    render_figure(cases, FIGURE)
    print(f"rewrote three tables in {DOCUMENT.relative_to(ROOT)} and {FIGURE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())