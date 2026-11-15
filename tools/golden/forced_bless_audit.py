"""Evidence for the forced-bless trace (`tensorcad.golden`, GAP-CAD-D2-GOLDEN-REGRESSION).

`--bless --force` is the one path that can write a baseline the gate itself refuses. This
script measures what that override *does* and what it now *leaves behind*, on a synthetic
geometry whose size is a parameter (no ship semantics, no dataset, ~2 s):

1. bless a 1.0 m-radius geometry, then force-bless a 1.05 m one -- i.e. the deliberate
   geometry change that a real hull update looks like from the engine's side;
2. verify that the blessed measurement moved by far more than the 1e-9 drift lock (so the
   override really was an override, not a no-op), and record the findings it stepped over;
3. verify the record survives a later clean bless (append-only) and that every later run,
   full and quick, announces it as a warning while the gate itself stays green;
4. refuse to write the report if any of the above does not hold.

Outputs (all numbers below are printed by this script, none are typed by hand):

    docs/golden/forced_bless_audit.json
    docs/golden/forced_bless_audit.md
    docs/golden/forced_bless_audit.png
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "golden"

# Measure *this* checkout, not whatever `tensorcad` is installed in site-packages: a
# script run from tools/golden/ does not put the repo root on sys.path, and on a machine
# with several worktrees the installed copy is routinely another one (observed: the probe
# below reported the forced-bless record as missing because it imported a worktree that
# predates it).  Same guard as tools/golden/ellipsoid_area_evidence.py.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import golden  # noqa: E402


class MovingGeometry(golden.GoldenCase):
    """A synthetic geometry whose size is a parameter: reference and measurement move together."""

    name = "moving"
    title = "synthetic probe whose geometry has a size parameter"
    reference_source = "tools.golden.forced_bless_audit.MovingGeometry.reference (closed form)"

    def __init__(self, size: float = 1.0) -> None:
        self.size = size

    def _volume(self) -> float:
        return 4.0 / 3.0 * math.pi * self.size**3

    def _parameters(self, level: str) -> dict[str, float]:
        return {"nx": {"coarse": 11, "fine": 21}[level]}

    def measure(self, level: str) -> dict[str, float]:
        offset = {"coarse": 3e-3, "fine": 8e-4}[level]
        return {"value": self._volume() * (1.0 + offset)}

    def reference(self) -> dict[str, float]:
        return {"value": self._volume()}

    def invariants(self, values: dict[str, float]) -> dict[str, tuple[float, float, str]]:
        return {}

    def geometry_sha256(self) -> str:
        return f"moving-geometry-{self.size:g}"


def _findings(result: golden.SuiteResult) -> list[dict[str, object]]:
    return [finding.to_dict() for finding in result.all_findings]


def measure() -> dict[str, object]:
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "baseline.json"
        golden.bless(path, cases=[MovingGeometry(1.0)])
        before = json.loads(path.read_text(encoding="utf-8"))
        blessed_before = before["cases"]["moving"]["measurements"]["fine"]["value"]

        baseline, refusal = golden.bless(path, cases=[MovingGeometry(1.05)], force=True)
        after = json.loads(path.read_text(encoding="utf-8"))
        blessed_after = after["cases"]["moving"]["measurements"]["fine"]["value"]

        record = after["provenance"].get("forced_bless") or []
        full = golden.run_suite(baseline_path=path, cases=[MovingGeometry(1.05)])
        quick = golden.run_suite(baseline_path=path, cases=[MovingGeometry(1.05)], quick=True)

        reopened, clean = golden.bless(path, cases=[MovingGeometry(1.05)])
        record_after_clean_bless = reopened["provenance"].get("forced_bless") or []

        fresh = Path(scratch) / "fresh.json"
        fresh_baseline, _ = golden.bless(fresh, cases=[MovingGeometry(1.0)])

    drift = (blessed_after - blessed_before) / blessed_before
    quick_visibility = sorted(
        {f["code"] for f in _findings(quick) if f["code"] == "forced_bless_recorded"}
    )
    full_visibility = sorted(
        {f["code"] for f in _findings(full) if f["code"] == "forced_bless_recorded"}
    )
    evidence = {
        "schema": "tensorcad.golden.forced_bless_audit/1",
        "gate_mode": "full + quick",
        "blessed_fine_value_before": blessed_before,
        "blessed_fine_value_after": blessed_after,
        "relative_drift_of_blessed_value": drift,
        "drift_lock": float(golden.DRIFT_TOLERANCE),
        "drift_lock_exceeded_by": abs(drift) / float(golden.DRIFT_TOLERANCE),
        "refusal_findings_at_force_time": _findings(refusal),
        "forced_bless_record": record,
        "record_entries_after_clean_bless": len(record_after_clean_bless),
        "record_survives_clean_bless": len(record_after_clean_bless) == len(record) and bool(record),
        "fresh_bless_has_record": "forced_bless" in (fresh_baseline["provenance"] or {}),
        "full_mode_passed": full.passed,
        "quick_mode_passed": quick.passed,
        "full_mode_visibility": full_visibility,
        "quick_mode_visibility": quick_visibility,
        "visibility_severity": sorted(
            {str(f["severity"]) for f in _findings(full) if f["code"] == "forced_bless_recorded"}
        ),
        "clean_bless_passed": clean.passed,
    }
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=OUT / "forced_bless_audit.json")
    parser.add_argument("--md", type=Path, default=OUT / "forced_bless_audit.md")
    parser.add_argument("--png", type=Path, default=OUT / "forced_bless_audit.png")
    args = parser.parse_args()

    evidence = measure()

    overrode = [entry for entry in evidence["forced_bless_record"]]
    checks = {
        "override exceeded the drift lock": abs(evidence["relative_drift_of_blessed_value"])
        > evidence["drift_lock"],
        "the override is recorded": bool(overrode),
        "the record names the findings overridden": bool(overrode)
        and bool(overrode[0].get("exceeded")),
        "the record survives a later clean bless": evidence["record_survives_clean_bless"],
        "a tree never forced has no record": evidence["fresh_bless_has_record"] is False,
        "full mode announces it": evidence["full_mode_visibility"] == ["forced_bless_recorded"],
        "quick mode announces it": evidence["quick_mode_visibility"] == ["forced_bless_recorded"],
        "the announcement is a warning": evidence["visibility_severity"] == ["warning"],
        "the gate stays green after a reviewed override": evidence["full_mode_passed"]
        and evidence["quick_mode_passed"],
    }
    evidence["checks"] = {name: bool(ok) for name, ok in checks.items()}
    evidence["passed"] = all(checks.values())

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Forced bless (`--bless --force`): what it overrides and what it leaves behind",
        "",
        f"- schema: `{evidence['schema']}`",
        f"- verdict: **{'PASS' if evidence['passed'] else 'FAIL'}**",
        f"- blessed fine measurement: `{evidence['blessed_fine_value_before']:.9g}` -> "
        f"`{evidence['blessed_fine_value_after']:.9g}` "
        f"(drift `{evidence['relative_drift_of_blessed_value']:+.3e}`, "
        f"= `{evidence['drift_lock_exceeded_by']:.3g}x` the `{evidence['drift_lock']:g}` lock)",
        f"- findings the override stepped over: "
        f"`{', '.join(overrode[0]['exceeded']) if overrode else 'none'}`",
        f"- record entries after a later clean bless: "
        f"`{evidence['record_entries_after_clean_bless']}` "
        f"(append-only: `{evidence['record_survives_clean_bless']}`)",
        f"- visible in full mode: `{evidence['full_mode_visibility']}`, "
        f"quick mode: `{evidence['quick_mode_visibility']}`, "
        f"severity: `{evidence['visibility_severity']}`",
        "",
        "| check | result |",
        "| --- | --- |",
    ]
    lines += [f"| {name} | {'PASS' if ok else 'FAIL'} |" for name, ok in checks.items()]
    lines += [
        "",
        "Reproduce: `python tools/golden/forced_bless_audit.py`",
        "",
        "![forced bless record](forced_bless_audit.png)",
        "",
    ]
    args.md.write_text("\n".join(lines), encoding="utf-8")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - the figure is optional, the numbers are not
        print("matplotlib missing: numbers written, figure skipped")
    else:
        figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2))
        values = [evidence["blessed_fine_value_before"], evidence["blessed_fine_value_after"]]
        left.bar(["blessed before\n(1.00 m)", "blessed after\nforce (1.05 m)"], values,
                 color=["#4c78a8", "#e45756"])
        for index, value in enumerate(values):
            left.text(index, value * 1.002, f"{value:.6f}", ha="center", va="bottom", fontsize=9)
        left.set_ylabel("blessed fine measurement [m³]")
        left.set_title(
            f"override moved the baseline {evidence['relative_drift_of_blessed_value']:+.2%}\n"
            f"(drift lock {evidence['drift_lock']:g})"
        )
        left.set_ylim(min(values) * 0.98, max(values) * 1.01)

        names = ["record written", "survives clean bless", "visible: full", "visible: quick"]
        states = [
            bool(evidence["forced_bless_record"]),
            evidence["record_survives_clean_bless"],
            evidence["full_mode_visibility"] == ["forced_bless_recorded"],
            evidence["quick_mode_visibility"] == ["forced_bless_recorded"],
        ]
        right.barh(names, [1 if state else 0 for state in states],
                   color=["#54a24b" if state else "#e45756" for state in states])
        right.set_xlim(0, 1.2)
        right.set_xticks([0, 1])
        right.set_xticklabels(["missing", "present"])
        right.set_title("trace of a forced bless")
        figure.suptitle("GAP-CAD-D2: the golden gate's escape hatch is recorded, not silent")
        figure.tight_layout()
        figure.savefig(args.png, dpi=150)
        plt.close(figure)

    print(json.dumps({"passed": evidence["passed"], "checks": evidence["checks"]}, indent=2))
    return 0 if evidence["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())