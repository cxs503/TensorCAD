"""Generate the evidence table for the neural-format verification audit.

Reads nothing by hand: every number in the audit document is written here from a
live measurement, so the doc can be regenerated and diffed instead of trusted.

    python examples/neural_format_audit/run_audit.py --output runs/2026-11-02-neural-format

Produces ``evidence.json`` (machine readable) and ``figure.png`` (the two-panel
plot embedded in the platform report). Standard library only for the audit
itself; matplotlib is needed for the figure only.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.neural_geometry import (
    Lineage,
    NeuralGeometryFormatError,
    ObjectStore,
    PayloadRef,
    TransactionRecord,
    build_manifest,
    dumps_canonical,
    parse_manifest,
    resolve_in_bundle,
    sha256_bytes,
)

CONTRACT = REPO_ROOT / "tests" / "test_neural_format_contract.py"


def load_contract_module() -> Any:
    """Reuse the test module's hostile corpus instead of maintaining a second one."""
    spec = importlib.util.spec_from_file_location("_tc_contract", CONTRACT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hostile_corpus(contract: Any) -> list[dict[str, Any]]:
    rows = []
    for case, overrides, expected in contract.HOSTILE_DERIVED:
        manifest = contract.make_manifest()
        if case == "derived_is_not_a_list" or case == "record_is_not_an_object":
            manifest["derived"] = overrides
        else:
            manifest["derived"] = [dict(contract.mesh_record(), **overrides)]
        try:
            parse_manifest(manifest)
        except NeuralGeometryFormatError as exc:
            rows.append({"case": case, "rejected": True, "message": str(exc),
                         "raises_other": False, "expected": expected})
        except Exception as exc:
            rows.append({"case": case, "rejected": False, "message": repr(exc),
                         "raises_other": True, "expected": expected})
        else:
            rows.append({"case": case, "rejected": False, "message": "",
                         "raises_other": False, "expected": expected})
    return rows


def fuzz(contract: Any, cases: int) -> dict[str, Any]:
    manifest = contract.make_manifest()
    manifest["derived"] = [contract.mesh_record()]
    rejected = accepted = 0
    undocumented: list[str] = []
    for candidate in contract._mutations(manifest, cases):
        try:
            parse_manifest(candidate)
        except NeuralGeometryFormatError:
            rejected += 1
        except Exception as exc:
            undocumented.append(f"{type(exc).__name__}: {exc}")
        else:
            accepted += 1
    return {"cases": cases, "rejected": rejected, "accepted": accepted,
            "undocumented_exceptions": undocumented}


def subprocess_checks(contract: Any, workdir: Path) -> dict[str, Any]:
    bundle = contract.publish(workdir, derived=[contract.mesh_record()])
    results: dict[str, Any] = {}
    for name, argv in (
        ("info", ("info", "--object", str(bundle))),
        ("verify", ("verify", "--json", "--object", str(bundle))),
    ):
        run = contract.run_cli(*argv)
        entry = {"exit_code": run.returncode, "stdout_head": run.stdout.strip().splitlines()[:6]}
        if name == "verify":
            payload = json.loads(run.stdout)
            entry["status"] = payload["status"]
            entry["checks"] = {check["name"]: {"passed": check["passed"],
                                              "severity": check["severity"]}
                               for check in payload["checks"]}
            entry["warnings"] = payload["warnings"]
        results[name] = entry
    return results


def containment_demo(workdir: Path, contract: Any) -> dict[str, Any]:
    """Show what a *naive* reader would do with a symlinked derived artifact.

    The naive branch is a faithful reimplementation of the pre-hardening read
    (``store.root / path`` with no resolution), kept here so the audit can quote
    the failure it fixed rather than describe it from memory.
    """
    outside = workdir / "outside.obj"
    outside.write_bytes(b"o secret mesh")
    digest = sha256_bytes(outside.read_bytes())
    root = workdir / "symlink.ngo"
    store = ObjectStore(root)
    manifest = build_manifest(
        object_id="urn:tensorcad:audit:symlink",
        revision=1,
        weights=PayloadRef("weights.safetensors", sha256_bytes(b"weights-v1"), 10),
        architecture_config={"coordinate_dim": 3, "output_dim": 1, "state_dim": 1,
                             "hidden_width": 8, "hidden_layers": 2},
        field_spec={"quantity": "truncated_signed_distance", "inside_sign": "negative",
                    "surface_isovalue": 0.0, "length_unit": "m", "truncation_distance": 0.15},
        coordinates={"dimension": 3, "frame": "object_local"},
        state={"schema_version": "0.2.0"},
        training={"schema_version": "0.2.0"},
        validation={"schema_version": "0.2.0", "status": "passed"},
        lineage=Lineage(parent_revision=None, root_revision=1, operation="fit"),
        provenance={"producer": "audit"},
    )
    manifest["derived"] = [{
        "kind": "surface_mesh", "format": "obj", "path": "derived/surface.obj",
        "sha256": digest, "bytes": outside.stat().st_size, "metrics": {},
    }]
    store.publish(
        files={"weights.safetensors": b"weights-v1", "derived/surface.obj": b"o mesh",
               "state.json": dumps_canonical({"schema_version": "0.2.0"}),
               "training.json": dumps_canonical({"schema_version": "0.2.0"}),
               "validation.json": dumps_canonical({"schema_version": "0.2.0"})},
        manifest=manifest,
        expected_parent=None,
        record=TransactionRecord(revision=1, parent_revision=None, status="committed",
                                 operation="fit", weights_sha256=sha256_bytes(b"weights-v1"),
                                 recorded_utc="2026-11-02T00:00:00+00:00"),
    )
    (root / "derived" / "surface.obj").unlink()
    os.symlink(outside, root / "derived" / "surface.obj")

    naive_target = root / "derived/surface.obj"          # pre-hardening read
    naive_hash = sha256_bytes(naive_target.read_bytes())
    try:
        resolve_in_bundle(root, "derived/surface.obj", "derived artifact")
        hardened = "accepted"
    except NeuralGeometryFormatError as exc:
        hardened = str(exc)
    verify = contract.run_cli("verify", "--json", "--object", str(root))
    payload = json.loads(verify.stdout) if verify.stdout.strip() else {"checks": []}
    checks = {check["name"]: check for check in payload.get("checks", [])}
    return {
        "outside_file": str(outside),
        "recorded_sha256": digest,
        "naive_read_sha256": naive_hash,
        "naive_read_would_match": naive_hash == digest,
        "hardened_verdict": hardened,
        "verify_exit_code": verify.returncode,
        "verify_containment_check": checks.get("derived_containment:derived/surface.obj"),
        "verify_hashed_the_escaped_payload": "derived_hash:derived/surface.obj" in checks,
    }


def figure(evidence: dict[str, Any], target: Path) -> str | None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        return f"matplotlib unavailable: {exc}"
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.4), dpi=160)

    corpus = evidence["hostile_corpus"]
    labels = [row["case"] for row in corpus]
    colours = ["#2e7d32" if row["rejected"] else "#c62828" for row in corpus]
    left.barh(range(len(corpus)), [1] * len(corpus), color=colours, height=0.72)
    left.set_yticks(range(len(corpus)), labels, fontsize=6.5)
    left.invert_yaxis()
    left.set_xlim(0, 1.15)
    left.set_xticks([])
    rejected = sum(1 for row in corpus if row["rejected"])
    left.set_title(f"Hostile derived[] manifests: {rejected}/{len(corpus)} rejected\n"
                   f"(all with NeuralGeometryFormatError)", fontsize=9)
    left.set_xlabel("rejected = green bar", fontsize=7.5)

    demo = evidence["containment_demo"]
    names = ["naive read\n(pre-hardening)", "resolve_in_bundle\n(hardened)"]
    verdicts = [1.0 if demo["naive_read_would_match"] else 0.0,
                1.0 if demo["verify_exit_code"] != 0 else 0.0]
    bars = right.bar(names, verdicts, color=["#c62828", "#2e7d32"], width=0.55)
    right.set_ylim(0, 1.45)
    right.set_yticks([])
    right.set_title("Symlink escape of derived/surface.obj", fontsize=9)
    notes = [
        f"hash matches outside payload: {demo['naive_read_would_match']}",
        f"verify exit {demo['verify_exit_code']}, payload unopened: "
        f"{not demo['verify_hashed_the_escaped_payload']}",
    ]
    for bar, note in zip(bars, notes, strict=False):
        right.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.06, note,
                   ha="center", va="bottom", fontsize=7.5, wrap=True)
    right.set_xlabel("caught / not caught", fontsize=7.5)

    fig.suptitle("TensorCAD NGO reader contract -- measured verification debt closure "
                 "(2026-11-02)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(target)
    plt.close(fig)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="directory for evidence.json and figure.png")
    args = parser.parse_args()
    out = Path(args.output).resolve()  # symlink demos need absolute targets
    out.mkdir(parents=True, exist_ok=True)

    contract = load_contract_module()
    evidence = {
        "generated_utc": subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                                        capture_output=True, text=True).stdout.strip(),
        "python": sys.version.split()[0],
        "banned_modules": list(contract.BANNED_MODULES),
        "installed_safetensors": importlib.util.find_spec("safetensors") is not None,
        "hostile_corpus": hostile_corpus(contract),
        "fuzz": fuzz(contract, contract.FUZZ_CASES),
        "subprocess": subprocess_checks(contract, out / "cli"),
        "containment_demo": containment_demo(out, contract),
    }
    evidence["figure_error"] = figure(evidence, out / "figure.png")
    (out / "evidence.json").write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: evidence[key] for key in ("fuzz", "containment_demo", "figure_error")},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
