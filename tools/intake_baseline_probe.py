"""Before/after probe for foreign mesh intake (anchor GAP-CAD-B2-MESH-FORMAT-MATRIX).

Runs the same five foreign fixtures through the **committed** (``git show <ref>:
tensorcad/mesh_formats.py``) reader set and through the **working-tree** reader set,
and records what each one does.  ``mesh_formats`` has no package-relative imports,
so the baseline copy can be imported as a standalone module — the comparison is
therefore against the real previous code, not a paraphrase of it.

Usage::

    python3 tools/intake_baseline_probe.py [--ref HEAD] [--out artifacts/.../before_after.json]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import numpy as np  # noqa: E402

import mesh_intake_matrix as fixtures  # noqa: E402
from tensorcad import mesh_formats as current  # noqa: E402


def load_baseline(ref: str):
    blob = subprocess.run(
        ["git", "show", f"{ref}:tensorcad/mesh_formats.py"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout
    tmp = Path(tempfile.mkdtemp()) / "mesh_formats_baseline.py"
    tmp.write_text(blob, encoding="utf-8")
    spec = importlib.util.spec_from_file_location("mesh_formats_baseline", tmp)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses needs the module registered
    spec.loader.exec_module(module)
    return module


CASES = (
    ("quad_cage", "gmsh", fixtures.write_quad_cage),
    ("quad_cage", "su2", fixtures.write_quad_cage),
    ("quad_cage", "cgns", fixtures.write_quad_cage),
    ("planar_ndime2", "su2", fixtures.write_planar_su2),
    ("marker_only_quads", "su2", fixtures.write_marker_only_su2),
    ("multiblock_cgns", "cgns", fixtures.write_multiblock_cgns),
)


def probe(module, outdir: Path) -> list[dict]:
    rows = []
    for case, fmt, writer in CASES:
        path = outdir / "files" / f"in_{case}.{'msh' if fmt == 'gmsh' else fmt}"
        writer(path, fmt)
        row = {"case": case, "format": fmt}
        try:
            got = module.read(path, fmt)
            row.update(
                {
                    "outcome": "ok",
                    "faces": int(len(got.faces)),
                    "volume": float(got.volume()),
                    "euler": int(got.euler()),
                }
            )
        except Exception as exc:
            row.update({"outcome": "raised", "error": f"{type(exc).__name__}: {exc}"})
        rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD")
    ap.add_argument("--out", default=str(ROOT / "artifacts" / "m20943_foreign_intake" / "before_after.json"))
    args = ap.parse_args()

    outdir = Path(args.out).parent
    (outdir / "files").mkdir(parents=True, exist_ok=True)
    baseline = load_baseline(args.ref)
    before = probe(baseline, outdir)
    after = probe(current, outdir)
    payload = {
        "baseline_ref": args.ref,
        "baseline_commit": subprocess.run(
            ["git", "rev-parse", args.ref], cwd=ROOT, capture_output=True, text=True
        ).stdout.strip(),
        "before": before,
        "after": after,
        "newly_accepted": [
            f"{b['format']}/{b['case']}"
            for b, a in zip(before, after)
            if b["outcome"] == "raised" and a["outcome"] == "ok"
        ],
        "regressions": [
            f"{b['format']}/{b['case']}"
            for b, a in zip(before, after)
            if b["outcome"] == "ok" and a["outcome"] != "ok"
        ],
    }
    Path(args.out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    for b, a in zip(before, after):
        print(f"{a['format']:5s} {a['case']:20s} before={b['outcome']:7s} after={a['outcome']:7s} "
              f"faces={a.get('faces', '-')} volume={a.get('volume', float('nan')):.6f}"
              if a["outcome"] == "ok"
              else f"{a['format']:5s} {a['case']:20s} before={b['outcome']:7s} after={a['outcome']:7s}")
    print("newly accepted:", payload["newly_accepted"])
    print("regressions:", payload["regressions"])
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())