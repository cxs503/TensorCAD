"""Re-render or verify the large-geometry evidence tables in the document.

The tables in ``docs/perf_scale_10m.md`` are generated from committed
measurement JSON (``benchmarks/perf_scale/artifacts/*.json``), never typed by
hand.  ``--check`` is what CI/tests run; ``--write`` updates the document after
a fresh measurement run.

Usage::

    python scripts/perf_scale_doc.py --check
    python scripts/perf_scale_doc.py --write \
        --mesh benchmarks/perf_scale/artifacts/ladder_mesh.json \
        --voxel benchmarks/perf_scale/artifacts/ladder_voxel.json \
        --fleet3pt benchmarks/perf_scale/ladder_fleet_10m_3pt.json
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_evidence import (  # noqa: E402
    check_document,
    load_artifact,
    render_document,
)

DEFAULT_DOC = ROOT / "docs" / "perf_scale_10m.md"
DEFAULT_MESH = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_mesh.json"
DEFAULT_VOXEL = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_voxel.json"
DEFAULT_FLEET3PT = ROOT / "benchmarks" / "perf_scale" / "ladder_fleet_10m_3pt.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--doc", default=str(DEFAULT_DOC))
    parser.add_argument("--mesh", default=str(DEFAULT_MESH))
    parser.add_argument("--voxel", default=str(DEFAULT_VOXEL))
    parser.add_argument("--fleet3pt", default=str(DEFAULT_FLEET3PT),
                        help="multi-body 10^6-10^7 ladder JSON (scripts/perf_scale_10m.py --bases)")
    parser.add_argument("--write", action="store_true", help="rewrite the document")
    parser.add_argument("--check", action="store_true", help="verify (default)")
    args = parser.parse_args(argv)

    doc_path = Path(args.doc)
    artifacts = {
        "mesh": load_artifact(args.mesh),
        "voxel": load_artifact(args.voxel),
        "fleet3pt": load_artifact(args.fleet3pt),
    }
    text = doc_path.read_text(encoding="utf-8")

    if args.write:
        doc_path.write_text(render_document(text, artifacts), encoding="utf-8")
        print(f"rewrote {doc_path}")
        return 0

    problems = check_document(text, artifacts)
    for p in problems:
        print(f"DRIFT: {p}")
    print("ok" if not problems else f"{len(problems)} problem(s)")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())