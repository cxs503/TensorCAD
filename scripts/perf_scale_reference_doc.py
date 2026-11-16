"""Re-render or verify the evidence tables of ``docs/perf_scale_reference.md``.

The reference document is the anchor's synthesis of **every** committed perf
ladder (mesh, voxel, fleet, single-body subdivision, file import).  Its tables
are rendered from committed measurement JSON, never typed by hand; ``--check`` is
what CI and the test suite run, ``--write`` regenerates the document after a new
measurement run.

Why this script exists: the four ``perf-evidence`` markers in that document were
previously bound to nothing -- no renderer was registered under those names and no
script or test read the file, so the tables were hand-maintained while carrying
markers that promise programmatic injection.  Three of the four had drifted from
their artifacts (the single-body table stopped at 3.3e6 triangles while the
artifact reaches 1.3e7).  Binding the names makes the promise checkable.

Usage::

    python scripts/perf_scale_reference_doc.py --check
    python scripts/perf_scale_reference_doc.py --write
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_evidence import (  # noqa: E402
    REFERENCE_TABLE_NAMES,
    check_document,
    load_artifact,
    render_document,
)

DEFAULT_DOC = ROOT / "docs" / "perf_scale_reference.md"
PERF = ROOT / "benchmarks" / "perf_scale"
DEFAULT_ARTIFACTS = {
    "mesh": PERF / "artifacts" / "ladder_mesh.json",
    "voxel": PERF / "artifacts" / "ladder_voxel.json",
    "fleet3pt": PERF / "ladder_fleet_10m_3pt.json",
    "single10m": PERF / "artifacts" / "ladder_single_10m.json",
    "import_single10m": PERF / "artifacts" / "ladder_import_single_10m.json",
    "import_fleet": PERF / "artifacts" / "ladder_import_fleet_10m.json",
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--doc", default=str(DEFAULT_DOC))
    parser.add_argument("--write", action="store_true", help="rewrite the document")
    parser.add_argument("--check", action="store_true", help="verify (default)")
    for key, path in DEFAULT_ARTIFACTS.items():
        parser.add_argument(f"--{key.replace('_', '-')}", dest=key, default=str(path))
    args = parser.parse_args(argv)

    artifacts = {k: load_artifact(getattr(args, k)) for k in DEFAULT_ARTIFACTS}
    doc_path = Path(args.doc)
    text = doc_path.read_text(encoding="utf-8")

    if args.write:
        doc_path.write_text(
            render_document(text, artifacts, names=REFERENCE_TABLE_NAMES),
            encoding="utf-8",
        )
        print(f"rewrote {doc_path}")
        return 0

    problems = check_document(text, artifacts, names=REFERENCE_TABLE_NAMES)
    for p in problems:
        print(f"DRIFT: {p}")
    print("ok" if not problems else f"{len(problems)} problem(s)")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())