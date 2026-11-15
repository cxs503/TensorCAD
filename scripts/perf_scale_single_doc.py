"""Re-render or verify the single-body 10^6-10^7 evidence tables.

The tables in ``docs/perf_scale_single_body_10m.md`` are generated from the
committed measurement JSON
(``benchmarks/perf_scale/artifacts/ladder_single_10m.json``), never typed by
hand.  ``--check`` is what tests/CI run; ``--write`` updates the document after a
fresh measurement run.

Usage::

    python scripts/perf_scale_single_doc.py --check
    python scripts/perf_scale_single_doc.py --write
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_evidence import (  # noqa: E402
    SINGLE10M_TABLE_NAMES,
    check_document,
    load_artifact,
    render_document,
)

DEFAULT_DOC = ROOT / "docs" / "perf_scale_single_body_10m.md"
DEFAULT_ARTIFACT = (
    ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_single_10m.json"
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--doc", default=str(DEFAULT_DOC))
    parser.add_argument("--artifact", default=str(DEFAULT_ARTIFACT))
    parser.add_argument("--write", action="store_true", help="rewrite the document")
    parser.add_argument("--check", action="store_true", help="verify (default)")
    args = parser.parse_args(argv)

    doc_path = Path(args.doc)
    artifacts = {"single10m": load_artifact(args.artifact)}
    text = doc_path.read_text(encoding="utf-8")

    if args.write:
        doc_path.write_text(
            render_document(text, artifacts, names=SINGLE10M_TABLE_NAMES),
            encoding="utf-8",
        )
        print(f"rewrote {doc_path}")
        return 0

    problems = check_document(text, artifacts, names=SINGLE10M_TABLE_NAMES)
    for problem in problems:
        print(f"DRIFT: {problem}")
    print("ok" if not problems else f"{len(problems)} problem(s)")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())