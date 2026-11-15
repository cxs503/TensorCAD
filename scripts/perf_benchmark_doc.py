"""Re-render or verify the tables of ``docs/perf_benchmark.md``.

The gate-conventions manual embeds two tables, both rendered from committed
sources, never typed by hand:

* ``baseline-gate-table`` <- ``benchmarks/perf_scale/baseline.json``
* ``code-facts-table``    <- the constants of ``tensorcad.perf_benchmark``

``--check`` is what tests/CI run; ``--write`` refreshes the document after a
baseline or threshold change.  This is the companion of
``scripts/perf_scale_doc.py`` (which covers the measurement report
``docs/perf_scale_10m.md``): different document, different table bundle.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_evidence import (  # noqa: E402
    BENCHMARK_DOC_TABLE_NAMES,
    check_document,
    load_artifact,
    render_document,
)

DEFAULT_DOC = ROOT / "docs" / "perf_benchmark.md"
DEFAULT_BASELINE = ROOT / "benchmarks" / "perf_scale" / "baseline.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--doc", default=str(DEFAULT_DOC))
    parser.add_argument("--baseline", default=str(DEFAULT_BASELINE))
    parser.add_argument("--write", action="store_true", help="rewrite the document")
    parser.add_argument("--check", action="store_true", help="verify (default)")
    args = parser.parse_args(argv)

    doc_path = Path(args.doc)
    artifacts = {"baseline": load_artifact(args.baseline)}
    text = doc_path.read_text(encoding="utf-8")

    if args.write:
        doc_path.write_text(
            render_document(text, artifacts, names=BENCHMARK_DOC_TABLE_NAMES),
            encoding="utf-8",
        )
        print(f"rewrote {doc_path}")
        return 0

    problems = check_document(text, artifacts, names=BENCHMARK_DOC_TABLE_NAMES)
    for p in problems:
        print(f"DRIFT: {p}")
    print("ok" if not problems else f"{len(problems)} problem(s)")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())