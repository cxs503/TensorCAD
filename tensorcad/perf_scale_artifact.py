"""Loss-free transport of large-scale benchmark results out of a cluster job.

Why this module exists
----------------------
The 10^6--10^7 triangle band of ``GAP-CAD-D4-PERF-SCALE`` can only be measured
on the compute nodes, and the first two attempts to do so (``cj2967``,
``cj2976``) produced **no data**, not because the physics or the code failed
but because the artifact did not survive the trip: the jobs wrote
``artifacts/perf_scale/ladder_fleet_10m.json`` inside the node worktree, and the
fetch step only returned files that were already tracked at submit time.  A
multi-minute measurement was silently thrown away by a packaging convention.

This module makes the result travel on several independent channels, so a
single packaging quirk cannot destroy a run:

1. ``write_result``  - JSON file at a *tracked* path inside the worktree.
2. ``encode_payload``/``emit`` - the whole result, gzipped and base64-encoded,
   printed to stdout between explicit begin/end markers.  Job logs are always
   returned, so this is the channel that cannot be lost.
3. ``emit`` also prints a compact human table, so a truncated log still shows
   the numbers.

``parse_emitted`` is the inverse of ``emit``: it reconstructs the result dict
from a captured log and **fails loudly** on a truncated payload (a missing
chunk is an error, never a partial result).  ``python -m
tensorcad.perf_scale_artifact <log> -o <json>`` does that recovery from the
shell.

The encoding is stdlib only (``gzip``, ``base64``, ``json``) and stable across
Python versions.
"""
from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path
from typing import Any, Iterable

#: Marker lines delimiting the machine-readable payload in a job log.
BEGIN_MARKER = "PERF10M_JSON_BEGIN"
END_MARKER = "PERF10M_JSON_END"

#: Prefix of one base64 chunk line.  Chunks are emitted so that no single line
#: exceeds common log-wrapping widths (some collectors hard-wrap at 1 KiB).
CHUNK_PREFIX = "PERF10M_B64"

#: Characters per base64 chunk (must be a multiple of 4).
CHUNK_CHARS = 900


class PayloadError(ValueError):
    """Raised when an emitted payload is missing, truncated or corrupt."""


def encode_payload(payload: Any) -> str:
    """Encode a JSON-serialisable payload as base64(gzip(utf-8 json))."""
    raw = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")


def decode_payload(encoded: str) -> Any:
    """Inverse of :func:`encode_payload`."""
    try:
        raw = gzip.decompress(base64.b64decode(encoded.encode("ascii"), validate=True))
    except Exception as exc:  # noqa: BLE001 - normalised into one error type
        raise PayloadError(f"payload is not valid base64+gzip: {exc}") from exc
    try:
        return json.loads(raw.decode("utf-8"))
    except Exception as exc:  # noqa: BLE001
        raise PayloadError(f"payload is not valid JSON: {exc}") from exc


def chunk_payload(encoded: str, chunk_chars: int = CHUNK_CHARS) -> list[str]:
    """Split an encoded payload into fixed-size chunks."""
    if chunk_chars <= 0:
        raise ValueError("chunk_chars must be positive")
    return [encoded[i:i + chunk_chars] for i in range(0, len(encoded), chunk_chars)] or [""]


def format_emitted(payload: Any, chunk_chars: int = CHUNK_CHARS) -> str:
    """Render the marker block (without the human table)."""
    chunks = chunk_payload(encode_payload(payload), chunk_chars=chunk_chars)
    total = len(chunks)
    lines = [BEGIN_MARKER, f"# chunks={total} chunk_chars={chunk_chars}"]
    lines += [f"{CHUNK_PREFIX} {i + 1}/{total} {chunk}" for i, chunk in enumerate(chunks)]
    lines.append(END_MARKER)
    return "\n".join(lines) + "\n"


def parse_emitted(text: str) -> Any:
    """Recover the payload emitted by :func:`format_emitted` from a log text.

    Raises :class:`PayloadError` when the begin/end markers are absent or when
    a chunk is missing - a partially delivered measurement must not be mistaken
    for a complete one.
    """
    lines = text.splitlines()
    try:
        start = lines.index(BEGIN_MARKER)
    except ValueError:
        raise PayloadError(f"begin marker {BEGIN_MARKER!r} not found in log") from None
    try:
        end = lines.index(END_MARKER, start + 1)
    except ValueError:
        raise PayloadError(f"end marker {END_MARKER!r} not found after begin marker") from None

    expected: int | None = None
    chunks: dict[int, str] = {}
    for line in lines[start + 1:end]:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            for token in line[1:].split():
                if token.startswith("chunks="):
                    expected = int(token.split("=", 1)[1])
            continue
        if not line.startswith(CHUNK_PREFIX):
            continue
        parts = line.split()
        if len(parts) != 3:
            raise PayloadError(f"malformed chunk line: {line[:60]!r}")
        seq_txt, total_txt = parts[1].split("/")
        data = parts[2]
        seq, total = int(seq_txt), int(total_txt)
        if expected is None:
            expected = total
        if total != expected:
            raise PayloadError(f"chunk count mismatch: expected {expected}, chunk says {total}")
        chunks[seq] = data

    if expected is None:
        raise PayloadError("no payload chunks found between markers")
    missing = [i for i in range(1, expected + 1) if i not in chunks]
    if missing:
        raise PayloadError(
            f"payload truncated: {len(missing)} of {expected} chunks missing "
            f"(first missing: {missing[0]})"
        )
    return decode_payload("".join(chunks[i] for i in range(1, expected + 1)))


def _format_table(payload: dict[str, Any]) -> str:
    """Compact human table: one line per measurement."""
    rows = payload.get("results", [])
    out = ["subtask             bodies  triangles    grid_cells      seconds   build_s   index_MB   rss_MB"]
    for r in rows:
        secs = r.get("seconds")
        build = r.get("build_seconds")
        out.append(
            "{:<19} {:>6} {:>10} {:>13} {:>12} {:>9} {:>10} {:>8}".format(
                str(r.get("subtask"))[:19],
                r.get("copies", 1),
                f"{r.get('triangles', 0):,}",
                f"{r.get('grid_cells', 0):,}",
                f"{secs:.3f}" if isinstance(secs, (int, float)) else str(r.get("status")),
                f"{build:.2f}" if isinstance(build, (int, float)) else "-",
                f"{(r.get('index_bytes') or 0) / 1e6:.1f}",
                f"{r.get('peak_rss_mb'):.0f}" if r.get("peak_rss_mb") else "-",
            )
        )
    for r in rows:
        if r.get("subtask") == "volume" and r.get("note"):
            out.append(f"  volume check: {r['note']}")
    for gate in payload.get("gates", []):
        out.append(f"  gate {gate.get('name')}: passed={gate.get('passed')} {gate.get('detail', '')}")
    for warn in payload.get("cluster_warnings", []):
        out.append(f"  cluster: {warn}")
    return "\n".join(out) + "\n"


def emit(payload: Any, stream=None, chunk_chars: int = CHUNK_CHARS) -> str:
    """Write the marker block and the human table; return the rendered block."""
    import sys

    stream = stream if stream is not None else sys.stdout
    block = format_emitted(payload, chunk_chars=chunk_chars)
    stream.write(block)
    stream.write(_format_table(payload))
    stream.flush()
    return block


def write_result(payload: Any, path: str | Path) -> Path:
    """Write the result JSON to ``path`` (created with its parents)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def recover_log(log_path: str | Path, out_path: str | Path) -> dict[str, Any]:
    """Parse a captured job log and write the recovered result JSON."""
    payload = parse_emitted(Path(log_path).read_text(encoding="utf-8", errors="replace"))
    write_result(payload, out_path)
    return payload


def _main(argv: Iterable[str] | None = None) -> int:
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Recover a perf_scale result payload from a captured job log"
    )
    parser.add_argument("log", help="captured job log (stdout+stderr)")
    parser.add_argument("-o", "--output", required=True, help="destination JSON path")
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        payload = recover_log(args.log, args.output)
    except PayloadError as exc:
        print(f"PERF10M_RECOVER_FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"recovered {len(payload.get('results', []))} measurements -> {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())