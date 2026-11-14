"""Tests for the cluster-artifact transport used by the 10^7 perf band.

The point of these tests is the failure mode that actually lost two cluster
runs (``cj2967``, ``cj2976``): the measurement was fine, the packaging was not.
So the assertions are about *invariants of transport* - roundtrip fidelity,
truncation detection, and chunk accounting - not about re-running a benchmark.
"""
from __future__ import annotations

import io
import json

import pytest

from tensorcad import perf_scale_artifact as art


def _payload(n_results: int = 3) -> dict:
    return {
        "schema": "tensorcad.perf_benchmark/2",
        "passed": True,
        "results": [
            {
                "subtask": "point_locate_accel",
                "copies": i + 1,
                "triangles": 4_000_000 * (i + 1),
                "grid_cells": 0,
                "seconds": 1.234 + i,
                "build_seconds": 31.1,
                "index_bytes": 319_400_000,
                "peak_rss_mb": 6151.0 + i,
                "status": "measured",
                "note": f"volume=7.5 m^3, analytic=7.5 m^3, rel_err=+0.001%, bodies={i + 1}",
            }
            for i in range(n_results)
        ],
        "gates": [{"name": "scaling_volume", "passed": True, "detail": "exp=1.01<=1.5"}],
        "cluster_warnings": ["10^7 triangles"],
    }


def test_roundtrip_is_lossless():
    payload = _payload()
    assert art.decode_payload(art.encode_payload(payload)) == payload


def test_parse_emitted_recovers_payload_from_mixed_log():
    payload = _payload(5)
    buf = io.StringIO()
    art.emit(payload, stream=buf, chunk_chars=120)
    log = "some noisy job header\n" + buf.getvalue() + "trailing noise\n"
    assert art.parse_emitted(log) == payload


def test_payload_is_chunked_for_long_results():
    payload = _payload(40)
    block = art.format_emitted(payload, chunk_chars=120)
    chunk_lines = [ln for ln in block.splitlines() if ln.startswith(art.CHUNK_PREFIX)]
    assert len(chunk_lines) > 1, "large payloads must be split into several lines"
    assert all(len(ln) < 400 for ln in chunk_lines)
    assert art.parse_emitted(block) == payload


def test_truncated_payload_is_an_error_not_a_partial_result():
    payload = _payload(40)
    block = art.format_emitted(payload, chunk_chars=120)
    lines = block.splitlines()
    # Drop one chunk line in the middle (what a log truncation/collector limit does).
    drop = next(i for i, ln in enumerate(lines) if ln.startswith(art.CHUNK_PREFIX) and i % 2)
    del lines[drop]
    with pytest.raises(art.PayloadError, match="truncated"):
        art.parse_emitted("\n".join(lines))


def test_missing_markers_raise():
    with pytest.raises(art.PayloadError, match="begin marker"):
        art.parse_emitted("no payload here\n")
    block = art.format_emitted(_payload(), chunk_chars=120)
    with pytest.raises(art.PayloadError, match="end marker"):
        art.parse_emitted(block.replace(art.END_MARKER, ""))


def test_corrupt_base64_raises_payload_error():
    with pytest.raises(art.PayloadError):
        art.decode_payload("!!!not base64!!!")


def test_recover_log_writes_json(tmp_path):
    payload = _payload(2)
    log = tmp_path / "job.log"
    log.write_text("header\n" + art.format_emitted(payload), encoding="utf-8")
    out = tmp_path / "sub" / "recovered.json"
    recovered = art.recover_log(log, out)
    assert recovered == payload
    assert json.loads(out.read_text(encoding="utf-8")) == payload


def test_recover_log_exit_code_on_missing_payload(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("nothing useful\n", encoding="utf-8")
    assert art._main([str(log), "-o", str(tmp_path / "x.json")]) == 1
    assert not (tmp_path / "x.json").exists()