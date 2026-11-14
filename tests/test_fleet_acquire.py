"""Tests for digest-pinned acquisition (:mod:`tensorcad.fleet.acquire`).

Two kinds of assertion here, and the difference matters:

* **transport** — the bytes arrived and hash to the pinned digest (file://, http:// on a
  loopback server, and the failure modes: mismatch, no pin, oversize, bad scheme, dead
  port);
* **honesty invariant** — ``verified = True`` appears if and only if a pin was matched,
  an admitted block's recorded digest equals the bytes on disk, and every withheld block
  names a reason. A gate that only executes code proves nothing; these assert the
  property the anchor's credibility rests on.
"""

from __future__ import annotations

import dataclasses
import hashlib
import http.server
import json
import socket
import threading
from pathlib import Path

import pytest

from tensorcad.fleet.acquire import (
    PLAN_SCHEMA,
    AcquisitionPlan,
    DigestMismatch,
    PlanError,
    SchemeRefused,
    SizeExceeded,
    SourcePlan,
    SourceUnavailable,
    UnpinnedDigest,
    acquire_fleet,
    acquire_ship,
    fetch,
    load_plan,
    read_custodian,
    self_audit,
    selftest,
)
from tensorcad.fleet.build import build_ship_geometry, write_binary_stl
from tensorcad.fleet.manifest import load_manifest


def _digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@pytest.fixture(scope="module")
def surrogate_stl(tmp_path_factory) -> tuple[Path, str]:
    """A hull-shaped STL built by this repository's route, written and digested once."""
    directory = tmp_path_factory.mktemp("surrogate")
    manifest = load_manifest()
    record = manifest.by_name("KCS")
    built = build_ship_geometry(record, nx=241, nz=49)
    path = directory / "kcs_surrogate.stl"
    digest = write_binary_stl(built.mesh, path, header="acquisition test surrogate")
    return path, digest


# ------------------------------------------------------------------ transport layer


def test_file_url_pins_the_digest(tmp_path):
    payload = b"iges-ish bytes for the digest test\n"
    source = tmp_path / "hull.igs"
    source.write_bytes(payload)
    fetched = fetch(source.as_uri(), expected_sha256=_digest(payload))
    assert fetched.payload == payload
    assert fetched.verified is True
    assert fetched.provenance["sha256"] == _digest(payload)
    assert fetched.provenance["pinned_sha256"] == _digest(payload)
    assert fetched.provenance["scheme"] == "file"
    assert fetched.provenance["bytes"] == len(payload)
    assert fetched.provenance["retrieved_at"].endswith("+00:00")


def test_http_roundtrip_records_transport_trace(tmp_path):
    payload = b"\x00\x01official hull bytes\x02" * 64

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802 - stdlib naming
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("ETag", '"abc123"')
            self.end_headers()
            self.wfile.write(payload)

        def log_message(self, *args):  # silence the test output
            return

    server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        url = f"http://127.0.0.1:{server.server_port}/hull.bin"
        fetched = fetch(url, expected_sha256=_digest(payload), timeout=10.0)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    assert fetched.verified is True
    assert fetched.payload == payload
    assert fetched.provenance["content_type"] == "application/octet-stream"
    assert fetched.provenance["etag"] == '"abc123"'
    assert fetched.provenance["status"] == 200


def test_one_flipped_byte_is_a_refusal(tmp_path):
    payload = b"the official file" * 10
    source = tmp_path / "hull.igs"
    source.write_bytes(payload)
    with pytest.raises(DigestMismatch) as excinfo:
        fetch(source.as_uri(), expected_sha256=_digest(b"the official file" * 10 + b"!"))
    assert excinfo.value.expected != excinfo.value.observed
    assert excinfo.value.expected in str(excinfo.value)
    assert excinfo.value.observed in str(excinfo.value)
    assert excinfo.value.reason == "digest_mismatch"


def test_unpinned_payload_is_refused_and_reports_what_to_pin():
    payload = b"bytes without a pin"
    # No pin, so the digest of the bytes is the only thing that can be reported.
    path = Path("/tmp/_tcad_unpinned_probe.bin")
    path.write_bytes(payload)
    try:
        with pytest.raises(UnpinnedDigest) as excinfo:
            read_custodian(path, expected_sha256=None)
    finally:
        path.unlink()
    assert excinfo.value.observed == _digest(payload)
    assert excinfo.value.observed in str(excinfo.value)
    assert excinfo.value.reason == "unpinned_digest"


def test_size_cap_refuses_a_truncated_read(tmp_path):
    source = tmp_path / "big.igs"
    source.write_bytes(b"x" * 4096)
    with pytest.raises(SizeExceeded):
        read_custodian(source, expected_sha256=None, max_bytes=1024)


def test_scheme_allowlist():
    with pytest.raises(SchemeRefused):
        fetch("ftp://example.invalid/hull.igs", expected_sha256="0" * 64)


def test_dead_port_is_source_unavailable_not_a_crash():
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    with pytest.raises(SourceUnavailable) as excinfo:
        fetch(f"http://127.0.0.1:{port}/never.igs", expected_sha256="0" * 64, timeout=2.0)
    assert excinfo.value.reason == "source_unavailable"


# --------------------------------------------------------------------------- plan


def test_packaged_plan_names_every_ship_and_its_open_requirements():
    plan = load_plan()
    assert plan.ships == ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"]
    assert plan.sha256 and len(plan.sha256) == 64
    for entry in plan.entries:
        problems = entry.problems()
        # An order is either pinned (url + digest retrieved, so nothing is missing) or it
        # says out loud what is missing. A plan entry may never be silently partial: that
        # is how 'official geometry is unavailable' turns into an unexamined assumption.
        assert entry.license.strip()
        if entry.url and entry.sha256:
            assert not [item for item in problems if item.startswith(("url:", "sha256:"))], entry.ship
        else:
            assert any(item.startswith("url:") for item in problems), entry.ship
            assert any(item.startswith("sha256:") for item in problems), entry.ship


def test_plan_rejects_an_unknown_format(tmp_path):
    document = {
        "schema": PLAN_SCHEMA,
        "ships": [{"ship": "KCS", "file_format": "dwg", "license": "x", "locator": "y"}],
    }
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(document))
    with pytest.raises(PlanError):
        load_plan(path)


def test_plan_rejects_duplicate_ships(tmp_path):
    entry = {"ship": "KCS", "file_format": "stl", "license": "x", "locator": "y"}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps({"schema": PLAN_SCHEMA, "ships": [entry, entry]}))
    with pytest.raises(PlanError):
        load_plan(path)


def test_overrides_for_ships_outside_the_plan_are_rejected(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text(
        json.dumps(
            {"schema": PLAN_SCHEMA, "ships": [{"ship": "KCS", "file_format": "stl",
                                               "license": "x", "locator": "y"}]}
        )
    )
    with pytest.raises(PlanError):
        acquire_fleet(plan_path=path, out_dir=tmp_path / "out", pins=["KVLCC2=" + "0" * 64, "KCS=" + "0" * 64])


# --------------------------------------------------- acquire_ship / honesty invariant


def _kcs_entry(**overrides) -> SourcePlan:
    base = dict(
        ship="KCS",
        anchor_id="test-kcs",
        file_format="stl",
        url=None,
        sha256=None,
        license="test fixture, internal use",
        locator="surrogate STL built by the series60_form route (test fixture)",
    )
    base.update(overrides)
    return SourcePlan(**base)


@pytest.fixture(scope="module")
def fleet():
    return load_manifest()


def test_no_url_and_no_custody_is_withheld_with_a_reason(fleet, tmp_path):
    block = acquire_ship(_kcs_entry(), fleet.by_name("KCS"), manifest=fleet, geometry_dir=tmp_path / "geo")
    assert block["status"] == "withheld"
    assert block["reason"] == "plan_incomplete"
    assert block["passed"] is False
    assert block["provenance"] is None


def test_custodian_file_without_a_pin_is_withheld(fleet, tmp_path, surrogate_stl):
    path, _ = surrogate_stl
    block = acquire_ship(_kcs_entry(), fleet.by_name("KCS"), manifest=fleet,
                         geometry_dir=tmp_path / "geo", custodian=path)
    assert block["status"] == "withheld"
    assert block["reason"] == "unpinned_digest"
    assert block["acquisition"] is None  # nothing was written, nothing is claimed


def test_tampered_custodian_file_is_withheld(fleet, tmp_path, surrogate_stl):
    path, digest = surrogate_stl
    tampered = tmp_path / "tampered.stl"
    payload = bytearray(path.read_bytes())
    payload[-100] ^= 0xFF
    tampered.write_bytes(bytes(payload))
    block = acquire_ship(_kcs_entry(), fleet.by_name("KCS"), manifest=fleet,
                         geometry_dir=tmp_path / "geo", custodian=tampered, pin=digest)
    assert block["status"] == "withheld"
    assert block["reason"] == "digest_mismatch"
    assert not (tmp_path / "geo" / "KCS.stl").exists()


def test_pinned_custodian_file_is_admitted_and_passes_the_gates(fleet, tmp_path, surrogate_stl):
    path, digest = surrogate_stl
    block = acquire_ship(_kcs_entry(), fleet.by_name("KCS"), manifest=fleet,
                         geometry_dir=tmp_path / "geo", custodian=path, pin=digest)
    assert block["status"] == "admitted", block.get("reason_detail")
    assert block["provenance"]["verified"] is True
    assert block["provenance"]["retrieval"] == "custodian_file"
    written = Path(block["geometry"]["path"])
    assert written.exists()
    assert _digest(written.read_bytes()) == digest == block["geometry"]["sha256"]
    names = {gate["name"] for gate in block["gates"]}
    assert {"G0_provenance", "G1_not_placeholder", "G3_watertight", "G7_volume_vs_reference",
            "G8_principal_dimensions"} <= names
    assert block["failed_gates"] == []
    assert abs(block["volume_relative_error"]) < 0.02


def test_custodian_format_overrides_the_ordered_format(fleet, tmp_path, surrogate_stl):
    path, digest = surrogate_stl
    entry = _kcs_entry(file_format="iges")  # ordered IGES, custodian delivered STL
    block = acquire_ship(entry, fleet.by_name("KCS"), manifest=fleet,
                         geometry_dir=tmp_path / "geo", custodian=path, pin=digest)
    assert block["geometry"]["route"] == "imported_stl"
    assert block["geometry"]["planned_route"] == "iges_field"
    assert "format_deviation" in block
    assert block["status"] == "admitted"


def test_self_audit_catches_a_dishonest_block(fleet, tmp_path, surrogate_stl):
    path, digest = surrogate_stl
    block = acquire_ship(_kcs_entry(), fleet.by_name("KCS"), manifest=fleet,
                         geometry_dir=tmp_path / "geo", custodian=path, pin=digest)
    honest = {"ships": [block]}
    assert self_audit(honest)["ok"] is True
    assert self_audit(honest)["digests_rechecked"] == 1

    lied = {**block, "provenance": {**block["provenance"], "verified": False}}
    audit = self_audit({"ships": [lied]})
    assert audit["ok"] is False
    assert any("verified" in problem for problem in audit["problems"])

    withheld = {"ships": [{**block, "status": "withheld", "reason": None}]}
    assert self_audit(withheld)["ok"] is False


def test_acquire_fleet_counts_admitted_and_withheld(tmp_path, surrogate_stl):
    path, digest = surrogate_stl
    plan_document = {
        "schema": PLAN_SCHEMA,
        "ships": [
            {"ship": "KCS", "file_format": "stl", "license": "fixture", "locator": "fixture"},
            {"ship": "KVLCC2", "file_format": "iges", "license": "fixture", "locator": "fixture"},
        ],
    }
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan_document))
    evidence = acquire_fleet(plan_path=plan_path, out_dir=tmp_path / "out",
                             custodians=[f"KCS={path}"], pins=[f"KCS={digest}"],
                             nx=241, nz=49)
    assert evidence["summary"]["ships_ordered"] == 2
    assert evidence["summary"]["ships_admitted"] == 1
    assert evidence["summary"]["withheld_reasons"] == {"plan_incomplete": 1}
    assert evidence["self_audit"]["ok"] is True
    assert (tmp_path / "out" / "acquisition_evidence.json").exists()


def test_selftest_controls_all_behave(tmp_path):
    report = selftest(tmp_path, nx=241, nz=49)
    assert report["all_controls_correct"] is True
    observed = {control["control"]: control["observed"] for control in report["controls"]}
    assert observed == {
        "pinned_payload": "admitted",
        "tampered_payload": "digest_mismatch",
        "unpinned_payload": "unpinned_digest",
        "unreachable_source": "source_unavailable",
    }
    assert "SURROGATE" in report["surrogate"]["path"].upper() or report["surrogate"]["disclaimer"]


def test_packaged_plan_is_a_real_plan_object():
    plan = load_plan()
    assert isinstance(plan, AcquisitionPlan)
    assert plan.for_ship("KCS") is not None
    assert plan.missing_for(["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"]) == []
    assert plan.for_ship("KVLCC1") is None