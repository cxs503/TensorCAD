"""Digest-pinned acquisition of official hull geometry (GAP-CAD-T1-BASELINE-FLEET).

The intake gate in :mod:`tensorcad.fleet.gates` decides *whether* a geometry file may
join the baseline fleet, and :mod:`tensorcad.fleet.sources` records *where* each
literal in the manifest comes from. Between those two there was a hole: nothing in the
package could take an official IGES/STL/offset file — from an URL or from a file a
custodian dropped on disk — check it against a pinned digest, write the provenance that
the gate demands, and hand it to the geometry route. That hole is why the acquisition
half of the anchor is reported as WITHHELD: the *code path* to receive the file did not
exist, so it could not be exercised even if the file arrived.

This module closes it, under one rule that is enforced in code rather than in prose:

    **``verified = True`` requires a digest that was pinned before the bytes were read
    and matched the bytes that were read.**

Consequences, all of them refusals instead of warnings:

* no pinned digest → :class:`UnpinnedDigest`. The digest of the bytes is reported so a
  custodian can pin it, but the run cannot call the file verified.
* pinned digest ≠ observed digest → :class:`DigestMismatch` (both digests in the
  message). A one-byte difference is a refusal, not a rounding note.
* transport failure → :class:`SourceUnavailable`; the ship is recorded ``withheld`` with
  the reason, never ``admitted`` with a surrogate substituted behind it.
* payload larger than the entry's ``max_bytes`` → :class:`SizeExceeded`; a truncated
  read can never be mistaken for the file.

An admitted file then goes through exactly the same path as a generated hull
(:func:`tensorcad.fleet.build.build_ship_geometry` →
:func:`tensorcad.fleet.gates.evaluate_gates`), so a vendor file is judged by the
placeholders / watertight / manifold / volume / L-B-T gates and not by its provenance
alone.

Entry point::

    python -m tensorcad.fleet.acquire --out-dir runs/m19948_acquire
    python -m tensorcad.fleet.acquire --custodian KCS=/data/kcs.igs --pin KCS=<sha256>

The ``--custodian``/``--url``/``--pin`` overrides are the "one command" property: when
the official file finally lands, nothing needs editing — the plan stays honest
(``url: null``, ``sha256: null``) and the custody facts enter on the command line and
land in the evidence.
"""

from __future__ import annotations

import dataclasses
import datetime as _dt
import hashlib
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from tensorcad.fleet.build import RouteRefusal, build_ship_geometry
from tensorcad.fleet.gates import VOLUME_TOLERANCE, evaluate_gates
from tensorcad.fleet.manifest import (
    FleetManifest,
    Provenance,
    ShipRecord,
    load_manifest,
    sha256_bytes,
)
from tensorcad.fleet.sources import SourceRegistry, load_registry

PLAN_SCHEMA = "tensorcad.fleet.acquisition_plan.v1"
EVIDENCE_SCHEMA = "tensorcad.fleet.acquisition.v1"
PLAN_PATH = Path(__file__).with_name("data") / "acquisition_plan.json"

#: Formats a plan entry may name, and the geometry route that consumes each one.
ROUTE_FOR_FORMAT = {
    "iges": "iges_field",
    "stl": "imported_stl",
    "offset_csv": "offset_table",
}
SUFFIX_FOR_FORMAT = {"iges": ".igs", "stl": ".stl", "offset_csv": ".csv"}
#: What a custodian's file extension implies; used only when it disagrees with the plan.
FORMAT_FOR_SUFFIX = {".igs": "iges", ".iges": "iges", ".stl": "stl", ".csv": "offset_csv"}
#: Extensions that count as the payload for each format (inside an archive as well).
EXTENSIONS_FOR_FORMAT = {"iges": (".igs", ".iges"), "stl": (".stl",), "offset_csv": (".csv",)}
#: Inverse of :data:`ROUTE_FOR_FORMAT`, so the file on disk carries the suffix its route expects.
FORMAT_FOR_ROUTE = {route: fmt for fmt, route in ROUTE_FOR_FORMAT.items()}
ALLOWED_SCHEMES = ("file", "http", "https")
DEFAULT_MAX_BYTES = 64 * 1024 * 1024
DEFAULT_TIMEOUT = 20.0


class AcquisitionError(RuntimeError):
    """A payload could not be turned into a verified geometry file."""

    reason = "acquisition_error"


class PlanError(AcquisitionError):
    reason = "plan_error"


class SourceUnavailable(AcquisitionError):
    """Transport failed: no bytes, so no claim. Recorded as ``withheld``."""

    reason = "source_unavailable"


class SchemeRefused(AcquisitionError):
    reason = "scheme_refused"


class SizeExceeded(AcquisitionError):
    reason = "size_exceeded"


class UnpinnedDigest(AcquisitionError):
    reason = "unpinned_digest"

    def __init__(self, observed: str, message: str = "") -> None:
        self.observed = observed
        super().__init__(
            message
            or f"no digest was pinned before reading: observed sha256 = {observed}. "
            "Pin it (--pin SHIP=<sha256>) to make the file verifiable; until then it is "
            "evidence of transport, not of identity"
        )


class DigestMismatch(AcquisitionError):
    reason = "digest_mismatch"

    def __init__(self, expected: str, observed: str) -> None:
        self.expected = expected
        self.observed = observed
        super().__init__(
            f"pinned sha256 = {expected} but the bytes read hash to {observed}: the payload is "
            "not the file the pin was taken from"
        )


# --------------------------------------------------------------------------- plan


@dataclass(frozen=True)
class SourcePlan:
    """One ship's acquisition order: what file, from where, pinned to what digest."""

    ship: str
    anchor_id: str
    file_format: str
    url: str | None
    sha256: str | None
    license: str
    locator: str
    expected_file: str = ""
    max_bytes: int = DEFAULT_MAX_BYTES
    notes: str = ""

    @property
    def route(self) -> str:
        return ROUTE_FOR_FORMAT[self.file_format]

    @property
    def suffix(self) -> str:
        return SUFFIX_FOR_FORMAT[self.file_format]

    def problems(self) -> list[str]:
        """Everything that stops this order from producing a *verified* file."""
        issues: list[str] = []
        if self.file_format not in ROUTE_FOR_FORMAT:
            issues.append(f"unknown file_format {self.file_format!r} (expected one of {sorted(ROUTE_FOR_FORMAT)})")
        if not self.url:
            issues.append("url: the source portal publishes the hull package but no single-file URL was "
                          "confirmed in this environment, so none is invented")
        elif urllib.parse.urlparse(self.url).scheme.lower() not in ALLOWED_SCHEMES:
            issues.append(f"url scheme {urllib.parse.urlparse(self.url).scheme!r} is not one of {ALLOWED_SCHEMES}")
        if not self.sha256:
            issues.append("sha256: no digest pinned, so a payload cannot be called verified")
        if not self.license.strip():
            issues.append("license: a file with no stated licence is not admissible")
        if self.max_bytes <= 0:
            issues.append("max_bytes must be positive")
        return issues

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "anchor_id": self.anchor_id,
            "file_format": self.file_format,
            "route": self.route,
            "url": self.url,
            "sha256": self.sha256,
            "license": self.license,
            "locator": self.locator,
            "expected_file": self.expected_file,
            "max_bytes": self.max_bytes,
            "notes": self.notes,
            "open_requirements": self.problems(),
        }


@dataclass(frozen=True)
class AcquisitionPlan:
    entries: tuple[SourcePlan, ...]
    sha256: str | None = None
    path: Path | None = None

    def for_ship(self, ship: str) -> SourcePlan | None:
        for entry in self.entries:
            if entry.ship == ship:
                return entry
        return None

    @property
    def ships(self) -> list[str]:
        return [entry.ship for entry in self.entries]

    def missing_for(self, ships: Iterable[str]) -> list[str]:
        return [name for name in ships if self.for_ship(name) is None]


def load_plan(path: str | Path | None = None) -> AcquisitionPlan:
    target = Path(path) if path is not None else PLAN_PATH
    payload = target.read_bytes()
    try:
        document = json.loads(payload)
    except json.JSONDecodeError as exc:  # pragma: no cover - malformed file
        raise PlanError(f"{target}: not valid JSON ({exc})") from exc
    if document.get("schema") != PLAN_SCHEMA:
        raise PlanError(f"{target}: schema {document.get('schema')!r} is not {PLAN_SCHEMA!r}")
    entries: list[SourcePlan] = []
    for raw in document.get("ships", []):
        try:
            entry = SourcePlan(
                ship=str(raw["ship"]),
                anchor_id=str(raw.get("anchor_id", "")),
                file_format=str(raw["file_format"]),
                url=None if raw.get("url") in (None, "") else str(raw["url"]),
                sha256=None if raw.get("sha256") in (None, "") else str(raw["sha256"]),
                license=str(raw.get("license", "")),
                locator=str(raw.get("locator", "")),
                expected_file=str(raw.get("expected_file", "")),
                max_bytes=int(raw.get("max_bytes", DEFAULT_MAX_BYTES)),
                notes=str(raw.get("notes", "")),
            )
        except KeyError as exc:
            raise PlanError(f"{target}: entry missing key {exc.args[0]!r}") from exc
        if entry.file_format not in ROUTE_FOR_FORMAT:
            raise PlanError(f"{target}: {entry.ship} file_format {entry.file_format!r} has no geometry route")
        entries.append(entry)
    if not entries:
        raise PlanError(f"{target}: no ship entries")
    names = [entry.ship for entry in entries]
    if len(set(names)) != len(names):
        raise PlanError(f"{target}: duplicate ship entries in {names}")
    return AcquisitionPlan(entries=tuple(entries), sha256=sha256_bytes(payload), path=target)


# -------------------------------------------------------------------------- fetch


@dataclass(frozen=True)
class Fetch:
    """The bytes that arrived plus the trace that lets a third party refetch them."""

    payload: bytes
    sha256: str
    provenance: dict[str, Any]

    @property
    def verified(self) -> bool:
        return bool(self.provenance.get("verified"))


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def fetch(
    url: str,
    *,
    expected_sha256: str | None,
    max_bytes: int = DEFAULT_MAX_BYTES,
    timeout: float = DEFAULT_TIMEOUT,
    retrieved_at: str | None = None,
) -> Fetch:
    """Read ``url`` and return the bytes only if they match ``expected_sha256``.

    ``file://`` is accepted as well as ``http(s)://`` so a custodian-supplied file and a
    downloaded file follow the same code path; the scheme is recorded in the provenance
    exactly because the two are *not* equally strong evidence.
    """
    parsed = urllib.parse.urlparse(url)
    scheme = parsed.scheme.lower()
    if scheme not in ALLOWED_SCHEMES:
        raise SchemeRefused(f"{url}: scheme {scheme!r} is not one of {ALLOWED_SCHEMES}")
    request = urllib.request.Request(url, headers={"User-Agent": "tensorcad-fleet-acquire/1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - scheme allowlisted
            payload = response.read(max_bytes + 1)
            headers = getattr(response, "headers", None)
            final_url = response.geturl() if hasattr(response, "geturl") else url
            status = getattr(response, "status", None)
    except urllib.error.HTTPError as exc:
        raise SourceUnavailable(f"{url}: HTTP {exc.code} {exc.reason}") from exc
    except (urllib.error.URLError, OSError, ValueError) as exc:
        raise SourceUnavailable(f"{url}: transport failed ({type(exc).__name__}: {exc})") from exc
    if len(payload) > max_bytes:
        raise SizeExceeded(f"{url}: payload exceeds max_bytes = {max_bytes} (read {len(payload)} bytes, truncated)")
    observed = sha256_bytes(payload)
    provenance: dict[str, Any] = {
        "url": url,
        "final_url": final_url,
        "scheme": scheme,
        "status": status,
        "bytes": len(payload),
        "sha256": observed,
        "content_type": headers.get("Content-Type") if headers else None,
        "etag": headers.get("ETag") if headers else None,
        "last_modified": headers.get("Last-Modified") if headers else None,
        "retrieved_at": retrieved_at or _utc_now(),
        "verified": False,
    }
    if not expected_sha256:
        raise UnpinnedDigest(observed)
    if expected_sha256.strip().lower() != observed:
        raise DigestMismatch(expected_sha256.strip().lower(), observed)
    provenance["verified"] = True
    provenance["pinned_sha256"] = expected_sha256.strip().lower()
    return Fetch(payload=payload, sha256=observed, provenance=provenance)


def _materialise(payload: bytes, *, fmt: str, source: str) -> tuple[bytes, dict[str, Any]]:
    """Return the geometry bytes, unwrapping a published ZIP that holds exactly one hull file.

    The *archive* is what the publisher ships and therefore what the digest pins; the
    member is extracted deterministically from it, and its own digest is recorded so a
    reader can tell which bytes became the mesh.
    """
    if not payload[:2] == b"PK":
        return payload, {"archive": None}
    import io
    import zipfile

    wanted_suffixes = EXTENSIONS_FOR_FORMAT[fmt]
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        members = [name for name in archive.namelist() if name.lower().endswith(wanted_suffixes)]
        if not members:
            raise AcquisitionError(
                f"{source}: archive holds no member with extension {wanted_suffixes} "
                f"(members: {archive.namelist()})"
            )
        if len(members) > 1:
            raise AcquisitionError(
                f"{source}: archive holds {len(members)} candidate members {members}; which one is the "
                "hull is a human decision, not a glob"
            )
        inner = archive.read(members[0])
    return inner, {
        "archive": "zip",
        "members": len(zipfile.ZipFile(io.BytesIO(payload)).namelist()),
        "member": members[0],
        "member_bytes": len(inner),
        "member_sha256": sha256_bytes(inner),
    }


def read_custodian(path: str | Path, *, expected_sha256: str | None,
                   max_bytes: int = DEFAULT_MAX_BYTES, retrieved_at: str | None = None) -> Fetch:
    """A file handed over by a custodian: same digest rule, no transport."""
    target = Path(path).resolve()
    if not target.is_file():
        raise SourceUnavailable(f"{target}: custodian file does not exist")
    size = target.stat().st_size
    if size > max_bytes:
        raise SizeExceeded(f"{target}: {size} bytes exceeds max_bytes = {max_bytes}")
    payload = target.read_bytes()
    observed = sha256_bytes(payload)
    provenance = {
        "url": target.as_uri(),
        "final_url": str(target),
        "scheme": "file",
        "status": None,
        "bytes": len(payload),
        "sha256": observed,
        "content_type": None,
        "etag": None,
        "last_modified": _dt.datetime.fromtimestamp(
            target.stat().st_mtime, _dt.timezone.utc
        ).replace(microsecond=0).isoformat(),
        "retrieved_at": retrieved_at or _utc_now(),
        "verified": False,
        "custody": "custodian_file",
    }
    if not expected_sha256:
        raise UnpinnedDigest(observed)
    if expected_sha256.strip().lower() != observed:
        raise DigestMismatch(expected_sha256.strip().lower(), observed)
    provenance["verified"] = True
    provenance["pinned_sha256"] = expected_sha256.strip().lower()
    return Fetch(payload=payload, sha256=observed, provenance=provenance)


# ------------------------------------------------------------------------ acquire


def _admitted_block(
    record: ShipRecord,
    entry: SourcePlan,
    fetched: Fetch,
    *,
    geometry_dir: Path,
    manifest: FleetManifest,
    registry: SourceRegistry | None,
    volume_tolerance: float,
    nx: int,
    nz: int,
    route: str | None = None,
    build_kwargs: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    route = route or entry.route
    geometry_path = geometry_dir / f"{record.name}{SUFFIX_FOR_FORMAT[FORMAT_FOR_ROUTE[route]]}"
    geometry_path.parent.mkdir(parents=True, exist_ok=True)
    inner, archive_info = _materialise(
        fetched.payload, fmt=FORMAT_FOR_ROUTE[route], source=str(fetched.provenance.get("url") or record.name)
    )
    geometry_path.write_bytes(inner)
    on_disk = sha256_bytes(geometry_path.read_bytes())
    # For an archived payload the geometry on disk is the extracted member, so that is what
    # has to match: the archive digest pins the transport, the member digest pins the mesh.
    expected_geometry = archive_info.get("member_sha256") or fetched.sha256
    if on_disk != expected_geometry:  # pragma: no cover - filesystem fault
        raise AcquisitionError(f"{geometry_path}: bytes on disk do not hash to the verified payload")
    provenance = Provenance(
        source=entry.locator or entry.anchor_id or entry.ship,
        url=str(fetched.provenance["final_url"] or fetched.provenance["url"]),
        license=entry.license,
        retrieval=("custodian_file" if fetched.provenance.get("custody") == "custodian_file" else "network_fetch"),
        verified=True,
    )
    admitted = dataclasses.replace(
        record,
        geometry_route=route,
        geometry_path=str(geometry_path),
        geometry_sha256=on_disk,
        provenance=provenance,
    )
    route_kwargs = {k: v for k, v in dict(build_kwargs or {}).items()
                    if route == "iges_field" and k in ("dx", "dz", "ring_points", "samples_per_metre")}
    built = build_ship_geometry(admitted, nx=nx, nz=nz, stl_path=str(geometry_path),
                                iges_path=str(geometry_path), **route_kwargs)
    verdict = evaluate_gates(admitted, built.mesh, manifest, volume_tolerance=volume_tolerance, registry=registry)
    return {
        "ship": record.name,
        "status": "admitted" if verdict["passed"] else "refused_by_gates",
        "reason": None if verdict["passed"] else "failed_gates: " + ",".join(verdict["failed_gates"]),
        "provenance": provenance.as_dict(),
        "acquisition": fetched.provenance,
        "geometry": {
            "path": str(geometry_path),
            "sha256": on_disk,
            "route": built.route,
            "planned_route": entry.route,
            "details": built.details,
            "payload": archive_info,
        },
        "gates": verdict["gates"],
        "failed_gates": verdict["failed_gates"],
        "passed": verdict["passed"],
        "topology": verdict["topology"],
        "volume_relative_error": verdict["volume_relative_error"],
        "reference_basis": verdict["reference_basis"],
    }


def acquire_ship(
    entry: SourcePlan,
    record: ShipRecord,
    *,
    manifest: FleetManifest,
    registry: SourceRegistry | None = None,
    geometry_dir: str | Path = "runs/acquire/geometry",
    url: str | None = None,
    pin: str | None = None,
    custodian: str | Path | None = None,
    timeout: float = DEFAULT_TIMEOUT,
    volume_tolerance: float = VOLUME_TOLERANCE,
    nx: int = 241,
    nz: int = 49,
    build_kwargs: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Try to acquire one ship. Returns a block; never raises for a transport/identity failure."""
    target_url = url or entry.url
    target_pin = pin or entry.sha256
    base = {
        "ship": record.name,
        "plan": entry.as_dict(),
        "open_requirements": entry.problems(),
    }
    try:
        if custodian is not None:
            fetched = read_custodian(custodian, expected_sha256=target_pin, max_bytes=entry.max_bytes)
        elif target_url:
            fetched = fetch(target_url, expected_sha256=target_pin, max_bytes=entry.max_bytes, timeout=timeout)
        else:
            return {**base, "status": "withheld", "reason": "plan_incomplete",
                    "reason_detail": "; ".join(entry.problems()),
                    "provenance": None, "acquisition": None, "geometry": None, "gates": None,
                    "passed": False, "topology": None, "volume_relative_error": None,
                    "reference_basis": record.reference_basis(),
                    "failed_gates": ["acquisition"]}
    except AcquisitionError as exc:
        return {**base, "status": "withheld", "reason": exc.reason, "reason_detail": str(exc),
                "provenance": None, "acquisition": None, "geometry": None, "gates": None,
                "passed": False, "topology": None, "volume_relative_error": None,
                "reference_basis": record.reference_basis(),
                "failed_gates": ["acquisition"]}
    route = entry.route
    if custodian is not None:
        implied = FORMAT_FOR_SUFFIX.get(Path(custodian).suffix.lower())
        if implied is not None and implied != entry.file_format:
            route = ROUTE_FOR_FORMAT[implied]
            base["format_deviation"] = (
                f"plan ordered {entry.file_format} ({entry.route}) but the custodian file is "
                f"{Path(custodian).suffix} -> {implied} ({route}); the payload is judged on {route}"
            )
    try:
        return {**base, **_admitted_block(
            record, entry, fetched,
            geometry_dir=Path(geometry_dir), manifest=manifest, registry=registry,
            volume_tolerance=volume_tolerance, nx=nx, nz=nz, route=route,
            build_kwargs=build_kwargs,
        )}
    except (RouteRefusal, AcquisitionError, ValueError, OSError) as exc:
        reason = getattr(exc, "reason", "geometry_route_refused")
        return {**base, "status": "refused_geometry", "reason": reason,
                "reason_detail": f"{type(exc).__name__}: {exc}",
                "provenance": None, "acquisition": fetched.provenance, "geometry": None,
                "gates": None, "passed": False, "topology": None,
                "volume_relative_error": None, "reference_basis": record.reference_basis(),
                "failed_gates": ["geometry_route"]}


def _split_overrides(values: Sequence[str] | None, what: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in values or ():
        if "=" not in raw:
            raise PlanError(f"--{what} expects SHIP=VALUE, got {raw!r}")
        ship, value = raw.split("=", 1)
        out[ship.strip()] = value.strip()
    return out


def self_audit(evidence: dict[str, Any]) -> dict[str, Any]:
    """Machine check of the honesty invariant, run on the evidence it is about.

    * every ``admitted`` block carries ``provenance.verified = True``;
    * every admitted block's recorded digest equals the digest of the file it points at
      *when that file still exists* (judged as ``not_run`` otherwise, not as a pass);
    * every non-admitted block names a reason.
    """
    problems: list[str] = []
    digest_checked = 0
    for block in evidence["ships"]:
        status = block["status"]
        if status == "admitted":
            provenance = block.get("provenance") or {}
            if not provenance.get("verified"):
                problems.append(f"{block['ship']}: admitted with provenance.verified != True")
            path = (block.get("geometry") or {}).get("path")
            if path and Path(path).is_file():
                observed = sha256_bytes(Path(path).read_bytes())
                digest_checked += 1
                if observed != (block.get("geometry") or {}).get("sha256"):
                    problems.append(f"{block['ship']}: file digest {observed} != recorded digest")
        elif status in ("withheld", "refused_by_gates", "refused_geometry"):
            if not block.get("reason"):
                problems.append(f"{block['ship']}: status {status} with no reason")
        else:
            problems.append(f"{block['ship']}: unknown status {status!r}")
    return {
        "ok": not problems,
        "problems": problems,
        "digests_rechecked": digest_checked,
        "digests_not_rechecked": sum(
            1 for block in evidence["ships"]
            if block["status"] == "admitted"
            and not Path((block.get("geometry") or {}).get("path", "")).is_file()
        ),
        "invariant": "verified=true iff a pinned digest matched the bytes actually read",
    }


def acquire_fleet(
    *,
    plan_path: str | Path | None = None,
    manifest_path: str | Path | None = None,
    registry_path: str | Path | None = None,
    out_dir: str | Path = "runs/m19948_acquire",
    geometry_dir: str | Path | None = None,
    urls: Sequence[str] | None = None,
    pins: Sequence[str] | None = None,
    custodians: Sequence[str] | None = None,
    timeout: float = DEFAULT_TIMEOUT,
    volume_tolerance: float = VOLUME_TOLERANCE,
    nx: int = 241,
    nz: int = 49,
    build_kwargs: Mapping[str, Any] | None = None,
    write: bool = True,
) -> dict[str, Any]:
    """Attempt every order in the plan and record the verdict for each."""
    plan = load_plan(plan_path)
    manifest = load_manifest(manifest_path)
    try:
        registry: SourceRegistry | None = load_registry(registry_path)
    except (OSError, ValueError):
        registry = None
    url_overrides = _split_overrides(urls, "url")
    pin_overrides = _split_overrides(pins, "pin")
    custodian_overrides = _split_overrides(custodians, "custodian")
    planned = set(plan.ships)
    unknown = (set(url_overrides) | set(pin_overrides) | set(custodian_overrides)) - planned
    if unknown:
        raise PlanError(f"overrides name ships that are not in the plan: {sorted(unknown)}")
    root = Path(out_dir)
    geo = Path(geometry_dir) if geometry_dir is not None else root / "geometry"
    if write:
        geo.mkdir(parents=True, exist_ok=True)
    ships = [
        acquire_ship(
            entry,
            manifest.by_name(entry.ship),
            manifest=manifest,
            registry=registry,
            geometry_dir=geo,
            url=url_overrides.get(entry.ship),
            pin=pin_overrides.get(entry.ship),
            custodian=custodian_overrides.get(entry.ship),
            timeout=timeout,
            volume_tolerance=volume_tolerance,
            nx=nx,
            nz=nz,
            build_kwargs=build_kwargs,
        )
        for entry in plan.entries
    ]
    admitted = [block["ship"] for block in ships if block["status"] == "admitted"]
    withheld = [block["ship"] for block in ships if block["status"] == "withheld"]
    refused = [block["ship"] for block in ships if block["status"].startswith("refused")]
    reasons: dict[str, int] = {}
    for block in ships:
        if block["status"] != "admitted":
            reasons[block.get("reason") or "unknown"] = reasons.get(block.get("reason") or "unknown", 0) + 1
    evidence: dict[str, Any] = {
        "schema": EVIDENCE_SCHEMA,
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "settings": {
            "plan_path": None if plan.path is None else str(plan.path),
            "plan_sha256": plan.sha256,
            "manifest_path": None if manifest.path is None else str(manifest.path),
            "geometry_dir": str(geo),
            "timeout": timeout,
            "nx": nx,
            "nz": nz,
            "volume_tolerance": volume_tolerance,
            "overrides": {
                "url": url_overrides,
                "pin": pin_overrides,
                "custodian": {k: str(v) for k, v in custodian_overrides.items()},
            },
        },
        "plan": [entry.as_dict() for entry in plan.entries],
        "ships": ships,
        "summary": {
            "ships_ordered": len(ships),
            "ships_admitted": len(admitted),
            "ships_withheld": len(withheld),
            "ships_refused": len(refused),
            "admitted": admitted,
            "withheld": withheld,
            "refused": refused,
            "withheld_reasons": reasons,
            "ships_verified_from_source": len(admitted),
            "honest_statement": (
                "verified_from_source counts ships whose bytes matched a digest pinned before reading. "
                "A withheld ship is not a failure of the fleet, it is the fleet telling you the official "
                "file has not arrived; no surrogate is substituted behind it."
            ),
        },
    }
    evidence["self_audit"] = self_audit(evidence)
    if write:
        root.mkdir(parents=True, exist_ok=True)
        (root / "acquisition_evidence.json").write_text(json.dumps(evidence, indent=2, default=str))
    return evidence

# ----------------------------------------------------------------------- selftest


def selftest(
    out_dir: str | Path = "runs/m19948_acquire",
    *,
    nx: int = 241,
    nz: int = 49,
    volume_tolerance: float = VOLUME_TOLERANCE,
) -> dict[str, Any]:
    """Exercise the intake path end to end without any network.

    A surrogate hull (this repository's Series-60 route at the KCS particulars) is written
    as a binary STL and then pushed through the *same* functions an official file would
    use. Three controls decide whether the path behaves:

    * ``pinned_payload`` — correct digest + custodian file → **admitted** and every gate
      re-evaluated on the re-read file;
    * ``tampered_payload`` — one byte flipped, original digest → **withheld
      (digest_mismatch)**;
    * ``unpinned_payload`` — correct bytes, no digest → **withheld (unpinned_digest)**.

    The surrogate is labelled as a surrogate in the STL header and in the custody record:
    passing this self-test proves the *pipe* is real, and says nothing about whether the
    official hulls have arrived.
    """
    from tensorcad.fleet.build import write_binary_stl

    root = Path(out_dir) / "selftest"
    geometry_dir = root / "geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest()
    registry = load_registry()
    record = manifest.by_name("KCS")
    built = build_ship_geometry(record, nx=nx, nz=nz)
    surrogate = root / "kcs_surrogate.stl"
    digest = write_binary_stl(
        built.mesh,
        surrogate,
        header="tensorcad fleet acquisition selftest SURROGATE (not the official KCS file)",
    )
    entry = SourcePlan(
        ship="KCS",
        anchor_id="selftest-surrogate",
        file_format="stl",
        url=None,
        sha256=digest,
        license="internal self-test surrogate, no third-party rights",
        locator="surrogate KCS built by the series60_form route for the acquisition self-test",
        expected_file="surrogate STL written by this self-test",
    )
    tampered = root / "kcs_surrogate_tampered.stl"
    payload = bytearray(surrogate.read_bytes())
    payload[len(payload) // 2] ^= 0x01
    tampered.write_bytes(bytes(payload))
    controls: list[dict[str, Any]] = []

    def record_control(name: str, expected: str, block: dict[str, Any]) -> None:
        observed = block["status"] if block["status"] == "admitted" else block.get("reason")
        controls.append(
            {
                "control": name,
                "expected": expected,
                "observed": observed,
                "correct": observed == expected,
                "detail": block.get("reason_detail"),
                "pipelines": str(block.get("geometry", {}).get("route")) if block.get("geometry") else None,
                "gates_failed": block.get("failed_gates"),
            }
        )

    pinned = acquire_ship(entry, record, manifest=manifest, registry=registry,
                          geometry_dir=geometry_dir, pin=digest, custodian=surrogate,
                          nx=nx, nz=nz, volume_tolerance=volume_tolerance)
    record_control("pinned_payload", "admitted", pinned)
    achieved = acquire_ship(entry, record, manifest=manifest, registry=registry,
                            geometry_dir=geometry_dir, pin=digest, custodian=tampered,
                            nx=nx, nz=nz, volume_tolerance=volume_tolerance)
    record_control("tampered_payload", "digest_mismatch", achieved)
    # Both pins must be absent for the control to mean anything: dropping only the
    # command-line pin leaves the plan's pin in force, which is the intended fallback.
    unpinned_entry = dataclasses.replace(entry, sha256=None)
    unpinned = acquire_ship(unpinned_entry, record, manifest=manifest, registry=registry,
                            geometry_dir=geometry_dir, pin=None, custodian=surrogate,
                            nx=nx, nz=nz, volume_tolerance=volume_tolerance)
    record_control("unpinned_payload", "unpinned_digest", unpinned)
    blocked = acquire_ship(entry, record, manifest=manifest, registry=registry,
                           geometry_dir=geometry_dir, pin=digest, custodian=None,
                           url="http://127.0.0.1:9/never-arrives.stl", timeout=2.0,
                           nx=nx, nz=nz, volume_tolerance=volume_tolerance)
    record_control("unreachable_source", "source_unavailable", blocked)
    return {
        "surrogate": {
            "path": str(surrogate),
            "sha256": digest,
            "route": built.route,
            "nx": nx,
            "nz": nz,
            "details": built.details,
            "disclaimer": "surrogate built by the series60_form route; not the official KCS geometry",
        },
        "admitted": pinned,
        "controls": controls,
        "all_controls_correct": all(control["correct"] for control in controls),
    }


# ---------------------------------------------------------------------------- CLI


def _plot(evidence: dict[str, Any], path: Path) -> str | None:
    """Two panels: what arrived per ship, and whether the intake pipe itself behaves."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is optional
        return f"matplotlib unavailable: {exc}"
    ships = evidence["ships"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    names = [block["ship"] for block in ships]
    codes = [0 if block["status"] == "admitted" else 1 for block in ships]
    colours = ["#2f855a" if code == 0 else "#b7791f" for code in codes]
    left.bar(names, codes, color=colours)
    left.set_ylim(0, 1.45)
    left.set_yticks([0, 1])
    left.set_yticklabels(["admitted (digest verified)", "withheld / refused"])
    left.set_title(f"Official-geometry acquisition, {len(ships)} ships ordered\n"
                   f"{evidence['summary']['ships_admitted']} admitted, "
                   f"{evidence['summary']['ships_withheld']} withheld")
    for bar, block in zip(left.bars if hasattr(left, "bars") else left.patches, ships):
        label = block["status"] if block["status"] == "admitted" else str(block.get("reason"))
        left.annotate(label, (bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.03),
                      ha="center", va="bottom", fontsize=7, rotation=90 if len(label) > 12 else 0)
    left.tick_params(axis="x", labelsize=8)

    controls = evidence.get("selftest", {}).get("controls", [])
    if controls:
        labels = [control["control"].replace("_", "\n") for control in controls]
        values = [0 if control["correct"] else 1 for control in controls]
        right.bar(labels, values, color=["#2f855a" if v == 0 else "#c53030" for v in values])
        right.set_ylim(0, 1.45)
        right.set_yticks([0, 1])
        right.set_yticklabels(["intake behaved", "INTAKE HOLE"])
        right.set_title(f"Intake pipeline controls "
                        f"({sum(1 for v in values if v == 0)}/{len(values)} correct)")
        for bar, control in zip(right.patches, controls):
            right.annotate(control["observed"], (bar.get_x() + bar.get_width() / 2, 0.04),
                           ha="center", va="bottom", fontsize=7, rotation=90)
    else:
        right.axis("off")
        right.text(0.5, 0.5, "self-test not run", ha="center", va="center")
    for axis in (left, right):
        axis.tick_params(axis="x", labelsize=8)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return None


def build_parser():  # pragma: no cover - argparse wiring
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.fleet.acquire",
        description="Digest-pinned acquisition of the six standard hulls "
                    "(GAP-CAD-T1-BASELINE-FLEET). A file is admitted only when its sha256 "
                    "matches a digest pinned before the bytes were read.",
    )
    parser.add_argument("--plan", default=None, help="acquisition plan JSON (default: packaged plan)")
    parser.add_argument("--manifest", default=None, help="fleet manifest JSON (default: packaged one)")
    parser.add_argument("--registry", default=None, help="source registry JSON (default: packaged one)")
    parser.add_argument("--out-dir", default="runs/m19948_acquire", help="directory for evidence + figure")
    parser.add_argument("--geometry-dir", default=None, help="where admitted payloads are written")
    parser.add_argument("--custodian", action="append", default=None, metavar="SHIP=PATH",
                        help="admit a file already on disk (retrieval=custodian_file)")
    parser.add_argument("--url", action="append", default=None, metavar="SHIP=URL",
                        help="override the plan URL for one ship")
    parser.add_argument("--pin", action="append", default=None, metavar="SHIP=SHA256",
                        help="pin the expected digest for one ship")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, help="per-request timeout [s]")
    parser.add_argument("--volume-tolerance", type=float, default=VOLUME_TOLERANCE)
    parser.add_argument("--nx", type=int, default=241, help="stations for the intake rebuild")
    parser.add_argument("--nz", type=int, default=49, help="levels for the intake rebuild")
    parser.add_argument("--iges-dx", type=float, default=1.0, help="half-breadth field spacing along x [m]")
    parser.add_argument("--iges-dz", type=float, default=0.2, help="half-breadth field spacing along z [m]")
    parser.add_argument("--samples-per-metre", type=float, default=3.0,
                        help="patch sampling density for the IGES route [1/m]")
    parser.add_argument("--ring-points", type=int, default=200, help="waterline ring resolution")
    parser.add_argument("--selftest", action="store_true",
                        help="run the offline intake-pipeline controls (surrogate STL, not official geometry)")
    parser.add_argument("--no-plot", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:  # pragma: no cover - thin wrapper
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    root = Path(args.out_dir)
    build_kwargs = {
        "dx": args.iges_dx,
        "dz": args.iges_dz,
        "samples_per_metre": args.samples_per_metre,
        "ring_points": args.ring_points,
    }
    try:
        evidence = acquire_fleet(
            plan_path=args.plan,
            manifest_path=args.manifest,
            registry_path=args.registry,
            out_dir=root,
            geometry_dir=args.geometry_dir,
            urls=args.url,
            pins=args.pin,
            custodians=args.custodian,
            timeout=args.timeout,
            volume_tolerance=args.volume_tolerance,
            nx=args.nx,
            nz=args.nz,
            build_kwargs=build_kwargs,
        )
    except PlanError as exc:
        print(f"plan error: {exc}", file=sys.stderr)
        return 2
    if args.selftest:
        evidence["selftest"] = selftest(root, nx=args.nx, nz=args.nz,
                                        volume_tolerance=args.volume_tolerance)
        (root / "acquisition_evidence.json").write_text(json.dumps(evidence, indent=2, default=str))
    if not args.no_plot:
        note = _plot(evidence, root / "acquisition_summary.png")
        if note:
            print(f"figure skipped: {note}", file=sys.stderr)
    print(json.dumps(evidence["summary"], indent=2))
    print(f"self_audit: ok={evidence['self_audit']['ok']} problems={evidence['self_audit']['problems']}")
    if args.selftest:
        print(f"pipeline controls correct: {evidence['selftest']['all_controls_correct']}")
    print(f"artefacts: {root / 'acquisition_evidence.json'}")
    ok = evidence["self_audit"]["ok"] and (not args.selftest or evidence["selftest"]["all_controls_correct"])
    return 0 if ok else 1


if __name__ == "__main__":  # pragma: no cover
    import sys

    raise SystemExit(main())
