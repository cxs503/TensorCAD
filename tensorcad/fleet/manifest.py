"""Baseline-fleet registry: provenance-carrying ship records and their gates.

The registry is the answer to a concrete failure mode recorded by the owner on
2026-11-13: the ship-platform geometry library ships eight STL files of which six
are 304-triangle placeholders, and two of the "different" ships are byte-identical.
Neither defect is detectable by looking at a geometry in isolation, so the fleet is
the unit of analysis: every record states where its numbers come from, every input
file is pinned by sha256, and the gates run over the whole set at once.

The module is deliberately literate about what it cannot prove. A record whose
provenance was entered offline reports ``retrieval="offline_recollection"`` and
``verified=False``, and a displacement volume that was *derived* from
``Cb*L*B*T`` rather than quoted declares ``basis="closure"``. The gate prints that
label next to the number, so a closure check can never be mistaken for an
independent validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

MANIFEST_PATH = Path(__file__).with_name("data") / "fleet_manifest.json"
SCHEMA = "tensorcad.fleet.baseline.v1"

#: Displacement-volume bases that are *not* an external measurement.
DERIVED_BASES = frozenset({"closure"})
#: Bases that require a provenance source to be stated.
ALL_BASES = frozenset({"literature_value_recalled", "literature_value_verified", "closure"})


class ManifestError(ValueError):
    """Raised when a fleet manifest is malformed or self-contradictory."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


@dataclass(frozen=True)
class Provenance:
    """Where a record's numbers come from, and how far that has been checked."""

    source: str
    url: str
    license: str
    retrieval: str
    verified: bool

    @property
    def complete(self) -> bool:
        return bool(self.source.strip()) and bool(self.url.strip()) and bool(self.license.strip())

    def as_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "url": self.url,
            "license": self.license,
            "retrieval": self.retrieval,
            "verified": self.verified,
        }


@dataclass(frozen=True)
class ShipRecord:
    """One hull of the baseline fleet plus everything needed to audit it."""

    name: str
    ship_class: str
    length: float
    beam: float
    draft: float
    block_coefficient: float
    displacement_basis: str
    provenance: Provenance
    displacement_volume: float | None = None
    lcb_fraction: float | None = None
    geometry_route: str = "series60_form"
    notes: str = ""
    midship_coefficient: float | None = None
    transom_ratio: float = 0.0
    geometry_sha256: str | None = None
    geometry_path: str | None = None

    # -- derived quantities -------------------------------------------------
    def closed_form_volume(self) -> float:
        """``Cb * L * B * T``: the volume implied by the declared particulars."""
        return self.block_coefficient * self.length * self.beam * self.draft

    def reference_volume(self) -> float:
        """Volume the mesh is compared against: published when present, else closure."""
        if self.displacement_volume is not None:
            return float(self.displacement_volume)
        return self.closed_form_volume()

    def reference_basis(self) -> str:
        if self.displacement_volume is not None:
            return self.displacement_basis
        return "closure"

    def is_closure_reference(self) -> bool:
        return self.reference_basis() in DERIVED_BASES

    def reference_is_external(self) -> bool:
        return self.displacement_volume is not None and self.provenance.verified

    def ratios(self) -> dict[str, float]:
        return {"length_to_beam": self.length / self.beam, "beam_to_draft": self.beam / self.draft}

    def problems(self) -> list[str]:
        issues: list[str] = []
        for name in ("length", "beam", "draft", "block_coefficient"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                issues.append(f"{name} must be finite and positive, got {value!r}")
        if not 0.2 <= self.block_coefficient <= 0.95:
            issues.append(f"block_coefficient {self.block_coefficient!r} outside [0.2, 0.95]")
        if self.displacement_basis not in ALL_BASES:
            issues.append(f"unknown displacement_basis {self.displacement_basis!r}")
        if self.displacement_volume is not None and self.displacement_volume <= 0.0:
            issues.append("displacement_volume must be positive when it is given")
        if self.midship_coefficient is not None and not (
            self.block_coefficient <= self.midship_coefficient <= 1.0
        ):
            issues.append(
                f"midship_coefficient {self.midship_coefficient!r} must lie in "
                f"[Cb, 1] = [{self.block_coefficient!r}, 1]"
            )
        if not 0.0 <= self.transom_ratio < 0.5:
            issues.append(f"transom_ratio {self.transom_ratio!r} outside [0, 0.5)")
        if not self.provenance.complete:
            issues.append("provenance must state source, url and license")
        if not self.name.strip():
            issues.append("name must not be empty")
        return issues

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "ship_class": self.ship_class,
            "length": self.length,
            "beam": self.beam,
            "draft": self.draft,
            "block_coefficient": self.block_coefficient,
            "displacement_volume": self.displacement_volume,
            "displacement_basis": self.displacement_basis,
            "reference_volume": self.reference_volume(),
            "reference_basis": self.reference_basis(),
            "reference_is_external": self.reference_is_external(),
            "lcb_fraction": self.lcb_fraction,
            "midship_coefficient": self.midship_coefficient,
            "transom_ratio": self.transom_ratio,
            "geometry_route": self.geometry_route,
            "geometry_path": self.geometry_path,
            "geometry_sha256": self.geometry_sha256,
            "provenance": self.provenance.as_dict(),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class FleetManifest:
    ships: tuple[ShipRecord, ...]
    placeholder: Mapping[str, Any] = field(default_factory=dict)
    path: Path | None = None
    sha256: str | None = None

    def by_name(self, name: str) -> ShipRecord:
        for ship in self.ships:
            if ship.name == name:
                return ship
        raise KeyError(name)

    def names(self) -> tuple[str, ...]:
        return tuple(ship.name for ship in self.ships)

    def duplicate_geometry(self) -> dict[str, tuple[str, ...]]:
        """Group ship names whose geometry sha256 collide (the KVLCC1 == KVLCC2 defect)."""
        buckets: dict[str, list[str]] = {}
        for ship in self.ships:
            if ship.geometry_sha256:
                buckets.setdefault(ship.geometry_sha256, []).append(ship.name)
        return {sha: tuple(names) for sha, names in buckets.items() if len(names) > 1}

    def problems(self) -> list[str]:
        issues: list[str] = []
        seen: set[str] = set()
        for ship in self.ships:
            if ship.name in seen:
                issues.append(f"duplicate ship name {ship.name!r}")
            seen.add(ship.name)
            issues.extend(f"{ship.name}: {problem}" for problem in ship.problems())
        if len(self.ships) < 1:
            issues.append("manifest contains no ships")
        for sha, names in self.duplicate_geometry().items():
            issues.append(f"byte-identical geometry {sha[:12]} shared by {names}")
        return issues


def _require(mapping: Mapping[str, Any], key: str, where: str) -> Any:
    if key not in mapping:
        raise ManifestError(f"{where}: missing required field {key!r}")
    return mapping[key]


def record_from_dict(payload: Mapping[str, Any], *, where: str = "ship") -> ShipRecord:
    provenance_payload = _require(payload, "provenance", where)
    provenance = Provenance(
        source=str(provenance_payload.get("source", "")),
        url=str(provenance_payload.get("url", "")),
        license=str(provenance_payload.get("license", "")),
        retrieval=str(provenance_payload.get("retrieval", "unspecified")),
        verified=bool(provenance_payload.get("verified", False)),
    )
    volume = payload.get("displacement_volume")
    lcb = payload.get("lcb_fraction")
    return ShipRecord(
        name=str(_require(payload, "name", where)),
        ship_class=str(payload.get("ship_class", "")),
        length=float(_require(payload, "length", where)),
        beam=float(_require(payload, "beam", where)),
        draft=float(_require(payload, "draft", where)),
        block_coefficient=float(_require(payload, "block_coefficient", where)),
        displacement_basis=str(payload.get("displacement_basis", "closure")),
        provenance=provenance,
        displacement_volume=None if volume is None else float(volume),
        lcb_fraction=None if lcb is None else float(lcb),
        midship_coefficient=None if payload.get("midship_coefficient") is None
        else float(payload["midship_coefficient"]),
        transom_ratio=float(payload.get("transom_ratio", 0.0)),
        geometry_route=str(payload.get("geometry_route", "series60_form")),
        notes=str(payload.get("notes", "")),
        geometry_sha256=payload.get("geometry_sha256"),
        geometry_path=payload.get("geometry_path"),
    )


def load_manifest(path: str | Path | None = None) -> FleetManifest:
    """Load and structurally validate the fleet manifest."""
    target = Path(path) if path is not None else MANIFEST_PATH
    raw = target.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    if payload.get("schema") != SCHEMA:
        raise ManifestError(f"{target}: schema {payload.get('schema')!r} != {SCHEMA!r}")
    ships_payload = _require(payload, "ships", str(target))
    if not isinstance(ships_payload, list) or not ships_payload:
        raise ManifestError(f"{target}: ships must be a non-empty list")
    ships = tuple(
        record_from_dict(entry, where=f"{target}:ships[{index}]")
        for index, entry in enumerate(ships_payload)
    )
    manifest = FleetManifest(
        ships=ships,
        placeholder=dict(payload.get("placeholder_signature", {})),
        path=target,
        sha256=sha256_bytes(raw),
    )
    issues = manifest.problems()
    if issues:
        raise ManifestError("; ".join(issues))
    return manifest


def fleet_sha_index(manifest: FleetManifest) -> dict[str, Any]:
    """Provenance block for a run: the manifest identity plus per-ship sources."""
    return {
        "manifest_path": None if manifest.path is None else str(manifest.path),
        "manifest_sha256": manifest.sha256,
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "ships": {
            ship.name: {
                "geometry_sha256": ship.geometry_sha256,
                "reference_basis": ship.reference_basis(),
                "reference_is_external": ship.reference_is_external(),
                "provenance": ship.provenance.as_dict(),
            }
            for ship in manifest.ships
        },
    }


def iter_names(manifest: FleetManifest) -> Iterable[str]:
    return manifest.names()