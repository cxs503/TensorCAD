"""One-shot patch: wire the source registry and the artefact re-ingest into the fleet run.

Three changes, each tied to an acceptance criterion of GAP-CAD-T1-BASELINE-FLEET:

1. ``G0b_source_anchor`` (gates) — a manifest digit with no named source fails, so the
   "source URL + licence per ship" requirement is gated instead of narrated.
2. ``write_binary_stl`` + re-ingest (run) — every admitted hull is written to
   ``manifests/<ship>.stl``, the file's sha256 is recorded, and the file is read back
   and re-gated, so the manifest's sha256 is the digest of the geometry actually
   delivered rather than of a mesh that only exists in memory.
3. ``verification_conclusion`` (run) — one machine-composed sentence per ship stating
   the gates passed, the volume error, the reference basis and its provenance state.
"""

from pathlib import Path

GATES = Path("tensorcad/fleet/gates.py")
text = GATES.read_text()

if "G0b_source_anchor" not in text:
    text = text.replace(
        "from tensorcad.fleet.manifest import DERIVED_BASES, FleetManifest, ShipRecord",
        "from tensorcad.fleet.manifest import DERIVED_BASES, FleetManifest, ShipRecord\n"
        "from tensorcad.fleet.sources import SourceRegistry",
        1,
    )
    GATE = '''
def gate_source_anchor(record: ShipRecord, registry: SourceRegistry) -> GateResult:
    """The record must name a source anchor: citation, URL, licence and locator.

    ``G0_provenance`` checks that the record states *a* URL and licence. This gate checks
    that the fleet's source registry holds the citation the digits would be re-derived
    from, and reports the registry's own sha256 so the anchor set is itself auditable.
    """
    anchor = registry.for_ship(record.name)
    if anchor is None:
        return GateResult(
            "G0b_source_anchor",
            False,
            f"no source anchor registered for {record.name!r}: a manifest digit with no named source "
            f"cannot be re-derived by a reader",
            {"registry_sha256": registry.sha256, "anchor_id": None},
        )
    covered = sorted(set(anchor.anchors))
    return GateResult(
        "G0b_source_anchor",
        True,
        f"anchor {anchor.anchor_id!r}: {anchor.citation[:64]!r}; locator {anchor.locator[:48]!r}; "
        f"licence {anchor.license!r}; retrieval={anchor.retrieval}; verified={anchor.verified}; "
        f"anchors={covered}",
        {
            "registry_sha256": registry.sha256,
            "anchor_id": anchor.anchor_id,
            "citation": anchor.citation,
            "url": anchor.url,
            "license": anchor.license,
            "locator": anchor.locator,
            "anchors": covered,
            "retrieval": anchor.retrieval,
            "verified": anchor.verified,
        },
    )


'''
    text = text.replace("def gate_reference_consistency(", GATE.lstrip("\n") + "def gate_reference_consistency(", 1)
    text = text.replace(
        """def evaluate_gates(
    record: ShipRecord,
    mesh: Mesh,
    manifest: FleetManifest,
    *,
    volume_tolerance: float = VOLUME_TOLERANCE,
) -> dict[str, Any]:""",
        """def evaluate_gates(
    record: ShipRecord,
    mesh: Mesh,
    manifest: FleetManifest,
    *,
    volume_tolerance: float = VOLUME_TOLERANCE,
    registry: SourceRegistry | None = None,
) -> dict[str, Any]:""",
        1,
    )
    text = text.replace(
        "        gate_principal_dimensions(record, facts, volume_tolerance),\n    ]",
        "        gate_principal_dimensions(record, facts, volume_tolerance),\n    ]\n"
        "    if registry is not None:\n"
        "        gates.insert(1, gate_source_anchor(record, registry))",
        1,
    )
GATES.write_text(text)

RUN = Path("tensorcad/fleet/run.py")
text = RUN.read_text()

text = text.replace(
    "from tensorcad.fleet.manifest import FleetManifest, fleet_sha_index, load_manifest, sha256_file",
    "from tensorcad.fleet.manifest import FleetManifest, fleet_sha_index, load_manifest, sha256_file\n"
    "from tensorcad.fleet.sources import load_registry",
    1,
)
text = text.replace(
    """from tensorcad.fleet.build import (
    RouteRefusal,
    build_ship_geometry,
    loft_station_table,
    sample_offset_table,
)""",
    """from tensorcad.fleet.build import (
    RouteRefusal,
    build_ship_geometry,
    loft_station_table,
    read_stl,
    sample_offset_table,
    write_binary_stl,
)""",
    1,
)
text = text.replace(
    """from tensorcad.fleet.gates import (
    VOLUME_TOLERANCE,
    evaluate_gates,
    gate_particulars_consistent,""",
    """from tensorcad.fleet.gates import (
    VOLUME_TOLERANCE,
    evaluate_gates,
    gate_particulars_consistent,
    gate_reference_consistency,""",
    1,
)

text = text.replace(
    """def _ship_block(record, manifest: FleetManifest, *, nx: int, nz: int,
                volume_tolerance: float, convergence: bool = True) -> dict[str, Any]:""",
    """def _verification_conclusion(verdict: dict[str, Any], artifact: dict[str, Any] | None) -> str:
    \"\"\"One sentence a reviewer can check the manifest against, composed from the gates.\"\"\"
    topology = verdict.get("topology") or {}
    record = verdict.get("record", {})
    gate_names = [gate["name"] for gate in verdict.get("gates", [])]
    failed = verdict.get("failed_gates") or []
    if verdict.get("passed"):
        head = f"ADMITTED ({len(gate_names) - len(failed)}/{len(gate_names)} gates passed)"
    else:
        head = f"REFUSED/FAILED (failed: {', '.join(failed)})"
    volume = verdict.get("volume_relative_error")
    reference = record.get("reference_volume")
    volume_text = (
        f"V_mesh {topology.get('signed_volume', float('nan')):.6g} m^3 vs V_ref {reference:.6g} m^3 "
        f"= {volume * 100:.4f}% (basis {verdict.get('reference_basis')}, "
        f"external={'yes' if verdict.get('reference_is_external') else 'no'})"
        if volume is not None and reference is not None
        else "no volume comparison"
    )
    topology_text = (
        f"watertight (boundary edges {topology.get('boundary_edge_count')}), "
        f"manifold+oriented (non-manifold {topology.get('nonmanifold_edge_count')}, "
        f"inconsistent {topology.get('inconsistently_oriented_edge_count')}), "
        f"Euler {topology.get('euler_characteristic')}"
        if topology
        else "no geometry built"
    )
    artifact_text = (
        f"delivered artefact {artifact['file']} sha256 {artifact['sha256'][:16]}... "
        f"({artifact['bytes']} bytes, re-ingest {'agrees' if artifact.get('reingest_ok') else 'DISAGREES'})"
        if artifact
        else "no artefact written"
    )
    return f"{head}; {topology_text}; {volume_text}; {artifact_text}."


def _artifact_block(record, mesh, verdict: dict[str, Any], artifact_dir: Path | None) -> dict[str, Any] | None:
    \"\"\"Write the admitted mesh to disk, then read it back and re-measure it.

    The digest has to belong to the file a consumer receives, so the file is written
    first and gated second: the re-ingest repeats the topology and volume gates on the
    bytes that were actually stored, which also catches a writer that silently loses
    or merges triangles.
    \"\"\"
    if artifact_dir is None:
        return None
    path = artifact_dir / f"{record.name}.stl"
    written_sha = write_binary_stl(mesh, path)
    stored_sha = sha256_file(path)
    reread = read_stl(path)
    facts = measure(reread)
    reference = record.reference_volume()
    volume_rel = abs(facts.signed_volume - reference) / reference
    reingest_ok = (
        stored_sha == written_sha
        and facts.boundary_edge_count == 0
        and facts.nonmanifold_edge_count == 0
        and facts.inconsistently_oriented_edge_count == 0
        and abs(facts.euler_characteristic - 2) == 0
    )
    return {
        "file": path.name,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": stored_sha,
        "sha256_on_write": written_sha,
        "triangles": facts.triangle_count,
        "reingest_volume": facts.signed_volume,
        "reingest_volume_relative_error": volume_rel,
        "reingest_volume_agrees_with_memory": volume_rel < 1e-6,
        "reingest_ok": reingest_ok and volume_rel < 1e-6,
    }


def _ship_block(record, manifest: FleetManifest, *, nx: int, nz: int,
                volume_tolerance: float, convergence: bool = True,
                registry=None, artifact_dir: Path | None = None) -> dict[str, Any]:""",
    1,
)

text = text.replace(
    "    verdict = evaluate_gates(record, built.mesh, manifest, volume_tolerance=volume_tolerance)",
    "    verdict = evaluate_gates(record, built.mesh, manifest, volume_tolerance=volume_tolerance,\n"
    "                             registry=registry)",
    1,
)
text = text.replace(
    """    if convergence:
        verdict["convergence"] = _convergence(record, nx=nx, nz=nz)
    return verdict""",
    """    artifact = _artifact_block(record, built.mesh, verdict, artifact_dir)
    verdict["artifact"] = artifact
    verdict["source_anchor"] = None if registry is None else (
        registry.for_ship(record.name).as_dict() if registry.for_ship(record.name) else None
    )
    verdict["verification_conclusion"] = _verification_conclusion(verdict, artifact)
    if convergence:
        verdict["convergence"] = _convergence(record, nx=nx, nz=nz)
    return verdict""",
    1,
)

text = text.replace(
    """    manifest = load_manifest(manifest_path)
    root = Path(out_dir)
    ships = [
        _ship_block(record, manifest, nx=nx, nz=nz, volume_tolerance=volume_tolerance)
        for record in manifest.ships
    ]""",
    """    manifest = load_manifest(manifest_path)
    registry = load_registry()
    root = Path(out_dir)
    artifact_dir = root / "manifests" if write else None
    if artifact_dir is not None:
        artifact_dir.mkdir(parents=True, exist_ok=True)
    ships = [
        _ship_block(record, manifest, nx=nx, nz=nz, volume_tolerance=volume_tolerance,
                    registry=registry, artifact_dir=artifact_dir)
        for record in manifest.ships
    ]""",
    1,
)

text = text.replace(
    '        "provenance": fleet_sha_index(manifest),',
    '        "provenance": fleet_sha_index(manifest),\n'
    '        "source_registry": {\n'
    '            "path": str(registry.path) if registry.path else None,\n'
    '            "sha256": registry.sha256,\n'
    '            "anchors_registered": [anchor.ship for anchor in registry.sources],\n'
    '            "anchors_missing": registry.missing_for(ship.name for ship in manifest.ships),\n'
    '        },',
    1,
)

text = text.replace(
    '            "hydrostatics_solver_calls": len(measured),',
    '            "hydrostatics_solver_calls": len(measured),\n'
    '            "artifacts_written": sum(1 for block in measured if block.get("artifact")),\n'
    '            "artifacts_reingest_ok": all(\n'
    '                block["artifact"]["reingest_ok"] for block in measured if block.get("artifact")\n'
    '            ),\n'
    '            "source_anchors_complete": not registry.missing_for(ship.name for ship in manifest.ships),\n'
    '            "worst_reference_consistency_relative_error": max(\n'
    '                (\n'
    '                    abs(gate["values"].get("relative_error") or 0.0)\n'
    '                    for block in measured\n'
    '                    for gate in block["gates"]\n'
    '                    if gate["name"] == "G6b_reference_consistency"\n'
    '                ),\n'
    '                default=None,\n'
    '            ),',
    1,
)

text = text.replace(
    '                        "topology": block["topology"],\n                        "gates": block["gates"],',
    '                        "topology": block["topology"],\n'
    '                        "gates": block["gates"],\n'
    '                        "source_anchor": block.get("source_anchor"),\n'
    '                        "artifact": block.get("artifact"),\n'
    '                        "verification_conclusion": block.get("verification_conclusion"),',
    1,
)
RUN.write_text(text)
print("run patched:", "verification_conclusion" in RUN.read_text())