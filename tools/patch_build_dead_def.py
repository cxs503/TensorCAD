"""One-shot patch: drop the shadowed duplicate build_ship_geometry and add an STL writer.

Run once from the repository root; it is kept for reproducibility of the change.
"""

from pathlib import Path

TARGET = Path("tensorcad/fleet/build.py")
text = TARGET.read_text()

marker_first = "def build_ship_geometry(\n    record: ShipRecord,\n    *,\n    nx: int = 121,\n    nz: int = 25,\n    stl_path: str | Path | None = None,\n) -> BuiltGeometry:"
marker_iges = "# --------------------------------------------------------------------------- #\n# IGES field route (fourth route)"

start = text.index(marker_first)
end = text.index(marker_iges)
removed = text[start:end]
assert "unknown geometry route" in removed, "unexpected block boundary"
assert removed.count("def build_ship_geometry") == 1
text = text[:start] + text[end:]

WRITER = '''
def write_binary_stl(mesh: Mesh, path: str | Path, *, header: str = "tensorcad fleet artifact") -> str:
    """Write ``mesh`` as a binary STL and return the sha256 of the bytes written.

    The delivered artefact, not the in-memory mesh, is what downstream consumers
    (and the manifest) have to be able to audit, so the writer pins the file's
    digest at the moment of writing. Coordinates are stored as float32, which is
    the format's own precision; :func:`read_stl` therefore returns a mesh whose
    volume differs from the input by float32 rounding only, and the fleet run
    re-reads every written file and re-measures it rather than trusting this copy.
    """
    payload = bytearray(header.encode("ascii", "replace")[:80].ljust(80, b"\\0"))
    payload += struct.pack("<I", len(mesh.faces))
    for a, b, c in mesh.faces:
        pa, pb, pc = mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        nx_, ny_, nz_ = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        norm = math.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
        if norm > 0.0:
            nx_, ny_, nz_ = nx_ / norm, ny_ / norm, nz_ / norm
        payload += struct.pack("<3f", nx_, ny_, nz_)
        for point in (pa, pb, pc):
            payload += struct.pack("<3f", point[0], point[1], point[2])
        payload += b"\\0\\0"
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(bytes(payload))
    return hashlib.sha256(bytes(payload)).hexdigest()


'''

anchor = "def build_ship_geometry(\n    record: ShipRecord,"
text = text.replace(anchor, WRITER.lstrip("\n") + anchor, 1)
if "import hashlib" not in text:
    text = text.replace("import math\nimport struct", "import hashlib\nimport math\nimport struct", 1)
TARGET.write_text(text)
print("patched", TARGET, "removed", removed.count("\\n"), "lines")