"""Diagnose how one 142 record decomposes: refs, in-domain flags, deviations."""
import sys

sys.path.insert(0, "/root/tp-agent/work/m22464")
from tensorcad.geomquality.iges_trim import parse_iges, _Builder, _model_extent  # noqa: E402
from tensorcad.geomquality import iges_trim as T  # noqa: E402


def bbox(points):
    return [(round(min(p[i] for p in points), 3), round(max(p[i] for p in points), 3))
            for i in range(3)]


path = sys.argv[1] if len(sys.argv) > 1 else "work/assets/kvlcc2.igs"
seq = int(sys.argv[2]) if len(sys.argv) > 2 else 35
document = parse_iges(path)
extent = _model_extent(document)
builder = _Builder(document, 1e-4 * extent, extent)
entry = document.entry(seq)
print("142 record", seq, "tokens", entry.tokens[:10], "extent", extent)
surface, kind = builder.surface(entry.integer(1))
print("surface", kind, "domain", surface.domain, "extent-tol", 1e-4 * extent)
refs = []
for index in range(len(entry.tokens)):
    token = entry.tokens[index].strip()
    if not token or any(c in token for c in ".eEdD"):
        continue
    value = entry.integer(index)
    if value <= 0 or value not in document.entries or value == entry.sequence:
        continue
    if document.entries[value].entity_type in (100, 102, 106, 110, 126):
        refs.append(value)
print("curve refs", refs)
for ref in refs:
    loop = T._loop_from_curve(builder, ref)
    print(f"  ref {ref} type {document.entries[ref].entity_type} points {len(loop)} bbox {bbox(loop)}")
    print(f"      in_domain={T._uv_in_domain(surface, loop)}")
    if T._uv_in_domain(surface, loop):
        uv = [(p[0], p[1]) for p in loop]
        print(f"      uv bbox {bbox([(u, v, 0.0) for u, v in uv])}")
        for other in refs:
            if other == ref:
                continue
            oloop = T._loop_from_curve(builder, other)
            print(f"      deviation vs {other}: {T._uv_deviation(surface, uv, oloop):.3e}")