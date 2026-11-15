"""Read-only probe of the official KCS IGES (KCS_hull_SVA.igs) structure."""
import sys, json
from tensorcad.iges import read_iges

path = sys.argv[1]
f = read_iges(path)
print("unit scale (m/unit):", f.metres_per_unit)
print("stats:", json.dumps(f.stats(), indent=1))
print("census:", json.dumps(f.surface_census(), indent=1)[:1500])

surfs = f.surfaces()
print("n surfaces:", len(surfs))
xs = []
for i, s in enumerate(surfs):
    bb = s.bounding_box()
    xs.append((i, s.domain_u, s.domain_v, bb))
    if i < 8 or i >= len(surfs) - 3:
        print(i, "domU", tuple(round(v, 4) for v in s.domain_u),
              "domV", tuple(round(v, 4) for v in s.domain_v),
              "bbox", tuple(tuple(round(c, 4) for c in p) for p in bb),
              "nu,nv", (len(s.control_u) if hasattr(s, "control_u") else "?"))
allb = [b for _, _, _, b in xs]
lo = [min(p[0] for _, _, _, b in xs for p in b), min(p[1] for _, _, _, b in xs for p in b), min(p[2] for _, _, _, b in xs for p in b)]
hi = [max(p[0] for _, _, _, b in xs for p in b), max(p[1] for _, _, _, b in xs for p in b), max(p[2] for _, _, _, b in xs for p in b)]
print("global bbox", lo, hi)
print("dir entries:", [(e.sequence, e.entity_type, e.form, e.blanked()) for e in f.of_type(128)][:10] if hasattr(f, "of_type") else "n/a")
print("all entity types:", json.dumps(f.stats()))