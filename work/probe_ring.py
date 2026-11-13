import math
from tensorcad import generate_hull, generate_wigley
from tensorcad.trim import hydrostatic_particulars

for name, hull in (("hull", generate_hull()), ("wigley", generate_wigley())):
    zs = sorted({round(v[2], 12) for v in hull.mesh.vertices})
    keel = min(zs)
    print("==", name, "keel", keel, "levels", len(zs))
    for z in zs:
        if z <= keel + 1e-12:
            continue
        row = hydrostatic_particulars(hull.mesh, z, station_count=21)
        print(f"z={z:+.6f} A={row.waterplane_area:.10f} V={row.volume:.8f} tpc={row.tpc:.6f} lcb={row.lcb:.6f} nudge={row.waterline_nudged}")
        for dz in (-1e-7, 1e-7):
            r = hydrostatic_particulars(hull.mesh, z + dz, station_count=21)
            dA = abs(r.waterplane_area - row.waterplane_area) / max(row.waterplane_area, 1e-30)
            dV = abs(r.volume - row.volume) / max(row.volume, 1e-30)
            dt = abs(r.tpc - row.tpc) / max(row.tpc, 1e-30)
            print(f"    dz={dz:+.0e} dA/A={dA:.2e} dV/V={dV:.2e} dTPC={dt:.2e} nudge={r.waterline_nudged}")
