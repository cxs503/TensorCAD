import math
from tensorcad import generate_hull
from tensorcad.trim import (hydrostatic_particulars, waterline_for_displacement,
                            trimmed_hydrostatics, equilibrium_trim, TrimError)

hull = generate_hull()
T = hull.parameters.draft
for frac in (0.25, 0.5, 0.75):
    z = -frac * T
    row = hydrostatic_particulars(hull.mesh, z, station_count=41)
    target = 0.9 * row.volume
    print(f"target from ring z={z:+.4f}: V={row.volume:.6f} A={row.waterplane_area:.6f}")
    try:
        h = waterline_for_displacement(hull.mesh, trim_angle=0.0, displaced_volume=target, station_count=41)
        back = hydrostatic_particulars(hull.mesh, h, station_count=41)
        print(f"   solved h={h:+.10f} V={back.volume:.10f} relerr={abs(back.volume-target)/target:.2e} nudge={back.waterline_nudged}")
    except TrimError as exc:
        print("   TrimError:", exc)

# trimmed hydrostatics preserving a ring-level displacement
z = -0.5 * T
up = hydrostatic_particulars(hull.mesh, z, station_count=41)
for deg in (0.1, 0.5, 1.0):
    try:
        t = trimmed_hydrostatics(hull.mesh, trim_angle=math.radians(deg), preserved=up.volume, station_count=41)
        print(f"trim {deg}deg: V={t.volume:.8f} rel={abs(t.volume-up.volume)/up.volume:.2e} A={t.waterplane_area:.6f} nudge={t.waterline_nudged} lcb={t.lcb:.5f}")
    except TrimError as exc:
        print(f"trim {deg}deg TrimError:", exc)

sol = equilibrium_trim(hull.mesh, weight=300.0, shift=1.0, waterline_z=z, station_count=41)
print("equilibrium at ring waterline:", sol)
