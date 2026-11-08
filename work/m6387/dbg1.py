import math, sys
sys.path.insert(0, ".")
from tensorcad.planing import *
p = PlaningHullParameters(length=10.0, beam=3.0, draft=0.05, deadrise_deg=30.0, planform_exponent=0.05, nx=21, nz=9)
p.validate()
h = analytic_planing_hydrostatics(p)
print({k: h[k] for k in ("displacement_volume","waterplane_area","wetted_surface_area","wetted_length","transom_section_area","vertical_centre_buoyancy")})
print("x0 =", wetted_length(p))
print("section_area(L)=", section_area(p, p.length))
print("Awp=", analytic_waterplane_area(p))
print("topside=", _wetted_topside_area(p, wetted_length(p)))
for s in (400, 4000, 40000, 400000):
    print("samples", s, _wetted_topside_area(p, wetted_length(p), s))
