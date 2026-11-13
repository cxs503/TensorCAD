import time, math, json
import tensorcad
from tensorcad import wigley, series60, geometry_field as gf
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.geometry_field import SphereField, EllipsoidField, to_mesh
import inspect

print("SphereField sig:", inspect.signature(SphereField.__init__))
print("EllipsoidField sig:", inspect.signature(EllipsoidField.__init__))
print("to_mesh sig:", inspect.signature(to_mesh))
print("generate_wigley sig:", inspect.signature(wigley.generate_wigley))
print("generate_series60 sig:", inspect.signature(series60.generate_series60))

t0=time.time()
w = wigley.generate_wigley()
print("wigley gen %.2fs" % (time.time()-t0), w.hydrostatics.displacement_volume, len(w.faces))
h=w.hydrostatics
print("wigley hydro keys:", h.to_dict().keys())
print("wigley analytic:", json.dumps(w.metadata().get("analytic"), indent=1)[:1200])

t0=time.time()
s = series60.generate_series60()
print("series60 gen %.2fs" % (time.time()-t0), s.hydrostatics.displacement_volume, len(s.faces))
print("series60 attrs:", [a for a in dir(s) if not a.startswith('_')][:30])
print("series60 analytic keys:", list(s.metadata().get("analytic", {}).keys()) if hasattr(s,'metadata') else "no meta")

t0=time.time()
for res in (24, 48, 96):
    m = to_mesh(SphereField(center=(0,0,0), radius=1.0), resolution=res)
    hh = compute_hydrostatics(m.mesh, waterline_z=-10.0, closed=True)
    print("sphere res",res,"V",hh.displacement_volume,"S",hh.wetted_surface_area,"%.2fs"%(time.time()-t0))
t0=time.time()
m = to_mesh(EllipsoidField(center=(0,0,0), radii=(2.0,1.0,0.5)), resolution=48)
hh = compute_hydrostatics(m.mesh, waterline_z=-10.0, closed=True)
print("ellipsoid", hh.displacement_volume, hh.wetted_surface_area, "%.2fs"%(time.time()-t0))

from tensorcad.dtmb5415 import DTMB5415Space, DRAFT_M, WATERLINE_Z
t0=time.time()
space = DTMB5415Space()
mesh = space.mesh()
print("5415 mesh %.2fs" % (time.time()-t0), len(mesh.faces))
h5 = compute_hydrostatics(mesh, waterline_z=WATERLINE_Z, closed=True)
print("5415 hydro:", h5.to_dict())
print("5415 space volume:", space.volume(), "hydro:", space.hydrostatics())