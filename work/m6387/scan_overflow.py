import math, sys, itertools
sys.path.insert(0, ".")
from tensorcad.planing import PlaningHullParameters, analytic_planing_hydrostatics, analytic_wetted_area, generate_planing_hull
bad=[]
n=0
for beta,q,T in itertools.product((1.0,6.0,12.0,18.0,30.0,45.0,59.0),
                                  (0.05,0.1,0.2,0.4,0.6,1.0,1.5,2.0,3.0,4.0),
                                  (0.05,0.2,0.35,0.55,1.0,2.0)):
    p=PlaningHullParameters(length=10.0,beam=3.0,draft=T,deadrise_deg=beta,planform_exponent=q,nx=21,nz=9)
    try: p.validate()
    except ValueError: continue
    n+=1
    h=analytic_planing_hydrostatics(p)
    nonfinite=[k for k,v in h.items() if isinstance(v,float) and not math.isfinite(v)]
    if nonfinite: bad.append((beta,q,T,nonfinite,h["wetted_surface_area"]))
print("cases",n,"non-finite",len(bad))
for b in bad[:20]: print(b)
