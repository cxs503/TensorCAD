import math, sys
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 30
from tensorcad.planing import *
def route_mp(p):
    L,T,tau = mp.mpf(p.length), mp.mpf(p.draft), mp.mpf(p.tan_beta)
    q = mp.mpf(p.planform_exponent); c = mp.mpf(p.half_beam_transom)/L**q
    def b(x): return c*x**q if x>0 else mp.mpf(0)
    def area(x):
        bb=b(x); h=bb*tau
        return T*T/tau if h>=T else bb*h + 2*bb*(T-h)
    def mom(x):
        bb=b(x); h=bb*tau
        return -T**3/(3*tau) if h>=T else -T*T*bb + T*tau*bb*bb - tau*tau*bb**3/3
    x0 = mp.mpf(wetted_length(p))
    vol = mp.quad(area, [0, L])
    awp = mp.quad(lambda x: 2*min(b(x), T/tau), [0, L])
    lcb = mp.quad(lambda x: x*area(x), [0, L])/vol
    lcf = mp.quad(lambda x: x*min(b(x), T/tau), [0, L])/(awp/2)
    vcb = mp.quad(mom, [0, L])/vol
    bottom = awp*mp.sqrt(1+tau*tau)
    # topside in t substitution x = x0 t**(2/q) : independent variable change,
    # but evaluate with mpmath tanh-sinh in x on [0,x0] as well
    def topx(x):
        bb=b(x); h=T-bb*tau
        if h<=0: return mp.mpf(0)
        db = c*q*x**(q-1) if x>0 else mp.mpf(0)
        return 2*h*mp.sqrt(1+db*db)
    top = mp.quad(topx, [0, x0]) if x0>0 else mp.mpf(0)
    return dict(displacement_volume=vol, waterplane_area=awp, longitudinal_centre_buoyancy=lcb,
                longitudinal_centre_flotation=lcf, vertical_centre_buoyancy=vcb,
                wetted_surface_area=bottom+top+area(L), _top=top)
import itertools
worst={}
for beta,q,T in itertools.product((1.0,12.0,18.0,30.0,45.0,59.0),(0.05,0.2,0.6,1.0,2.0,4.0),(0.05,0.30,0.55,2.0)):
    p=PlaningHullParameters(length=10.0,beam=3.0,draft=T,deadrise_deg=beta,planform_exponent=q,nx=21,nz=9)
    try: p.validate()
    except ValueError: continue
    a=analytic_planing_hydrostatics(p); r=route_mp(p)
    for k in ("displacement_volume","waterplane_area","longitudinal_centre_buoyancy","longitudinal_centre_flotation","vertical_centre_buoyancy","wetted_surface_area"):
        d=float((mp.mpf(a[k])-r[k])/r[k])
        if k not in worst or abs(d)>abs(worst[k][0]): worst[k]=(d,beta,q,T,float(r[k]))
for k,v in sorted(worst.items()):
    print(f"{k:32s} {v[0]:+.3e} at beta={v[1]} q={v[2]} T={v[3]} ref={v[4]:.6g}")
