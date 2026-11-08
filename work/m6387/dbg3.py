import math, sys
sys.path.insert(0,".")
import mpmath as mp
from tensorcad.planing import *
p = PlaningHullParameters(length=10.0, beam=3.0, draft=2.0, deadrise_deg=45.0, planform_exponent=0.05, nx=21, nz=9)
p.validate()
x0 = wetted_length(p)
tau = p.tan_beta
bottom = analytic_waterplane_area(p)*math.sqrt(1+tau*tau)
transom = section_area(p,p.length)
ship_total = analytic_wetted_area(p)
print("x0", x0, "bottom", bottom, "transom", transom, "total", ship_total, "=> top", ship_total-bottom-transom)

def top_of(x):
    b = chine_half_beam(p,x); h = p.draft - b*tau
    if h<=0: return 0.0
    slope = (p.half_beam_transom/p.length**p.planform_exponent)*p.planform_exponent*x**(p.planform_exponent-1.0)
    return 2*h*math.sqrt(1+slope*slope)
import scipy.integrate as si
val,err = si.quad(top_of, 0.0, x0, limit=500, points=[1e-8,1e-4,1e-2,0.1,1.0])
print("scipy top", val, "err", err)

mp.mp.dps=40
def mptop(x):
    L=mp.mpf(p.length); T=mp.mpf(p.draft); t=mp.mpf(tau); q=mp.mpf(p.planform_exponent)
    c=mp.mpf(p.half_beam_transom)/L**q
    b=c*x**q if x>0 else mp.mpf(0); h=T-b*t
    if h<=0: return mp.mpf(0)
    slope=c*q*x**(q-1) if x>0 else mp.mpf(0)
    return 2*h*mp.sqrt(1+slope*slope)
for pts in ([0,x0], [0,1e-12,1e-8,1e-4,1e-2,0.5,2.0,5.0,9.0,x0]):
    print("mp pts", len(pts), mp.nstr(mp.quad(mptop, pts), 20))
# in t variable with e=40 substitution (analytic identity, cross-check of the same math)
def mptop_t(tt):
    L=mp.mpf(p.length); T=mp.mpf(p.draft); t=mp.mpf(tau); q=mp.mpf(p.planform_exponent)
    c=mp.mpf(p.half_beam_transom)/L**q
    X = x0*tt**40
    b = c*X**q; h=T-b*t
    dx = x0*40*tt**39
    db = c*q*X**(q-1)*dx
    return 2*h*mp.sqrt(dx*dx+db*db)
print("mp t-var", mp.nstr(mp.quad(mptop_t, [0,1]), 20))
