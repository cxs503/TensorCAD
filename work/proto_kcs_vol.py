"""Prototype: tessellate the official KCS IGES (55 untrimmed type-128 patches),
close at the waterplane, measure displacement volume vs published 52,030 m3."""
import math, sys, time
from tensorcad.iges import read_iges

PATH = '/root/sim-assets/m21253-379b12515423/artifacts/baseline_fleet/raw/closure/kcs/KCS_hull_SVA.igs'
NU = NV = int(sys.argv[1]) if len(sys.argv) > 1 else 16
T = 10.8


def tess(f, nu, nv):
    tris = []
    for s in f.surfaces():
        grid = s.sample(nu, nv)          # row-major, v outer, u inner
        for j in range(nv - 1):
            for i in range(nu - 1):
                a = grid[j * nu + i]
                b = grid[j * nu + i + 1]
                c = grid[(j + 1) * nu + i + 1]
                d = grid[(j + 1) * nu + i]
                for t in ((a, b, c), (a, c, d)):
                    if area(t) > 1e-12:
                        tris.append(t)
    return tris


def area(t):
    a, b, c = t
    ux, uy, uz = b[0]-a[0], b[1]-a[1], b[2]-a[2]
    vx, vy, vz = c[0]-a[0], c[1]-a[1], c[2]-a[2]
    cx, cy, cz = uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx
    return 0.5*math.sqrt(cx*cx+cy*cy+cz*cz)


def tet(t):
    a, b, c = t
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) + a[1]*(b[2]*c[0]-b[0]*c[2]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0


def mirror(t):
    out = []
    for tri in t:
        m = tuple((p[0], -p[1], p[2]) for p in tri)
        out.append((m[0], m[2], m[1]))
    return out


def clip(tris, zmax):
    """Keep z <= zmax, return clipped triangles."""
    out = []
    for tri in tris:
        poly = list(tri)
        new = []
        for k in range(3):
            p, q = poly[k], poly[(k+1) % 3]
            inp, inq = p[2] <= zmax, q[2] <= zmax
            if inp:
                new.append(p)
            if inp != inq:
                s = (zmax - p[2]) / (q[2] - p[2])
                new.append((p[0]+s*(q[0]-p[0]), p[1]+s*(q[1]-p[1]), zmax))
        for k in range(1, len(new)-1):
            out.append((new[0], new[k], new[k+1]))
    return out


f = read_iges(PATH)
t0 = time.time()
tris = tess(f, NU, NV)
print(f"patches=55 grid={NU}x{NV} tris={len(tris)} in {time.time()-t0:.1f}s")
full = tris + mirror(tris)
print("mirrored tris:", len(full))
v_signed = sum(tet(t) for t in full)
print("signed volume of full shell (open at deck):", round(v_signed, 2))
if v_signed < 0:
    full = [(t[0], t[2], t[1]) for t in full]
    print("flipped to outward")

clo = clip(full, T)
print("clipped tris:", len(clo))
# boundary edges
from collections import Counter
def key(p): return (round(p[0], 6), round(p[1], 6), round(p[2], 6))
edges = Counter()
for t in clo:
    for k in range(3):
        a, b = key(t[k]), key(t[(k+1) % 3])
        edges[(a, b)] += 1
bnd = [(a, b) for (a, b), n in edges.items() if n == 1]
print("boundary edges:", len(bnd), "max |z-T| on boundary:",
      max(abs(a[2]-T) for a, b in bnd) if bnd else None)
# loops
adj = {}
for a, b in bnd:
    adj.setdefault(a, []).append(b)
loops, seen = [], set()
for start in adj:
    if start in seen: continue
    loop = [start]; seen.add(start); cur = start
    while True:
        nxt = [n for n in adj[cur] if n not in seen]
        if not nxt: break
        cur = nxt[0]; loop.append(cur); seen.add(cur)
    loops.append(loop)
print("open loops:", [len(l) for l in loops])
# waterplane area by shoelace on the closed loop (x,y)
def shoelace(loop):
    s = 0.0
    for i in range(len(loop)):
        x1, y1 = loop[i][0], loop[i][1]
        x2, y2 = loop[(i+1) % len(loop)][0], loop[(i+1) % len(loop)][1]
        s += x1*y2 - x2*y1
    return s/2.0
Awp_signed = sum(shoelace(l) for l in loops if len(l) > 2)
Awp = abs(Awp_signed)
print("waterplane area:", round(Awp, 2))
V = sum(tet(t) for t in clo) + Awp*T/3.0
print("V below T:", round(V, 1), "published 52030 -> err %", round(100*(V-52030)/52030, 3))
xs = [p[0] for t in clo for p in t]; ys = [p[1] for t in clo for p in t]
print("bbox clipped: x", round(min(xs), 2), round(max(xs), 2), "y", round(min(ys), 2), round(max(ys), 2))
print("wetted surface:", round(sum(area(t) for t in clo), 2))
# centroid (LCB) via divergence: need full closed body; use half+mirror already
lcb_num = sum((t[0][0]+t[1][0]+t[2][0])/4.0 * tet(t) for t in clo)
print("sum x*tet:", round(lcb_num, 1))