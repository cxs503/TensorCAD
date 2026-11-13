import numpy as np, scipy.io as sio, trimesh, time
from scipy.interpolate import CubicSpline

d = sio.loadmat('/root/ship-data-official/dtmb5415/DTMB5415-database.mat',
                variable_names=['xMatrix', 'yMatrix', 'zMatrix'])
P = np.stack([d['xMatrix'][:, 0], d['yMatrix'][:, 0], d['zMatrix'][:, 0]], 1)
blk = P.reshape(50, 45, 3)
fore = blk[0::2]   # (25 girth, 45 station, 3)
aft = blk[1::2]
NG, NS = 25, 45


def refine(grid, ns_out, ng_out):
    """Bicubic (natural) interpolation in (station, girth) index space -> (ng_out, ns_out, 3)."""
    ng, ns, _ = grid.shape
    si = np.arange(ns)
    sj = np.arange(ng)
    fs = np.arange(ns_out) * (ns - 1) / (ns_out - 1)
    fj = np.arange(ng_out) * (ng - 1) / (ng_out - 1)
    tmp = np.empty((ng, ns_out, 3))
    for j in range(ng):
        cs = CubicSpline(si, grid[j], axis=0, bc_type='natural')
        tmp[j] = cs(fs)
    out = np.empty((ng_out, ns_out, 3))
    for i in range(ns_out):
        cg = CubicSpline(sj, tmp[:, i, :], axis=0, bc_type='natural')
        out[:, i] = cg(fj)
    return out


def build(grids, ns, ng):
    """Two patches + bridge band + mirror -> (V, F)."""
    fg, ag = grids
    fore_idx = np.arange(ns * ng).reshape(ng, ns)
    aft_idx = np.arange(ns * ng, 2 * ns * ng).reshape(ng, ns)
    verts = np.vstack([fg.reshape(-1, 3), ag.reshape(-1, 3)])
    quads = []
    for patch in (fore_idx, aft_idx):
        for j in range(ng - 1):
            for i in range(ns - 1):
                quads.append((patch[j, i], patch[j, i + 1], patch[j + 1, i + 1], patch[j + 1, i]))
    for j in range(ng - 1):
        quads.append((fore_idx[j, -1], fore_idx[j + 1, -1], aft_idx[j + 1, 0], aft_idx[j, 0]))
    nv = len(verts)
    on_centre = np.abs(verts[:, 1]) < 1e-9
    mirror_map = np.full(nv, -1, dtype=int)
    new_verts = []
    for k in range(nv):
        if on_centre[k]:
            mirror_map[k] = k
        else:
            mirror_map[k] = nv + len(new_verts)
            new_verts.append([verts[k, 0], -verts[k, 1], verts[k, 2]])
    allv = np.vstack([verts, np.array(new_verts)]) if new_verts else verts
    tris = []
    for q in quads:
        tris.append((q[0], q[1], q[2]))
        tris.append((q[0], q[2], q[3]))
    allt = np.array(tris)
    mirrored = np.array([[mirror_map[a], mirror_map[c], mirror_map[b]] for a, b, c in tris])
    return allv, np.vstack([allt, mirrored])


t0 = time.time()
p = '/root/ship-performance-platform-incoming/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull_exact.stl'
raw = open(p, 'rb').read(); n = np.frombuffer(raw[80:84], dtype='<u4')[0]
rec = np.frombuffer(raw[84:84 + 50 * n], dtype=np.dtype([('n', '<3f4'), ('v', '<9f4'), ('a', '<u2')]))
tris_stl = rec['v'].reshape(n, 3, 3).astype(np.float64)
exact = trimesh.Trimesh(vertices=tris_stl.reshape(-1, 3), faces=np.arange(n * 3).reshape(n, 3), process=True)
print('exact stl watertight', exact.is_watertight, 'verts', len(exact.vertices), 'bbox', exact.bounds.tolist())

for ns_out, ng_out, tag in ((NS, NG, 'base'), (133, 49, 'refine3x2')):
    fg = refine(fore, ns_out, ng_out) if ns_out != NS else fore.copy()
    ag = refine(aft, ns_out, ng_out) if ns_out != NS else aft.copy()
    V, F = build((fg, ag), ns_out, ng_out)
    print(tag, 'verts', len(V), 'faces', len(F), 'seam x: fore', round(fg[0, 0, 0], 2), round(fg[0, -1, 0], 2),
          'aft', round(ag[0, 0, 0], 2), round(ag[0, -1, 0], 2))
    m = trimesh.Trimesh(vertices=V, faces=F, process=False)
    el = m.edges_unique_length
    print('  edge median %.3f max %.3f; watertight %s; euler %d; boundary edges %d' %
          (np.median(el), el.max(), m.is_watertight, m.euler_number, (np.bincount(np.asarray(m.edges_unique_inverse).ravel()) == 1).sum() if False else -1))
    # point-to-mesh distance: KD-tree over triangle centroids + exact point-triangle
    from scipy.spatial import cKDTree
    tri = V[F]
    cen = tri.mean(1)
    tree = cKDTree(cen)
    kq = min(12, len(cen))
    _, idxs = tree.query(exact.vertices, k=kq)

    def pt_tri(p, a, b, c):
        ab, ac, ap = b - a, c - a, p - a
        d1 = (ab * ap).sum(-1); d2 = (ac * ap).sum(-1)
        bp = p - b; d3 = (ab * bp).sum(-1); d4 = (ac * bp).sum(-1)
        cp = p - c; d5 = (ab * cp).sum(-1); d6 = (ac * cp).sum(-1)
        vc = d1 * d4 - d3 * d2; vb = d5 * d2 - d1 * d6; va = d3 * d6 - d5 * d4
        denom = va + vb + vc
        out = np.full(len(p), np.inf)
        with np.errstate(divide='ignore', invalid='ignore'):
            v = vb / denom; w = vc / denom
            out = np.where(np.isfinite(out), out, np.inf)
            cand = a + v[..., None] * ab + w[..., None] * ac
        # regions
        e = np.where((d1 <= 0) & (d2 <= 0), np.linalg.norm(ap, axis=-1), np.inf)
        e = np.minimum(e, np.where((d3 >= 0) & (d4 <= d3), np.linalg.norm(bp, axis=-1), np.inf))
        e = np.minimum(e, np.where((d6 >= 0) & (d5 <= d6), np.linalg.norm(cp, axis=-1), np.inf))
        inside = (va >= 0) & (vb >= 0) & (vc >= 0)
        e = np.minimum(e, np.where(inside, np.linalg.norm(p - cand, axis=-1), np.inf))
        # edge AB and AC and BC
        def seg(p, a, b):
            ab = b - a; t = ((p - a) * ab).sum(-1) / (ab * ab).sum(-1)
            t = np.clip(t, 0, 1)
            return np.linalg.norm(p - (a + t[..., None] * ab), axis=-1)
        e = np.minimum(e, seg(p, a, b)); e = np.minimum(e, seg(p, a, c)); e = np.minimum(e, seg(p, b, c))
        return e

    pts = exact.vertices
    d_all = np.full(len(pts), np.inf)
    for kk in range(kq):
        fi = idxs[:, kk]
        a, b, c = tri[fi, 0], tri[fi, 1], tri[fi, 2]
        d_all = np.minimum(d_all, pt_tri(pts, a, b, c))
    dist = d_all
    res = float(np.median(el))
    print('  dist exact->mesh: median %.4f p95 %.4f max %.4f ; frac<res %.4f  frac<res/2 %.4f' %
          (np.median(dist), np.percentile(dist, 95), dist.max(), (dist < res).mean(), (dist < res / 2).mean()))
print('elapsed', round(time.time() - t0, 1))