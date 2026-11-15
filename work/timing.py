import sys, time
sys.path.insert(0, '/srv/tensor/agent/work/m22464')
from tensorcad.geomquality.trim_mesh import import_cad
for f in sys.argv[1:]:
    t = time.time()
    r = import_cad(f, weld_tolerance=1e-3, seam_tolerance=1e-3)
    rep = r.report
    g = {k: v['passed'] for k, v in rep['gates'].items()}
    print(f"{f}: faces={len(r.mesh.faces)} verts={len(r.mesh.vertices)} t={time.time()-t:.1f}s refused={len(rep['import'].get('refused', {}))} gates={g}", flush=True)
