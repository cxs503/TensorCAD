import cProfile, pstats, signal, sys, io, time
sys.path.insert(0, '/srv/tensor/agent/work/m22464')
from tensorcad.geomquality.trim_mesh import import_cad

class TO(Exception):
    pass

def handler(s, f):
    raise TO()

path = sys.argv[1]
budget = int(sys.argv[2]) if len(sys.argv) > 2 else 180
signal.signal(signal.SIGALRM, handler)
signal.alarm(budget)
prof = cProfile.Profile()
t0 = time.time()
prof.enable()
status = "done"
try:
    r = import_cad(path, weld_tolerance=1e-3, seam_tolerance=1e-3)
    print("done faces=%d verts=%d t=%.1fs" % (len(r.mesh.faces), len(r.mesh.vertices), time.time() - t0))
except TO:
    status = "interrupted"
    print("interrupted at %.1fs" % (time.time() - t0))
finally:
    prof.disable()
    signal.alarm(0)
    s = io.StringIO()
    pstats.Stats(prof, stream=s).sort_stats('cumulative').print_stats(25)
    print(s.getvalue())