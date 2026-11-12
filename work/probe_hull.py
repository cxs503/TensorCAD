"""Probe: hull plating solid for mid-surface extraction (throw-away exploration)."""
from __future__ import annotations

import math
import time

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.midsurface import extract_midsurface, mid_surface_between, defeature_planar_regions
from tensorcad.series60 import Series60Parameters, generate_series60, analytic_series60_hydrostatics
from tensorcad.shell import offset_inward

p = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.60, nx=41, nz=13)
t0 = time.time()
model = generate_series60(p)
print("hull build", round(time.time() - t0, 2), "faces", len(model.faces))
hull = model.mesh
print("watertight", meshlib.is_watertight(hull), "volume", meshlib.total_signed_volume(hull))
ref = analytic_series60_hydrostatics(p)
cb_mesh = meshlib.total_signed_volume(hull) / (p.length * p.beam * p.draft)
print("Cb mesh", cb_mesh, "analytic", ref["block_coefficient"],
      "rel", cb_mesh / ref["block_coefficient"] - 1)
box = hull.bounding_box()
print("bbox", box.minimum, box.maximum)

# deck cap = faces whose vertices are all at z ~ 0
def cap_face(face):
    return all(abs(hull.vertices[i][2]) < 1e-12 for i in face)

cap = [i for i, f in enumerate(hull.faces) if cap_face(f)]
print("cap faces", len(cap), "of", len(hull.faces))

t = 0.02
t0 = time.time()
inner_full = offset_inward(hull, t)
print("offset", round(time.time() - t0, 2))

wetted = [i for i in range(len(hull.faces)) if i not in set(cap)]
outer = Mesh(vertices=hull.vertices, faces=tuple(hull.faces[i] for i in wetted))
print("wetted faces", len(outer.faces), "boundary edges", len(meshlib.boundary_edges(outer)))
inner_inner = Mesh(vertices=inner_full.vertices, faces=tuple(hull.faces[i] for i in wetted))

vertices = list(outer.vertices) + list(inner_inner.vertices)
n_out = len(outer.vertices)
faces = list(outer.faces)
faces.extend(tuple(i + n_out for i in f)[::-1] for f in inner_inner.faces)
solid = Mesh(vertices=tuple(vertices), faces=tuple(faces))
print("plating faces", len(solid.faces), "watertight", meshlib.is_watertight(solid),
      "area", meshlib.surface_area(solid))

t0 = time.time()
result = extract_midsurface(solid, max_thickness=0.5)
print("extract", round(time.time() - t0, 2), "mid faces", len(result.mesh.faces),
      "area", result.area, "hull area", meshlib.surface_area(outer),
      "rel", result.area / meshlib.surface_area(outer) - 1)
print("paired frac", result.paired_fraction, "unpaired", len(result.unpaired),
      "thk mean", result.mean_thickness, "min", result.min_measured, "max", result.max_measured)

# independent route: index-aligned averaging of the two skins
mid_ref = mid_surface_between(outer, inner_inner)
print("index-aligned mid area", meshlib.surface_area(mid_ref),
      "rel vs extracted", result.area / meshlib.surface_area(mid_ref) - 1)
# max distance from extracted mid vertices to the reference mid surface (nearest vertex distance)
import itertools
best = 0.0
ref_points = mid_ref.vertices
for v in result.mesh.vertices:
    d = min(math.dist(v, w) for w in ref_points)
    best = max(best, d)
print("max vertex-to-reference distance", best, "fraction of L", best / p.length)

t0 = time.time()
d = defeature_planar_regions(result.mesh)
print("defeature", round(time.time() - t0, 2), d.face_count_before, d.face_count_after,
      d.face_reduction, "skipped", len(d.skipped_regions), "area err", d.area_relative_error)