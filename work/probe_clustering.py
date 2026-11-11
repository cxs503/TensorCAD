"""Probe 2: can u-clustering alone fix the SUBOFF cap fans?"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.parametric_hull import (  # noqa: E402
    loft_nurbs_surface,
    suboff_cross_sections,
    suboff_station_fractions,
    tessellate_surface,
    mesh_quality_report,
)

stations = suboff_station_fractions(41)
sections = suboff_cross_sections(n_points=32, station_fractions=stations)
surface = loft_nurbs_surface(sections)

for spacing in ("uniform_x", "cosine"):
    for n_u in (160, 320, 640):
        m = tessellate_surface(surface, n_u=n_u, n_v=16, u_spacing=spacing)
        q = mesh_quality_report(m.mesh, n_faces=m.hull_face_count)
        print("%-10s n_u=%-4d faces=%-6d minJ=%.4f maxAsp=%6.2f meanAsp=%.2f watertight=%s"
              % (spacing, n_u, m.hull_face_count, q["min_scaled_jacobian"],
                 q["max_aspect_ratio"], q["mean_aspect_ratio"], q["watertight"]))