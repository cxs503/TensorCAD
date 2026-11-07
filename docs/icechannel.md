# The ice channel corridor: what the plane says and what is true

`tensorcad.icechannel` turns the campaign route from the de Long Islands to Pevek
into a corridor a DEM ice field or a FEM ice belt can be built in. It has shipped
inside the PC6 contract (`docs/pc6.md`) with a geometry test but without an
evidence page of its own, and its length test only says "plan length agrees with
the great-circle total to 0.2 %" without saying *why* there is anything to agree
at all. This page closes both gaps and records one defect found on the way.

## 1. The case

| item | value |
|---|---|
| route | de Long (Bennett Island) -> Pevek roadstead, East Siberian Sea |
| waypoints | 8 (0.1 deg input precision) |
| great-circle length | 1061.791 km (573.3 n mile) |
| corridor | 250 m wide x 11.0 m deep |
| frame | local azimuthal-equidistant plane, origin at waypoint 0 |
| use | DEM ice-field / FEM ice-belt domain around an icebreaker escort |

The numbers that matter downstream are the corridor's *true* length and *true*
width: a DEM ice field built on a corridor that is 0.5 % short or 0.5 % narrow
puts the hull and the ice in slightly the wrong place, silently. Both are
therefore measured here against spherical trigonometry, not against themselves.

## 2. Method

**Projection.** Waypoints go to a plane by
``x = d sin(azimuth)``, ``y = d cos(azimuth)`` with ``d`` the great-circle
distance from the origin (haversine on a sphere of R = 6371.0088 km) and the
azimuth the *initial* bearing. Distances from the origin are therefore exact by
construction; nothing else is.

**Closed form for the rest.** In central-angle polar coordinates about the
origin the plane metric is ``ds_plane^2 = R^2 dc^2 + R^2 c^2 dpsi^2`` while the
sphere's is ``ds_true^2 = R^2 dc^2 + R^2 sin^2(c) dpsi^2``. The map is diagonal
with radial scale exactly 1 and tangential scale ``k = c/sin c``, so a plane
vector of length ``L`` at angle ``beta`` to the radial direction is worth

    L_true = L * sqrt(cos^2(beta) + (sin(beta)/k)^2)

in true metres. That is `tensorcad.icechannel.plane_to_true_factor`, and it is
first order in ``L/R``: over a 250 m corridor the dropped term is 1.5e-9.

**Verification routes.** (i) The closed form is checked against the *direct
geodesic problem* — a 1 m true step from a point at 8 distances x 6 headings x 24
bearings, 1152 cases, plane displacement compared with the great-circle
length of the same step. (ii) The corridor solid is checked by topology
(watertight, edge-manifold, Euler characteristic, self-intersection) and by
volume against the analytic prism ``plan_length x width x depth``. (iii) The
width is re-expressed per station with the same factor, so the plane chord and
the true width are the same measurement in two metrics.

## 3. Results

### 3.1 The distortion law against spherical truth

| quantity | value |
|---|---|
| cases | 1152 |
| worst relative error, closed form vs direct geodesic | 1.252e-09 |
| worst case at (distance km, approach deg, bearing deg) | (1500.0, 330.0, 0) |

The closed form is not an approximation fitted to these points: it reproduces
them to 1.3e-09, which is the second-order remainder of a
1 m step at 1500 km.

### 3.2 Route length, leg by leg

| leg | geodesic [km] | plane [km] | plane - true [%] | cumulative start [km] | angle to radial [deg] | true/plane |
|---|---:|---:|---:|---:|---:|---:|
| de Long (Bennett Island) -> East Siberian Sea north | 187.387 | 187.387 | -0.0000 | 0.0 | 0.0 | 1.000000 |
| East Siberian Sea north -> Medvezhyi Islands north | 216.283 | 216.283 | +0.0000 | 187.4 | 0.0 | 1.000000 |
| Medvezhyi Islands north -> Kolyma drift | 210.884 | 210.884 | +0.0000 | 403.7 | 0.1 | 1.000000 |
| Kolyma drift -> Chaun approach | 185.318 | 185.320 | +0.0007 | 614.6 | 3.9 | 0.999993 |
| Chaun approach -> Cape Billings south | 138.062 | 138.126 | +0.0460 | 799.9 | 24.8 | 0.999540 |
| Cape Billings south -> Pevek roadstead | 77.925 | 78.006 | +0.1039 | 937.9 | 32.9 | 0.998962 |
| Pevek roadstead -> Pevek port | 45.932 | 46.034 | +0.2221 | 1015.9 | 47.8 | 0.997783 |

Totals: plane 1062.039 km vs geodesic
1061.791 km = +0.0233 %,
which is the number the old test bounded at 0.2 %. The per-leg error is not
uniform and not small at the far end (+0.222 % on
the last leg): it grows with distance from the origin and with the leg's angle to
the radial direction, exactly as ``k = c/sin c`` says. Applying the factor per leg
recovers the geodesic total to -4.206e-06 %.

### 3.3 The corridor as a solid

| quantity | value | reference | error |
|---|---:|---:|---:|
| volume | 2920607332.311 m3 | 2920607332.310 m3 (plan x width x depth) | 2.75e-13 |
| vertices / faces / edges | 32 / 60 / 90 | — | — |
| Euler characteristic | 2 | 2 | exact |
| watertight / edge-manifold | True / True | True | — |
| self-intersection pairs | 0 | 0 | — |
| surface area | 554.390 km2 | — | — |

### 3.4 Width in true metres

| quantity | minimum [m] | maximum [m] | mean [m] |
|---|---:|---:|---:|
| plane station chord | 250.000 | 254.599 | 251.118 |
| true width | 249.479 | 253.970 | 250.735 |

The nominal 250 m corridor is, in true metres, between
249.479 m and 253.970 m. The
minimum is the projection narrowing (0.279 %); the
maximum is the mitre at waypoint 4, where the station chord
across a turn is longer than the perpendicular width by construction — a widened
diagonal, not a widened corridor. Perpendicular width at every station is the
nominal 250 m to within the projection factor.

### 3.5 The two defects this page found

Both are the same class: route guards that were not there, so the module happily
built a corridor that closed topologically and pierced itself.

**(a) A repeated waypoint.** A duplicated coordinate is real route data. It used to
build a solid that passed every topological test and folded through itself:
watertight True, edge-manifold True, Euler
2, and 10
self-intersecting face pairs, because a zero-length leg has no tangent and the wall
was laid on top of itself. Now refused by name (`MIN_SEGMENT_KM` = 1 m).

**(b) A leg shorter than the mitre excursion.** Even with a legal 2.2 m leg, a 250 m
corridor cannot turn: the mitred wall vertex reaches ``half_width * tan(turn/2)``
along each leg, so a leg shorter than that has its inner wall laid back over the
previous strip. Measured (guards lifted on purpose, so the geometry answers)
against the criterion now enforced:

| turn [deg] | excursion half*tan(turn/2) [m] | legs measured [m] | folded below [m] | simple from [m] | criterion |
|---:|---:|---|---:|---:|---|
| 20 | 22.0 | 6 - 55 | 21 | 23 | fold region ends at the excursion |
| 45 | 51.8 | 13 - 129 | 49 | 54 | fold region ends at the excursion |
| 90 | 125.0 | 31 - 312 | 119 | 131 | fold region ends at the excursion |

The criterion, not a fitted constant: for a 20 deg turn the fold disappears between
0.95 and 1.05 excursions, for 45 deg and 90 deg likewise, and `ice_channel` now
refuses exactly the region that folds, quoting both the turn and the offending leg.
The regression test lifts both guards and shows the folded solid the guards prevent.

## 4. Error verdict

| key quantity | reference | error | < 3 % ? |
|---|---|---:|---|
| closed-form scale factor | direct geodesic problem, 1152 cases | 1.3e-09 | yes, by 8 orders |
| corridor volume | analytic prism, 2920607332.310 m3 | 2.7e-13 | yes |
| corridor length (plane) | geodesic total, 1061.791 km | +0.0233 % | yes |
| corridor length (de-biased) | geodesic total | -4.2e-06 % | yes |
| corridor width, true metres | nominal 250 m | -0.279 % worst narrowing | yes |
| Euler characteristic | 2 for a genus-0 closed solid | exact | yes |

No quantity on this page is above 3 %, and the two that are nonzero are now
explained by a closed form rather than tolerated by a bound.

## 5. Traps and scope

* **The plane is a co-registration frame, not a survey.** Hulls, the corridor and
  every other object of the campaign are placed in the same plane, so the frame
  must not be silently changed to a per-leg one: that would remove the distortion
  and break the alignment between hull and channel. The correction is published as
  a factor instead.
* **`measured_width_m > width` is the mitre, not an error.** The mitre cap is
  documented as clipping sharp turns; the same cap widens the station chord.
* **Spherical earth.** Haversine on R = 6371.0088 km ignores the ellipsoid
  (~0.3 % in the meridian radius at 70 deg N). Since the corridor is 250 m wide,
  the ellipsoid matters more than the projection: a WGS84 geodesic would shift the
  route length by ~0.05 % and every station position by up to tens of metres. Out
  of scope here, stated rather than hidden.
* **No inverse projection.** The walls are defined in the plane; their true
  latitude/longitude would need the inverse map, which this module does not
  provide, so widths are corrected by the local factor rather than by
  reconstructing wall points on the sphere.
* **The route is a campaign polyline at 0.1 deg precision**, not a surveyed
  track: the waypoint coordinates themselves are ±5 km.

Figures: `docs/icechannel/scale_law.png` (distortion vs distance, closed form and
this route), `docs/icechannel/corridor.png` (plan view, colour = narrowing).
