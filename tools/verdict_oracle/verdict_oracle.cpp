// Verdict oracle: a thin C ABI over the Verdict element-quality library.
//
// Why this file exists
// --------------------
// ``docs/mesh_quality.md`` §5 records the meshmagick comparison as *not run*: meshmagick's
// ``print_quality()`` delegates to ``vtk.vtkMeshQuality``, and no VTK wheel exists for
// loongarch64.  But ``vtkMeshQuality`` is a wrapper around Verdict -- the element metrics
// live in VTK's ``ThirdParty/verdict``.  Compiling those same sources here gives a runtime
// oracle with the *same code path* as meshmagick, so the quality numbers this tree reports
// can be compared against the library meshmagick would have called, without VTK and
// without re-implementing the definitions a third time.
//
// This wrapper deliberately contains no metric algebra: every number comes from Verdict.
// It only maps a cell-kind name to a node count, a metric name to a Verdict function, and
// a flat coordinate buffer to the ``double[n][3]`` layout Verdict expects.
//
// Build: ``python tools/verdict_oracle/build.py`` (see that file for the pinned VTK
// reference).  The resulting shared library is not committed; the recipe is.

#include <cstddef>
#include <cstring>

#include "verdict.h"

namespace {

// Every Verdict quality function used here has the signature
// ``double f(int num_nodes, const double coordinates[][3])``.
using MetricFn = double (*)(int, const double[][3]);

struct Metric {
  const char *name;
  MetricFn fn;
};

constexpr int kTri = 0;
constexpr int kQuad = 1;
constexpr int kTet = 2;
constexpr int kHex = 3;

const Metric kTriMetrics[] = {
    {"area", &verdict::tri_area},
    {"scaled_jacobian", &verdict::tri_scaled_jacobian},
    {"aspect_ratio", &verdict::tri_aspect_ratio},
    {"condition", &verdict::tri_condition},
    {"radius_ratio", &verdict::tri_radius_ratio},
    {"min_angle", &verdict::tri_minimum_angle},
    {"max_angle", &verdict::tri_maximum_angle},
    {"edge_ratio", &verdict::tri_edge_ratio},
    {"distortion", &verdict::tri_distortion},
    {"equiangle_skew", &verdict::tri_equiangle_skew},
};

const Metric kQuadMetrics[] = {
    {"area", &verdict::quad_area},
    {"scaled_jacobian", &verdict::quad_scaled_jacobian},
    {"aspect_ratio", &verdict::quad_aspect_ratio},
    {"condition", &verdict::quad_condition},
    {"radius_ratio", &verdict::quad_radius_ratio},
    {"min_angle", &verdict::quad_minimum_angle},
    {"max_angle", &verdict::quad_maximum_angle},
    {"edge_ratio", &verdict::quad_edge_ratio},
    {"skew", &verdict::quad_skew},
    {"taper", &verdict::quad_taper},
    {"warpage", &verdict::quad_warpage},
    {"equiangle_skew", &verdict::quad_equiangle_skew},
};

const Metric kTetMetrics[] = {
    {"volume", &verdict::tet_volume},
    {"scaled_jacobian", &verdict::tet_scaled_jacobian},
    {"aspect_ratio", &verdict::tet_aspect_ratio},
    {"condition", &verdict::tet_condition},
    {"radius_ratio", &verdict::tet_radius_ratio},
    {"min_angle", &verdict::tet_minimum_angle},
    {"edge_ratio", &verdict::tet_edge_ratio},
    {"distortion", &verdict::tet_distortion},
    {"equiangle_skew", &verdict::tet_equiangle_skew},
};

const Metric kHexMetrics[] = {
    {"volume", &verdict::hex_volume},
    {"scaled_jacobian", &verdict::hex_scaled_jacobian},
    {"condition", &verdict::hex_condition},
    {"edge_ratio", &verdict::hex_edge_ratio},
    {"skew", &verdict::hex_skew},
    {"taper", &verdict::hex_taper},
    {"oddy", &verdict::hex_oddy},
    {"distortion", &verdict::hex_distortion},
    {"equiangle_skew", &verdict::hex_equiangle_skew},
};

struct Kind {
  const char *name;
  int nodes;
  const Metric *metrics;
  int count;
};

const Kind kKinds[] = {
    {"tri", 3, kTriMetrics, static_cast<int>(sizeof(kTriMetrics) / sizeof(Metric))},
    {"quad", 4, kQuadMetrics, static_cast<int>(sizeof(kQuadMetrics) / sizeof(Metric))},
    {"tet", 4, kTetMetrics, static_cast<int>(sizeof(kTetMetrics) / sizeof(Metric))},
    {"hex", 8, kHexMetrics, static_cast<int>(sizeof(kHexMetrics) / sizeof(Metric))},
};

constexpr int kKindCount = static_cast<int>(sizeof(kKinds) / sizeof(Kind));

const Kind *kind_at(int kind) {
  if (kind < 0 || kind >= kKindCount) return nullptr;
  return &kKinds[kind];
}

}  // namespace

extern "C" {

// Handshake used by the Python wrapper to reject a stale or foreign build.
const char *verdict_oracle_info(void) {
  return "verdict-oracle/1 (VTK ThirdParty/verdict; no VTK, no meshmagick)";
}

// Number of cell kinds the build exposes.
int verdict_oracle_kind_count(void) { return kKindCount; }

const char *verdict_oracle_kind_name(int kind) {
  const Kind *k = kind_at(kind);
  return k ? k->name : nullptr;
}

// Nodes per element for a kind, or -1 for an unknown kind.
int verdict_oracle_nodes_per_element(int kind) {
  const Kind *k = kind_at(kind);
  return k ? k->nodes : -1;
}

int verdict_oracle_metric_count(int kind) {
  const Kind *k = kind_at(kind);
  return k ? k->count : -1;
}

// The i-th metric name of a kind, or nullptr when the index is out of range.
const char *verdict_oracle_metric_name(int kind, int index) {
  const Kind *k = kind_at(kind);
  if (!k || index < 0 || index >= k->count) return nullptr;
  return k->metrics[index].name;
}

// Evaluate every metric of one element.  ``coords`` holds nodes*3 doubles in the element's
// own node order (the same order the metric module uses); ``out`` must hold at least
// ``verdict_oracle_metric_count(kind)`` doubles.  Returns 0 on success, -1 on a bad kind.
int verdict_oracle_evaluate(int kind, const double *coords, double *out) {
  const Kind *k = kind_at(kind);
  if (!k || coords == nullptr || out == nullptr) return -1;
  // Verdict indexes coordinates[num_nodes][3]; the buffer we are handed is already in that
  // layout, so a reinterpret_cast is enough -- but keep it explicit so the alignment
  // assumption is visible in one place.
  const double(*xyz)[3] = reinterpret_cast<const double(*)[3]>(coords);
  for (int i = 0; i < k->count; ++i) out[i] = k->metrics[i].fn(k->nodes, xyz);
  return 0;
}

}  // extern "C"