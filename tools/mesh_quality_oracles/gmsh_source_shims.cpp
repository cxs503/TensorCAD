// Two symbols the Gmsh element classes reference but the quality measures never
// call. Gmsh's own definitions of both (`GmshMessage.cpp`, `Context.cpp`) pull the
// Options / ColorTable / graphics graph into the build; that graph is graphing, not
// geometry, so it is stubbed here instead of vendored.
//
// * `Msg::Error` / `Msg::Warning` -- reached only from Gmsh's own error paths
//   (malformed meshes). The oracle validates its input before it gets there.
// * `Msg::Info` -- reached from `_getQualityFunctionSpace` (`mesh/qualityMeasuresJacobian.cpp`)
//   when it has to fall back to a full-space basis, i.e. on the `minSICN`/`minSIGE` path
//   this harness now exercises. Stubbed as a real message (printf), so the fallback is
//   visible on stderr instead of swallowed.
// * `CTX::instance()` -- reached from the *display* methods behind the element
//   vtables (`getEdgeRep`, `getFaceRep`, `getNumEdgesRep`). The pointer is never
//   dereferenced by a quality measure; storage is supplied so the linker is
//   satisfied without linking the Options struct.
#include <cstdarg>
#include <cstdio>

#include "Context.h"
#include "GModel.h"
#include "GmshMessage.h"
#include "MElement.h"

void Msg::Error(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
}

void Msg::Info(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
}

void Msg::Warning(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
}

static char ctxStorage[1 << 16] __attribute__((aligned(64)));

// * `GModel::current(int)` -- `MVertex`'s constructor (reached by the oracle when it
//   builds the vertices of a corpus mesh) asks the current model for its vertex
//   numbering (`setMaxVertexNumber` / `incrementAndGetMaxVertexNumber`, both inline in
//   `GModel.h`). Gmsh's own definition lives in `GModel.cpp`, which drags the full
//   model/Geo/Mesh graph; the numbering is bookkeeping the quality measures never read,
//   so a zero-initialised model image is supplied instead. Only inline member access
//   happens on it; no virtual call, no allocation.
static char modelStorage[sizeof(GModel)] __attribute__((aligned(64)));

GModel *GModel::current(int) { return reinterpret_cast<GModel *>(modelStorage); }

// * `GModel::getMeshElementByTag(n, entityTag)` -- the MSub* element constructors
//   (`MSubElement.cpp`) resolve their parent element through the model. Those
//   constructors are only reached from the element factory, which nothing on
//   the quality path calls (see the retention note in `build_gmsh_source_oracle.sh`);
//   they are linked only because the factory is. Returning `nullptr` is therefore never
//   observed: the only code that would read it is unreachable.
MElement *GModel::getMeshElementByTag(std::size_t, int &entityTag)
{
  entityTag = 0;
  return nullptr;
}

// * `tetcircumcenter` -- declared in `mesh/meshGRegionDelaunayInsertion.h` and defined in
//   `mesh/meshGRegionDelaunayInsertion.cpp`, which is the Delaunay insertion driver and
//   drags the whole mesh/Geo graph. The only caller reachable from this harness is
//   `MTetrahedron::circumcenter()` (and `getOuterRadius`), which no quality measure here
//   uses: `gamma`, `eta` and `getVolume` are edge-length/determinant arithmetic and
//   `minSICN`/`minSIGE` go through `JacobianBasis`. The stub answers NaN -- not 0 -- so
//   that an accidental future use produces a visible NaN rather than a plausible number.
double tetcircumcenter(double a[3], double b[3], double c[3], double d[3],
                       double circumcenter[3], double *xi, double *eta,
                       double *zeta)
{
  (void)a; (void)b; (void)c; (void)d;
  (void)xi; (void)eta; (void)zeta;
  for(int i = 0; i < 3; ++i) circumcenter[i] = __builtin_nan("");
  return __builtin_nan("");
}
