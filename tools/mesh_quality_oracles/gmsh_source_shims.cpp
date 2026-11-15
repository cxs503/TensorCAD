// Two symbols the Gmsh element classes reference but the quality measures never
// call. Gmsh's own definitions of both (`GmshMessage.cpp`, `Context.cpp`) pull the
// Options / ColorTable / graphics graph into the build; that graph is graphing, not
// geometry, so it is stubbed here instead of vendored.
//
// * `Msg::Error` / `Msg::Warning` -- reached only from Gmsh's own error paths
//   (malformed meshes). The oracle validates its input before it gets there.
// * `CTX::instance()` -- reached from the *display* methods behind the element
//   vtables (`getEdgeRep`, `getFaceRep`, `getNumEdgesRep`). The pointer is never
//   dereferenced by a quality measure; storage is supplied so the linker is
//   satisfied without linking the Options struct.
#include <cstdarg>
#include <cstdio>

#include "Context.h"
#include "GmshMessage.h"

void Msg::Error(const char *fmt, ...)
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

CTX *CTX::instance() { return reinterpret_cast<CTX *>(ctxStorage); }