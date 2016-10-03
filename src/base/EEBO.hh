#ifndef EEBO_HH
#define EEBO_HH

#include <memory>

// libMesh includes
#include "libmesh/perf_log.h"
#include "libmesh/parallel.h"
#include "libmesh/libmesh_common.h"

using namespace libMesh;

namespace EEBO {
using libMesh::out;
using libMesh::err;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&& ... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

#endif
