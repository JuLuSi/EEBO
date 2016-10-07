#ifndef EEBOINIT_HH
#define EEBOINIT_HH

#include "EEBO.hpp"
// libMesh
#include "libmesh/libmesh.h"

namespace EEBO {
class Init : public libMesh::LibMeshInit
{
public:
  Init(int argc, char* argv[], MPI_Comm COMM_WORLD_IN = MPI_COMM_WORLD);
  ~Init() override = default;
};
} // namespace EEBO

#endif
