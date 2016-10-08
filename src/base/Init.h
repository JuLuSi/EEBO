#ifndef EEBO_INIT_HPP
#define EEBO_INIT_HPP

#include "EEBO.h"
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
