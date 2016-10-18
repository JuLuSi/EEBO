#ifndef EEBO_INIT_HPP
#define EEBO_INIT_HPP

#include "EEBO.h"
// libMesh
#include "libmesh/libmesh.h"

namespace EEBO {

/// Init
class Init : public libMesh::LibMeshInit {
 public:
  /// Constructor
  /// \param argc
  /// \param argv
  /// \param COMM_WORLD_IN
  /// \return
  Init(int argc, char* argv[], MPI_Comm COMM_WORLD_IN = MPI_COMM_WORLD);

  /// Destructor
  ~Init() override = default;
};
} // namespace EEBO

#endif
