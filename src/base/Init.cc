#include "Init.hh"

// PETSc
#ifdef LIBMESH_HAVE_PETSC
#include "petscsys.h"
#endif

using namespace libMesh;
using namespace EEBO;

Init::Init(int argc, char** argv, MPI_Comm COMM_WORLD_IN) :
    LibMeshInit(argc, argv, COMM_WORLD_IN)
{
// Get rid of PETSc error handler.
#ifdef LIBMESH_HAVE_PETSC
  PetscPopSignalHandler();
#endif
}
