#include "EEBO.hh"
#include "Init.hh"
#include "App.hh"
#include "FEProblem.hh"
#include "HeatTransfer.hh"
#include "libmesh/mesh_generation.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"

using namespace libMesh;
using namespace EEBO;

int main(int argc, char** argv)
{
  Init init(argc, argv);
  auto app = std::make_unique<App>(nullptr, init.comm());

  Mesh mesh(init.comm());
  MeshTools::Generation::build_square(mesh,
                                      10, 10,
                                      -1., 1.,
                                      -1., 1.,
                                      QUAD9);
  mesh.print_info();
  
  FEProblem<HeatTransfer> model(mesh);
  model.init();
  model.solve();

  auto exoio = std::make_unique<ExodusII_IO>(mesh);
  exoio->write_equation_systems("test.e", model);

  return 0;
}
