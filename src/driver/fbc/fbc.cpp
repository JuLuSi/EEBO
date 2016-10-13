#include "EEBO.h"
#include "Init.h"
#include "App.h"
#include "FEProblem.h"
#include "HeatTransfer.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"

using namespace libMesh;
using namespace EEBO;

int main(int argc, char** argv) {
  auto init = std::make_unique<Init>(argc, argv);
  auto app = std::make_unique<App>(nullptr, init->comm());

  FEProblem<HeatTransfer> model(*app->mesh());
  model.init();
  model.solve();

  auto exoio = std::make_unique<ExodusII_IO>(*app->mesh());
  exoio->write_equation_systems("test.e", model);

  return 0;
}
