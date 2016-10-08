#include "gtest/gtest.h"

#include "test_init.h"
#include "App.h"
#include "FEProblem.h"
#include "HeatTransfer.h"

#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

using namespace EEBO;
using namespace libMesh;

TEST(DriverTest, CreateSimpleDriver)
{
  auto app = std::make_unique<App>(nullptr, init->comm());
  Mesh mesh(init->comm());
  MeshTools::Generation::build_square(mesh,
                                      10, 10,
                                      -1., 1.,
                                      -1., 1.,
                                      QUAD9);

  auto model = std::make_unique<FEProblem<HeatTransfer> >(mesh);

  EXPECT_NO_THROW(model->init());
}