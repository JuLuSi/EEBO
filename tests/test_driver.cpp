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

  auto model = std::make_unique<FEProblem<HeatTransfer> >(*app->mesh());

  EXPECT_NO_THROW(model->init());
}