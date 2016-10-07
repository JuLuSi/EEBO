#include "gtest/gtest.h"

#include "Init.hpp"
#include "App.hpp"
#include "FEProblem.hpp"
#include "HeatTransfer.hpp"

#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

using namespace EEBO;
using namespace libMesh;

class DriverTest : public ::testing::Test
{
 protected:
  virtual void SetUp() {
    int argc = 0;
    char program_name[] = "InitGetComm";
    char program_args[] = "";

    char *argv[] = {program_name, program_args};
    init_ = std::make_unique<Init>(argc, argv);
  }

  std::unique_ptr<Init> init_;
};

TEST_F(DriverTest, CreateSimpleDriver) {
  auto app = std::make_unique<App>(nullptr, init_->comm());
  Mesh mesh(init_->comm());
  MeshTools::Generation::build_square(mesh,
                                      10, 10,
                                      -1., 1.,
                                      -1., 1.,
                                      QUAD9);

  auto model = std::make_unique<FEProblem<HeatTransfer> >(mesh);

  EXPECT_NO_THROW(model->init());
}