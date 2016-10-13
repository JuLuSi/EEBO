#include "gtest/gtest.h"

#include "test_init.h"
#include "App.h"

using namespace EEBO;
using namespace libMesh;

TEST(AppTest, CreateAppWithOptionfile)
{
  auto app = std::make_unique<App>("../tests/data/CreateAppWithOptionfile.in", init->comm());
}

TEST(AppTest, CreateAppWithOutOptionfile)
{
  auto app = std::make_unique<App>(nullptr, init->comm());
}