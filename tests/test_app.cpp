#include "gtest/gtest.h"

#include "test_init.h"
#include "App.h"

using namespace EEBO;
using namespace libMesh;

TEST(App, CreateAppWithOptionfile) {
  auto app = std::make_unique<App>("../tests/data/CreateAppWithOptionfile.in", init->comm());
}

TEST(App, CreateAppWithOutOptionfile) {
  auto app = std::make_unique<App>(nullptr, init->comm());
}