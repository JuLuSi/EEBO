#include "test_init.h"

#include "gtest/gtest.h"

using namespace EEBO;

Init* init;

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);

  init = new Init(argc, argv);
  auto return_code = RUN_ALL_TESTS();
  delete init;

  return return_code;
}