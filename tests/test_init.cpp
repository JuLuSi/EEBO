#include "gtest/gtest.h"

#include "Init.hpp"

using namespace EEBO;

TEST(Init, InitGetComm) {
  int argc = 0;
  char program_name[] = "InitGetComm";
  char program_args[] = "";

  char *argv[] = {program_name, program_args};
  EXPECT_NO_THROW(Init init(argc, argv));
}