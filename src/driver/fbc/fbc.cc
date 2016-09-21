#include "EEBOApp.hh"

using namespace EEBO;

int main(int argc, const char** argv)
{
  auto fbcCtx = std::unique_ptr<App>(new App(nullptr));
  fbcCtx->init(argc, argv);

  return 0;
}