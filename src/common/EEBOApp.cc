#include "EEBOApp.hh"

EEBO::App::App(const char* inputFile)
{
  if (inputFile != nullptr) {
    _inputFile = new GetPot("inputFile");
  }
}

EEBO::App::~App()
{

}

void EEBO::App::init(int argc, const char** argv)
{
  LibMeshInit init(argc, argv);
  _comm = &(init.comm());
}