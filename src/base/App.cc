#include "App.hh"

EEBO::App::App(const char* infile, Parallel::Communicator& comm) :
    ParallelObject(comm)
{
  if (infile != nullptr) {
    _infile = GetPot(infile);
  }
  displayInfoBanner();
}

EEBO::App::~App()
{

}

void EEBO::App::displayInfoBanner()
{
  out
      << "--------------------------------------------------------------------------" << std::endl
      << "                   Energy Efficient Building Operation DFG Project        " << std::endl
      << std::endl
      << " Version: " << "0.0" << " Date: " << __DATE__ << " " << __TIME__ << std::endl
      << "--------------------------------------------------------------------------" << std::endl
      << std::endl;
}