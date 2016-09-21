#ifndef EEBOAPP_HH
#define EEBOAPP_HH

#include "EEBO.hh"

#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

namespace EEBO {
class App
{

public:
  App(const char* inputFile);

  ~App();

  void init(int argc, const char** argv);

  Parallel::Communicator* comm()
  { return _comm; }

private:
  GetPot* _inputFile;
  Parallel::Communicator* _comm;
};
}

#endif
