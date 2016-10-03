#ifndef EEBOAPP_HH
#define EEBOAPP_HH

#include "EEBO.hh"
// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parallel_object.h"

namespace EEBO {
class App : public ParallelObject
{

public:
  App(const char* infile, Parallel::Communicator& comm);

  ~App();

  GetPot infile()
  { return _infile; }

private:
  void displayInfoBanner();
  GetPot _infile;
};
}

#endif
