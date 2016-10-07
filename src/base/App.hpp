#ifndef EEBOAPP_HH
#define EEBOAPP_HH

#include "EEBO.hpp"
// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parallel_object.h"

namespace EEBO {
class App : public libMesh::ParallelObject
{

public:
  App(const char* infile, libMesh::Parallel::Communicator& comm);

  ~App() override;

  GetPot infile()
  { return _infile; }

private:
  void displayInfoBanner();
  GetPot _infile;
};
}  // namespace EEBO

#endif
