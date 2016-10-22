#include "App.h"
#include "libmesh/mesh_generation.h"

using namespace libMesh;
using namespace EEBO;

App::App(const char* infile, Parallel::Communicator& comm) :
    ParallelObject(comm) {

  if (infile != nullptr) {
    infile_ = GetPot(infile);
  }
  auto mesh_type = infile_("Mesh/type", "generated");

  if (strcmp(mesh_type, "generated") == 0) {
    generateMesh();
  } else if (strcmp(mesh_type, "file") == 0) {
    auto mesh_filename = infile_("Mesh/filename", "");
    mesh_ = std::make_shared<Mesh>(this->comm());
    mesh_->read(mesh_filename);
  }

  displayInfoBanner();
  parseInputFile();
}

App::~App() {

}

void App::displayInfoBanner() {
  out << std::endl
      << "-------------------------------------------------------------------------- \n"
      << "                   Energy Efficient Building Operation DFG Project         \n"
      << "\n"
      << " Version: " << "0.0" << " Date: " << __DATE__ << " " << __TIME__ << "\n"
      << "-------------------------------------------------------------------------- \n"
      << std::endl;
}

void App::parseInputFile() {

}

void App::generateMesh() {
  mesh_ = std::make_shared<Mesh>(this->comm());
  MeshTools::Generation::build_square(*mesh_,
                                      10, 10,
                                      0., 1.,
                                      0., 1.,
                                      QUAD9);
}
