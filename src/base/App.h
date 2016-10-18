#ifndef EEBOAPP_HPP
#define EEBOAPP_HPP

#include "EEBO.h"
// libMesh
#include "libmesh/mesh.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel_object.h"

namespace EEBO {

/// Application singleton
class App : public libMesh::ParallelObject {
 public:
  /// Constructor.
  /// \param infile Path to file with input options
  /// \param comm Communicator
  /// \return
  App(const char* infile, libMesh::Parallel::Communicator& comm);

  ~App() override;

  /// Parse the input file and publicate valid options to the application.
  void parseInputFile();

  /// Generate a mesh.
  void generateMesh();

  GetPot infile() { return infile_; }

  std::shared_ptr<libMesh::Mesh> mesh() { return mesh_; };

 private:
  /// Displays information about the initialised application.
  void displayInfoBanner();

  /// Input file object.
  GetPot infile_;

  /// Generated or loaded mesh.
  std::shared_ptr<libMesh::Mesh> mesh_;
};
}  // namespace EEBO

#endif
