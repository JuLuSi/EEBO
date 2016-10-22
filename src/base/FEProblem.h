#ifndef FEPROBLEM_HH
#define FEPROBLEM_HH

#include "EEBO.h"
#include "SystemBase.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"

#include <type_traits>

namespace EEBO {

template<typename T>
class FEProblem : public libMesh::EquationSystems {
 public:
  explicit FEProblem(libMesh::Mesh& mesh);
  ~FEProblem() override = default;

  void init() override;

  T& sys();

 protected:
  T* sys_;
};

template<typename T>
FEProblem<T>::FEProblem(libMesh::Mesh& mesh) :
    EquationSystems(mesh) {
  static_assert(std::is_base_of<SystemBase, T>::value,
                "FEProblem needs to be instantiated with a type which is inherited from SystemBase");

  sys_ = &(add_system<T>("sys0"));
}

template<typename T>
void FEProblem<T>::init() {
  libMesh::EquationSystems::init();
  libMesh::EquationSystems::print_info(out);
}

template<typename T>
T& FEProblem<T>::sys() {
  return get_system<T>("sys0");
}

} // namespace EEBO
#endif
