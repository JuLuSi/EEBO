#ifndef FEPROBLEM_HH
#define FEPROBLEM_HH

#include "EEBO.hh"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"

#include <type_traits>

namespace EEBO {

template<typename T>
class FEProblem : public libMesh::EquationSystems
{
public:
  explicit FEProblem(libMesh::Mesh& mesh);
  ~FEProblem() override = default;

  void init() override;
  void step();

  T& sys();

 protected:
  T* _sys;
  double _dt = 0;
};

template<typename T>
FEProblem<T>::FEProblem(libMesh::Mesh& mesh) :
    EquationSystems(mesh)
{
  static_assert(std::is_base_of<libMesh::TransientNonlinearImplicitSystem, T>::value,
                "FEProblem needs to be instantiated with a type which is inherited from SystemBase");

  _sys = &(add_system<T>("sys0"));
}

template<typename T>
void FEProblem<T>::init()
{
  libMesh::EquationSystems::init();
  libMesh::EquationSystems::print_info(out);
}

template<typename T>
void FEProblem<T>::step()
{
  *(_sys->old_local_solution) = *(_sys->current_local_solution);
  _sys->solve();
}

template<typename T>
T& FEProblem<T>::sys() {
  return *_sys;
}

} // namespace EEBO
#endif
