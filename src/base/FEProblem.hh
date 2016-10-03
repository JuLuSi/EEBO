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
class FEProblem : public EquationSystems
{
public:
  FEProblem(Mesh& mesh);
  virtual ~FEProblem()
  {};

  void init();
  void step();

protected:
  T* _sys;
  double _dt = 0;
};

template<typename T>
FEProblem<T>::FEProblem(Mesh& mesh) :
    EquationSystems(mesh)
{
  static_assert(std::is_base_of<TransientNonlinearImplicitSystem, T>::value,
                "FEProblem needs to be instantiated with a type which is inherited from SystemBase");

  _sys = &(add_system<T>("sys0"));
}

template<typename T>
void FEProblem<T>::init()
{
  EquationSystems::init();
  EquationSystems::print_info(out);
}

template<typename T>
void FEProblem<T>::step()
{
  *(_sys->old_local_solution) = *(_sys->current_local_solution);
  _sys->solve();
}

}
#endif
