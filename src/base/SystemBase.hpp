#ifndef SYSTEMBASE_HPP
#define SYSTEMBASE_HPP

#include "EEBO.hpp"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"

namespace EEBO {

class SystemBase : public libMesh::TransientNonlinearImplicitSystem,
                   public libMesh::NonlinearImplicitSystem::ComputeJacobian,
                   public libMesh::NonlinearImplicitSystem::ComputeResidual,
                   public libMesh::System::Initialization
{
 public:
  virtual ~SystemBase();

 protected:
  SystemBase(libMesh::EquationSystems& eqs, const std::string& name, const unsigned int number);

  bool _verbose = true;
  unsigned int _dim;
};

} // namespace EEBO

#endif
