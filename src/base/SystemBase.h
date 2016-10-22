#ifndef SYSTEMBASE_HPP
#define SYSTEMBASE_HPP

#include "EEBO.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"

namespace EEBO {

/// Base class for every simulation model.
/// Describes a basis for steady or transient linear or nonlinear systems.
/// The base classes demand certain functions to be implemented, which ensure a common functionality of the models.
class SystemBase : public libMesh::TransientNonlinearImplicitSystem,
                   public libMesh::NonlinearImplicitSystem::ComputeJacobian,
                   public libMesh::NonlinearImplicitSystem::ComputeResidual,
                   public libMesh::System::Initialization {
 public:
  virtual ~SystemBase();

  /// Get a pointer to the modifiable Jacobian matrix.
  libMesh::SparseMatrix<libMesh::Number>* jacobianMatrix() { return matrix; };

  /// Get a pointer to the modifiable mass matrix.
  libMesh::SparseMatrix<libMesh::Number>* massMatrix() { return this->request_matrix("mass"); };

 protected:
  /// Protected constructor to ensure this base class is not created without subclassing first.
  /// \param eqs
  /// \param name Name of the system. Will be set by FEProblem.
  /// \param number Number of the system in the list of equation systems. Will be set by EquationSystems.
  /// \return
  SystemBase(libMesh::EquationSystems& eqs, const std::string& name, const unsigned int number);

  virtual void timeDerivative() = 0;

  virtual void assembleMass() = 0;

  /// Verbose output of member functions.
  bool verbose_ = true;

  /// Physical dimension of the Model.
  unsigned int dim_;

  /// Timestep for transient models.
  unsigned int dt_ = 0;
};

} // namespace EEBO

#endif
