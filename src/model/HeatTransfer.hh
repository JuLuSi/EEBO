#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include "EEBO.hh"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"

namespace EEBO {

class HeatTransfer : public libMesh::TransientNonlinearImplicitSystem,
                     public libMesh::NonlinearImplicitSystem::ComputeJacobian,
                     public libMesh::NonlinearImplicitSystem::ComputeResidual,
                     public libMesh::System::Initialization
{
public:
  HeatTransfer(libMesh::EquationSystems& eqs, const std::string& name, const unsigned int number);

  ~HeatTransfer() override;

  void initialize() override;

  void jacobian(const libMesh::NumericVector<libMesh::Number>& X,
                libMesh::SparseMatrix<libMesh::Number>& J,
                libMesh::NonlinearImplicitSystem& S) override;

  void residual(const libMesh::NumericVector<libMesh::Number>& X,
                libMesh::NumericVector<libMesh::Number>& F,
                libMesh::NonlinearImplicitSystem& S) override;

  static libMesh::Number initialSolution(const libMesh::Point& p,
                                const libMesh::Parameters& parameters,
                                const std::string& sys_name,
                                const std::string& unknown_name);

private:
  bool _verbose = true;
  unsigned int _dim;
  unsigned int _temperature_varnum;
};

} // namespace EEBO
#endif
