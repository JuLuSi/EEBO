#ifndef HEATTRANSFER_HPP
#define HEATTRANSFER_HPP

#include "EEBO.h"
#include "SystemBase.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"

namespace EEBO {

class HeatTransfer : public SystemBase {
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
  unsigned int temperature_varnum_;
};

} // namespace EEBO
#endif
