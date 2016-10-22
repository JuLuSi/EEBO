#ifndef HEATTRANSFER_HPP
#define HEATTRANSFER_HPP

#include "EEBO.h"
#include "SystemBase.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"

namespace EEBO {

///
/// Heat transfer model
/// A simple single variable system representing the heat transfer in a domain.
///
/// F(u) = \Delta u + f = 0
///
/// Default forcing term f is 0.0
///
class HeatTransfer final : public SystemBase {
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

  static libMesh::Number initialState(const libMesh::Point& p,
                                      const libMesh::Parameters& parameters,
                                      const std::string& sys_name,
                                      const std::string& unknown_name);

  void setForcing(const double forcing_term) { forcing_term_ = forcing_term; };

 private:
  void timeDerivative() override;

  void assembleMass() override;

  unsigned int temperature_varnum_; ///< Variable number of the temperature.

  double forcing_term_ = 0.0; ///< Forcing term in the residual (without test function).
};

} // namespace EEBO
#endif
