#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include "EEBO.hh"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_solver.h"

namespace EEBO {

class HeatTransfer : public TransientNonlinearImplicitSystem,
                     public NonlinearImplicitSystem::ComputeJacobian,
                     public NonlinearImplicitSystem::ComputeResidual,
                     public System::Initialization
{
public:
  HeatTransfer(EquationSystems& eqs, const std::string& name, const unsigned int number);

  ~HeatTransfer();

  virtual void initialize() override;

  virtual void jacobian(const NumericVector<Number>& X,
                        SparseMatrix<Number>& J,
                        NonlinearImplicitSystem& S) override;

  virtual void residual(const NumericVector<Number>& X,
                        NumericVector<Number>& F,
                        NonlinearImplicitSystem& S) override;

  static Number initialSolution(const Point& p,
                                const Parameters& parameters,
                                const std::string& sys_name,
                                const std::string& unknown_name);

private:
  bool _verbose = true;
  unsigned int _dim;
  unsigned int _temperature_varnum;
};

}
#endif
