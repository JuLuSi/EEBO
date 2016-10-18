#include "HeatTransfer.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/parameters.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/equation_systems.h"
#include "libmesh/analytic_function.h"

using namespace libMesh;
using namespace EEBO;

HeatTransfer::HeatTransfer(EquationSystems& eqs, const std::string& name, const unsigned int number) :
    SystemBase(eqs, name, number) {

  temperature_varnum_ = add_variable("Temperature", SECOND, LAGRANGE);

  dim_ = 2;

  attach_init_object(*this);
  nonlinear_solver->jacobian_object = this;
  nonlinear_solver->residual_object = this;
}

HeatTransfer::~HeatTransfer() {

}

void HeatTransfer::initialize() {
  project_solution(HeatTransfer::initialState, nullptr, get_equation_systems().parameters);

  *old_local_solution = *current_local_solution;

  if (verbose_)
    out << "<<< Initializing HeatTransfer" << std::endl;
}

void HeatTransfer::jacobian(const NumericVector<Number>& X,
                            SparseMatrix<Number>& J,
                            NonlinearImplicitSystem& /* S */) {
  const DofMap& dof_map = get_dof_map();
  FEType fe_type = dof_map.variable_type(temperature_varnum_);

  UniquePtr<FEBase> fe(FEBase::build(dim_, fe_type));
  QGauss qrule(dim_, FIFTH);
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim_, fe_type));
  QGauss qface(dim_ - 1, FIFTH);
  fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real>>& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

  DenseMatrix<Number> Je;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator el = get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    fe->reinit(elem);

    Je.resize(dof_indices.size(), dof_indices.size());

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
      Number u = 0;
      Gradient grad_u;

      for (unsigned int j = 0; j < phi.size(); j++) {
        u += phi[j][qp] * X(dof_indices[j]);
        grad_u += dphi[j][qp] * X(dof_indices[j]);
      }

      for (unsigned int i = 0; i < phi.size(); i++) {
        for (unsigned int j = 0; j < phi.size(); j++) {
          Je(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
        }
      }
    }
    dof_map.constrain_element_matrix(Je, dof_indices);
    J.add_matrix(Je, dof_indices);
  }
}

void HeatTransfer::residual(const NumericVector<Number>& X,
                            NumericVector<Number>& F,
                            NonlinearImplicitSystem& /* S */) {
  const DofMap& dof_map = get_dof_map();
  FEType fe_type = dof_map.variable_type(temperature_varnum_);

  UniquePtr<FEBase> fe(FEBase::build(dim_, fe_type));
  QGauss qrule(dim_, FIFTH);
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim_, fe_type));
  QGauss qface(dim_ - 1, FIFTH);
  fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real>>& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  F.zero();

  MeshBase::const_element_iterator el = get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    fe->reinit(elem);
    Fe.resize(dof_indices.size());

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
      Number u = 0;
      Gradient grad_u;

      for (unsigned int j = 0; j < phi.size(); j++) {
        u += phi[j][qp] * X(dof_indices[j]);
        grad_u += dphi[j][qp] * X(dof_indices[j]);
      }

      for (unsigned int i = 0; i < phi.size(); i++) {
        Fe(i) += JxW[qp] * (grad_u * dphi[i][qp] - (forcing_term_ * phi[i][qp]));
      }
    }
    dof_map.constrain_element_vector(Fe, dof_indices);
    F.add_vector(Fe, dof_indices);
  }
}

Number HeatTransfer::initialState(const Point& /* p */,
                                  const Parameters& /* parameters */,
                                  const std::string& /* sys_name */,
                                  const std::string& /* unknown_name */) {
  return 0.0;
}

void HeatTransfer::timeDerivative() {

}
