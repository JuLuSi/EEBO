#include "gtest/gtest.h"

#include "test_init.h"
#include "App.h"
#include "FEProblem.h"
#include "HeatTransfer.h"
#include "libmesh/dof_map.h"
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/exact_solution.h"
#include "libmesh/exodusII_io.h"
#include "petsc.h"

using namespace EEBO;
using namespace libMesh;

TEST(Driver, CreateSimpleDriver) {
  auto app = std::make_unique<App>(nullptr, init->comm());

  auto model = std::make_unique<FEProblem<HeatTransfer> >(*app->mesh());

  EXPECT_NO_THROW(model->init());
}

// Build a driver with the heat transfer model and calculate for the exact solution
// with pure Dirichlet boundary conditions
// u = x^2 + y^2
// f = -4
// s.t.
// \Delta u + f = 4 - 4 = 0
Real exactDirichletAnalytic(const Real x,
                            const Real y) {
  return x*x + y*y;
}

void exactDirichletAnalyticWrapper(libMesh::DenseVector<libMesh::Number>& value,
                                   const libMesh::Point& p,
                                   const libMesh::Real) {
  value(0) = exactDirichletAnalytic(p(0), p(1));
}

TEST(Driver, RunHeatTransferModelWithExactDirichlet) {
  auto app = std::make_unique<App>(nullptr, init->comm());
  auto model = std::make_unique<FEProblem<HeatTransfer> >(*app->mesh());

  HeatTransfer& sys = model->sys();

  // Apply analytic function to Dirichlet boundary conditions.
  const boundary_id_type all_ids[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> all_bdys(all_ids, all_ids + (2 * 2));
  std::vector<unsigned int> vars(1, 0);
  AnalyticFunction<> exact_solution_object(exactDirichletAnalyticWrapper);
  sys.get_dof_map().add_dirichlet_boundary(DirichletBoundary(all_bdys, vars, exact_solution_object));

  sys.setForcing(-4.0);

  model->init();

  PetscOptionsSetValue(PETSC_NULL, "-ksp_type", "preonly");
  PetscOptionsSetValue(PETSC_NULL, "-pc_type", "lu");

  sys.solve();

  // Compute L2 error against the analytic solution
  ExactSolution exact_sol(*model);
  exact_sol.attach_exact_value(0, &exact_solution_object);
  exact_sol.compute_error("sys0", "Temperature");
  Real l2err = exact_sol.l2_error("sys0", "Temperature");

  out << "L2 error: " << l2err << std::endl;
  EXPECT_LE(l2err, 1E-10);
}