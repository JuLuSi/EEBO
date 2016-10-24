#include "gtest/gtest.h"

#include "test_init.h"
#include "App.h"
#include "FEProblem.h"
#include "HeatTransfer.h"
#include "libmesh/dof_map.h"
#include "libmesh/analytic_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/petsc_matrix.h"

using namespace EEBO;
using namespace libMesh;

TEST(Eigenproblem, AllDirichlet) {
  auto app = std::make_unique<App>(nullptr, init->comm());

  auto model = std::make_unique<FEProblem<HeatTransfer> >(*app->mesh());

  HeatTransfer& sys = model->sys();

  // Apply analytic function to Dirichlet boundary conditions.
  const boundary_id_type all_ids[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> all_bdys(all_ids, all_ids + (2 * 2));
  std::vector<unsigned int> vars(1, 0);
  ZeroFunction<Number> zf;
  sys.get_dof_map().add_dirichlet_boundary(DirichletBoundary(all_bdys, vars, zf));

  model->init();

  // Extract non-essential dofs
  const DofMap& dof_map = sys.get_dof_map();
  std::vector<dof_id_type> local_non_condensed_dofs_set;
  for (dof_id_type i = dof_map.first_dof(); i < dof_map.end_dof(); i++) {
    if (!dof_map.is_constrained_dof(i)) {
      local_non_condensed_dofs_set.push_back(i);
    }
  }

  SparseMatrix<Number>* J = sys.jacobianMatrix();
  SparseMatrix<Number>* M = sys.massMatrix();

  auto Jc = std::make_unique<PetscMatrix<Number> >(init->comm());
  J->create_submatrix(*Jc.get(), local_non_condensed_dofs_set, local_non_condensed_dofs_set);
  auto Mc = std::make_unique<PetscMatrix<Number> >(init->comm());
  M->create_submatrix(*Mc.get(), local_non_condensed_dofs_set, local_non_condensed_dofs_set);

  auto eps = EigenSolver<Number>::build(init->comm());
  eps->set_eigenproblem_type(GHEP);
  eps->set_position_of_spectrum(SMALLEST_MAGNITUDE);
  eps->solve_generalized(*Jc, *Mc, 5, 10, 10e-8, 1000);

  EXPECT_NEAR(eps->get_eigenvalue(0).first, static_cast<double>(pi*pi + pi*pi), 10e-4);
  EXPECT_NEAR(eps->get_eigenvalue(1).first, static_cast<double>(2*pi*pi + 3*pi*pi), 10e-3);
  EXPECT_NEAR(eps->get_eigenvalue(2).first, static_cast<double>(3*pi*pi + 2*pi*pi), 10e-3);
  EXPECT_NEAR(eps->get_eigenvalue(3).first, static_cast<double>(4*pi*pi + 4*pi*pi), 10e-2);
}
