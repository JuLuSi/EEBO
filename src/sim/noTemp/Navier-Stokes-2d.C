// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Systems Example 2 - Unsteady Nonlinear Navier-Stokes</h1>
// \author John W. Peterson
// \date 2004
//
// This example shows how a simple, unsteady, nonlinear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Navier-Stokes equations for low-speed incompressible fluid flow.  This
// example introduces the concept of the inner nonlinear loop for each
// timestep, and requires a good deal of linear algebra number-crunching
// at each step.  If you have a ExodusII viewer such as ParaView installed,
// the script movie.sh in this directory will also take appropriate screen
// shots of each of the solution files in the time sequence.  These rgb files
// can then be animated with the "animate" utility of ImageMagick if it is
// installed on your system.  On a PIII 1GHz machine in debug mode, this
// example takes a little over a minute to run.  If you would like to see
// a more detailed time history, or compute more timesteps, that is certainly
// possible by changing the n_timesteps and dt variables below.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes (EquationSystems & es,
                      const std::string & system_name);
                      
// Perform explicit Newton step as starting guess for Newtons method

void perform_implicit_euler (EquationSystems & es, const std::string & system_name);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // This example NaNs with the Eigen sparse linear solvers and
  // Trilinos solvers, but should work OK with either PETSc or
  // Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements in 2D.  Building these
  // higher-order elements allows us to use higher-order
  // approximation, as in example 3.
  MeshTools::Generation::build_square (mesh,
                                       10, 10,
                                       0., 1.,
                                       0., 1.,
                                       QUAD4);

  mesh.all_second_order();

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Creates a transient system named "Navier-Stokes"
  TransientLinearImplicitSystem & system =
    equation_systems.add_system<TransientLinearImplicitSystem> ("Navier-Stokes");

  // Add the variables "u" & "v" to "Navier-Stokes".  They
  // will be approximated using second-order approximation.
  system.add_variable ("u", SECOND);
  system.add_variable ("v", SECOND);

  // Add the variable "p" to "Navier-Stokes". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  system.add_variable ("p", FIRST);

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_stokes);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Create a performance-logging object for this example
  PerfLog perf_log("Systems Example 2");

  // Get a reference to the Stokes system to use later.
  TransientLinearImplicitSystem & navier_stokes_system =
    equation_systems.get_system<TransientLinearImplicitSystem>("Navier-Stokes");

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  const Real dt = 0.01;
  navier_stokes_system.time     = 0.0;
  const unsigned int n_timesteps = 15;

  // The number of steps and the stopping criterion are also required
  // for the nonlinear iterations.
  const unsigned int n_nonlinear_steps = 100;
  const Real nonlinear_tolerance       = 1.e-3;

  // We also set a standard linear solver flag in the EquationSystems object
  // which controls the maxiumum number of linear solver iterations allowed.
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;

  // problem specifix parameters
  equation_systems.parameters.set<double>("rho") = 1;					// first line, density of the fluid
  equation_systems.parameters.set<double>("nu") = 1;					// first line, kinematic viscosity
  equation_systems.parameters.set<double>("eps") = 0; 				// for second line, pressure-velocity coupling

  // Tell the system of equations what the timestep is by using
  // the set_parameter function.  The matrix assembly routine can
  // then reference this parameter.make
  equation_systems.parameters.set<Real> ("dt")   = dt;

  // The first thing to do is to get a copy of the solution at
  // the current nonlinear iteration.  This value will be used to
  // determine if we can exit the nonlinear loop.
  UniquePtr<NumericVector<Number> >
    last_nonlinear_soln (navier_stokes_system.solution->clone());
    
std::cout<<"Norm of the initial value: "<<navier_stokes_system.solution->l2_norm()<<std::endl;

  for (unsigned int t_step=0; t_step<n_timesteps; ++t_step)
    {
  //~ equation_systems.parameters.set<double>("T_Dir") = 100*t_step; 				// heating on the Dirichlet boundary
      // Incremenet the time counter, set the time step size as
      // a parameter in the EquationSystem.
      navier_stokes_system.time += dt;

      // A pretty update message
      libMesh::out << "\n\n*** Solving time step "
                   << t_step
                   << ", time = "
                   << navier_stokes_system.time
                   << " ***"
                   << std::endl;

      // Now we need to update the solution vector from the
      // previous time step.  This is done directly through
      // the reference to the Stokes system.
      *navier_stokes_system.old_local_solution = *navier_stokes_system.current_local_solution;

      // At the beginning of each solve, reset the linear solver tolerance
      // to a "reasonable" starting value.
      const Real initial_linear_solver_tol = 1.e-6;
      equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;
	  std::cout<<"Norm of the Newton starting guess "<<navier_stokes_system.solution->l2_norm()<<std::endl;
      // Now we begin the nonlinear loop
      for (unsigned int l=0; l<n_nonlinear_steps; ++l)
        {
          // This is the solution of the last Newton loop
          last_nonlinear_soln->zero();
          last_nonlinear_soln->add(*navier_stokes_system.solution);

          // Assemble & solve the linear system to get the new Newton 
          // iterate, NOT THE UPDATE!!
          perf_log.push("linear solve");
          equation_systems.get_system("Navier-Stokes").solve();
          perf_log.pop("linear solve");

          // Compute the difference between this solution and the last
          // nonlinear iterate.
          last_nonlinear_soln->add (-1., *navier_stokes_system.solution);
		  //~ std::cout<<"norm of the solution:"<<navier_stokes_system.solution->l2_norm()<<std::endl;
          // Close the vector before computing its norm
          last_nonlinear_soln->close();

          // Compute the l2 norm of the differencenavier_stokes_system
          const Real norm_delta = last_nonlinear_soln->l2_norm();

          // How many iterations were required to solve the linear system?
          const unsigned int n_linear_iterations = navier_stokes_system.n_linear_iterations();

          // What was the final residual of the linear system?
          const Real final_linear_residual = navier_stokes_system.final_linear_residual();

          // Print out convergence information for the linear and
          // nonlinear iterations.
          libMesh::out << "Linear solver converged at step: "
                       << n_linear_iterations
                       << ", final residual: "
                       << final_linear_residual
                       <<"( "<< equation_systems.parameters.get<Real>("linear solver tolerance")<<" )"
                       << "  Nonlinear convergence: ||u - u_old|| = "
                       << norm_delta
                       << std::endl;

          // Terminate the solution iteration if the difference between
          // this nonlinear iterate and the last is sufficiently small, AND
          // if the most recent linear system was solved to a sufficient tolerance.
          if ((norm_delta < nonlinear_tolerance) &&
              (navier_stokes_system.final_linear_residual() < nonlinear_tolerance))
            {
              libMesh::out << " Nonlinear solver converged at step "
                           << l
                           << std::endl;
              break;
            }

          // Otherwise, decrease the linear system tolerance.  For the inexact Newton
          // method, the linear solver tolerance needs to decrease as we get closer to
          // the solution to ensure quadratic convergence.  The new linear solver tolerance
          // is chosen (heuristically) as the square of the previous linear system residual norm.

          equation_systems.parameters.set<Real> ("linear solver tolerance") =
            std::min(final_linear_residual*final_linear_residual,initial_linear_solver_tol);
            //std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol);
            //~ std::cout<<"linear solver tolerance set to "<<std::max(1e-16,std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol))<<std::endl;
        } // end nonlinear loop

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out every nth timestep to file.
      const unsigned int write_interval = 1;

      if ((t_step+1)%write_interval == 0)
        {
          std::ostringstream file_name;

          //~ // We write the file in the ExodusII format.
          file_name << "out_"
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step + 1
                    << ".pvtu";

          //~ ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                    //~ equation_systems);
          // output as vtk files
          VTKIO(mesh).write_equation_systems (file_name.str(),
                                                    equation_systems);
          //.gmv file
           //~ std::ostringstream file_name2;

          //~ // We write the file in the ExodusII format.
          //~ file_name2 << "out_"
                    //~ << std::setw(3)
                    //~ << std::setfill('0')
                    //~ << std::right
                    //~ << t_step + 1
                    //~ << ".gmv";         
           //~ GMVIO(mesh).write_equation_systems (file_name2.str(), equation_systems);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    } // end timestep loop.

  // All done.
  return 0;
}

// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void assemble_stokes (EquationSystems & es,
                      const std::string & system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Navier-Stokes");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("Navier-Stokes");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);
  
  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  // Build a Finite Element object of the specified type for
  // the pressure variables.
  UniquePtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  // Note that the QRule of the velocity (second order) is also used for 
  // the first order FE-Base of the pressure and temperature.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe_vel->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dphi = fe_vel->get_dphi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & psi = fe_pres->get_phi();

  // The value of the linear shape function gradients at the quadrature points
  const std::vector<std::vector<RealGradient> > & dpsi = fe_pres->get_dphi();
  
  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = navier_stokes_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;

  // Find out what the timestep size parameter is from the system, and
  // the value of theta for the theta method.  We use implicit Euler (theta=1)
  // for this simulation even though it is only first-order accurate in time.
  // The reason for this decision is that the second-order Crank-Nicolson
  // method is notoriously oscillatory for problems with discontinuous
  // initial data such as the lid-driven cavity.  Therefore,
  // we sacrifice accuracy in time for stability, but since the solution
  // reaches steady state relatively quickly we can afford to take small
  // timesteps.  If you monitor the initial nonlinear residual for this
  // simulation, you should see that it is monotonically decreasing in time.
  const Real dt    = es.parameters.get<Real>("dt");
  
  // We will just use implicit Euler here
  // const Real theta = 1.;
  
  // get the parameters of the system
  const auto rho = es.parameters.get<double>("rho");
  const auto nu = es.parameters.get<double>("nu");
  const auto eps = es.parameters.get<double>("eps");


  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Reposition the submatrices...  The idea is this:
      //
      //         -           -          -  -
      //        | Kuu Kuv Kup |        | Fu |
      //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
      //        | Kpu Kpv Kpp |        | Fp |
      //         -           -          -  -
      //
      // The DenseSubMatrix.repostition () member takes the
      // (row_offset, column_offset, row_size, column_size).
      //
      // Similarly, the DenseSubVector.reposition () member
      // takes the (row_offset, row_size)

      Kuu.reposition (0, 0, n_u_dofs, n_u_dofs);
      Kuv.reposition (0, n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (0, n_u_dofs + n_v_dofs, n_u_dofs, n_p_dofs);
      
      Kvu.reposition (n_u_dofs, 0, n_v_dofs, n_u_dofs);
      Kvv.reposition (n_u_dofs, n_u_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (n_u_dofs, n_u_dofs + n_v_dofs,  n_v_dofs, n_p_dofs);

      Kpu.reposition (n_u_dofs + n_v_dofs, 0, n_p_dofs, n_u_dofs);
      Kpv.reposition (n_u_dofs + n_v_dofs, n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (n_u_dofs + n_v_dofs, n_u_dofs + n_v_dofs, n_p_dofs, n_p_dofs);

      Fu.reposition (0, n_u_dofs);
      Fv.reposition (n_u_dofs, n_v_dofs);
      Fp.reposition (n_u_dofs + n_v_dofs, n_p_dofs);

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This must be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.
          Number u = 0., u_old = 0.;
          Number v = 0., v_old = 0.;
          Number p_old = 0.;

          Gradient grad_u, grad_u_old;
          Gradient grad_v, grad_v_old;

          // Compute the velocity & its gradient 
          // from the previous timestep
          // and the old Newton iterate.
          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              // Velocity from the old timestep:
              u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              grad_u_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_u[l]));
              grad_v_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_v[l]));

              // Velocity from the previous Newton iterate:
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
			}

            
          // Compute the old pressure value at this quadrature point.
          for (unsigned int l=0; l<n_p_dofs; l++)
            p_old += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);

          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (u_old, v_old);
          const NumberVectorValue U     (u,     v);
          const Number u_x = grad_u(0);
          const Number u_y = grad_u(1);
          const Number v_x = grad_v(0);
          const Number v_y = grad_v(1);

          // First, an i-loop over the velocity degrees of freedom.
          // We know that n_u_dofs == n_v_dofs so we can compute contributions
          // for both at the same time.
          // Note that, contrary to the systems_of_equations_ex2.C example
          // , we assume theta = 1 here, meaning we use implicit Euler.
          for (unsigned int i=0; i<n_u_dofs; i++)
            {
              Fu(i) += JxW[qp]*(u_old*phi[i][qp]                        // -C constant term from old timestep
                                +dt*( (U*grad_u)*phi[i][qp]            // -N(x_i) , nonlinear part at old Newton iterate.
								     +(u*u_x + u*u_x + v*u_y + v*u_y)*phi[i][qp]// N'(x_i)x_i, from lhs of Newton update formula
								    )
                                );              


              Fv(i) += JxW[qp]*(v_old*phi[i][qp]                        // constant term from old timestep
                                +dt*( (U*grad_v)*phi[i][qp]            // -N(x_i) , nonlinear part at old Newton iterate.
								     +(v*v_y + v*v_y + u*v_x+ u*v_x)*phi[i][qp]// N'(x_i)x_i, from lhs of Newton update formula
								     )		// constant term in F
                                ); 
                                

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j<n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp]*(phi[j][qp]*phi[i][qp] +                         // u^(i+1) * v_11
								    dt*( (phi[j][qp]*u_x + U*dphi[j][qp]) * phi[i][qp] // nonlinear part
										+ nu*dphi[j][qp]*dphi[i][qp]));                  // linear part
										
                  Kuv(i,j) += JxW[qp]*dt*u_y*phi[i][qp]*phi[j][qp];    				   // linear part

                  Kvv(i,j) += JxW[qp]*(phi[j][qp]*phi[i][qp] +                         // v^(i+1) * v12
								    dt*( (phi[j][qp]*v_y + U*dphi[j][qp]) * phi[i][qp] // nonlinear part
										+ nu*dphi[j][qp]*dphi[i][qp]));            	   // linear part                                    
           

                  Kvu(i,j) += JxW[qp]*dt*v_x*phi[i][qp]*phi[j][qp];                    // linear part
                }

              // Matrix contributions for the up and vp couplings.
              for (unsigned int j=0; j<n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*(dt*dpsi[j][qp](0)*phi[i][qp]) / rho; // (dF_1.1/dp) (p^(i+1))
                 Kvp(i,j) += JxW[qp]*(dt*dpsi[j][qp](1)*phi[i][qp]) / rho; // (dF_1.2/dp) (p^(i+1))
                 //   Kup(i,j) += JxW[qp]*(-dt*psi[j][qp]*dphi[i][qp](0)); //libmesh example
                   // Kvp(i,j) += JxW[qp]*(-dt*psi[j][qp]*dphi[i][qp](1));

                }
            }

          // Now an i-loop over the pressure degrees of freedom. 
          for (unsigned int i=0; i<n_p_dofs; i++)
            {
			  // Note that the Fp block is identically zero as F has no
              // nonlinearity and no constant, and G = F (no timestepping as no temporal derivative)
              
			  for (unsigned int j=0; j<n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*(-psi[i][qp]*dphi[j][qp](0));  // (dF_2/du) (u^(i+1))
                  Kpv(i,j) += JxW[qp]*(-psi[i][qp]*dphi[j][qp](1));  // (dF_2/dv) (v^(i+1))
                   // Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0); //libmesh example
                 //   Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);

                }
              for (unsigned int j= 0;j< n_p_dofs;j++)
                {
				  Kpp(i,j) += eps * psi[j][qp] * psi[i][qp]; // some regularisation term.
			    }
		    }
        } // end of the quadrature point qp-loop


      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method. The penalty method used here
      // is equivalent (for Lagrange basis functions) to lumpingqface
      // the matrix resulting from the L2 projection penalty
      // approach introduced in example 3.
      {
        // The penalty value.  \f$ \frac{1}{\epsilon} \f$
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == libmesh_nullptr)
            {
              UniquePtr<Elem> side (elem->build_side(s));
			
              // Loop over the nodes on the side.
              for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                {
                  // Boundary ids are set internally by
                  // build_square().
                  // 0=bottom
                  // 1=right
                  // 2=top
                  // 3=left

                  // Set u = 0  everywhere
               //   const Real u_value = 0.;
                    const Real u_value =
                            (mesh.get_boundary_info().has_boundary_id(elem, s, 2))
                            ? 1. : 0.;
                  // Set v = 0 everywhere
                  const Real v_value = 0.;
                    
                  // Find the node on the element matching this node on
                  // the side.  That defined where in the element matrix
                  // the boundary condition will be applied.
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                    if (elem->node_id(n) == side->node_id(ns))
                      {
                        // Matrix contribution.
                        Kuu(n,n) += penalty;
                        Kvv(n,n) += penalty;
						
                        // Right-hand-side contribution.
                        Fu(n) += penalty*u_value;
                        Fv(n) += penalty*v_value;
                      }                
                } // end face node loop
			  
            } // end if (elem->neighbor(side) == libmesh_nullptr)

//          const unsigned int pressure_node = 0;
//          const Real p_value               = 0.0;
//          for (unsigned int c=0; c<elem->n_nodes(); c++)
//              if (elem->node_id(c) == pressure_node)
//              {
//                  Kpp(c,c) += penalty;
//                  Fp(c)    += penalty*p_value;
//              }


      } // end boundary condition section

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);
       //~ //check for NaNs
	  //~ std::cout<<Ke.max()<<std::endl;
	  //~ std::cout<<" K : u-line "<<std::endl;
	  //~ Kuu.print();
	  //~ Kuv.print();
	  //~ Kup.print();
	  //~ KuT.print();
	  
	  //~ std::cout<<" K : v-line "<<std::endl;	  
	  //~ Kvu.print();
	  //~ Kvv.print();
	  //~ Kvp.print();
	  //~ KvT.print();
	  
	  //~ std::cout<<" K : p-line "<<std::endl;	  
	  //~ Kpu.print();
	  //~ Kpv.print();
	  //~ Kpp.print();
	  //~ KpT.print();

	  //~ std::cout<<" K : T-line "<<std::endl;	  
	  //~ KTu.print();
	  //~ KTv.print();
	  //~ KTp.print();
	  //~ KTT.print();
	  
	  //~ std::cout<<KTT.max()<<std::endl;
	  //~ std::cout<<Kpp.max()<<std::endl;
    } // end of element loop
    
}

void perform_implicit_euler (EquationSystems & es, const std::string & system_name)
{
  //Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("Navier-Stokes");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");
  const unsigned int T_var = navier_stokes_system.variable_number ("T");
  

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);
  
  // Get the Finite Element type for "T".
  FEType fe_temp_type = navier_stokes_system.variable_type(T_var);
  
  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  // Build a Finite Element object of the specified type for
  // the pressure variables.
  UniquePtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

  // Build a Finite Element object of the specified type for
  // the temperature variables.
  UniquePtr<FEBase> fe_temp (FEBase::build(dim, fe_temp_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  // Note that the QRule of the velocity (second order) is also used for 
  // the first order FE-Base of the pressure and temperature.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  fe_temp->attach_quadrature_rule (&qrule);

//########################### Boundary Integration Stuff ###########################
  // Finite element object for boundary integration
  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_temp_type));
  
  // Build a Quadrature Rule with one dimension less
  QGauss qface (dim-1, fe_temp_type.default_quadrature_order());
  
  // attach it
  fe_face->attach_quadrature_rule (&qface);
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe_vel->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dphi = fe_vel->get_dphi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & psi = fe_pres->get_phi();

  // The value of the linear shape function gradients at the quadrature points
  const std::vector<std::vector<RealGradient> > & dpsi = fe_pres->get_dphi();
  
  // The element shape functions for the temperature 
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & tau = fe_temp->get_phi();

  // The element shape function gradients for the temperature
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dtau = fe_temp->get_dphi();
  
  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = navier_stokes_system.get_dof_map();
  
  // the mass matrix for the lhs
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke), KuT(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke), KvT(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke), KpT(Ke),
    KTu(Ke), KTv(Ke), KTp(Ke), KTT(Ke);
    
  // the rhs
  DenseVector<Number> Fe;
  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe),
    FT(Fe); 

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_T;

  const Real dt    = es.parameters.get<Real>("dt");
  
  // get the parameters of the system
  const auto alpha = es.parameters.get<double>("alpha");
  const auto rho = es.parameters.get<double>("rho");
  const auto g_1 = es.parameters.get<double>("g_1");
  const auto g_2 = es.parameters.get<double>("g_2");
  const auto nu = es.parameters.get<double>("nu");
  const auto eps = es.parameters.get<double>("eps");
  const auto kappa = es.parameters.get<double>("kappa");
  const auto gamma = es.parameters.get<double>("gamma");
  const auto T_0 = es.parameters.get<double>("T_0");
  const auto T_Dir = es.parameters.get<double>("T_Dir");
  const auto T_out = es.parameters.get<double>("T_out");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_T, T_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();
      const unsigned int n_T_dofs = dof_indices_T.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);
	  fe_temp->reinit (elem);
	  fe_face->reinit (elem);
	  
	  // resize mass matrix
	  Ke.resize(n_dofs,n_dofs);
      Kuu.reposition (0, 0, n_u_dofs, n_u_dofs);
      Kuv.reposition (0, n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (0, n_u_dofs + n_v_dofs, n_u_dofs, n_p_dofs);
      KuT.reposition (0, n_u_dofs + n_v_dofs + n_p_dofs, n_u_dofs, n_T_dofs);
      
      Kvu.reposition (n_u_dofs, 0, n_v_dofs, n_u_dofs);
      Kvv.reposition (n_u_dofs, n_u_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (n_u_dofs, n_u_dofs + n_v_dofs,  n_v_dofs, n_p_dofs);
      KvT.reposition (n_u_dofs, n_u_dofs + n_v_dofs + n_p_dofs, n_v_dofs, n_T_dofs);
      
      Kpu.reposition (n_u_dofs + n_v_dofs, 0, n_p_dofs, n_u_dofs);
      Kpv.reposition (n_u_dofs + n_v_dofs, n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (n_u_dofs + n_v_dofs, n_u_dofs + n_v_dofs, n_p_dofs, n_p_dofs);
      KpT.reposition (n_u_dofs + n_v_dofs, n_u_dofs + n_v_dofs + n_p_dofs, n_p_dofs, n_T_dofs);
      
      KTu.reposition (n_u_dofs + n_v_dofs + n_p_dofs, 0, n_T_dofs, n_u_dofs);
      KTv.reposition (n_u_dofs + n_v_dofs + n_p_dofs, n_u_dofs, n_T_dofs, n_v_dofs);
      KTp.reposition (n_u_dofs + n_v_dofs + n_p_dofs, n_u_dofs + n_v_dofs, n_T_dofs, n_p_dofs);
      KTT.reposition (n_u_dofs + n_v_dofs + n_p_dofs, n_u_dofs + n_v_dofs + n_p_dofs, n_T_dofs, n_T_dofs);	
        
	  // Zero the vector and resize (element shape might have changed)
	  Fe.resize (n_dofs);
      Fu.reposition (0, n_u_dofs);
      Fv.reposition (n_u_dofs, n_v_dofs);
      Fp.reposition (n_u_dofs + n_v_dofs, n_p_dofs);
      FT.reposition (n_u_dofs + n_v_dofs + n_p_dofs, n_T_dofs);	  //~ std::cout<<" K : u-line "<<std::endl;
	  //~ Kuu.print();
	  //~ Kuv.print();
	  //~ Kup.print();
	  //~ KuT.print();
	  
	  //~ std::cout<<" K : v-line "<<std::endl;	  
	  //~ Kvu.print();
	  //~ Kvv.print();
	  //~ Kvp.print();
	  //~ KvT.print();
	  
	  //~ std::cout<<" K : p-line "<<std::endl;	  
	  //~ Kpu.print();
	  //~ Kpv.print();
	  //~ Kpp.print();
	  //~ KpT.print();

	  //~ std::cout<<" K : T-line "<<std::endl;	  
	  //~ KTu.print();
	  //~ KTv.print();
	  //~ KTp.print();
	  //~ KTT.print();
	  //~ std::cout<<" F "<<std::endl;
	  //~ Fu.print(libMesh::out);
	  //~ Fv.print(libMesh::out);
	  //~ Fp.print(libMesh::out);
	  //~ FT.print(libMesh::out);
      
      // Set up the Fe vector
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
		  Number u = 0., v = 0., p = 0., T = 0.;  // x_i (at "starting" time of Euler)
		  Gradient grad_u, grad_v, grad_p, grad_T;        // grad x_i (at "starting" time of Euler)
          const NumberVectorValue U (u,v);		  
          // Compute the velocity & its gradient at time i
          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));	
		    }	
		      
          // Compute the pressure value at this quadrature point at time i
          for (unsigned int l=0; l<n_p_dofs; l++)
            {
			  p += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);
			  grad_p.add_scaled (dpsi[l][qp],navier_stokes_system.current_solution (dof_indices_p[l]));
			} 
			
          // Compute the temperature & its gradient at time i
          for (unsigned int l=0; l<n_T_dofs; l++)
            {      
              T += tau[l][qp]*navier_stokes_system.current_solution (dof_indices_T[l]);
              grad_T.add_scaled (dtau[l][qp],navier_stokes_system.current_solution (dof_indices_T[l]));
            }
            
          // First, an i-loop over the velocity degrees of freedom.
          // We know that n_u_dofs == n_v_dofs so we can compute contributions
          // for both at the same time.
          // Note that, contrary to the systems_of_equations_ex2.C example
          // , we assume theta = 1 here, meaning we use implicit Euler.
          for (unsigned int i=0; i<n_u_dofs; i++)
            {
              Fu(i) += JxW[qp]*( u*phi[i][qp] +
								dt*( (-U*grad_u - grad_p(0)/rho - g_1*alpha*(T-T_0))* phi[i][qp]
									-nu*grad_u*dphi[i][qp]) );
              
              
              Fv(i) += JxW[qp]*( v*phi[i][qp] +
								dt*( (-U*grad_v - grad_p(1)/rho - g_2*alpha*(T-T_0))* phi[i][qp]
									-nu*grad_v*dphi[i][qp]) );
              
              // mass matrix for the lhs
              for (unsigned int j=0; j<n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp]*(phi[j][qp]*phi[i][qp]);
										
                  Kvv(i,j) += JxW[qp]*(phi[j][qp]*phi[i][qp]);                                           
                }
			}
			
		  // loop over the pressure dofs
          for (unsigned int i=0; i<n_p_dofs; i++)
            {
              Fp(i) += JxW[qp]*( grad_u(0) + grad_v(1))*psi[i][qp];
              
              // mass matrix for the lhs
              for (unsigned int j=0; j<n_p_dofs; j++)
                  Kpp(i,j) += JxW[qp]*(eps*psi[j][qp]*psi[i][qp]);               
		    }
		    
          // loop over the temperature dofs
          for (unsigned int i=0; i<n_T_dofs; i++)
            {
              FT(i) += JxW[qp]*( T*tau[i][qp] +
								dt*( -U*grad_T*tau[i][qp] - kappa*grad_T*dtau[i][qp]) );
              
              // mass matrix for the lhs
              for (unsigned int j=0; j<n_T_dofs; j++)
                  KTT(i,j) += JxW[qp]*(tau[j][qp]*tau[i][qp]);                                   

			}
		}// end qp loop
		
      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method. The penalty method used here
      // is equivalent (for Lagrange basis functions) to lumpingqface
      // the matrix resulting from the L2 projection penalty
      // approach introduced in example 3.
      {
        // The penalty value.  \f$ \frac{1}{\epsilon} \f$
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == libmesh_nullptr)
            {
              UniquePtr<Elem> side (elem->build_side(s));
			
              // Loop over the nodes on the side.
              for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                {
                  // Boundary ids are set internally by
                  // build_square().
                  // 0=bottom
                  // 1=right
                  // 2=top
                  // 3=left

                  // Set u = 0  everywhere 
                  const Real u_value = 0.;
                  // Set v = 0 everywhere
                  const Real v_value = 0.;
                    
                  // Find the node on the element matching this node on
                  // the side.  That defined where in the element matrix
                  // the boundary condition will be applied.
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                    if (elem->node_id(n) == side->node_id(ns))
                      {
                        // Matrix contribution.
                        Kuu(n,n) += penalty;
                        Kvv(n,n) += penalty;
						
                        // Right-hand-side contribution.
                        Fu(n) += penalty*u_value;
                        Fv(n) += penalty*v_value;
                        // Set T = T_Dir on the bottom boundary
                        if(mesh.get_boundary_info().has_boundary_id(elem, s, 0))
						  {
						    KTT(n,n) += penalty;
						    FT(n) += penalty*T_Dir;
						  }
                        
                      }                
                } // end face node loop
                
              // the boundary integrals  
              const std::vector<std::vector<Real> > & tau_face = fe_face->get_phi();
              const std::vector<std::vector<RealGradient> > & dtau_face = fe_face->get_dphi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();
			  //~ const std::vector<libMesh::Point>& normal_face = fe_face->get_normals();
			  //~ for(auto ele:*normal_face) std::cout<<ele<<std::endl;
              fe_face->reinit(elem,s);
              Number T = 0;
              Gradient grad_T;
              if(mesh.get_boundary_info().has_boundary_id(elem, s, 0)) // bottom Gamma_1
              {
			    for(unsigned int qp = 0; qp<qface.n_points(); qp++)
                  {
				    // Compute the temperature & its gradient at time i
				    for (unsigned int l=0; l<n_T_dofs; l++)
				  	  {      
					    T += tau[l][qp]*navier_stokes_system.current_solution (dof_indices_T[l]);
					    grad_T.add_scaled (dtau[l][qp],navier_stokes_system.current_solution (dof_indices_T[l]));
					  }				
				    for(unsigned int i = 0; i<n_T_dofs; i++)
					    FT(i) += -JxW[qp]*dt*kappa*grad_T(1)*tau_face[i][qp];
				  }
			  }
			  else // Gamma / Gamma_1, here we hav <grad T, eta > = f = gamma*(T-T_out)
			  {
			    for(unsigned int qp = 0; qp<qface.n_points(); qp++)
                  {
				    // Compute the temperature & its gradient at time i
				    for (unsigned int l=0; l<n_T_dofs; l++)
				  	  {      
					    T += tau[l][qp]*navier_stokes_system.current_solution (dof_indices_T[l]);
					    grad_T.add_scaled (dtau[l][qp],navier_stokes_system.current_solution (dof_indices_T[l]));
					  }	
				    for(unsigned int i = 0; i<n_T_dofs; i++)
					    FT(i) += JxW[qp]*dt*kappa*gamma*(T-T_out)*tau_face[i][qp];
				  }
			  }
			  
            } // end if (elem->neighbor(side) == libmesh_nullptr)

      } // end boundary condition section
      
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);
  }// end element loop
  
       //~ //check for NaNs
	  //~ std::cout<<" K : u-line "<<std::endl;
	  //~ Kuu.print();
	  //~ Kuv.print();
	  //~ Kup.print();
	  //~ KuT.print();
	  
	  //~ std::cout<<" K : v-line "<<std::endl;	  
	  //~ Kvu.print();
	  //~ Kvv.print();
	  //~ Kvp.print();
	  //~ KvT.print();
	  
	  //~ std::cout<<" K : p-line "<<std::endl;	  
	  //~ Kpu.print();
	  //~ Kpv.print();
	  //~ Kpp.print();
	  //~ KpT.print();

	  //~ std::cout<<" K : T-line "<<std::endl;	  
	  //~ KTu.print();
	  //~ KTv.print();
	  //~ KTp.print();
	  //~ KTT.print();
	  //~ std::cout<<" F "<<std::endl;
	  //~ Fu.print(libMesh::out);
	  //~ Fv.print(libMesh::out);
	  //~ Fp.print(libMesh::out);
	  //~ FT.print(libMesh::out);

}// end of function




