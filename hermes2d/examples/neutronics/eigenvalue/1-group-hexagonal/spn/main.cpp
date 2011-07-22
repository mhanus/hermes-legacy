#define HERMES_REPORT_ALL
#include "definitions.h"
#include "../problem_data.h"

using namespace RefinementSelectors;
using namespace WeakFormsNeutronics::Multigroup::MaterialProperties::SPN;

const unsigned int N_GROUPS = 1;  // Monoenergetic (single group) problem.
const unsigned int SPN_ORDER = 3; // SP3 approximation

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

const int INIT_REF_NUM[N_EQUATIONS] = {  // Initial uniform mesh refinement for the individual solution components.
  1, 1//, 2, 2                              
};
const int P_INIT[N_EQUATIONS] = {        // Initial polynomial orders for the individual solution components. 
  1, 1//, 2, 2                             
};      
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  
const bool HERMES_VISUALIZATION = true;  // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;     // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = true;       // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool INTERMEDIATE_VISUALIZATION = true; // Set to "true" to display coarse mesh solutions during adaptivity.

// Power iteration control.
double k_eff = 1.0;         // Initial eigenvalue approximation.
double TOL_PIT_CM = 1e-7;   // Tolerance for eigenvalue convergence on the coarse mesh.
double TOL_PIT_FM = 1e-8;   // Tolerance for eigenvalue convergence on the fine mesh.
bool PROJECT_INITIAL_SOLUTION = true;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
   
  MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, rm_map);
  matprop.set_nuSigma_f(nSf);
  matprop.set_nu(nu);
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
  
  matprop.validate();
  
  cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  H2DReader mloader;
  mloader.load((std::string("../") + mesh_file).c_str(), meshes[0]);
  
  // Convert the mesh so that it has one type of elements (optional). 
  //meshes[0]->convert_quads_to_triangles();
  meshes[0]->convert_triangles_to_quads();
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  
  SupportClasses::SPN::Views views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);

  // Create pointers to solutions on coarse and fine meshes and from the latest power iteration, respectively.
  Hermes::vector<Solution*> coarse_solutions, power_iterates;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    coarse_solutions.push_back(new Solution());
    power_iterates.push_back(new Solution(meshes[i], 1.0));   
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space *> spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces.push_back(new H1Space(meshes[i], P_INIT[i]));
  
  // Initialize the eigenvalue iterator.
  SupportClasses::SourceIteration si(NEUTRONICS_SPN, matprop, hermes2d, fission_regions);
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, SPN_ORDER, power_iterates, fission_regions, k_eff, bdy_vacuum);
  
  // Initialize the discrete algebraic representation of the problem.
  DiscreteProblem dp(&wf, spaces);
      
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  report_num_dof("Coarse mesh power iteration, NDOF: ", spaces);
  si.eigenvalue_iteration(power_iterates, dp, TOL_PIT_CM, matrix_solver);
  
  if (STRATEGY >= 0)
  {
    // DOF and CPU convergence graphs
    GnuplotGraph graph_dof("Error convergence", "NDOF", "log(error)");
    graph_dof.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_dof.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Error convergence", "CPU time [s]", "log(error)");
    graph_cpu.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_cpu.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    // Create pointers to coarse mesh solutions used for error estimation.
    Hermes::vector<Solution*> coarse_solutions;
    
    // Initialize all the new solution variables.
    for (unsigned int i = 0; i < N_EQUATIONS; i++) 
      coarse_solutions.push_back(new Solution());
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
    
    // Adaptivity loop:
    int as = 1; bool done = false;
    Hermes::vector<Space *>* fine_spaces;
    do 
    {
      info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      info("Solving on fine meshes.");
      
      fine_spaces = Space::construct_refined_spaces(spaces);
      
      if (PROJECT_INITIAL_SOLUTION)
      {
        if (as == 1)
          for (unsigned int i = 0; i < N_EQUATIONS; i++)  
            coarse_solutions[i]->copy(power_iterates[i]);
          
        report_num_dof("Projecting previous solutions onto new meshes, NDOF: ", *fine_spaces);
        OGProjection::project_global(*fine_spaces, coarse_solutions, power_iterates, matrix_solver);
      }
      
      // Solve the fine mesh problem.
      DiscreteProblem dp_ref(&wf, *fine_spaces);
      report_num_dof("Fine mesh power iteration, NDOF: ", *fine_spaces);
      si.eigenvalue_iteration(power_iterates, dp_ref, TOL_PIT_FM, matrix_solver);
            
      report_num_dof("Projecting fine mesh solutions on coarse meshes, NDOF: ", spaces);
      OGProjection::project_global(spaces, power_iterates, coarse_solutions, matrix_solver);
      
      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION && INTERMEDIATE_VISUALIZATION)
      {
        cpu_time.tick();
        info("Visualizing.");
        views.show_solutions(coarse_solutions);
        views.show_orders(spaces);
        cpu_time.tick(HERMES_SKIP);
      }
      
      Adapt adaptivity(spaces);
      
      info("Calculating errors.");
      Hermes::vector<double> h1_moment_errors;
      double h1_err_est = adaptivity.calc_err_est(coarse_solutions, power_iterates, &h1_moment_errors) * 100;
      
      // Time measurement.
      cpu_time.tick();
      double cta = cpu_time.accumulated();
      
      // Report results.
      
      // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
      double keff_err = 1e5*fabs(wf.get_keff() - REF_K_EFF)/REF_K_EFF;
      
      report_errors("odd moment err_est_coarse (H1): ", h1_moment_errors);
      info("total err_est_coarse (H1): %g%%", h1_err_est);
      info("k_eff err: %g milli-percent", keff_err);
      
      // Add entry to DOF convergence graph.
      int ndof_coarse = Space::get_num_dofs(spaces);
      graph_dof.add_values(0, ndof_coarse, h1_err_est);
      graph_dof.add_values(1, ndof_coarse, keff_err);
      
      // Add entry to CPU convergence graph.
      graph_cpu.add_values(0, cta, h1_err_est);
      graph_cpu.add_values(1, cta, keff_err);
            
      cpu_time.tick(HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (h1_err_est < ERR_STOP || as == MAX_ADAPT_NUM) 
        done = true;
      else 
      {
        info("Adapting the coarse meshes.");
        done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);        
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) 
          done = true;
      }
      
      if (!done)
      {
        for(unsigned int i = 0; i < N_EQUATIONS; i++)
          delete (*fine_spaces)[i]->get_mesh();
        delete fine_spaces;
        
        // Increase counter.
        as++;
      }
    }
    while (!done);
  }
  else
  {
    if (HERMES_VISUALIZATION)
    {
      views.show_solutions(power_iterates);
      views.show_orders(spaces);
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", power_iterates);
      views.save_orders_vtk("space", spaces);
    }
    
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double keff_err = 1e5*fabs(wf.get_keff() - REF_K_EFF)/REF_K_EFF;
    info("K_eff error = %g pcm", keff_err);
  }
  
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  View::wait(HERMES_WAIT_KEYPRESS);
  
  // Test the unit source normalization.
  SupportClasses::PostProcessor pp(NEUTRONICS_SPN);
  pp.normalize_to_unit_fission_source(&power_iterates, matprop);
  views.show_solutions(power_iterates);
  
  SupportClasses::SPN::SourceFilter sf(power_iterates, matprop, fission_regions);
  info("Total fission source by normalized flux: %g.", sf.integrate());
    
  // Wait for the view to be closed.  
  View::wait();
  
  return 0;
}
