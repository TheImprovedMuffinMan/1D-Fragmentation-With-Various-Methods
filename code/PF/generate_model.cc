#include "generate_model.hh"

#include <non_linear_solver_tao.hh>
#include <solver_vector_petsc.hh>

using namespace akantu;

// -------------------------------------------------------------------------- //
void applyInitialVelocity(SolidMechanicsModel& solid, Real strain_rate) {
  auto& velocities = solid.getVelocity();
  const auto& positions = solid.getMesh().getNodes();

  for (Int n = 0; n < solid.getMesh().getNbNodes(); ++n) {
    velocities(n, 0) = strain_rate * positions(n, 0);
  }
}

// -------------------------------------------------------------------------- //
void applyInitialDisplacement(SolidMechanicsModel& solid,
                              Real displacement,
                              Real bar_length) {
  auto& u = solid.getDisplacement();
  const auto& positions = solid.getMesh().getNodes();

  for (Int n = 0; n < solid.getMesh().getNbNodes(); ++n) {
    u(n, 0) = 2. * displacement * positions(n, 0) / bar_length;
  }

  solid.applyBC(BC::Dirichlet::FixedValue(displacement, _x), "right");
  solid.applyBC(BC::Dirichlet::FixedValue(-displacement, _x), "left");
}

// -------------------------------------------------------------------------- //
void initialiseModel(CouplerSolidPhaseField& model,
                     const PhaseFieldParameters& params) {
  auto& solid = model.getSolidMechanicsModel();
  auto& phase = model.getPhaseFieldModel();

  // Solid-mechanics model
  solid.initFull(_analysis_method = _static);

  auto& solver = solid.getNonLinearSolver("static");

  solver.set("max_iterations", 2000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  solid.initNewSolver(_explicit_lumped_mass);

  // Phase-field model
  PetscOptionsSetValue(NULL, "-tao_type", "gpcg");

  phase.initDOFManager("petsc");

  phase.initFull(_analysis_method = _static,
                 _solver_options = ModelSolverOptions{
                     .non_linear_solver_type = NonLinearSolverType::_petsc_tao,
                     .sparse_solver_type = SparseSolverType::_petsc});

  solid.setBaseNameToDumper("solid_mechanics_model", "phasefield_solid");
  phase.setBaseNameToDumper("phase_field", "phasefield_damage");
  phase.addDumpFieldToDumper("phase_field", "damage");
  phase.addDumpFieldToDumper("phase_field", "blocked_dofs");

  if (model.getMesh().isDistributed()) {
    phase.addDumpFieldToDumper("phase_field", "partitions");
  }

  phase.dump("phase_field");

  // Boundary conditions
  solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "Yblocked");

  phase.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");

  phase.applyBC(BC::Dirichlet::FixedValue(0., _x), "right");

  // Optional initial displacement
  if (params.increment_displacement) {
    applyInitialDisplacement(solid, params.displacement, params.bar_length);

    solid.solveStep("static");
  }
}