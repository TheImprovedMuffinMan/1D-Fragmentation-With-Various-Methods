#include <aka_common.hh>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <mesh.hh>
#include <node_synchronizer.hh>
#include <non_linear_solver_tao.hh>
#include <solver_vector_petsc.hh>
#include <stdexcept>
#include <string>
#include <vector>

#include "diagnostic_writer.hh"
#include "generate_model.hh"
#include "model_functions.hh"
#include "phasefield_fragments.hh"

using namespace akantu;

// ========================================================================== //

PhaseFieldParameters readSimulationParameters() {
  const auto& user = getUserParser();

  PhaseFieldParameters params;

  params.bar_length = user.getParameterValue<Real>("bar_length");
  params.number_elements = user.getParameterValue<UInt>("number_elements");
  params.strain_rate = user.getParameterValue<Real>("strain_rate");

  params.displacement = user.getParameterValue<Real>("displacement");
  params.increment_displacement =
      user.getParameterValue<bool>("increment_displacement");

  params.damage_limit = user.getParameterValue<Real>("damage_limit");
  params.time_step_factor = user.getParameterValue<Real>("time_step_factor");

  params.dump_freq = user.getParameterValue<UInt>("dump_freq");

  params.output_dir = user.getParameterValue<std::string>("output_dir");
  params.diagnostics_file =
      user.getParameterValue<std::string>("diagnostics_file");
  params.fragment_file = user.getParameterValue<std::string>("fragment_file");

  params.phase_energy_tolerance =
      user.getParameterValue<Real>("phase_energy_tolerance");
  params.stable_steps = user.getParameterValue<UInt>("stable_steps");
  params.max_time = user.getParameterValue<Real>("max_time");

  return params;
}

// ========================================================================== //
void writeSimulationDump(CouplerSolidPhaseField& model,
                         const PhaseFieldParameters& params,
                         Int step,
                         Real time,
                         Real dt,
                         Real E_tot_0) {
  auto& solid = model.getSolidMechanicsModel();
  auto& phase = model.getPhaseFieldModel();
  auto& mesh = solid.getMesh();

  const std::string diagnostics_file =
      params.output_dir + "/" + params.diagnostics_file;

  const std::string fragments_file =
      params.output_dir + "/" + params.fragment_file;

  // Energies
  const Real E_pot = solid.getEnergy("potential");
  const Real E_kin = solid.getEnergy("kinetic");
  const Real E_phase = phase.getEnergy();

  const Real E_tot = E_pot + E_kin + E_phase;

  const Real W_ext = solid.getEnergy("external_work");
  const Real balance_error = W_ext - (E_tot - E_tot_0);

  // Global displacement and strain
  const Real displacement = params.strain_rate * params.bar_length * time;
  const Real strain = displacement / params.bar_length;

  // Fragmentation analysis
  PhaseFieldFragmentManager fragment_data(
      phase, params.damage_limit, "fragments");

  const Int n_fragments = fragment_data.countFragments("fragment");

  const auto& coordinates = mesh.getNodes();
  const auto& mass = solid.getMass();
  const auto& velocity = solid.getVelocity();

  std::vector<FragmentRecord> fragments;
  fragments.reserve(n_fragments);

  Real mean_fragment_size = 0.0;
  Real min_fragment_size = 0.0;
  Real max_fragment_size = 0.0;

  if (n_fragments > 0) {
    min_fragment_size = params.bar_length;
  }

  for (Int i = 0; i < n_fragments; ++i) {
    auto& fragment_group = fragment_data.getFragment(i);
    const auto& nodes = fragment_group.getNodeGroup().getNodes();

    Real x_min = std::numeric_limits<Real>::max();
    Real x_max = std::numeric_limits<Real>::lowest();

    Real fragment_mass = 0.0;
    Real weighted_centre = 0.0;
    Real weighted_velocity = 0.0;

    for (auto node : nodes) {
      const Real x = coordinates(node, 0);
      const Real node_mass = mass(node, 0);

      x_min = std::min(x_min, x);
      x_max = std::max(x_max, x);

      fragment_mass += node_mass;
      weighted_centre += node_mass * x;
      weighted_velocity += node_mass * velocity(node, 0);
    }

    const Real fragment_size = x_max - x_min;

    Real fragment_centre = 0.0;
    Real fragment_velocity = 0.0;

    if (fragment_mass > 0.0) {
      fragment_centre = weighted_centre / fragment_mass;
      fragment_velocity = weighted_velocity / fragment_mass;
    }

    FragmentRecord fragment{};

    fragment.fragment_id = static_cast<std::size_t>(i);
    fragment.size = fragment_size;
    fragment.mass = fragment_mass;
    fragment.centre = fragment_centre;
    fragment.mean_velocity = fragment_velocity;

    fragments.push_back(fragment);

    mean_fragment_size += fragment_size;

    if (fragment_size < min_fragment_size) {
      min_fragment_size = fragment_size;
    }

    if (fragment_size > max_fragment_size) {
      max_fragment_size = fragment_size;
    }
  }

  if (n_fragments > 0) {
    mean_fragment_size /= static_cast<Real>(n_fragments);
  }

  // Momentum
  Real momentum = 0.0;

  for (Int node = 0; node < mesh.getNbNodes(); ++node) {
    momentum += mass(node, 0) * velocity(node, 0);
  }

  // Stress and strain
  const auto& bulk_material = solid.getMaterial(0);

  const auto& stress = bulk_material.getStress(_triangle_3);
  const auto& grad_u = bulk_material.getGradU(_triangle_3);

  const auto& element_filter = bulk_material.getElementFilter(_triangle_3);

  Array<Real> quadrature_coordinates(0, mesh.getSpatialDimension());

  solid.getFEEngine().computeIntegrationPointsCoordinates(
      quadrature_coordinates, _triangle_3, _not_ghost, element_filter);

  std::vector<MaterialFieldRecord> material_fields;
  material_fields.reserve(stress.size());

  Real average_stress = 0.0;
  Real max_stress = std::numeric_limits<Real>::lowest();

  Real average_strain = 0.0;
  Real max_strain = std::numeric_limits<Real>::lowest();

  for (Int q = 0; q < stress.size(); ++q) {
    const Real stress_xx = stress(q, 0);
    const Real strain_xx = grad_u(q, 0);

    MaterialFieldRecord field{};

    field.x = quadrature_coordinates(q, 0);
    field.strain = strain_xx;
    field.stress = stress_xx;

    material_fields.push_back(field);

    average_stress += stress_xx;
    average_strain += strain_xx;

    if (stress_xx > max_stress) {
      max_stress = stress_xx;
    }

    if (strain_xx > max_strain) {
      max_strain = strain_xx;
    }
  }

  if (stress.size() > 0) {
    average_stress /= static_cast<Real>(stress.size());
    average_strain /= static_cast<Real>(stress.size());
  } else {
    average_stress = 0.0;
    max_stress = 0.0;
    average_strain = 0.0;
    max_strain = 0.0;
  }

  // Nodal fields
  const auto& displacement_field = solid.getDisplacement();
  const auto& acceleration_field = solid.getAcceleration();
  const auto& damage_field = phase.getDamage();

  std::vector<NodalFieldRecord> nodal_fields;
  nodal_fields.reserve(mesh.getNbNodes());

  for (Int node = 0; node < mesh.getNbNodes(); ++node) {
    NodalFieldRecord field{};

    field.x = coordinates(node, 0);
    field.displacement = displacement_field(node, 0);
    field.velocity = velocity(node, 0);
    field.acceleration = acceleration_field(node, 0);
    field.damage = damage_field(node);

    nodal_fields.push_back(field);
  }

  // Boundary reactions
  const auto& internal_force = solid.getInternalForce();

  const auto& nodes_left =
      mesh.getElementGroup("left").getNodeGroup().getNodes();

  const auto& nodes_right =
      mesh.getElementGroup("right").getNodeGroup().getNodes();

  Real reaction_left = 0.0;
  Real reaction_right = 0.0;

  for (auto node : nodes_left) {
    reaction_left -= internal_force(node, 0);
  }

  for (auto node : nodes_right) {
    reaction_right -= internal_force(node, 0);
  }

  // Field output files
  const std::string nodal_file =
      params.output_dir + "/fields/nodal_" + std::to_string(step) + ".out";

  const std::string material_file =
      params.output_dir + "/fields/material_" + std::to_string(step) + ".out";

  // Global diagnostics
  GlobalDiagnostics data{};

  data.step = static_cast<std::size_t>(step);
  data.time = time;
  data.dt = dt;

  data.displacement = displacement;
  data.strain = strain;

  data.reaction_left = reaction_left;
  data.reaction_right = reaction_right;

  data.kinetic_energy = E_kin;
  data.elastic_energy = E_pot;
  data.reversible_cohesive_energy = 0.0;
  data.fracture_energy = E_phase;
  data.contact_energy = 0.0;
  data.external_work = W_ext;
  data.total_energy = E_tot;
  data.energy_error = balance_error;

  data.num_fragments = static_cast<std::size_t>(n_fragments);
  data.mean_fragment_size = mean_fragment_size;
  data.min_fragment_size = min_fragment_size;
  data.max_fragment_size = max_fragment_size;

  data.average_stress = average_stress;
  data.max_stress = max_stress;
  data.average_strain = average_strain;
  data.max_strain = max_strain;

  data.momentum = momentum;

  // Write output
  writeDiagnostics(diagnostics_file, data);
  writeFragments(fragments_file, time, fragments);
  writeNodalFieldSnapshot(nodal_file, nodal_fields);
  writeMaterialFieldSnapshot(material_file, material_fields);

  phase.dump("phase_field");
}

// ========================================================================== //
void runSimulation(CouplerSolidPhaseField& model,
                   const PhaseFieldParameters& params) {
  const auto wall_start = std::chrono::steady_clock::now();

  const auto& comm = Communicator::getStaticCommunicator();
  const Int prank = comm.whoAmI();

  auto& solid = model.getSolidMechanicsModel();
  auto& phase = model.getPhaseFieldModel();

  auto& DOFManager = aka::as_type<DOFManagerPETSc>(phase.getDOFManager());

  auto& tao_solver =
      aka::as_type<NonLinearSolverTAO>(phase.getNonLinearSolver());

  const Real dt = solid.getStableTimeStep() * params.time_step_factor;
  solid.setTimeStep(dt);

  applyInitialVelocity(solid, params.strain_rate);

  // Initial energy after imposing the initial velocity field
  const Real E_pot_0 = solid.getEnergy("potential");
  const Real E_kin_0 = solid.getEnergy("kinetic");
  const Real E_phase_0 = phase.getEnergy();

  const Real E_tot_0 = E_pot_0 + E_kin_0 + E_phase_0;

  // Phase-field bounds
  SolverVectorPETSc lower_bound(DOFManager, "lower_bound");
  lower_bound.resize();
  lower_bound.set(0.);

  SolverVectorPETSc upper_bound(DOFManager, "upper_bound");
  upper_bound.resize();
  upper_bound.set(1.);

  tao_solver.setBounds(lower_bound, upper_bound);

  if (prank == 0) {
    std::cout << "Time step: " << dt << "\n";
    std::cout << "Strain rate: " << params.strain_rate << "\n";
  }

  // Safety checks
  if (params.dump_freq == 0) {
    throw std::runtime_error("dump_freq must be greater than zero.");
  }

  if (params.stable_steps == 0) {
    throw std::runtime_error("stable_steps must be greater than zero.");
  }

  // Output files
  const std::string diagnostics_file =
      params.output_dir + "/" + params.diagnostics_file;

  const std::string fragments_file =
      params.output_dir + "/" + params.fragment_file;

  if (prank == 0) {
    std::filesystem::create_directories(params.output_dir);
    std::filesystem::create_directories(params.output_dir + "/fields");

    initialiseDiagnostics(diagnostics_file);
    initialiseFragments(fragments_file);
  }

  // Fragment-count convergence
  Int previous_nb_fragments = -1;
  UInt stable_fragment_checks = 0;

  const UInt check_frequency = params.number_elements;
  const UInt required_stable_checks = params.stable_steps;

  // Simulation variables
  Real time = 0.0;
  Int step = 0;

  double simulation_wall_time = 0.0;

  // Time integration
  std::string stopping_reason;
  while (stable_fragment_checks < required_stable_checks &&
         time < params.max_time) {
    const auto step_start = std::chrono::steady_clock::now();

    const Real current_du = 0.5 * params.bar_length * params.strain_rate * dt;

    solid.setTimeStep(dt);

    if (params.increment_displacement) {
      solid.applyBC(BC::Dirichlet::IncrementValue(current_du, _x), "right");

      solid.applyBC(BC::Dirichlet::IncrementValue(-current_du, _x), "left");
    }

    model.solve("explicit_lumped", "static");

    time += dt;
    ++step;

    // Update phase-field bounds
    lower_bound.set(0.);
    upper_bound.set(PETSC_INFINITY);
    tao_solver.setBounds(lower_bound, upper_bound);

    // Check whether the fragment count has stabilised
    if (step % check_frequency == 0) {
      PhaseFieldFragmentManager fragment_count(
          phase, params.damage_limit, "check");

      const Int nb_fragments = fragment_count.countFragments("check");

      if (nb_fragments == previous_nb_fragments && nb_fragments > 1) {
        ++stable_fragment_checks;
      } else {
        stable_fragment_checks = 0;
      }

      previous_nb_fragments = nb_fragments;

      if (prank == 0) {
        std::cout << "Step: " << step << " -- fragments = " << nb_fragments
                  << " -- stable checks = " << stable_fragment_checks << "/"
                  << required_stable_checks << "\n";
      }
    }

    const auto step_end = std::chrono::steady_clock::now();

    simulation_wall_time +=
        std::chrono::duration<double>(step_end - step_start).count();

    // Periodic diagnostics
    if (step % params.dump_freq == 0 && prank == 0) {
      writeSimulationDump(model, params, step, time, dt, E_tot_0);
    }
  }

  // if (prank == 0) {
  //   std::cout << "Fragment count stable for " << required_stable_checks
  //             << " consecutive checks.\n";

  //   std::cout << "Final time: " << time << "\n";
  // }
  if (stable_fragment_checks >= required_stable_checks) {
    stopping_reason = "fragment_count_converged";
  } else {
    stopping_reason = "maximum_time_reached";
  }

  if (prank == 0) {
    if (stopping_reason == "fragment_count_converged") {
      std::cout << "Fragment count stable for " << required_stable_checks
                << " consecutive checks.\n";
    } else {
      std::cout << "Maximum simulation time reached.\n";
    }

    std::cout << "Final time: " << time << "\n";
  }

  // Final fragment count
  PhaseFieldFragmentManager final_fragments(
      phase, params.damage_limit, "final_fragments");

  const Int n_fragments = final_fragments.countFragments("final_fragment");

  // Final energy balance
  const Real E_pot_final = solid.getEnergy("potential");
  const Real E_kin_final = solid.getEnergy("kinetic");
  const Real E_phase_final = phase.getEnergy();

  const Real E_tot_final = E_pot_final + E_kin_final + E_phase_final;

  const Real W_ext_final = solid.getEnergy("external_work");

  const Real final_energy_error = W_ext_final - (E_tot_final - E_tot_0);

  const auto wall_end = std::chrono::steady_clock::now();

  const double wall_time =
      std::chrono::duration<double>(wall_end - wall_start).count();

  // Run summary
  if (prank == 0) {
    const std::string summary_file = params.output_dir + "/summary.out";

    RunSummary summary{};

    summary.physical_time = time;
    summary.wall_time = wall_time;
    summary.simulation_wall_time = simulation_wall_time;
    summary.steps = static_cast<std::size_t>(step);
    summary.num_fragments = static_cast<std::size_t>(n_fragments);
    summary.energy_error = final_energy_error;
    summary.stopping_reason = stopping_reason;

    writeRunSummary(summary_file, summary);

    std::cout << "Number of fragments: " << n_fragments << "\n";
  }
}