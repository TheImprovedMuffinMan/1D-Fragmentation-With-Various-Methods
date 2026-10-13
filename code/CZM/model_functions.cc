#include <chrono>
#include <cmath>
#include <filesystem>
#include <fragment_manager.hh>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mesh.hh>
#include <solid_mechanics_model_cohesive.hh>
#include <stdexcept>
#include <string>
#include <vector>

#include "diagnostic_writer.hh"
#include "generate_model.hh"
#include "model_functions.hh"

using namespace akantu;

// ========================================================================== //

SimulationParameters readSimulationParameters() {
  const auto& user = getUserParser();

  SimulationParameters params;

  params.bar_length = user.getParameterValue<Real>("bar_length");
  params.number_elements = user.getParameterValue<UInt>("number_elements");
  params.strain_rate = user.getParameterValue<Real>("strain_rate");
  params.max_time = user.getParameterValue<Real>("max_time");
  params.dump_freq = user.getParameterValue<UInt>("dump_freq");
  params.output_dir = user.getParameterValue<std::string>("output_dir");
  params.output_file = user.getParameterValue<std::string>("output_file");
  params.fragment_file = user.getParameterValue<std::string>("fragment_file");
  params.cohesive_energy_tolerance =
      user.getParameterValue<Real>("cohesive_energy_tolerance");
  params.stable_steps = user.getParameterValue<UInt>("stable_steps");

  return params;
}

// ========================================================================== //
void writeSimulationDump(SolidMechanicsModelCohesive& model,
                         Mesh& mesh,
                         const SimulationParameters& params,
                         Int step,
                         Real dt,
                         Real applied_vel,
                         Real external_work,
                         Real E_tot_0,
                         Real reaction_left_current,
                         Real reaction_right_current) {
  const Real time = (step + 1) * dt;

  const std::string diagnostics_file =
      params.output_dir + "/" + params.output_file;

  const std::string fragments_file =
      params.output_dir + "/" + params.fragment_file;

  // Energies
  const Real E_pot = model.getEnergy("potential");
  const Real E_kin = model.getEnergy("kinetic");
  const Real E_rev = model.getEnergy("reversible");
  const Real E_dis = model.getEnergy("dissipated");
  const Real E_con = model.getEnergy("contact");

  const Real E_tot = E_pot + E_kin + E_rev + E_dis + E_con;

  const Real balance_error = external_work - (E_tot - E_tot_0);

  const Real displacement = 2.0 * applied_vel * time;
  const Real strain = displacement / params.bar_length;

  // Fragmentation analysis
  FragmentManager fragment_data(model);
  fragment_data.computeAllData();

  const UInt n_fragments = fragment_data.getNbFragment();

  const auto& nb_elements_per_fragment =
      fragment_data.getNbElementsPerFragment();

  const auto& fragment_mass = fragment_data.getMass();
  const auto& fragment_centre = fragment_data.getCenterOfMass();
  const auto& fragment_velocity = fragment_data.getVelocity();

  const Real h = params.bar_length / params.number_elements;

  std::vector<FragmentRecord> fragments;
  fragments.reserve(n_fragments);

  for (UInt i = 0; i < n_fragments; ++i) {
    FragmentRecord fragment{};

    fragment.fragment_id = static_cast<std::size_t>(i);
    fragment.size = nb_elements_per_fragment(i) * h;
    fragment.mass = fragment_mass(i);
    fragment.centre = fragment_centre(i, 0);
    fragment.mean_velocity = fragment_velocity(i, 0);

    fragments.push_back(fragment);
  }

  Real mean_fragment_size = 0.0;
  Real min_fragment_size = 0.0;
  Real max_fragment_size = 0.0;

  if (n_fragments > 0) {
    min_fragment_size = params.bar_length;

    for (UInt i = 0; i < n_fragments; ++i) {
      const Real fragment_size = nb_elements_per_fragment(i) * h;

      mean_fragment_size += fragment_size;

      if (fragment_size < min_fragment_size) {
        min_fragment_size = fragment_size;
      }

      if (fragment_size > max_fragment_size) {
        max_fragment_size = fragment_size;
      }
    }

    mean_fragment_size /= static_cast<Real>(n_fragments);
  }

  // Momentum
  const auto& mass = model.getMass();
  const auto& velocity = model.getVelocity();

  Real momentum = 0.0;

  for (UInt node = 0; node < mesh.getNbNodes(); ++node) {
    momentum += mass(node, 0) * velocity(node, 0);
  }

  // Stress and strain
  const auto& bulk_material = model.getMaterial(0);

  const auto& stress = bulk_material.getStress(_triangle_3);
  const auto& grad_u = bulk_material.getGradU(_triangle_3);

  const auto& element_filter = bulk_material.getElementFilter(_triangle_3);

  Array<Real> quadrature_coordinates(0, mesh.getSpatialDimension());

  model.getFEEngine().computeIntegrationPointsCoordinates(
      quadrature_coordinates, _triangle_3, _not_ghost, element_filter);

  std::vector<MaterialFieldRecord> material_fields;
  material_fields.reserve(stress.size());

  for (Int q = 0; q < stress.size(); ++q) {
    MaterialFieldRecord field{};

    field.x = quadrature_coordinates(q, 0);
    field.strain = grad_u(q, 0);
    field.stress = stress(q, 0);

    material_fields.push_back(field);
  }

  Real average_stress = 0.0;
  Real max_stress = std::numeric_limits<Real>::lowest();

  Real average_strain = 0.0;
  Real max_strain = std::numeric_limits<Real>::lowest();

  for (Int q = 0; q < stress.size(); ++q) {
    const Real stress_xx = stress(q, 0);
    const Real strain_xx = grad_u(q, 0);

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
  const auto& coordinates = mesh.getNodes();
  const auto& displacement_field = model.getDisplacement();
  const auto& velocity_field = model.getVelocity();
  const auto& acceleration_field = model.getAcceleration();

  std::vector<NodalFieldRecord> nodal_fields;
  nodal_fields.reserve(mesh.getNbNodes());

  for (UInt node = 0; node < mesh.getNbNodes(); ++node) {
    NodalFieldRecord field{};

    field.x = coordinates(node, 0);
    field.displacement = displacement_field(node, 0);
    field.velocity = velocity_field(node, 0);
    field.acceleration = acceleration_field(node, 0);
    field.damage = std::numeric_limits<double>::quiet_NaN();

    nodal_fields.push_back(field);
  }

  const std::string nodal_file =
      params.output_dir + "/fields/nodal_" + std::to_string(step + 1) + ".out";

  const std::string material_file = params.output_dir + "/fields/material_" +
                                    std::to_string(step + 1) + ".out";

  // Global diagnostics
  GlobalDiagnostics data{};

  data.step = static_cast<std::size_t>(step + 1);
  data.time = time;
  data.dt = dt;

  data.displacement = displacement;
  data.strain = strain;

  data.reaction_left = reaction_left_current;
  data.reaction_right = reaction_right_current;

  data.kinetic_energy = E_kin;
  data.elastic_energy = E_pot;
  data.reversible_cohesive_energy = E_rev;
  data.fracture_energy = E_dis;
  data.contact_energy = E_con;
  data.external_work = external_work;
  data.total_energy = E_tot;
  data.energy_error = balance_error;

  data.num_fragments = static_cast<std::size_t>(n_fragments);
  data.mean_fragment_size = mean_fragment_size;
  data.min_fragment_size = min_fragment_size;
  data.max_fragment_size = max_fragment_size;
  data.momentum = momentum;

  data.average_stress = average_stress;
  data.max_stress = max_stress;
  data.average_strain = average_strain;
  data.max_strain = max_strain;

  writeDiagnostics(diagnostics_file, data);
  writeFragments(fragments_file, time, fragments);
  writeNodalFieldSnapshot(nodal_file, nodal_fields);
  writeMaterialFieldSnapshot(material_file, material_fields);

  model.dump("cohesive elements");
}
// ========================================================================== //
void runSimulation(SolidMechanicsModelCohesive& model,
                   Mesh& mesh,
                   const SimulationParameters& params) {
  const auto wall_start = std::chrono::steady_clock::now();
  const Real dt_crit = model.getStableTimeStep();
  const Real dt = 0.1 * dt_crit;
  model.setTimeStep(dt);
  const Real applied_vel = params.strain_rate * params.bar_length * 0.5;

  // Initial energies
  const Real E_pot_0 = model.getEnergy("potential");
  const Real E_kin_0 = model.getEnergy("kinetic");
  const Real E_rev_0 = model.getEnergy("reversible");
  const Real E_dis_0 = model.getEnergy("dissipated");
  const Real E_con_0 = model.getEnergy("contact");

  const Real E_tot_0 = E_pot_0 + E_kin_0 + E_rev_0 + E_dis_0 + E_con_0;

  // Quantities needed to integrate the external work
  Real external_work = 0.0;
  Real reaction_left_previous = 0.0;
  Real reaction_right_previous = 0.0;

  // Boundary nodes
  const auto& nodes_left =
      mesh.getElementGroup("left").getNodeGroup().getNodes();

  const auto& nodes_right =
      mesh.getElementGroup("right").getNodeGroup().getNodes();

  // Safety checks for user input parameters
  if (params.dump_freq == 0) {
    throw std::runtime_error("dump_freq must be greater than zero.");
  }

  if (params.stable_steps == 0) {
    throw std::runtime_error("stable_steps must be greater than zero.");
  }

  // Output files
  const std::string diagnostics_file =
      params.output_dir + "/" + params.output_file;

  const std::string fragments_file =
      params.output_dir + "/" + params.fragment_file;

  initialiseDiagnostics(diagnostics_file);
  initialiseFragments(fragments_file);

  // Simulation loop variables
  // Real previous_dissipated_energy = 0.;
  // bool fragmentation_started = false;

  // UInt stable_steps = 0;
  // const UInt required_stable_steps = params.stable_steps;

  // Fragment-count convergence
  Int previous_nb_fragments = -1;
  UInt stable_fragment_checks = 0;

  const UInt check_frequency = 1000;
  const UInt required_stable_checks = params.stable_steps;

  Int step = 0;

  model.setBaseNameToDumper("cohesive elements", "cohesive");

  model.addDumpFieldToDumper("cohesive elements", "damage");

  model.addDumpFieldVectorToDumper("cohesive elements", "opening");

  model.addDumpFieldVectorToDumper("cohesive elements", "tractions");

  // simulation time
  double simulation_wall_time = 0.0;

  // max time check
  std::string stopping_reason;

  // Time integration
  // while (stable_steps < required_stable_steps && step * dt < params.max_time)
  // {
  while (stable_fragment_checks < required_stable_checks &&
         step * dt < params.max_time) {
    // Timing for the simulation step (excluding I/O)
    const auto step_start = std::chrono::steady_clock::now();
    // Apply boundary displacement corresponding to constant velocity
    applyVelocityBoundaryConditions(model, applied_vel, step, dt);

    // Check whether cohesive elements should be inserted
    model.checkCohesiveStress();

    // Advance one explicit time step
    model.solveStep("explicit_lumped");

    // Check convergence of dissipated cohesive energy
    // const Real dissipated_energy = model.getEnergy("dissipated");

    // if (dissipated_energy > 0.) {
    //   fragmentation_started = true;
    // }

    // if (fragmentation_started && previous_dissipated_energy > 0.) {
    //   const Real energy_change =
    //       std::abs(dissipated_energy - previous_dissipated_energy);

    //   const Real convergence_limit =
    //       params.cohesive_energy_tolerance * dissipated_energy;

    //   if (energy_change < convergence_limit) {
    //     ++stable_steps;
    //   } else {
    //     stable_steps = 0;
    //   }
    // }
    // previous_dissipated_energy = dissipated_energy;

    // Check whether the fragment count has stabilised
    if ((step + 1) % check_frequency == 0) {
      FragmentManager fragment_count(model);
      fragment_count.computeAllData();

      const Int nb_fragments = fragment_count.getNbFragment();

      if (nb_fragments == previous_nb_fragments && nb_fragments > 1) {
        ++stable_fragment_checks;
      } else {
        stable_fragment_checks = 0;
      }

      previous_nb_fragments = nb_fragments;

      std::cout << "Step: " << step + 1 << " -- fragments = " << nb_fragments
                << " -- stable checks = " << stable_fragment_checks << "/"
                << required_stable_checks << "\n";
    }

    // External work
    const auto& internal_force = model.getInternalForce();

    Real reaction_left_current = 0.0;
    Real reaction_right_current = 0.0;

    for (auto node : nodes_left) {
      reaction_left_current -= internal_force(node, 0);
    }

    for (auto node : nodes_right) {
      reaction_right_current -= internal_force(node, 0);
    }

    const Real reaction_left =
        0.5 * (reaction_left_previous + reaction_left_current);

    const Real reaction_right =
        0.5 * (reaction_right_previous + reaction_right_current);

    external_work +=
        (reaction_left * (-applied_vel) + reaction_right * applied_vel) * dt;

    reaction_left_previous = reaction_left_current;
    reaction_right_previous = reaction_right_current;

    // Timing for the simulation step (excluding I/O)
    const auto step_end = std::chrono::steady_clock::now();

    simulation_wall_time +=
        std::chrono::duration<double>(step_end - step_start).count();
    if (step % params.dump_freq == 0) {
      writeSimulationDump(model,
                          mesh,
                          params,
                          step,
                          dt,
                          applied_vel,
                          external_work,
                          E_tot_0,
                          reaction_left_current,
                          reaction_right_current);
    }

    ++step;
  }

  // if (stable_steps >= required_stable_steps) {
  //   stopping_reason = "dissipated_energy_converged";
  // } else {
  //   stopping_reason = "maximum_time_reached";
  // }

  if (stable_fragment_checks >= required_stable_checks) {
    stopping_reason = "fragment_count_converged";
  } else {
    stopping_reason = "maximum_time_reached";
  }

  // if (stopping_reason == "dissipated_energy_converged") {
  //   std::cout << "Dissipated cohesive energy converged for "
  //             << required_stable_steps << " consecutive steps.\n";
  // } else {
  //   std::cout << "Maximum simulation time reached.\n";
  // }
  // std::cout << "Final time: " << step * dt << "\n";

  if (stable_fragment_checks >= required_stable_checks) {
    stopping_reason = "fragment_count_converged";
  } else {
    stopping_reason = "maximum_time_reached";
  }

  if (stopping_reason == "fragment_count_converged") {
    std::cout << "Fragment count stable for " << required_stable_checks
              << " consecutive checks.\n";
  } else {
    std::cout << "Maximum simulation time reached.\n";
  }

  std::cout << "Final time: " << step * dt << "\n";

  FragmentManager fragment_data(model);
  fragment_data.computeAllData();
  const UInt n_fragments = fragment_data.getNbFragment();

  const Real E_pot_final = model.getEnergy("potential");
  const Real E_kin_final = model.getEnergy("kinetic");
  const Real E_rev_final = model.getEnergy("reversible");
  const Real E_dis_final = model.getEnergy("dissipated");
  const Real E_con_final = model.getEnergy("contact");

  const Real E_tot_final =
      E_pot_final + E_kin_final + E_rev_final + E_dis_final + E_con_final;

  const Real final_energy_error = external_work - (E_tot_final - E_tot_0);

  const auto wall_end = std::chrono::steady_clock::now();

  const double wall_time =
      std::chrono::duration<double>(wall_end - wall_start).count();

  const std::string summary_file = params.output_dir + "/summary.out";

  RunSummary summary{};

  summary.physical_time = step * dt;
  summary.wall_time = wall_time;
  summary.simulation_wall_time = simulation_wall_time;
  summary.steps = static_cast<std::size_t>(step);
  summary.num_fragments = static_cast<std::size_t>(n_fragments);
  summary.energy_error = final_energy_error;
  summary.stopping_reason = stopping_reason;

  writeRunSummary(summary_file, summary);

  std::cout << "Number of fragments: " << n_fragments << "\n";
}