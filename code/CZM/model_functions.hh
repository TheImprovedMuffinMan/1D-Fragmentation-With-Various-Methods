#pragma once

#include <aka_common.hh>
#include <mesh.hh>
#include <solid_mechanics_model_cohesive.hh>

#include <string>
#include <vector>

/**
 * @brief Structure to hold the simulation parameters.
 *
 * This structure contains the parameters that are read from the user input
 * and used to configure the simulation.
 */
struct SimulationParameters {
  akantu::Real bar_length;
  akantu::UInt number_elements;
  akantu::Real strain_rate;
  akantu::UInt dump_freq;
  akantu::Real cohesive_energy_tolerance;
  akantu::UInt stable_steps;
  std::string output_dir;
  std::string output_file;
  std::string fragment_file;
  akantu::Real max_time;
};

/**
 * @brief Reads the simulation parameters from the user parser and returns a
 * SimulationParameters structure.
 *
 * This function reads the simulation parameters from the user parser and
 * populates a SimulationParameters structure with the values. The structure is
 * then returned to the caller.
 *
 * @return SimulationParameters structure containing the parsed parameters.
 */
SimulationParameters readSimulationParameters();

/**
 * @brief Runs the simulation for the given model and mesh with the specified
 * parameters.
 *
 * This function runs the simulation for the given cohesive zone model and mesh
 * using the specified simulation parameters. It performs time integration,
 * computes diagnostics, and writes output files as specified in the parameters.
 *
 * @param model Corresponding to the cohesive zone model.
 * @param mesh Mesh of the model.
 * @param params Simulation parameters containing the bar length, strain rate,
 * total simulation time, dump frequency, output directory, and output file
 * name.
 */
void runSimulation(akantu::SolidMechanicsModelCohesive& model,
                   akantu::Mesh& mesh,
                   const SimulationParameters& params);

/**
 * @brief Writes a simulation dump for the given model and mesh at the specified
 * time step.
 *
 * This function writes a simulation dump for the given cohesive zone model and
 * mesh at the specified time step. The dump includes relevant data such as
 * nodal displacements, velocities, and other field variables.
 *
 * @param model Corresponding to the cohesive zone model.
 * @param mesh Mesh of the model.
 * @param params Simulation parameters containing the bar length, strain rate,
 * total simulation time, dump frequency, output directory, and output file
 * name.
 * @param step Current time step of the simulation.
 * @param dt Time step size.
 * @param applied_vel Applied velocity at the boundaries.
 * @param external_work External work done on the system.
 * @param E_tot_0 Total energy at the initial time step.
 * @param reaction_left_current Current reaction force at the left boundary.
 * @param reaction_right_current Current reaction force at the right boundary.
 */
void writeSimulationDump(akantu::SolidMechanicsModelCohesive& model,
                         akantu::Mesh& mesh,
                         const SimulationParameters& params,
                         akantu::Int step,
                         akantu::Real dt,
                         akantu::Real applied_vel,
                         akantu::Real external_work,
                         akantu::Real E_tot_0,
                         akantu::Real reaction_left_current,
                         akantu::Real reaction_right_current);