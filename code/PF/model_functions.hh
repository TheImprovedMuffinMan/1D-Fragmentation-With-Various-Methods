#pragma once
#include "coupler_solid_phasefield.hh"
#include <aka_common.hh>
#include <string>
#include <vector>

/**
 * @brief Structure to hold the parameters for the phase field simulation.
 */
struct PhaseFieldParameters {
  akantu::Real bar_length;
  akantu::UInt number_elements;
  akantu::Real strain_rate;
  akantu::Real displacement;
  bool increment_displacement;
  akantu::Real damage_limit;
  akantu::Real time_step_factor;
  akantu::Real phase_energy_tolerance;
  akantu::UInt stable_steps;
  akantu::UInt dump_freq;
  std::string output_dir;
  std::string diagnostics_file;
  std::string fragment_file;
};

/**
 * @brief Read the simulation parameters from the user input.
 *
 * @return PhaseFieldParameters The parameters for the phase field simulation.
 */
PhaseFieldParameters readSimulationParameters();

/**
 * @brief Write a snapshot of the simulation state to output files.
 *
 * @param model The coupled solid-phase field model.
 * @param params The parameters for the phase field simulation.
 * @param step The current simulation step.
 * @param time The current simulation time.
 * @param dt The current time step size.
 * @param E_tot_0 The initial total energy of the system.
 */
void writeSimulationDump(akantu::CouplerSolidPhaseField& model,
                         const PhaseFieldParameters& params,
                         akantu::Int step,
                         akantu::Real time,
                         akantu::Real dt,
                         akantu::Real E_tot_0);

/**
 * @brief Run the simulation for the given model and parameters.
 *
 * @param model The coupled solid-phase field model to run the simulation on.
 * @param params The parameters for the phase field simulation.
 */
void runSimulation(akantu::CouplerSolidPhaseField& model,
                   const PhaseFieldParameters& params);
