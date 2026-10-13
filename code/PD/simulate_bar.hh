#ifndef SIMULATE_BAR_HH
#define SIMULATE_BAR_HH

#include "bar.hh"

#include <cstddef>
#include <vector>

struct Parameters;

struct SimulationResults {
  std::vector<double> t;

  std::vector<std::vector<double>> displacement;
  std::vector<std::vector<double>> damage;

  std::vector<std::size_t> num_fragments;
  std::vector<std::vector<double>> cuts;
  std::vector<std::vector<double>> sizes;

  std::vector<double> kinetic;
  std::vector<double> strain_energy;
  std::vector<double> dissipated;
  std::vector<double> external_work;

  std::vector<double> x;
  std::vector<double> health;

  std::size_t steps;
  bool stopped_early;
  double run_time;
  double energy_error;
};

/**
 * @brief Simulates the dynamic fragmentation of the peridynamic bar.
 *
 * The bar begins unstrained with the velocity field
 *
 * @f[
 * v(x) = \dot{\varepsilon}(x-L/2),
 * @f]
 *
 * corresponding to a uniform strain rate. Boundary points continue to follow
 * this prescribed motion throughout the simulation.
 *
 * Time integration is performed using velocity-Verlet. At each time step the
 * bond stretches are computed, irreversible bond damage is updated, forces are
 * recalculated, and the point velocities are advanced.
 *
 * The simulation stops either when the prescribed maximum simulation time is
 * reached or when the dissipated fracture energy has changed by less than the
 * specified tolerance over several recorded snapshots.
 *
 * @param params Simulation and material parameters.
 * @param bar Initialised peridynamic bar.
 *
 * @return Complete simulation history.
 */
SimulationResults simulateBar(const Parameters& params, const Bar& bar);

#endif