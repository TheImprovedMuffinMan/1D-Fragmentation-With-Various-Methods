#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <cstddef>
#include <string>

struct Parameters {
  // Material properties
  double density;
  double E;
  double strength;
  double toughness;

  // Bar geometry and discretisation
  double length;
  double area;
  double dx;
  double horizon;

  // Material model
  std::string bond_law;
  std::string bond_shape;

  // Random field
  double scatter;
  unsigned int seed;

  // Loading
  double strain_rate;

  // Time integration
  double time_step_factor;
  double max_time;

  // Stopping and output
  std::size_t dump_freq;
  std::size_t check_frequency;
  std::size_t stable_steps;

  std::string output_dir;
  std::string diagnostics_file;
  std::string fragment_file;
};

/**
 * @brief Reads the simulation parameters from a text file.
 *
 * The parameter file contains one parameter per line in the form
 *
 *     parameter = value
 *
 * Lines beginning with '#' and text following '#' are ignored as comments.
 *
 * @param filename Path to the parameter file.
 *
 * @return Parameters structure containing the values read from the file.
 */
Parameters readParameters(const std::string& filename);

#endif