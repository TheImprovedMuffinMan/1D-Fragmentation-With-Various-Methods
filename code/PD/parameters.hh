#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <string>

struct Parameters {
  // Material properties
  double density;
  double E;
  double strength;
  double toughness;

  // Bar geometry
  double length;
  double dx;
  double m;

  // Material model
  std::string bond_law;
  std::string bond_shape;

  // Random field
  double scatter;
  unsigned int seed;

  // Loading
  double strain_rate;
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