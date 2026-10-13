#ifndef MAKE_BAR_HH
#define MAKE_BAR_HH

#include "bar.hh"

struct Parameters;

/**
 * @brief Builds the peridynamic bar, its points, bonds, stiffness and strength.
 *
 * The bar consists of equally spaced material points. Each point interacts
 * with every neighbouring point lying within the peridynamic horizon.
 *
 * Additional ghost points are placed beyond both ends of the bar and are used
 * to impose the boundary conditions.
 *
 * The function computes the bond geometry, partial-volume correction,
 * weighting, stiffness, strength, breaking stretch and stable time step.
 *
 * @param params Simulation and material parameters.
 *
 * @return Fully initialised Bar structure.
 */
Bar makeBar(const Parameters& params);
#endif