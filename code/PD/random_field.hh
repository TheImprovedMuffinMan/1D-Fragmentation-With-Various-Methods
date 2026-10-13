#ifndef RANDOM_FIELD_HH
#define RANDOM_FIELD_HH

#include <vector>

struct Parameters;

struct RandomFieldData {
  std::vector<double> strength;
  std::vector<double> toughness;
};

/**
 * @brief Generates spatially correlated random strength and toughness fields.
 *
 * Each material point receives a strength and toughness that vary randomly
 * around the mean values specified in @p params.
 *
 * The standard deviation of each field is controlled by @p params.scatter.
 * A scatter of zero produces a perfectly homogeneous material.
 *
 * The random variation is spatially smoothed over the correlation length
 *
 * @f[
 * l_c = \frac{E G_c}{\sigma_c^2}.
 * @f]
 *
 * The random values are first generated on a helper grid independent of the
 * discretisation spacing and are then interpolated onto the material points.
 * Therefore, changing @p dx does not redefine the underlying random field.
 *
 * @param params Material properties, scatter and random seed.
 * @param x Positions of all material and boundary points.
 * @param dx Material-point spacing.
 * @param horizon Peridynamic horizon.
 * @param L Bar length.
 *
 * @return Spatially varying strength and toughness fields.
 */
RandomFieldData randomField(const Parameters& params,
                            const std::vector<double>& x,
                            double dx,
                            double horizon,
                            double L);

#endif
