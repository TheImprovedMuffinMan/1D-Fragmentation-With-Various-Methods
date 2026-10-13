#include "random_field.hh"
#include "parameters.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

std::vector<double> smoothNoise(std::mt19937& generator,
                                const std::vector<double>& x,
                                double l_c,
                                double dx,
                                double horizon,
                                double L) {
  const std::size_t N = x.size();
  std::normal_distribution<double> normal(0.0, 1.0);

  // Smoothing below the grid scale is impossible.
  if (l_c <= 0.0 || l_c < dx) {
    std::vector<double> z(N);

    for (double& value : z) {
      value = normal(generator);
    }

    return z;
  }

  // helper grid

  const double spacing = l_c / 8.0;
  const double margin = std::max(8.0 * l_c, 4.0 * horizon);

  const double x_start = -margin;
  const double x_end = L + margin;

  const std::size_t M = static_cast<std::size_t>(
                            std::floor((x_end - x_start) / spacing + 1e-12)) +
                        1;

  std::vector<double> xr(M);
  std::vector<double> zr(M);

  for (std::size_t i = 0; i < M; ++i) {
    xr[i] = x_start + static_cast<double>(i) * spacing;
    zr[i] = normal(generator);
  }

  // ------------------------------------------------------------ smoothing

  const double width = l_c / spacing;

  const long h = static_cast<long>(std::ceil(3.0 * width));

  const std::size_t bell_size = static_cast<std::size_t>(2 * h + 1);

  std::vector<double> bell(bell_size);

  double bell_sum = 0.0;

  for (long q = -h; q <= h; ++q) {
    const std::size_t index = static_cast<std::size_t>(q + h);
    const double position = static_cast<double>(q) / width;
    bell[index] = std::exp(-0.5 * position * position);
    bell_sum += bell[index];
  }

  for (double& value : bell) {
    value /= bell_sum;
  }

  std::vector<double> zs(M, 0.0);

  for (long q = -h; q <= h; ++q) {
    const double weight = bell[static_cast<std::size_t>(q + h)];

    for (std::size_t k = 0; k < M; ++k) {
      long neighbour = static_cast<long>(k) + q;
      neighbour = std::max(0L, std::min(neighbour, static_cast<long>(M) - 1));
      zs[k] += weight * zr[static_cast<std::size_t>(neighbour)];
    }
  }

  double bell_squared_sum = 0.0;
  for (const double value : bell) {
    bell_squared_sum += value * value;
  }
  const double normalisation = std::sqrt(bell_squared_sum);
  for (double& value : zs) {
    value /= normalisation;
  }

  // interpolation

  std::vector<double> z(N);

  for (std::size_t p = 0; p < N; ++p) {
    const double xp = x[p];

    if (xp <= xr.front()) {
      const double slope = (zs[1] - zs[0]) / (xr[1] - xr[0]);
      z[p] = zs[0] + slope * (xp - xr[0]);

    } else if (xp >= xr.back()) {
      const std::size_t last = M - 1;
      const double slope =
          (zs[last] - zs[last - 1]) / (xr[last] - xr[last - 1]);

      z[p] = zs[last] + slope * (xp - xr[last]);

    } else {
      const auto upper = std::upper_bound(xr.begin(), xr.end(), xp);
      const std::size_t right = static_cast<std::size_t>(upper - xr.begin());
      const std::size_t left = right - 1;
      const double fraction = (xp - xr[left]) / (xr[right] - xr[left]);
      z[p] = zs[left] + fraction * (zs[right] - zs[left]);
    }
  }

  return z;
}

} // namespace

RandomFieldData randomField(const Parameters& params,
                            const std::vector<double>& x,
                            double dx,
                            double horizon,
                            double L) {
  if (params.scatter < 0.0) {
    throw std::runtime_error("scatter must be non-negative.");
  }
  const double E = params.E;
  const double mu_strength = params.strength;
  const double mu_toughness = params.toughness;

  RandomFieldData fields;
  const std::size_t N = x.size();
  fields.strength.resize(N);
  fields.toughness.resize(N);

  if (params.scatter == 0.0) {
    std::fill(fields.strength.begin(), fields.strength.end(), mu_strength);
    std::fill(fields.toughness.begin(), fields.toughness.end(), mu_toughness);
    return fields;
  }

  std::mt19937 generator(params.seed);
  const double l_c = E * mu_toughness / (mu_strength * mu_strength);

  const std::vector<double> z_strength =
      smoothNoise(generator, x, l_c, dx, horizon, L);

  const std::vector<double> z_toughness =
      smoothNoise(generator, x, l_c, dx, horizon, L);

  for (std::size_t i = 0; i < N; ++i) {
    fields.strength[i] =
        mu_strength + params.scatter * mu_strength * z_strength[i];
    fields.toughness[i] =
        mu_toughness + params.scatter * mu_toughness * z_toughness[i];
    fields.strength[i] = std::max(fields.strength[i], 0.05 * mu_strength);
    fields.toughness[i] = std::max(fields.toughness[i], 0.05 * mu_toughness);
  }

  return fields;
}