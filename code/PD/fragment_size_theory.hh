#ifndef FRAGMENT_SIZE_THEORY_HH
#define FRAGMENT_SIZE_THEORY_HH

struct Parameters;

struct FragmentSizeTheory {
  double s0;
  double rate0;
  double rate_ratio;
  double zmr;
  double grady;
};

/**
 * @brief Computes theoretical predictions for the average fragment size.
 *
 * The function evaluates two fragmentation theories for the material properties
 * and strain rate contained in @p params.
 *
 * The Zhou-Molinari-Ramesh prediction is based on their fitted relation for
 * one-dimensional dynamic fragmentation with cohesive cracks.
 *
 * The Grady prediction follows from an energy balance and generally predicts
 * larger fragments because it does not account for the finite time required
 * for cracks to open.
 *
 * The function also computes the natural length and strain-rate scales:
 *
 * @f[
 * s_0 = \frac{E G_c}{\sigma_c^2}
 * @f]
 *
 * and
 *
 * @f[
 * \dot{\varepsilon}_0 =
 * \frac{c \sigma_c^3}{E^2 G_c},
 * \qquad
 * c = \sqrt{\frac{E}{\rho}}.
 * @f]
 *
 * @param params Material properties and applied strain rate.
 *
 * @return FragmentSizeTheory containing the natural scales and theoretical
 * fragment-size predictions.
 */
FragmentSizeTheory fragmentSizeTheory(const Parameters& params);

#endif