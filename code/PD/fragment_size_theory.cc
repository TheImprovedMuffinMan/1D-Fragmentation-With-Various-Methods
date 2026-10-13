#include "fragment_size_theory.hh"

#include "parameters.hh"

#include <cmath>

FragmentSizeTheory fragmentSizeTheory(const Parameters& params) {
  const double E = params.E;
  const double rho = params.density;
  const double sigma = params.strength;
  const double G = params.toughness;
  const double rate = params.strain_rate;

  const double c = std::sqrt(E / rho);

  FragmentSizeTheory theory;

  theory.s0 = E * G / (sigma * sigma);
  theory.rate0 = c * sigma * sigma * sigma / (E * E * G);
  theory.rate_ratio = rate / theory.rate0;

  const double r = theory.rate_ratio;

  theory.zmr = theory.s0 * 4.5 /
               (1.0 + 0.77 * std::pow(r, 0.25) + 5.4 * std::pow(r, 0.75));

  theory.grady = std::cbrt(24.0 * G / (rho * rate * rate));

  return theory;
}