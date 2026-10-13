#include "bond_health.hh"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> bondHealth(const std::string& law,
                               const std::vector<double>& max_stretch,
                               const std::vector<double>& s_elastic,
                               const std::vector<double>& s_break) {
  if (max_stretch.size() != s_elastic.size() ||
      max_stretch.size() != s_break.size()) {
    throw std::invalid_argument(
        "Bond-health input vectors must have the same size.");
  }

  std::vector<double> health(max_stretch.size(), 1.0);

  std::string law_lower = law;

  std::transform(law_lower.begin(),
                 law_lower.end(),
                 law_lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (law_lower == "bilinear") {
    for (std::size_t i = 0; i < max_stretch.size(); ++i) {
      if (s_elastic[i] < max_stretch[i]) {
        health[i] = s_elastic[i] * (s_break[i] - max_stretch[i]) /
                    (max_stretch[i] * (s_break[i] - s_elastic[i]));
      }
    }

  } else if (law_lower == "pmb") {
    for (std::size_t i = 0; i < max_stretch.size(); ++i) {
      if (s_break[i] < max_stretch[i]) {
        health[i] = 0.0;
      }
    }

  } else {
    throw std::invalid_argument("Invalid bond law: " + law);
  }

  for (double& h : health) {
    h = std::clamp(h, 0.0, 1.0);
  }

  return health;
}