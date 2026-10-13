#ifndef COUNT_FRAGMENTS_HH
#define COUNT_FRAGMENTS_HH

#include <cstddef>
#include <vector>

struct Bar;

struct FragmentData {
  std::size_t num_fragments;
  std::vector<double> cut_x;
  std::vector<double> sizes;
};

/**
 * @brief Finds where the bar has split and computes the size of each fragment.
 *
 * The bar is split in the gap between two neighbouring points when no unbroken
 * bond spans that gap. A bond between points i and j spans all gaps from i to
 * j - 1.
 *
 * For every gap, the function counts the number of active breakable bonds that
 * cross it. Gaps crossed by no active bond are interpreted as cuts.
 *
 * Cuts separated by less than one horizon are merged. A fragment thinner than
 * the horizon is below the smallest length scale represented by the model and
 * is therefore merged with a neighbouring fragment.
 *
 * @param bar Peridynamic bar containing geometry and bond information.
 * @param health Current health of every bond, between 0 and 1.
 *
 * @return FragmentData containing the number of fragments, cut positions and
 * fragment sizes.
 */
FragmentData countFragments(const Bar& bar, const std::vector<double>& health);

#endif