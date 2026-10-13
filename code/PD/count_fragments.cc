#include "count_fragments.hh"
#include "bar.hh"

#include <vector>

FragmentData countFragments(const Bar& bar, const std::vector<double>& health) {
  const std::size_t N = bar.num_points;

  // Working bonds between two inside points, numbered from 0 to N-1 along bar.
  std::vector<long> change(N, 0);
  for (std::size_t b = 0; b < health.size(); ++b) {
    if (health[b] > 0.0 && bar.bond.can_break[b]) {
      const std::size_t i = bar.bond.i[b] - bar.num_boundary;
      const std::size_t j = bar.bond.j[b] - bar.num_boundary;
      ++change[i];
      --change[j];
    }
  }

  // each bond adds +1 at tis left and -1 at its right end; a running sum gives
  // the number of bonds crossing each gap
  std::vector<long> crossing(N, 0);
  if (N > 0) {
    crossing[0] = change[0];
  }
  for (std::size_t i = 1; i < N; ++i) {
    crossing[i] = crossing[i - 1] + change[i];
  }

  std::vector<double> cut_raw;
  for (std::size_t gap = 0; gap + 1 < N; ++gap) {
    if (crossing[gap] == 0) {
      const std::size_t left = bar.num_boundary + gap;
      const std::size_t right = left + 1;
      // cuts are midway in the gap
      cut_raw.push_back(0.5 * (bar.x[left] + bar.x[right]));
    }
  }

  // Merge cuts closer together than one horizon.
  const double min_size = bar.horizon;
  std::vector<double> cut_x;
  double last = 0.0;
  for (const double cut : cut_raw) {
    if (cut - last >= min_size) {
      cut_x.push_back(cut);
      last = cut;
    }
  }

  while (!cut_x.empty() && bar.L - cut_x.back() < min_size) {
    cut_x.pop_back();
  }

  // Store the results in a FragmentData struct.
  FragmentData result;

  result.num_fragments = cut_x.size() + 1;
  result.cut_x = cut_x;
  result.sizes.reserve(result.num_fragments);
  double previous = 0.0;
  for (const double cut : cut_x) {
    result.sizes.push_back(cut - previous);
    previous = cut;
  }
  result.sizes.push_back(bar.L - previous);

  return result;
}