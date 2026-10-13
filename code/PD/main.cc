#include "fragment_size_theory.hh"
#include "make_bar.hh"
#include "parameters.hh"
#include "simulate_bar.hh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

#include "fragment_size_theory.hh"
#include "make_bar.hh"
#include "parameters.hh"
#include "simulate_bar.hh"

#include <iostream>

int main() {
  const Parameters params = readParameters("../code/PD/parameters.dat");

  Bar bar = makeBar(params);

  std::cout << "Bar: " << bar.num_points << " points, " << bar.bond.count
            << " bonds, horizon " << bar.horizon * 1e6 << " um, time step "
            << bar.dt << " s\n";

  SimulationResults results = simulateBar(params, bar);

  const FragmentSizeTheory theory = fragmentSizeTheory(params);

  const std::vector<double>& final_sizes = results.sizes.back();
  const double average_size =
      std::accumulate(final_sizes.begin(), final_sizes.end(), 0.0) /
      static_cast<double>(final_sizes.size());

  std::cout
      << "\n------------------------------------------------------------\n";

  std::cout << "Fragments: " << results.num_fragments.back() << "\n";

  std::cout << "Average fragment size: " << average_size * 1e6 << " um\n";

  std::cout << "ZMR prediction: " << theory.zmr * 1e6 << " um\n";

  std::cout << "Grady prediction: " << theory.grady * 1e6 << " um\n";

  std::cout << "Energy error: " << results.energy_error << "\n";

  // --------------------------------------------------------------- output

  const std::string output_dir = "../data/PD";

  std::filesystem::create_directories(output_dir);

  std::ofstream output(output_dir + "/results.out");

  output << std::setprecision(16);

  output << "# snapshot time kinetic strain_energy dissipated "
            "external_work num_fragments mean_fragment_size\n";

  for (std::size_t k = 0; k < results.t.size(); ++k) {
    double mean_size = 0.0;

    if (!results.sizes[k].empty()) {
      mean_size = std::accumulate(
                      results.sizes[k].begin(), results.sizes[k].end(), 0.0) /
                  static_cast<double>(results.sizes[k].size());
    }

    output << k << " " << results.t[k] << " " << results.kinetic[k] << " "
           << results.strain_energy[k] << " " << results.dissipated[k] << " "
           << results.external_work[k] << " " << results.num_fragments[k] << " "
           << mean_size << "\n";
  }

  output.close();

  return 0;
}