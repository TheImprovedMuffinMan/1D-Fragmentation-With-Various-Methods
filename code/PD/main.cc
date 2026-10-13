#include "fragment_size_theory.hh"
#include "make_bar.hh"
#include "parameters.hh"
#include "simulate_bar.hh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  const Parameters params = readParameters("../code/PD/parameters.dat");

  Bar bar = makeBar(params);

  std::cout << "Bar: " << bar.num_points << " points, " << bar.bond.count
            << " bonds\n";

  std::cout << "dx: " << bar.dx * 1e6 << " um\n";

  std::cout << "Horizon: " << bar.horizon * 1e6 << " um\n";

  std::cout << "m = horizon/dx: " << bar.m << "\n";

  std::cout << "Time step: " << bar.dt << " s\n";

  SimulationResults results = simulateBar(params, bar);

  const FragmentSizeTheory theory = fragmentSizeTheory(params);

  const std::vector<double>& final_sizes = results.sizes.back();

  double average_size = 0.0;

  if (!final_sizes.empty()) {
    average_size =
        std::accumulate(final_sizes.begin(), final_sizes.end(), 0.0) /
        static_cast<double>(final_sizes.size());
  }

  std::cout
      << "\n------------------------------------------------------------\n";

  std::cout << "Fragments: " << results.num_fragments.back() << "\n";

  std::cout << "Average fragment size: " << average_size * 1e6 << " um\n";

  std::cout << "ZMR prediction: " << theory.zmr * 1e6 << " um\n";

  std::cout << "Grady prediction: " << theory.grady * 1e6 << " um\n";

  std::cout << "Energy error: " << results.energy_error << "\n";

  std::cout << "Steps: " << results.steps << "\n";

  std::cout << "Wall time: " << results.run_time << " s\n";

  std::cout << "Stopping reason: " << results.stopping_reason << "\n";

  // --------------------------------------------------------------- output

  std::filesystem::create_directories(params.output_dir);

  // ------------------------------------------------------ diagnostics output

  const std::string diagnostics_path =
      params.output_dir + "/" + params.diagnostics_file;

  std::ofstream output(diagnostics_path);

  if (!output.is_open()) {
    throw std::runtime_error("Could not open diagnostics file: " +
                             diagnostics_path);
  }

  output << std::setprecision(16);

  output << "# snapshot time kinetic strain_energy dissipated "
         << "external_work num_fragments mean_fragment_size\n";

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

  // ------------------------------------------------------ fragment output

  const std::string fragment_path =
      params.output_dir + "/" + params.fragment_file;

  std::ofstream fragment_output(fragment_path);

  if (!fragment_output.is_open()) {
    throw std::runtime_error("Could not open fragment file: " + fragment_path);
  }

  fragment_output << std::setprecision(16);

  fragment_output << "# snapshot time fragment_id "
                  << "x_left x_right centre size\n";

  for (std::size_t k = 0; k < results.t.size(); ++k) {
    const std::vector<double>& sizes = results.sizes[k];

    const std::vector<double>& cuts = results.cuts[k];

    for (std::size_t f = 0; f < sizes.size(); ++f) {
      const double x_left = (f == 0) ? 0.0 : cuts[f - 1];

      const double x_right = (f < cuts.size()) ? cuts[f] : bar.L;

      const double centre = 0.5 * (x_left + x_right);

      fragment_output << k << " " << results.t[k] << " " << f << " " << x_left
                      << " " << x_right << " " << centre << " " << sizes[f]
                      << "\n";
    }
  }

  return 0;
}