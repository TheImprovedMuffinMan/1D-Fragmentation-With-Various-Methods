// Disclaimer: This file was mainly written by ChatGPT
#include "diagnostic_writer.hh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

// ========================================================================== //
void initialiseDiagnostics(const std::string& filename) {
  const std::filesystem::path path(filename);

  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream file(filename);

  if (!file) {
    throw std::runtime_error("Could not create diagnostics file: " + filename);
  }

  file << "# step "
       << "time "
       << "dt "
       << "displacement "
       << "strain "
       << "reaction_left "
       << "reaction_right "
       << "kinetic_energy "
       << "elastic_energy "
       << "reversible_cohesive_energy "
       << "fracture_energy "
       << "contact_energy "
       << "external_work "
       << "total_energy "
       << "energy_error "
       << "num_fragments "
       << "mean_fragment_size "
       << "min_fragment_size "
       << "max_fragment_size "
       << "average_stress "
       << "max_stress "
       << "average_strain "
       << "max_strain "
       << "momentum\n";
}

// ========================================================================== //
void writeDiagnostics(const std::string& filename,
                      const GlobalDiagnostics& data) {
  std::ofstream file(filename, std::ios::app);

  if (!file) {
    throw std::runtime_error("Could not open diagnostics file: " + filename);
  }

  file << std::setprecision(16);

  file << data.step << " " << data.time << " " << data.dt << " "
       << data.displacement << " " << data.strain << " " << data.reaction_left
       << " " << data.reaction_right << " " << data.kinetic_energy << " "
       << data.elastic_energy << " " << data.reversible_cohesive_energy << " "
       << data.fracture_energy << " " << data.contact_energy << " "
       << data.external_work << " " << data.total_energy << " "
       << data.energy_error << " " << data.num_fragments << " "
       << data.mean_fragment_size << " " << data.min_fragment_size << " "
       << data.max_fragment_size << " " << data.average_stress << " "
       << data.max_stress << " " << data.average_strain << " "
       << data.max_strain << " " << data.momentum << "\n";
}

// ========================================================================== //
void initialiseFragments(const std::string& filename) {
  const std::filesystem::path path(filename);

  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream file(filename);

  if (!file) {
    throw std::runtime_error("Could not create fragment file: " + filename);
  }

  file << "# time "
       << "fragment_id "
       << "size "
       << "mass "
       << "centre "
       << "mean_velocity\n";
}

// ========================================================================== //
void writeFragments(const std::string& filename,
                    double time,
                    const std::vector<FragmentRecord>& fragments) {
  std::ofstream file(filename, std::ios::app);

  if (!file) {
    throw std::runtime_error("Could not open fragment file: " + filename);
  }

  file << std::setprecision(16);

  for (const auto& fragment : fragments) {
    file << time << " " << fragment.fragment_id << " " << fragment.size << " "
         << fragment.mass << " " << fragment.centre << " "
         << fragment.mean_velocity << "\n";
  }
}

// ========================================================================== //
void writeNodalFieldSnapshot(const std::string& filename,
                             const std::vector<NodalFieldRecord>& fields) {
  const std::filesystem::path path(filename);

  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream file(filename);

  if (!file) {
    throw std::runtime_error("Could not create nodal field file: " + filename);
  }

  file << std::setprecision(16);

  file << "# x "
       << "displacement "
       << "velocity "
       << "acceleration "
       << "damage\n";

  for (const auto& field : fields) {
    file << field.x << " " << field.displacement << " " << field.velocity << " "
         << field.acceleration << " " << field.damage << "\n";
  }
}

// ========================================================================== //
void writeMaterialFieldSnapshot(
    const std::string& filename,
    const std::vector<MaterialFieldRecord>& fields) {
  const std::filesystem::path path(filename);

  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream file(filename);

  if (!file) {
    throw std::runtime_error("Could not create material field file: " +
                             filename);
  }

  file << std::setprecision(16);

  file << "# x "
       << "strain "
       << "stress\n";

  for (const auto& field : fields) {
    file << field.x << " " << field.strain << " " << field.stress << "\n";
  }
}

// ========================================================================== //
void writeRunSummary(const std::string& filename, const RunSummary& summary) {
  const std::filesystem::path path(filename);

  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream file(filename);

  if (!file) {
    throw std::runtime_error("Could not create run summary file: " + filename);
  }

  file << std::setprecision(16);

  file << "physical_time " << summary.physical_time << "\n";
  file << "wall_time " << summary.wall_time << "\n";
  file << "simulation_wall_time " << summary.simulation_wall_time << "\n";
  file << "steps " << summary.steps << "\n";
  file << "num_fragments " << summary.num_fragments << "\n";
  file << "energy_error " << summary.energy_error << "\n";
  file << "stopping_reason " << summary.stopping_reason << "\n";
}