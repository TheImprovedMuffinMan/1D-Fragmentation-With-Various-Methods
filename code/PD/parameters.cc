#include "parameters.hh"

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

std::string trim(const std::string& text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");

  if (first == std::string::npos) {
    return "";
  }

  const std::size_t last = text.find_last_not_of(" \t\r\n");

  return text.substr(first, last - first + 1);
}

} // namespace

Parameters readParameters(const std::string& filename) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open parameter file: " + filename);
  }

  Parameters params{};

  std::string line;

  while (std::getline(file, line)) {
    const std::size_t comment = line.find('#');

    if (comment != std::string::npos) {
      line = line.substr(0, comment);
    }

    line = trim(line);

    if (line.empty()) {
      continue;
    }

    const std::size_t equal = line.find('=');

    if (equal == std::string::npos) {
      throw std::runtime_error("Invalid parameter line: " + line);
    }

    const std::string key = trim(line.substr(0, equal));
    const std::string value = trim(line.substr(equal + 1));

    if (key == "density") {
      params.density = std::stod(value);

    } else if (key == "E") {
      params.E = std::stod(value);

    } else if (key == "strength") {
      params.strength = std::stod(value);

    } else if (key == "toughness") {
      params.toughness = std::stod(value);

    } else if (key == "length") {
      params.length = std::stod(value);

    } else if (key == "area") {
      params.area = std::stod(value);

    } else if (key == "dx") {
      params.dx = std::stod(value);

    } else if (key == "horizon") {
      params.horizon = std::stod(value);

    } else if (key == "bond_law") {
      params.bond_law = value;

    } else if (key == "bond_shape") {
      params.bond_shape = value;

    } else if (key == "scatter") {
      params.scatter = std::stod(value);

    } else if (key == "seed") {
      params.seed = static_cast<unsigned int>(std::stoul(value));

    } else if (key == "strain_rate") {
      params.strain_rate = std::stod(value);

    } else if (key == "time_step_factor") {
      params.time_step_factor = std::stod(value);

    } else if (key == "max_time") {
      params.max_time = std::stod(value);

    } else if (key == "dump_freq") {
      params.dump_freq = static_cast<std::size_t>(std::stoull(value));

    } else if (key == "check_frequency") {
      params.check_frequency = static_cast<std::size_t>(std::stoull(value));

    } else if (key == "stable_steps") {
      params.stable_steps = static_cast<std::size_t>(std::stoull(value));

    } else if (key == "output_dir") {
      params.output_dir = value;

    } else if (key == "diagnostics_file") {
      params.diagnostics_file = value;

    } else if (key == "fragment_file") {
      params.fragment_file = value;

    } else {
      throw std::runtime_error("Unknown parameter: " + key);
    }
  }

  // ------------------------------------------------------- validation

  if (params.density <= 0.0) {
    throw std::runtime_error("density must be greater than zero.");
  }

  if (params.E <= 0.0) {
    throw std::runtime_error("E must be greater than zero.");
  }

  if (params.strength <= 0.0) {
    throw std::runtime_error("strength must be greater than zero.");
  }

  if (params.toughness <= 0.0) {
    throw std::runtime_error("toughness must be greater than zero.");
  }

  if (params.length <= 0.0) {
    throw std::runtime_error("length must be greater than zero.");
  }

  if (params.area <= 0.0) {
    throw std::runtime_error("area must be greater than zero.");
  }

  if (params.dx <= 0.0) {
    throw std::runtime_error("dx must be greater than zero.");
  }

  if (params.horizon <= 0.0) {
    throw std::runtime_error("horizon must be greater than zero.");
  }

  if (params.scatter < 0.0) {
    throw std::runtime_error("scatter must be non-negative.");
  }

  if (params.strain_rate <= 0.0) {
    throw std::runtime_error("strain_rate must be greater than zero.");
  }

  if (params.time_step_factor <= 0.0 || params.time_step_factor > 1.0) {
    throw std::runtime_error(
        "time_step_factor must be in the interval (0, 1].");
  }

  if (params.max_time <= 0.0) {
    throw std::runtime_error("max_time must be greater than zero.");
  }

  if (params.dump_freq == 0) {
    throw std::runtime_error("dump_freq must be greater than zero.");
  }

  if (params.check_frequency == 0) {
    throw std::runtime_error("check_frequency must be greater than zero.");
  }

  if (params.stable_steps == 0) {
    throw std::runtime_error("stable_steps must be greater than zero.");
  }

  if (params.bond_law.empty()) {
    throw std::runtime_error("bond_law must be specified.");
  }

  if (params.bond_shape.empty()) {
    throw std::runtime_error("bond_shape must be specified.");
  }

  if (params.output_dir.empty()) {
    throw std::runtime_error("output_dir must be specified.");
  }

  if (params.diagnostics_file.empty()) {
    throw std::runtime_error("diagnostics_file must be specified.");
  }

  if (params.fragment_file.empty()) {
    throw std::runtime_error("fragment_file must be specified.");
  }

  return params;
}