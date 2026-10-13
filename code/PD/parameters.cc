#include "parameters.hh"

#include <fstream>
#include <sstream>
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
    } else if (key == "dx") {
      params.dx = std::stod(value);
    } else if (key == "m") {
      params.m = std::stod(value);
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
    } else {
      throw std::runtime_error("Unknown parameter: " + key);
    }
  }

  return params;
}