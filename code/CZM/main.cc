#include <aka_common.hh>
#include <solid_mechanics_model_cohesive.hh>

#include <chrono>
#include <cstdlib>
#include <iostream>

#include "generate_mesh.hh"
#include "generate_model.hh"
#include "model_functions.hh"

using namespace akantu;

// -------------------------------------------------------------------------- //
int main(int argc, char* argv[]) {
  initialize("../code/CZM/material.dat", argc, argv);

  const SimulationParameters params = readSimulationParameters();

  const auto start_time = std::chrono::high_resolution_clock::now();

  const Real x0 = -0.5 * params.bar_length;
  const Real xf = 0.5 * params.bar_length;

  Mesh mesh(2);
  generateMesh(mesh, params.bar_length, x0, xf, params.number_elements);

  SolidMechanicsModelCohesive model(mesh);

  initialiseModel(model, mesh, params.strain_rate);

  runSimulation(model, mesh, params);

  const auto end_time = std::chrono::high_resolution_clock::now();

  const std::chrono::duration<double> total_time = end_time - start_time;

  std::cout << "Simulation time: " << total_time.count() << " seconds\n";

  finalize();

  return EXIT_SUCCESS;
}