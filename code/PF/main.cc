#include <aka_common.hh>
#include <mesh.hh>

#include "coupler_solid_phasefield.hh"
#include "generate_mesh.hh"
#include "generate_model.hh"
#include "model_functions.hh"

using namespace akantu;

// -------------------------------------------------------------------------- //
int main(int argc, char* argv[]) {
  initialize("../code/PF/material.dat", argc, argv);

  const PhaseFieldParameters params = readSimulationParameters();

  const auto& comm = Communicator::getStaticCommunicator();
  const Int prank = comm.whoAmI();

  Mesh mesh(2);

  const Real x0 = -0.5 * params.bar_length;
  const Real xf = 0.5 * params.bar_length;

  if (prank == 0) {
    generateMesh(mesh, params.bar_length, x0, xf, params.number_elements);
  }

  mesh.distribute();

  CouplerSolidPhaseField model(mesh);

  initialiseModel(model, params);
  runSimulation(model, params);

  finalize();

  return EXIT_SUCCESS;
}