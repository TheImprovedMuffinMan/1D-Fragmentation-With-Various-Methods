// Disclaimer: This file was mainly written by ChatGPT
#ifndef DIAGNOSTICS_WRITER_HH
#define DIAGNOSTICS_WRITER_HH

#include <cstddef>
#include <string>
#include <vector>

/**
 * @brief Global physical diagnostics shared by CZM, phase-field and
 * peridynamic simulations.
 *
 * This structure contains quantities that have the same physical meaning
 * across the three fragmentation methods. Each model is responsible for
 * computing its own values before passing them to the common writer.
 */
struct GlobalDiagnostics {
  /// Current simulation step
  std::size_t step;

  /// Current physical simulation time [s]
  double time;

  /// Current time step [s]
  double dt;

  /// Imposed or measured end displacement [m]
  double displacement;

  /// Global strain [-]
  double strain;

  /// Reaction force at the left boundary [N]
  double reaction_left;

  /// Reaction force at the right boundary [N]
  double reaction_right;

  /// Total kinetic energy [J]
  double kinetic_energy;

  /// Recoverable elastic or strain energy [J]
  double elastic_energy;

  /// Reversible cohesive energy [J]
  double reversible_cohesive_energy;

  /// Energy dissipated through fracture or damage [J]
  double fracture_energy;

  /// Contact energy [J]
  double contact_energy;

  /// Mechanical work supplied through the boundaries [J]
  double external_work;

  /// Total mechanical energy used for the energy-balance check [J]
  double total_energy;

  /// Absolute energy-balance residual [J]
  double energy_error;

  /// Number of fragments at the current time
  std::size_t num_fragments;

  /// Mean fragment size [m]
  double mean_fragment_size;

  /// Minimum fragment size [m]
  double min_fragment_size;

  /// Maximum fragment size [m]
  double max_fragment_size;

  /// Spatially averaged stress [Pa]
  double average_stress;

  /// Maximum stress in the body [Pa]
  double max_stress;

  /// Spatially averaged strain [-]
  double average_strain;

  /// Maximum strain in the body [-]
  double max_strain;

  /// Total linear momentum in the loading direction [kg m/s]
  double momentum;
};

/**
 * @brief Physical properties of one fragment at one recorded time.
 */
struct FragmentRecord {
  /// Integer identifier assigned to the fragment
  std::size_t fragment_id;

  /// Fragment length [m]
  double size;

  /// Position of the fragment centre [m]
  double centre;

  /// Mean velocity of the fragment in the loading direction [m/s]
  double mean_velocity;

  /// Mass of the fragment [kg]
  double mass;
};

/**
 * @brief Physical nodal field values recorded at one spatial position.
 *
 * These quantities naturally live at mesh nodes in the finite-element models
 * and at material points in the peridynamic model.
 */
struct NodalFieldRecord {
  /// Reference position [m].
  double x;

  /// Displacement in the loading direction [m].
  double displacement;

  /// Velocity in the loading direction [m/s].
  double velocity;

  /// Acceleration in the loading direction [m/s^2].
  double acceleration;

  /// Damage variable [-].
  double damage;
};

/**
 * @brief Material quantities recorded at an integration or material point.
 *
 * In Akantu, stress and the displacement gradient are material fields stored
 * at quadrature points rather than nodal fields. For the pseudo-1D problem,
 * strain and stress correspond to the axial components epsilon_xx and sigma_xx.
 */
struct MaterialFieldRecord {
  /// Position of the integration or material point [m].
  double x;

  /// Axial strain [-].
  double strain;

  /// Axial stress [Pa].
  double stress;
};

/**
 * @brief Summary information recorded once at the end of a simulation.
 */
struct RunSummary {
  /// Final physical simulation time reached [s].
  double physical_time;

  /// Wall-clock runtime of the simulation [s].
  double wall_time;

  /// Simulation runtime without I/O operations [s].
  double simulation_wall_time;

  /// Total number of completed time steps.
  std::size_t steps;

  /// Final number of fragments.
  std::size_t num_fragments;

  /// Maximum or final normalised energy-balance error [-].
  double energy_error;

  /// Human-readable description of why the simulation stopped.
  std::string stopping_reason;
};

/**
 * @brief Creates the global diagnostics file and writes its column header.
 *
 * Any existing file at the same path is overwritten.
 *
 * @param filename Path to the diagnostics output file.
 */
void initialiseDiagnostics(const std::string& filename);

/**
 * @brief Creates the fragment-history file and writes its column header.
 *
 * Any existing file at the same path is overwritten.
 *
 * @param filename Path to the fragment output file.
 */
void initialiseFragments(const std::string& filename);

/**
 * @brief Appends one line of global diagnostics to the output file.
 *
 * @param filename Path to the diagnostics output file.
 * @param data Global physical quantities for the current recorded time.
 */
void writeDiagnostics(const std::string& filename,
                      const GlobalDiagnostics& data);

/**
 * @brief Writes the fragment properties recorded at one simulation time.
 *
 * Each fragment is written on a separate line together with the supplied
 * physical time.
 *
 * @param filename Path to the fragment output file.
 * @param time Physical simulation time corresponding to the records [s].
 * @param fragments Fragment properties at the current time.
 */
void writeFragments(const std::string& filename,
                    double time,
                    const std::vector<FragmentRecord>& fragments);

/**
 * @brief Writes one complete nodal field snapshot.
 *
 * Any existing file at the same path is overwritten. A different filename
 * should therefore be supplied for every recorded snapshot.
 *
 * @param filename Path to the nodal field snapshot file.
 * @param fields Nodal field values at all recorded spatial positions.
 */
void writeNodalFieldSnapshot(const std::string& filename,
                             const std::vector<NodalFieldRecord>& fields);

/**
 * @brief Writes one complete material-field snapshot.
 *
 * For the Akantu models, the records correspond to quadrature points.
 * For peridynamics, equivalent material quantities may be evaluated at the
 * material points.
 *
 * @param filename Path to the material field snapshot file.
 * @param fields Stress and strain values at all recorded material positions.
 */
void writeMaterialFieldSnapshot(const std::string& filename,
                                const std::vector<MaterialFieldRecord>& fields);
/**
 * @brief Writes final summary information for one simulation run.
 *
 * Any existing file at the same path is overwritten.
 *
 * @param filename Path to the summary output file.
 * @param summary Final runtime and simulation information.
 */
void writeRunSummary(const std::string& filename, const RunSummary& summary);

#endif