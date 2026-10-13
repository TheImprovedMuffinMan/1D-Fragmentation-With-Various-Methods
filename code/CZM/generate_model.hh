#pragma once

#include <aka_common.hh>
#include <mesh.hh>
#include <solid_mechanics_model_cohesive.hh>

/// @brief Applies the initial affine velocity field v(x,0) = strain_rate * x.
void applyInitialVelocity(akantu::SolidMechanicsModelCohesive& model,
                          akantu::Real strain_rate);

/// @brief Restricts cohesive insertion to vertical facets.
void restrictCohesiveInsertion(akantu::SolidMechanicsModelCohesive& model,
                               akantu::Mesh& mesh);

/// @brief Initialises the extrinsic cohesive-zone model and solvers.
void initialiseModel(akantu::SolidMechanicsModelCohesive& model,
                     akantu::Mesh& mesh,
                     akantu::Real strain_rate);

/// @brief Applies the prescribed displacement at the two ends of the bar.
void applyVelocityBoundaryConditions(akantu::SolidMechanicsModelCohesive& model,
                                     akantu::Real applied_velocity,
                                     akantu::Int current_time_step,
                                     akantu::Real dt);