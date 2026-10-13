#pragma once

#include <aka_common.hh>
#include <solid_mechanics_model.hh>

#include "coupler_solid_phasefield.hh"
#include "model_functions.hh"

void applyInitialVelocity(akantu::SolidMechanicsModel& solid,
                          akantu::Real strain_rate);

void applyInitialDisplacement(akantu::SolidMechanicsModel& solid,
                              akantu::Real displacement,
                              akantu::Real bar_length);

void initialiseModel(akantu::CouplerSolidPhaseField& model,
                     const PhaseFieldParameters& params);