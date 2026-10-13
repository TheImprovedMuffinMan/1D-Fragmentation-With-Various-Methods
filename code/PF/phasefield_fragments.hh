#pragma once

#include <aka_common.hh>
#include <string>
#include <vector>

#include "coupler_solid_phasefield.hh"
#include "group_manager.hh"
#include "mesh_iterators.hh"
#include "phase_field_element_filter.hh"

using namespace akantu;

class PhaseFieldFragmentManager : public GroupManager {
public:
  PhaseFieldFragmentManager(PhaseFieldModel& model,
                            const Real& damage_limit,
                            const ID& id)
    : GroupManager(model.getMesh(), id), model(model),
      damage_limit(damage_limit) {}

  Int countFragments(const std::string& prefix = "fragment") {
    const Int dim = model.getSpatialDimension();

    nb_fragments = createClusters(
        dim, prefix, PhaseFieldElementFilter(model, damage_limit));

    cluster_prefix = prefix;

    return nb_fragments;
  }

  ElementGroup& getFragment(Int fragment_id) {
    return getElementGroup(cluster_prefix + "_" + std::to_string(fragment_id));
  }

  Int getNbFragments() const {
    return nb_fragments;
  }

  Real computeLowDamageEnergy(const Real& damage_limit) {
    PhaseFieldElementFilter element_filter(model, damage_limit);

    const std::string name = "low_damage_" + std::to_string(damage_limit);

    auto& group = createElementGroup(name);

    for_each_element(mesh, [&](auto&& element) {
      if (element_filter(element)) {
        group.add(element);
      }
    });

    return model.getEnergy(name);
  }

private:
  PhaseFieldModel& model;
  Real damage_limit;

  Int nb_fragments{0};
  std::string cluster_prefix{"fragment"};
};