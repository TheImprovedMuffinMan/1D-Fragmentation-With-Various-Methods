#include "make_bar.hh"

#include "parameters.hh"
#include "random_field.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

Bar makeBar(const Parameters& params) {
  Bar bar;

  // points inside the bar
  const std::size_t num_points =
      static_cast<std::size_t>(std::llround(params.length / params.dx));

  // bar length rounded to whole number of cells
  const double L = num_points * params.dx;

  // spacing between points, horizon and how far they reach
  const double dx = L / static_cast<double>(num_points);
  const double horizon = params.m * params.dx;

  // Each point owns a cell dx wide, so a neighbour with centre which is up to
  // half a cell beyond the horizon will still overlap. that sets how many
  // neighbours a point has on each side and how many points are needed in the
  // end
  const double radius = dx / 2.0;
  const std::size_t reach =
      static_cast<std::size_t>(std::floor((horizon + radius) / dx + 1e-12));

  const std::size_t num_boundary = reach;
  const std::size_t num_total = num_points + 2 * num_boundary;

  // coordinates of all points, including ghost points
  std::vector<double> x(num_total);
  for (std::size_t i = 0; i < num_total; ++i) {
    x[i] =
        (static_cast<double>(i) + 0.5 - static_cast<double>(num_boundary)) * dx;
  }

  // indices of points inside the bar, numbered from 0 to num_points - 1
  std::vector<std::size_t> inside;
  inside.reserve(num_points);

  for (std::size_t i = 0; i < num_points; ++i) {
    inside.push_back(num_boundary + i);
  }

  // boundary points are the first and last num_boundary points
  std::vector<std::size_t> boundary;
  boundary.reserve(2 * num_boundary);

  for (std::size_t i = 0; i < num_boundary; ++i) {
    boundary.push_back(i);
  }

  for (std::size_t i = num_boundary + num_points; i < num_total; ++i) {
    boundary.push_back(i);
  }

  const double area = dx * dx;
  const double volume = area * dx;

  // bonds
  std::size_t num_bonds = 0;
  for (std::size_t gap = 1; gap <= reach; ++gap) {
    num_bonds += num_total - gap;
  }

  std::vector<std::size_t> bi;
  std::vector<std::size_t> bj;
  bi.reserve(num_bonds);
  bj.reserve(num_bonds);

  for (std::size_t gap = 1; gap <= reach; ++gap) {
    const std::size_t n = num_total - gap;

    for (std::size_t i = 0; i < n; ++i) {
      bi.push_back(i);
      bj.push_back(i + gap);
    }
  }

  std::vector<double> bond_length(num_bonds);

  for (std::size_t b = 0; b < num_bonds; ++b) {
    bond_length[b] = x[bj[b]] - x[bi[b]];
  }

  // Partial-volume correction.

  std::vector<double> overlap(num_bonds);
  for (std::size_t b = 0; b < num_bonds; ++b) {
    overlap[b] = (horizon + radius - bond_length[b]) / (2.0 * radius);
    if (bond_length[b] <= horizon - radius) {
      overlap[b] = 1.0;
    }
    if (bond_length[b] > horizon + radius) {
      overlap[b] = 0.0;
    }
    overlap[b] = std::clamp(overlap[b], 0.0, 1.0);
  }

  // Bond weighting.
  std::string bond_shape = params.bond_shape;

  std::transform(bond_shape.begin(),
                 bond_shape.end(),
                 bond_shape.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  std::vector<double> weight(num_bonds);

  if (bond_shape == "constant") {
    std::fill(weight.begin(), weight.end(), 1.0);

  } else if (bond_shape == "conical") {
    for (std::size_t b = 0; b < num_bonds; ++b) {
      weight[b] = std::max(0.0, 1.0 - bond_length[b] / horizon);
      if (bond_length[b] >= horizon * (1.0 - 1e-12)) {
        weight[b] = 0.0;
      }
    }
  } else {
    throw std::invalid_argument("bond_shape must be 'constant' or 'conical'.");
  }

  // Remove bonds with zero overlap or zero weight.
  std::vector<std::size_t> bi_keep;
  std::vector<std::size_t> bj_keep;
  std::vector<double> length_keep;
  std::vector<double> overlap_keep;
  std::vector<double> weight_keep;

  bi_keep.reserve(num_bonds);
  bj_keep.reserve(num_bonds);
  length_keep.reserve(num_bonds);
  overlap_keep.reserve(num_bonds);
  weight_keep.reserve(num_bonds);

  // keep only bonds with non-zero overlap and weight
  for (std::size_t b = 0; b < num_bonds; ++b) {
    if (overlap[b] > 0.0 && weight[b] > 0.0) {
      bi_keep.push_back(bi[b]);
      bj_keep.push_back(bj[b]);
      length_keep.push_back(bond_length[b]);
      overlap_keep.push_back(overlap[b]);
      weight_keep.push_back(weight[b]);
    }
  }

  // move back to the original vectors, which are now smaller
  bi = std::move(bi_keep);
  bj = std::move(bj_keep);
  bond_length = std::move(length_keep);
  overlap = std::move(overlap_keep);
  weight = std::move(weight_keep);

  num_bonds = bi.size();

  // stiffness
  const double E = params.E;
  const double kc = E / (area * horizon * horizon);

  double c;

  if (bond_shape == "constant") {
    c = 2.0 * kc;
  } else {
    c = 6.0 * kc;
  }

  std::vector<double> k_bond(num_bonds);

  for (std::size_t b = 0; b < num_bonds; ++b) {
    k_bond[b] = c * weight[b] * volume * overlap[b];
  }

  const std::size_t mid =
      num_boundary + static_cast<std::size_t>(std::llround(num_points / 2.0)) -
      1;

  double E_as_built = 0.0;

  for (std::size_t b = 0; b < num_bonds; ++b) {
    if (bi[b] == mid || bj[b] == mid) {
      E_as_built += k_bond[b] * bond_length[b];
    }
  }

  E_as_built *= 0.5;

  const double stiffness_scale = E / E_as_built;

  for (double& k : k_bond) {
    k *= stiffness_scale;
  }

  std::vector<double> energy_coef(num_bonds);

  for (std::size_t b = 0; b < num_bonds; ++b) {
    energy_coef[b] = 0.5 * k_bond[b] * bond_length[b] * volume;
  }

  // bond strength
  const auto fields = randomField(params, x, dx, horizon, L);
  const std::vector<double>& strength_pt = fields.strength;
  const std::vector<double>& toughness_pt = fields.toughness;
  std::vector<double> strength(num_bonds);
  std::vector<double> toughness(num_bonds);

  for (std::size_t b = 0; b < num_bonds; ++b) {
    strength[b] = 0.5 * (strength_pt[bi[b]] + strength_pt[bj[b]]);
    toughness[b] = 0.5 * (toughness_pt[bi[b]] + toughness_pt[bj[b]]);
  }

  std::vector<double> s_elastic(num_bonds);
  for (std::size_t b = 0; b < num_bonds; ++b) {
    s_elastic[b] = strength[b] / E;
  }

  // Compute the energy required to produce one cut.
  double cut_sum = 0.0;
  long long previous_gap = -1;
  for (std::size_t b = 0; b < num_bonds; ++b) {
    const long long gap = std::llround(bond_length[b] / dx);

    if (gap != previous_gap) {
      cut_sum += static_cast<double>(gap) * energy_coef[b];
      previous_gap = gap;
    }
  }

  std::vector<double> target(num_bonds);
  for (std::size_t b = 0; b < num_bonds; ++b) {
    target[b] = toughness[b] * area / cut_sum;
  }

  std::string bond_law = params.bond_law;
  std::transform(bond_law.begin(),
                 bond_law.end(),
                 bond_law.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::vector<double> s_break(num_bonds);
  if (bond_law == "bilinear") {
    for (std::size_t b = 0; b < num_bonds; ++b) {
      s_break[b] = target[b] / s_elastic[b];
    }

  } else if (bond_law == "pmb") {
    for (std::size_t b = 0; b < num_bonds; ++b) {
      s_break[b] = std::sqrt(target[b]);
    }

  } else {
    throw std::invalid_argument("bond_law must be 'bilinear' or 'pmb'.");
  }

  double mean_s_break = 0.0;
  for (const double value : s_break) {
    mean_s_break += value;
  }
  mean_s_break /= static_cast<double>(num_bonds);

  // safety: a bond must break before it opens
  if (bond_law == "bilinear") {
    for (std::size_t b = 0; b < num_bonds; ++b) {
      if (s_break[b] <= s_elastic[b]) {
        s_break[b] = s_elastic[b] * (1.0 + 1e-6);
      }
    }
  }

  // only bonds between two inside points can break. Bonds to boundary points
  // are held together so the bar is always gripped at its ends
  std::vector<bool> both_inside(num_bonds, false);
  const std::size_t first_inside = num_boundary;
  const std::size_t last_inside = num_boundary + num_points;
  for (std::size_t b = 0; b < num_bonds; ++b) {
    both_inside[b] = bi[b] >= first_inside && bi[b] < last_inside &&
                     bj[b] >= first_inside && bj[b] < last_inside;
  }

  // The largest time step that keeps the calculation stable is set by the
  // stiffest point. We use half of it to be stable (can be changed to even
  // lower with contact).
  std::vector<double> total_spring(num_total, 0.0);
  for (std::size_t b = 0; b < num_bonds; ++b) {
    const double spring = k_bond[b] / bond_length[b];
    total_spring[bi[b]] += spring;
    total_spring[bj[b]] += spring;
  }

  double max_total_spring = 0.0;
  for (const std::size_t i : inside) {
    max_total_spring = std::max(max_total_spring, total_spring[i]);
  }
  const double dt_max = std::sqrt(2.0 * params.density / max_total_spring);
  const double dt = 0.5 * dt_max;

  // Save all the data in the Bar struct.
  bar.L = L;
  bar.dx = dx;
  bar.horizon = horizon;
  bar.num_points = num_points;
  bar.num_boundary = num_boundary;
  bar.num_total = num_total;
  bar.x = std::move(x);
  bar.inside = std::move(inside);
  bar.boundary = std::move(boundary);
  bar.area = area;
  bar.volume = volume;
  bar.wave_speed = std::sqrt(E / params.density);
  bar.dt = dt;
  bar.bond.i = std::move(bi);
  bar.bond.j = std::move(bj);
  bar.bond.length = std::move(bond_length);
  bar.bond.overlap = std::move(overlap);
  bar.bond.k = std::move(k_bond);
  bar.bond.energy_coef = std::move(energy_coef);
  bar.bond.s_elastic = std::move(s_elastic);
  bar.bond.s_break = std::move(s_break);
  bar.bond.can_break = std::move(both_inside);
  bar.bond.count = num_bonds;
  bar.mean_s_break = mean_s_break;
  bar.strength_pt = strength_pt;
  bar.toughness_pt = toughness_pt;

  return bar;
}