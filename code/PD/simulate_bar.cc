#include "simulate_bar.hh"

#include "bond_health.hh"
#include "count_fragments.hh"
#include "parameters.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace {

struct StretchData {
  std::vector<double> stretch;
  std::vector<double> direction;
};

StretchData bondStretch(const std::vector<double>& u,
                        const std::vector<double>& x,
                        const std::vector<std::size_t>& bi,
                        const std::vector<std::size_t>& bj,
                        const std::vector<double>& L0) {
  const std::size_t num_bonds = bi.size();

  StretchData data;

  data.stretch.resize(num_bonds);
  data.direction.resize(num_bonds);

  for (std::size_t b = 0; b < num_bonds; ++b) {
    const double yi = x[bi[b]] + u[bi[b]];
    const double yj = x[bj[b]] + u[bj[b]];
    const double d = yj - yi;

    data.stretch[b] = (std::abs(d) - L0[b]) / L0[b];

    if (d > 0.0) {
      data.direction[b] = 1.0;
    } else if (d < 0.0) {
      data.direction[b] = -1.0;
    } else {
      data.direction[b] = 0.0;
    }
  }

  return data;
}

std::vector<double> bondForces(const std::vector<double>& health,
                               const std::vector<double>& k_bond,
                               const std::vector<double>& stretch,
                               const std::vector<double>& direction,
                               const std::vector<std::size_t>& bi,
                               const std::vector<std::size_t>& bj,
                               std::size_t n_all) {
  std::vector<double> force(n_all, 0.0);

  for (std::size_t b = 0; b < health.size(); ++b) {
    const double fb = health[b] * k_bond[b] * stretch[b] * direction[b];
    force[bi[b]] += fb;
    force[bj[b]] -= fb;
  }

  return force;
}

void moveBoundary(std::vector<double>& u,
                  std::vector<double>& v,
                  double t,
                  double rate,
                  const std::vector<double>& x,
                  const std::vector<std::size_t>& boundary,
                  double L) {
  for (const std::size_t i : boundary) {
    const double relative_position = x[i] - L / 2.0;
    u[i] = rate * t * relative_position;
    v[i] = rate * relative_position;
  }
}

std::vector<std::size_t> makeSaveSteps(std::size_t n_steps,
                                       std::size_t n_snapshots) {
  std::vector<std::size_t> save_at;

  save_at.reserve(n_snapshots + 2);
  save_at.push_back(0);

  if (n_snapshots > 1) {
    for (std::size_t i = 0; i < n_snapshots; ++i) {
      const double fraction =
          static_cast<double>(i) / static_cast<double>(n_snapshots - 1);
      const std::size_t step = static_cast<std::size_t>(
          std::llround(fraction * static_cast<double>(n_steps)));
      save_at.push_back(step);
    }
  }

  save_at.push_back(n_steps);
  std::sort(save_at.begin(), save_at.end());
  save_at.erase(std::unique(save_at.begin(), save_at.end()), save_at.end());

  return save_at;
}

void recordSnapshot(SimulationResults& results,
                    double t,
                    const std::vector<double>& u,
                    const std::vector<double>& v,
                    const std::vector<double>& health,
                    const std::vector<double>& stretch,
                    double dissipated,
                    double work_in,
                    const Bar& bar,
                    double rho) {
  results.t.push_back(t);
  std::vector<double> displacement(bar.num_points);

  for (std::size_t p = 0; p < bar.num_points; ++p) {
    displacement[p] = u[bar.inside[p]];
  }

  results.displacement.push_back(std::move(displacement));

  // Damage at each point.
  std::vector<double> numerator(bar.num_total, 0.0);
  std::vector<double> denominator(bar.num_total, 0.0);

  for (std::size_t b = 0; b < bar.bond.count; ++b) {
    const double weight = bar.volume * bar.bond.overlap[b];
    numerator[bar.bond.i[b]] += health[b] * weight;
    numerator[bar.bond.j[b]] += health[b] * weight;
    denominator[bar.bond.i[b]] += weight;
    denominator[bar.bond.j[b]] += weight;
  }

  std::vector<double> damage(bar.num_points);
  const double eps = std::numeric_limits<double>::epsilon();
  for (std::size_t p = 0; p < bar.num_points; ++p) {
    const std::size_t i = bar.inside[p];
    damage[p] = 1.0 - numerator[i] / std::max(denominator[i], eps);
  }

  results.damage.push_back(std::move(damage));

  const FragmentData fragments = countFragments(bar, health);

  results.num_fragments.push_back(fragments.num_fragments);
  results.cuts.push_back(fragments.cut_x);
  results.sizes.push_back(fragments.sizes);

  double kinetic = 0.0;

  for (const std::size_t i : bar.inside) {
    kinetic += v[i] * v[i];
  }

  kinetic *= 0.5 * rho * bar.volume;

  double strain_energy = 0.0;

  for (std::size_t b = 0; b < bar.bond.count; ++b) {
    strain_energy +=
        bar.bond.energy_coef[b] * health[b] * stretch[b] * stretch[b];
  }

  results.kinetic.push_back(kinetic);
  results.strain_energy.push_back(strain_energy);
  results.dissipated.push_back(dissipated);
  results.external_work.push_back(work_in);
}

} // namespace

SimulationResults simulateBar(const Parameters& params, const Bar& bar) {
  const double rate = params.strain_rate;

  const std::vector<double>& x = bar.x;
  const std::vector<std::size_t>& inside = bar.inside;
  const std::vector<std::size_t>& boundary = bar.boundary;

  const double rho = params.density;
  const double volume = bar.volume;
  const double dt = bar.dt;
  const double L = bar.L;

  const std::vector<std::size_t>& bi = bar.bond.i;
  const std::vector<std::size_t>& bj = bar.bond.j;
  const std::vector<double>& L0 = bar.bond.length;
  const std::vector<double>& k_bond = bar.bond.k;
  const std::vector<double>& energy_coef = bar.bond.energy_coef;

  // Bonds allowed to break.

  std::vector<std::size_t> breakable;

  for (std::size_t b = 0; b < bar.bond.count; ++b) {
    if (bar.bond.can_break[b]) {
      breakable.push_back(b);
    }
  }

  std::vector<double> s_elastic(breakable.size());
  std::vector<double> s_break(breakable.size());
  std::vector<double> energy_coef_breakable(breakable.size());

  for (std::size_t q = 0; q < breakable.size(); ++q) {
    const std::size_t b = breakable[q];

    s_elastic[q] = bar.bond.s_elastic[b];
    s_break[q] = bar.bond.s_break[b];
    energy_coef_breakable[q] = energy_coef[b];
  }

  // ------------------------------------------------------------ run length

  const double t_fail = bar.mean_s_break / rate;
  const double t_total = 6.0 * t_fail;

  const std::size_t n_steps = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(t_total / dt)));

  const std::size_t n_snapshots = 120;
  const std::size_t n_quiet = 12;

  const double quiet_tol = 1e-3;

  const std::vector<std::size_t> save_at = makeSaveSteps(n_steps, n_snapshots);

  // --------------------------------------------------------- initial state

  std::vector<double> u(bar.num_total, 0.0);
  std::vector<double> v(bar.num_total);

  for (std::size_t i = 0; i < bar.num_total; ++i) {
    v[i] = rate * (x[i] - L / 2.0);
  }

  std::vector<double> health(bar.bond.count, 1.0);
  std::vector<double> health_breakable(breakable.size(), 1.0);
  std::vector<double> max_stretch(breakable.size(), 0.0);

  moveBoundary(u, v, 0.0, rate, x, boundary, L);

  StretchData bond_data = bondStretch(u, x, bi, bj, L0);

  std::vector<double> force = bondForces(health,
                                         k_bond,
                                         bond_data.stretch,
                                         bond_data.direction,
                                         bi,
                                         bj,
                                         bar.num_total);

  std::vector<double> acceleration(bar.num_points);

  for (std::size_t p = 0; p < bar.num_points; ++p) {
    acceleration[p] = force[inside[p]] / rho;
  }

  double dissipated = 0.0;
  double work_in = 0.0;

  double power_prev = 0.0;

  for (const std::size_t i : boundary) {
    power_prev -= force[i] * v[i] * volume;
  }

  SimulationResults results;

  recordSnapshot(results,
                 0.0,
                 u,
                 v,
                 health,
                 bond_data.stretch,
                 dissipated,
                 work_in,
                 bar,
                 rho);

  std::size_t save_index = 1;

  std::cout << "Running " << n_steps << " time steps of " << dt << " s ("
            << params.bond_law << " bonds, strain rate " << rate
            << " /s) ...\n";

  const std::size_t print_every = std::max<std::size_t>(1, n_steps / 10);
  bool stopped_early = false;
  std::size_t final_step = 0;
  const auto start_time = std::chrono::steady_clock::now();

  // =========================================================== time loop

  for (std::size_t step = 1; step <= n_steps; ++step) {
    final_step = step;

    const double t = static_cast<double>(step) * dt;

    // 1-2. Half velocity update and position update.
    std::vector<double> v_half(bar.num_points);

    for (std::size_t p = 0; p < bar.num_points; ++p) {
      const std::size_t i = inside[p];
      v_half[p] = v[i] + 0.5 * dt * acceleration[p];
      u[i] += dt * v_half[p];
    }

    // 3. Boundary points follow the prescribed motion.

    moveBoundary(u, v, t, rate, x, boundary, L);

    // 4. Bond stretch and irreversible damage.
    bond_data = bondStretch(u, x, bi, bj, L0);

    std::vector<double> stretch_breakable(breakable.size());

    for (std::size_t q = 0; q < breakable.size(); ++q) {
      const std::size_t b = breakable[q];
      stretch_breakable[q] = bond_data.stretch[b];
      max_stretch[q] = std::max(max_stretch[q], stretch_breakable[q]);
    }

    std::vector<double> new_health =
        bondHealth(params.bond_law, max_stretch, s_elastic, s_break);

    for (std::size_t q = 0; q < breakable.size(); ++q) {
      new_health[q] = std::min(health_breakable[q], new_health[q]);

      const double health_loss = health_breakable[q] - new_health[q];

      dissipated += energy_coef_breakable[q] * health_loss *
                    stretch_breakable[q] * stretch_breakable[q];
      health_breakable[q] = new_health[q];
      health[breakable[q]] = new_health[q];
    }

    // 5-6. New forces and final velocity update.

    force = bondForces(health,
                       k_bond,
                       bond_data.stretch,
                       bond_data.direction,
                       bi,
                       bj,
                       bar.num_total);

    for (std::size_t p = 0; p < bar.num_points; ++p) {
      const std::size_t i = inside[p];
      acceleration[p] = force[i] / rho;
      v[i] = v_half[p] + 0.5 * dt * acceleration[p];
    }

    // Work done by the prescribed boundary motion.
    double power_now = 0.0;
    for (const std::size_t i : boundary) {
      power_now -= force[i] * v[i] * volume;
    }
    work_in += 0.5 * dt * (power_prev + power_now);
    power_prev = power_now;
    // Record snapshots.
    if (save_index < save_at.size() && step == save_at[save_index]) {
      recordSnapshot(results,
                     t,
                     u,
                     v,
                     health,
                     bond_data.stretch,
                     dissipated,
                     work_in,
                     bar,
                     rho);
      ++save_index;
      if (dissipated > 0.0 && results.dissipated.size() > n_quiet) {
        const std::size_t old = results.dissipated.size() - 1 - n_quiet;

        if (results.dissipated[old] > 0.0 &&
            dissipated - results.dissipated[old] <= quiet_tol * dissipated) {
          stopped_early = true;
          break;
        }
      }
    }

    if (step % print_every == 0) {
      const double progress =
          100.0 * static_cast<double>(step) / static_cast<double>(n_steps);

      std::cout << progress << " %"
                << "   t = " << t * 1e6 << " us"
                << "   fragments = " << results.num_fragments.back() << "\n";
    }
  }

  const auto end_time = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = end_time - start_time;

  results.x.reserve(bar.num_points);

  for (const std::size_t i : inside) {
    results.x.push_back(x[i]);
  }

  results.health = health;
  results.steps = final_step;
  results.stopped_early = stopped_early;
  results.run_time = elapsed.count();

  // ---------------------------------------------------------- energy check

  double max_error = 0.0;
  double max_total = 0.0;

  const double initial_total = results.kinetic.front() +
                               results.strain_energy.front() +
                               results.dissipated.front();

  for (std::size_t k = 0; k < results.t.size(); ++k) {
    const double total =
        results.kinetic[k] + results.strain_energy[k] + results.dissipated[k];

    const double error = (total - initial_total) - results.external_work[k];
    max_error = std::max(max_error, std::abs(error));
    max_total = std::max(max_total, std::abs(total));
  }

  results.energy_error =
      max_error / std::max(max_total, std::numeric_limits<double>::epsilon());

  std::cout << "Done in " << results.run_time << " s. "
            << results.num_fragments.back() << " fragments.";

  if (stopped_early) {
    std::cout << " (Stopped early at step " << final_step << " of " << n_steps
              << ": nothing was breaking any more.)";
  }

  std::cout << "\n";

  return results;
}