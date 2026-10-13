function results = simulate_bar(params, bar)
%SIMULATE_BAR  Pull the bar apart at a constant strain rate and follow it.
%
%   results = SIMULATE_BAR(params, bar) runs the simulation for the bar
%   built by MAKE_BAR and returns its history.
%
%   LOADING
%   The bar starts unstrained but already moving: every point has velocity
%       v = strain_rate * (x - L/2),
%   so the two halves move apart and the whole bar stretches evenly, at the
%   same rate everywhere. The boundary points keep exactly this motion for
%   the whole run. Because the bar is stretched evenly from the first
%   instant, no shock wave is sent in from the ends.
%
%   TIME STEPPING (velocity-Verlet)
%   Each step:  1. half-update velocities from the current forces
%               2. move the points
%               3. move the boundary points to where they should be
%               4. measure every bond's stretch and update its damage
%               5. compute the new forces
%               6. finish the velocity update
%
%   HOW LONG TO RUN
%   Long enough for an evenly stretched bond to reach its breaking stretch
%   six times over -- far past the moment the bar breaks up. The run stops
%   if nothing more is breaking: when the energy used up in breaking
%   bonds has grown by less than 0.1 % last 12 time steps.

rate = params.strain_rate;

%% ------------------------------------------------------------- shortcuts
x   = bar.x;         n_all = bar.num_total;
in  = bar.inside;    bd    = bar.boundary;
rho = params.density;
vol = bar.volume;    dt    = bar.dt;
L   = bar.L;

bi = bar.bond.i;     bj = bar.bond.j;     L0 = bar.bond.length;
k_bond = bar.bond.k; ecoef = bar.bond.energy_coef;
both_ends = [bi; bj];                 % to add each bond force to both points

brk     = find(bar.bond.can_break);   % the bonds allowed to break
s_el    = bar.bond.s_elastic(brk);
s_br    = bar.bond.s_break(brk);
ecoef_b = ecoef(brk);

%% ------------------------------------------------------------ run length
t_fail  = bar.mean_s_break / rate;    % time for an even stretch to break bonds
t_total = 6 * t_fail;
n_steps = max(1, ceil(t_total / dt));

n_snapshots = 120;     % how many times to record the state
n_quiet     = 12;      % stop if nothing broke over this many snapshots ...
quiet_tol   = 1e-3;    % ... meaning broken-bond energy grew by < 0.1 %

%% --------------------------------------------------------- initial state
u = zeros(n_all, 1);                  % displacement
v = rate * (x - L/2);                 % velocity: even stretching
health      = ones(bar.bond.count, 1);   % every bond intact
health_b    = health(brk);
max_stretch = zeros(numel(brk), 1);
[u, v] = move_boundary(u, v, 0.0, rate, x, bd, L);

%% -------------------------------------------------------------- storage
save_at = unique([0, round(linspace(0, n_steps, n_snapshots)), n_steps]);
ns = numel(save_at);
results.t             = zeros(1, ns);
results.displacement  = zeros(bar.num_points, ns);
results.damage        = zeros(bar.num_points, ns);
results.num_fragments = zeros(1, ns);
results.cuts          = cell(1, ns);
results.sizes         = cell(1, ns);
results.kinetic       = zeros(1, ns);
results.strain_energy = zeros(1, ns);
results.dissipated    = zeros(1, ns);
results.external_work = zeros(1, ns);

%% ------------------------------------------- forces at the starting point
[s, direction] = bond_stretch(u, x, bi, bj, L0);
f = bond_forces(health, k_bond, s, direction, both_ends, n_all);
a = f(in) / rho;

dissipated = 0;     % energy used up in damaging bonds
work_in    = 0;     % work done by the boundary points
power_prev = -sum(f(bd) .* v(bd)) * vol;

k_save = 1;
[results, k_save] = record(results, k_save, 0.0, u, v, health, s, ...
                           dissipated, work_in, bar, rho, in);

fprintf('Running %d time steps of %.3g s (%s bonds, strain rate %.3g /s) ...\n', ...
        n_steps, dt, params.bond_law, rate);
print_every = max(1, floor(n_steps / 10));
stopped_early = false;
tic;

%% =========================================================== time loop
for step = 1:n_steps
    t = step * dt;

    % 1-2. half-update velocities, then move the points
    v_half = v(in) + 0.5 * dt * a;
    u(in)  = u(in) + dt * v_half;

    % 3. boundary points follow the prescribed motion
    [u, v] = move_boundary(u, v, t, rate, x, bd, L);

    % 4. stretch and damage
    [s, direction] = bond_stretch(u, x, bi, bj, L0);
    s_b = s(brk);
    max_stretch = max(max_stretch, s_b);
    new_health  = bond_health(params.bond_law, max_stretch, s_el, s_br);
    new_health  = min(health_b, new_health);         % damage never heals
    % Energy loss.
    dissipated  = dissipated + sum(ecoef_b .* (health_b - new_health) .* s_b.^2);
    health_b    = new_health;
    health(brk) = health_b;

    % 5-6. new forces, finish the velocity update
    f = bond_forces(health, k_bond, s, direction, both_ends, n_all);
    a = f(in) / rho;
    v(in) = v_half + 0.5 * dt * a;

    % work done on the bar by the boundary points (trapezoid rule)
    power_now  = -sum(f(bd) .* v(bd)) * vol;
    work_in    = work_in + 0.5 * dt * (power_prev + power_now);
    power_prev = power_now;

    % record a snapshot, and stop if nothing has broken for a while
    if k_save <= ns && step == save_at(k_save)
        [results, k_save] = record(results, k_save, t, u, v, health, s, ...
                                   dissipated, work_in, bar, rho, in);
        k_old = (k_save - 1) - n_quiet;
        if dissipated > 0 && k_old >= 1 && results.dissipated(k_old) > 0 && ...
           (dissipated - results.dissipated(k_old)) <= quiet_tol * dissipated
            stopped_early = true;
            break
        end
    end

    if mod(step, print_every) == 0
        fprintf('  %3.0f %%   t = %7.3f us   fragments = %d\n', ...
                100 * step / n_steps, t * 1e6, results.num_fragments(k_save - 1));
    end
end
run_time = toc;

%% ------------------------------------------------------------- 
n = k_save - 1;                      % snapshots actually taken
for name = {'t','num_fragments','kinetic','strain_energy','dissipated','external_work'}
    results.(name{1}) = results.(name{1})(1:n);
end
results.displacement = results.displacement(:, 1:n);
results.damage       = results.damage(:, 1:n);
results.cuts  = results.cuts(1:n);
results.sizes = results.sizes(1:n);

results.x            = x(in);
results.health       = health;
results.steps        = step;
results.stopped_early = stopped_early;
results.run_time     = run_time;

% Energy check: kinetic + stored + used-up energy must equal what went in.
total = results.kinetic + results.strain_energy + results.dissipated;
error_abs = (total - total(1)) - results.external_work;
results.energy_error = max(abs(error_abs)) / max(max(abs(total)), eps);

fprintf('Done in %.1f s. %d fragments.', run_time, results.num_fragments(end));
if stopped_early
    fprintf(' (Stopped early at step %d of %d: nothing was breaking any more.)', ...
            step, n_steps);
end
fprintf('\n');
end


%% =====================================================================
function [s, direction] = bond_stretch(u, x, bi, bj, L0)
% Stretch of every bond = (current length - original length) / original.
y = x + u;                       % current positions
d = y(bj) - y(bi);
s = (abs(d) - L0) ./ L0;
direction = sign(d);
end

function f = bond_forces(health, k_bond, s, direction, both_ends, n_all)
% Force in every bond, added to both of its points with opposite signs
fb = health .* k_bond .* s .* direction;
f  = accumarray(both_ends, [fb; -fb], [n_all 1]);
end

function [u, v] = move_boundary(u, v, t, rate, x, bd, L)
% The boundary points follow the even-stretching motion exactly.
xb    = x(bd);
u(bd) = (rate * t) * (xb - L/2);
v(bd) = rate * (xb - L/2);
end

function [res, k] = record(res, k, t, u, v, health, s, dissipated, work_in, bar, rho, in)
% Store one snapshot of the bar.
res.t(k) = t;
res.displacement(:, k) = u(in);

% Damage at a point = fraction of its bonds (weighted by overlap) that is
% gone: 0 = all bonds intact, 1 = all broken.
w   = bar.volume * bar.bond.overlap;
ends = [bar.bond.i; bar.bond.j];
num = accumarray(ends, [health .* w; health .* w], [bar.num_total 1]);
den = accumarray(ends, [w; w], [bar.num_total 1]);
dmg = 1 - num ./ max(den, eps);
res.damage(:, k) = dmg(in);

[res.num_fragments(k), res.cuts{k}, res.sizes{k}] = count_fragments(bar, health);

res.kinetic(k)       = 0.5 * rho * bar.volume * sum(v(in).^2);
res.strain_energy(k) = sum(bar.bond.energy_coef .* health .* s.^2);
res.dissipated(k)    = dissipated;
res.external_work(k) = work_in;
k = k + 1;
end
