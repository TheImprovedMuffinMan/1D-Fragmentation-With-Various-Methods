function [fig_size, fig_energy, fig_count] = plot_results(params, bar, results)
%PLOT_RESULTS  Three figure windows:
%   1. the average fragment size against Zhou, Molinari & Ramesh (2006)
%   2. where the energy goes during the run
%   3. the number of fragments over time

%%%%%% Change this as per required using results struct %%%%%%%%%%

%% ---------------------------------------- 1. fragment size against theory
% The average fragment size predicted by Zhou et al. (2006,
% Eq. 29a) for this material over a wide range of strain rates. T
avg = mean(results.sizes{end});           % average fragment size [m]
th  = fragment_size_theory(params);

lo = min(th.rate0 * 1e-4, params.strain_rate / 10);
hi = max(th.rate0 * 1e2,  params.strain_rate * 10);
curve = params;
curve.strain_rate = logspace(log10(lo), log10(hi), 300);
curve_th = fragment_size_theory(curve);

fig_size = figure('Color', 'w', 'Position', [80 120 720 540], ...
                  'Name', 'Fragment size');
loglog(curve.strain_rate, curve_th.zmr * 1e6, 'k-', 'LineWidth', 2); hold on;
loglog(params.strain_rate, avg * 1e6, 'o', 'MarkerSize', 11, ...
       'MarkerFaceColor', [0.85 0.2 0.1], 'MarkerEdgeColor', 'k');
grid on;
xlabel('strain rate  [1/s]');
ylabel('average fragment size  [\mum]');
legend('Zhou, Molinari & Ramesh (2006), Eq. 29a', ...
       sprintf('this simulation: %.0f \\mum, %d fragments', ...
               avg * 1e6, results.num_fragments(end)), ...
       'Location', 'southwest');
title(sprintf(params.bond_law, params.bond_shape, avg / th.zmr));

%% ------------------------------------------------------------ 2. energy
% Energy per square metre of the bar's cross-section. The bar is already
% moving at the start, so its kinetic energy is large; what matters is how
% much it CHANGES. The work done by the ends must equal the change in
% kinetic energy plus the energy stored in bonds plus the energy used up
% breaking them -- the dotted line checks that.
%
% Once the bar has broken up, the energy used up breaking bonds stops
% growing, but the work done and the stored energy keep rising: the grips
% go on pulling the two end pieces, and the bonds joining them to the grips
% are never allowed to break, so they keep stretching. That is the grips,
% not the fragmentation, and it does not change the fragment count.
t_us     = results.t * 1e6;
per_area = 1 / bar.area;
work   = results.external_work * per_area;
stored = results.strain_energy * per_area;
broken = results.dissipated    * per_area; % i think this is similar to cohesive energy
d_kin  = (results.kinetic - results.kinetic(1)) * per_area;
check  = d_kin + stored + broken - work;

fig_energy = figure('Color', 'w', 'Position', [820 120 720 540], ...
                    'Name', 'Energy');
plot(t_us, work,   'k-', 'LineWidth', 2); hold on;
plot(t_us, stored, 'g-', 'LineWidth', 1.8);
plot(t_us, broken, 'r-', 'LineWidth', 1.8);
plot(t_us, d_kin,  'b-', 'LineWidth', 1.8);
plot(t_us, check,  'k:', 'LineWidth', 1.5);
grid on;
xlabel('time  [\mus]');
ylabel('energy per unit area  [J/m^2]');
legend('work done by the ends', 'stored in bonds', 'used up breaking bonds', ...
       'change in kinetic energy', 'balance (should stay near 0)', ...
       'Location', 'best');
title('Energy');

%% ------------------------------------------------- 3. number of fragments
fig_count = figure('Color', 'w', 'Position', [450 60 720 540], ...
                   'Name', 'Number of fragments');
plot(t_us, results.num_fragments, 'LineWidth', 2);
grid on;
xlabel('time  [\mus]');
ylabel('number of fragments');
title("number of fragments");
end
