%% RUN_BAR_EXAMPLE  A bar pulled apart at a high strain rate.
%
% This is the file to run. It builds a 1D bar, pulls its two ends apart at a
% constant strain rate, and follows the bar as it breaks into fragments. It
% then prints how many fragments formed, compares their average size with
% theory, and plots what happened.
%
% Everything you might want to change is in section 1, with its unit.
% The default case takes roughly a minute on a laptop.
%
% main and functions required:
%   run_bar_example.m       this script (main)
%   make_bar.m              builds the points, bonds and their properties
%   random_field.m          small random variation of strength and toughness
%   simulate_bar.m          the time loop
%   bond_health.m           how a bond softens and breaks
%   count_fragments.m       finds the cuts and the fragment sizes
%   fragment_size_theory.m  formulas to compare against
%   plot_results.m          the figures 

clear; close all;

%% 1. Inputs ---------------------------------------------------------------

% Material props
params.density   = 2750;       % kg/m^3
params.E         = 275e9;      % Young's modulus [Pa]
params.strength  = 300e6;      % tensile strength [Pa]
params.toughness = 100;        % fracture energy Gc [J/m^2]: the energy
                               % needed to open one square metre of crack

% Bar geom
params.length = 0.050;         % bar length [m]  (50 mm)
params.dx     = 3e-6;          % spacing between points (structured)
params.m      = 5.015;         % horizon = m * dx: how many points each point
                               % reaches on either side. (can be changed for mesh convergence)

% material model (bilinear is better than PMB)
params.bond_law   = 'pmb';   % 'bilinear': bonds soften, then break (like CZM linear softening model)
                                  % 'pmb'     : bonds snap suddenly no softe
params.bond_shape = 'constant';   % 'constant': all bonds weighted equally (from what i found this has little effect)
                                  % 'conical' : shorter bonds weighted more

% Small random variation in strength and toughness 
params.scatter = 0.01;         % standard deviation as a fraction (0.01 = 1 %)
params.seed    = 1;           

% Loading
params.strain_rate = 1e5;     

%% 2. Build the bar and run ------------------------------------------------
bar = make_bar(params); % creates neighbours, ghost nodes for disp boundary conditions
fprintf('Bar: %d points, %d bonds, horizon %.1f um, time step %.3g s\n', ...
        bar.num_points, bar.bond.count, bar.horizon * 1e6, bar.dt);

results = simulate_bar(params, bar); % main kernel

%% 3. Results --------------------------------------------------------------
th    = fragment_size_theory(params);
sizes = results.sizes{end};

fprintf('\n------------------------------------------------------------\n');

plot_results(params, bar, results);
