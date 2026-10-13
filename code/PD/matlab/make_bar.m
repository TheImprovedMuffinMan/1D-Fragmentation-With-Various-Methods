function bar = make_bar(params)
%MAKE_BAR  Build the bar: its points, the bonds between them, and their
%   stiffness and strength.
%
%   bar = MAKE_BAR(params) uses the inputs set in RUN_BAR_EXAMPLE and returns
%   everything the time loop needs in one struct.
%
%   The bar is a row of equally spaced points. Each point is joined by a
%   spring (a "bond") to EVERY neighbour within a distance called the
%   horizon - not only to the next point. A bond carries force while it is
%   intact, softens once it is stretched past its elastic limit, and breaks
%   for good at its breaking stretch. Wherever every bond across a gap has
%   broken, the bar has split into two pieces.
%
%   A few extra points sit beyond each end ("ghost nodes"). They are not
%   free: the simulation moves them by hand to pull the bar apart.
%   basically to apply boundary conditions

%% ---------------------------------------------------------------- points
num_points = round(params.length / params.dx);   % points inside the bar
L  = num_points * params.dx;     % bar length, rounded to a whole number of cells
dx = L / num_points;             % spacing between points
horizon = params.m * params.dx;  % how far each point reaches

% Each point owns a cell dx wide, so a neighbour whose centre lies up to
% half a cell beyond the horizon still overlaps it. That sets how many
% neighbours a point has on each side, and so how many boundary points are
% needed at each end.
radius = dx / 2;
reach  = floor((horizon + radius) / dx + 1e-12);
num_boundary = reach;

num_total = num_points + 2 * num_boundary;
x = ((1:num_total)' - num_boundary - 0.5) * dx;     % positions, left to right

inside   = (num_boundary + 1 : num_boundary + num_points)';
boundary = [(1:num_boundary)'; (num_boundary + num_points + 1 : num_total)'];

% In 1D every point stands for a small cube of material, dx on each side.
area   = dx * dx;
volume = area * dx;

%% ----------------------------------------------------------------- bonds
% List every pair (i, j) with j - i = 1, 2, ..., reach, shortest bonds first.
num_bonds = sum(num_total - (1:reach));
bi = zeros(num_bonds, 1);
bj = zeros(num_bonds, 1);
filled = 0;
for gap = 1:reach
    n = num_total - gap;
    bi(filled + (1:n)) = (1:n)';
    bj(filled + (1:n)) = (1:n)' + gap;
    filled = filled + n;
end
bond_length = x(bj) - x(bi);

% Partial-volume correction: a neighbour near the edge of the horizon only
% partly overlaps it, so it only counts partly (a fraction between 0 and 1).
overlap = (horizon + radius - bond_length) / (2 * radius);
overlap(bond_length <= horizon - radius) = 1.0;
overlap(bond_length >  horizon + radius) = 0.0;
overlap = min(1.0, max(0.0, overlap));

% Weighting along the bond. 'constant': every bond in the horizon counts the
% same. 'conical': short bonds count more, falling to zero at the horizon.
switch lower(params.bond_shape)
    case 'constant'
        weight = ones(num_bonds, 1);
    case 'conical'
        weight = max(0, 1 - bond_length / horizon);
        % Bonds that land on the horizon edge itself must get exactly zero
        % (rounding can leave them at 1e-16).
        weight(bond_length >= horizon * (1 - 1e-12)) = 0;
    otherwise
        error('bond_shape must be ''constant'' or ''conical''.');
end

keep = (overlap > 0) & (weight > 0);
bi = bi(keep);  bj = bj(keep);  bond_length = bond_length(keep);
overlap = overlap(keep);  weight = weight(keep);

%% ---------------------------------------------------------- bond stiffness
% Start from the textbook value, chosen so the bonds store the same energy
% as an elastic bar, (1/2) E strain^2 ...
E  = params.E;
kc = E / (area * horizon^2);
switch lower(params.bond_shape)
    case 'constant', c = 2 * kc;
    case 'conical',  c = 6 * kc;
end
k_bond = c * weight * volume .* overlap;      % force = k_bond * stretch

% ... then scale it so that the bar AS BUILT, with a finite number of
% points, has a Young's modulus of exactly E.
mid       = num_boundary + round(num_points / 2);
touch     = (bi == mid) | (bj == mid);
E_as_built = 0.5 * sum(k_bond(touch) .* bond_length(touch));
k_bond    = k_bond * (E / E_as_built);

% Energy stored in a bond = energy_coef * health * stretch^2.
% This is the spring energy (1/2) K x^2 written with the stretch s instead
% of the extension x. Treat one bond as an ordinary spring:
%   force      F = k_bond * V * s      (k_bond*s is force per unit volume;
%                                       times the point volume V gives N)
%   extension  x = s * L0              (L0 = original bond length)
%   so         K = F / x = k_bond * V / L0
%   energy     (1/2) K x^2 = (1/2) * k_bond * L0 * V * s^2 = energy_coef * s^2
% A damaged bond has stiffness health * K, hence the health factor.
% energy_coef is in joules: the energy the bond would hold at a stretch of 1.
energy_coef = 0.5 * k_bond .* bond_length * volume;

%% --------------------------------------------------------- bond strength
% Strength and toughness are given point by point (with a little random
% scatter); each bond takes the average of its two end points.
[strength_pt, toughness_pt] = random_field(params, x, dx, horizon, L);
strength  = 0.5 * (strength_pt(bi)  + strength_pt(bj));
toughness = 0.5 * (toughness_pt(bi) + toughness_pt(bj));

% Elastic limit: the stretch at which the bond force peaks. Setting it to
% strength/E makes the bar's peak stress equal the input strength.
s_elastic = strength / E;

% Breaking stretch, chosen so that cutting the bar costs exactly Gc of
% energy per unit area. A cut must break every bond that crosses it; bonds
% that skip g points cross g possible cut positions, so the energy of one
% cut is the sum over bond lengths of (number crossing) x (energy each).
gap      = round(bond_length / dx);
[g, ia]  = unique(gap);
cut_sum  = sum(double(g) .* energy_coef(ia));
target   = toughness * area / cut_sum;

switch lower(params.bond_law)
    case 'bilinear', s_break = target ./ s_elastic;   % softening bond
    case 'pmb',      s_break = sqrt(target);          % brittle bond, snaps
    otherwise, error('bond_law must be ''bilinear'' or ''pmb''.');
end
mean_s_break = mean(s_break);      % sets the run length (SIMULATE_BAR)

% Safety: a softening bond must break after it starts softening.
if strcmpi(params.bond_law, 'bilinear')
    bad = s_break <= s_elastic;
    s_break(bad) = s_elastic(bad) * (1 + 1e-6);
end

% Only bonds between two inside points may break. Bonds to boundary points
% are held together so the bar is always gripped at its ends.
both_inside = (bi > num_boundary) & (bi <= num_boundary + num_points) & ...
              (bj > num_boundary) & (bj <= num_boundary + num_points);

%% -------------------------------------------------------------- time step
% The largest time step that keeps the calculation stable is set by the
% stiffest point. We use half of it to be stable (can be changed to even lower with contact).
spring = k_bond ./ bond_length;
total_spring = accumarray(bi, spring, [num_total 1]) + accumarray(bj, spring, [num_total 1]);
dt_max = sqrt(2 * params.density / max(total_spring(inside)));
dt     = 0.5 * dt_max;

%% ------------------------------------------------------------------ pack
bar.L = L;  bar.dx = dx;  bar.horizon = horizon;
bar.num_points = num_points;  bar.num_boundary = num_boundary;
bar.num_total  = num_total;
bar.x = x;  bar.inside = inside;  bar.boundary = boundary;
bar.area = area;  bar.volume = volume;
bar.wave_speed = sqrt(E / params.density);
bar.dt = dt;

bar.bond.i = bi;  bar.bond.j = bj;  bar.bond.length = bond_length;
bar.bond.overlap = overlap;
bar.bond.k = k_bond;  bar.bond.energy_coef = energy_coef;
bar.bond.s_elastic = s_elastic;  bar.bond.s_break = s_break;
bar.bond.can_break = both_inside;
bar.bond.count = numel(bi);
bar.mean_s_break = mean_s_break;
bar.strength_pt = strength_pt;  bar.toughness_pt = toughness_pt;
end
