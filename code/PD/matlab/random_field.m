function [strength, toughness] = random_field(params, x, dx, horizon, L)
%RANDOM_FIELD  Give every point its own strength and toughness.
%
%   Real materials are never perfectly uniform: grains and small flaws make
%   some spots a little weaker than others, and cracks start there. Here
%   each point gets a strength and a toughness that vary randomly around the
%   input values, with standard deviation params.scatter (0.01 = 1 %).
%   With params.scatter = 0 the bar is perfectly uniform.
%
%   The variation is SMOOTH: neighbouring points get similar values, and
%   values become independent only over a distance
%       l_c = E * Gc / strength^2      (about 0.3 mm for the ceramic),
%   the natural size of the zone where a crack forms. Without this, every
%   grid point would be an independent weak spot, and the number of
%   fragments would change when you change the grid.
%
%   The same params.seed always gives the same material -- even if you
%   change dx -- because the random numbers are drawn on their own fixed
%   grid and only then read off at the points.

E  = params.E;
mu_strength  = params.strength;
mu_toughness = params.toughness;

if params.scatter == 0
    strength  = mu_strength  * ones(numel(x), 1);
    toughness = mu_toughness * ones(numel(x), 1);
    return
end

% Our own random-number stream, so the result depends only on the seed.
stream = RandStream('mt19937ar', 'Seed', params.seed);

l_c = E * mu_toughness / mu_strength^2;          % correlation length [m]

z_strength  = smooth_noise(stream, x(:), l_c, dx, horizon, L);
z_toughness = smooth_noise(stream, x(:), l_c, dx, horizon, L);

strength  = mu_strength  + params.scatter * mu_strength  * z_strength;
toughness = mu_toughness + params.scatter * mu_toughness * z_toughness;

% Never let a point drop below 5 % of the average (keeps values positive).
strength  = max(strength,  0.05 * mu_strength);
toughness = max(toughness, 0.05 * mu_toughness);
end

% -------------------------------------------------------------------------
function z = smooth_noise(stream, x, l_c, dx, horizon, L)
%SMOOTH_NOISE  Random values with mean 0, standard deviation 1, that vary
%   smoothly over a distance l_c. Returned at the positions x.
N = numel(x);
if l_c <= 0 || l_c < dx
    z = randn(stream, N, 1);        % smoother than the grid is impossible
    return
end

% 1. Independent random numbers on a fine helper grid that covers the bar
%    plus a margin on each side. The helper grid is tied to the bar and to
%    l_c, never to dx.
spacing = l_c / 8;
margin  = max(8 * l_c, 4 * horizon);
xr = (-margin) : spacing : (L + margin);
zr = randn(stream, numel(xr), 1);

% 2. Smooth them: replace each value by a bell-shaped (Gaussian) weighted
%    average of its neighbours, about l_c wide.
width = l_c / spacing;
h = ceil(3 * width);
bell = exp(-0.5 * ((-h:h) / width).^2);
bell = bell(:) / sum(bell);

M  = numel(zr);
k  = (1:M)';
zs = zeros(M, 1);
for q = -h:h
    zs = zs + bell(q + h + 1) * zr(min(max(k + q, 1), M));
end
zs = zs / sqrt(sum(bell.^2));        % back to standard deviation 1

% 3. Read the smooth values off at the bar's points.
z = interp1(xr(:), zs, x, 'linear', 'extrap');
end
