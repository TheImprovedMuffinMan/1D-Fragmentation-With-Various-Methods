function th = fragment_size_theory(params)
%FRAGMENT_SIZE_THEORY  Predicted average fragment size, for comparison.
%
%   th = FRAGMENT_SIZE_THEORY(params) returns, for the material and strain
%   rate in params:
%
%     th.zmr    Zhou, Molinari & Ramesh (2006), their Eq. 29a. A fit to many
%               1D simulations of the same bar with cohesive cracks. This is
%               the main reference: our bar should land close to it.
%     th.grady  Grady (1982), from an energy balance. An upper limit: it
%               ignores how long a crack takes to open, so it predicts
%               fragments that are too big, especially at low strain rate.
%
%   and the two natural scales of the problem:
%     th.s0       = E*Gc/strength^2          a natural fragment length [m]
%     th.rate0    = c*strength^3/(E^2*Gc)    a natural strain rate [1/s]
%     th.rate_ratio = strain_rate / rate0   below ~0.01 is "slow", above ~1 "fast"
%
%   params.strain_rate may also be a list of rates; the results are then lists.

E = params.E;  rho = params.density;
sigma = params.strength;  G = params.toughness;
rate = params.strain_rate;

c = sqrt(E / rho);                          
th.s0    = E * G / sigma^2;
th.rate0 = c * sigma^3 / (E^2 * G);
th.rate_ratio = rate / th.rate0;

r = th.rate_ratio;
th.zmr   = th.s0 * 4.5 ./ (1 + 0.77 * r.^0.25 + 5.4 * r.^0.75);
th.grady = (24 * G ./ (rho * rate.^2)).^(1/3);
end
