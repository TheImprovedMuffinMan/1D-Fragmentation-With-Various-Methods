function [num_fragments, cut_x, sizes] = count_fragments(bar, health)
%COUNT_FRAGMENTS  Find where the bar has split, and the size of each piece.
%
%   [num_fragments, cut_x, sizes] = COUNT_FRAGMENTS(bar, health)
%
%   The bar is split in the gap between two neighbouring points when NO
%   unbroken bond spans that gap. (A bond from point i to point j spans the
%   gaps i, i+1, ..., j-1.) So we count, for every gap, how many working
%   bonds cross it, and the gaps with zero are the cuts.
%
%   Very small pieces are merged into their neighbour: two cuts closer than
%   one horizon count as a single cut. A piece thinner than the horizon is
%   below the smallest scale the model can describe, and in a region where
%   the bar is crushed every single gap would otherwise count as a cut.
%
%   Outputs
%     num_fragments   number of pieces
%     cut_x           positions of the cuts [m]
%     sizes           length of each piece [m]; they add up to the bar length

N = bar.num_points;

% Working bonds between two inside points, numbered 1..N along the bar.
use = (health > 0) & bar.bond.can_break;
i = bar.bond.i(use) - bar.num_boundary;
j = bar.bond.j(use) - bar.num_boundary;

% Each bond adds +1 at its left end and -1 at its right end; a running sum
% then gives the number of bonds crossing every gap.
if isempty(i)
    change = zeros(N, 1);
else
    one = ones(numel(i), 1);
    change = accumarray(i, one, [N 1]) - accumarray(j, one, [N 1]);
end
crossing = cumsum(change);

gaps   = find(crossing(1:N-1) == 0);            % gaps with nothing across
x_in   = bar.x(bar.num_boundary + (1:N));
cut_raw = 0.5 * (x_in(gaps) + x_in(gaps + 1));  % cut sits midway in the gap
cut_raw = cut_raw(:);

% Merge cuts closer together than one horizon.
min_size = bar.horizon;
keep = false(numel(cut_raw), 1);
last = 0;
for a = 1:numel(cut_raw)
    if cut_raw(a) - last >= min_size
        keep(a) = true;
        last = cut_raw(a);
    end
end
cut_x = cut_raw(keep);
% A last sliver thinner than the horizon belongs to the piece before it.
while ~isempty(cut_x) && (bar.L - cut_x(end)) < min_size
    cut_x(end) = [];
end

num_fragments = numel(cut_x) + 1;
sizes = diff([0; cut_x; bar.L]);
end
