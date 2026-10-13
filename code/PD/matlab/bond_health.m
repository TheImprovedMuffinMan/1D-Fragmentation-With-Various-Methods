function health = bond_health(law, max_stretch, s_elastic, s_break)
%BOND_HEALTH  How much of each bond still works: 1 = intact, 0 = broken.
%
%   health = BOND_HEALTH(law, max_stretch, s_elastic, s_break)
%
%   The force in a bond is always   health * k * stretch,   so health scales
%   the bond's stiffness down as it gets damaged. It depends on the LARGEST
%   stretch the bond has ever reached, so damage never heals.
%
%   'bilinear'  The force rises in a straight line up to the elastic limit
%               s_elastic, then falls in a straight line to zero at s_break.
%               In between, the bond is softening (partly damaged).
%
%                 force
%                   |      /\
%                   |     /  \
%                   |    /    \
%                   |   /      \
%                   |__/________\______ stretch
%                     0  s_elastic  s_break
%
%   'pmb'       The bond stays fully elastic and snaps at s_break
%               (the "prototype microelastic brittle" bond of Silling).

health = ones(size(max_stretch)); %beta from the paper

switch lower(law)
    case 'bilinear'
        % On the falling branch, force = k*s_elastic*(s_break - s)/(s_break - s_elastic).
        % Dividing by k*s gives the health factor.
        soft = max_stretch > s_elastic;
        health(soft) = s_elastic(soft) .* (s_break(soft) - max_stretch(soft)) ./ ...
                       (max_stretch(soft) .* (s_break(soft) - s_elastic(soft))); 
    case 'pmb'
        health(max_stretch > s_break) = 0;
end

health = min(1.0, max(0.0, health));
end
