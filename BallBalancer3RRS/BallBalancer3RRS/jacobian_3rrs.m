function J = jacobian_3rrs(q, P, th0)
% JACOBIAN_3RRS  Inverse-kinematic Jacobian  J = d theta / d q  (3x3).
%
%   theta_dot = J * q_dot   maps platform rates [alpha_dot; beta_dot; h_dot]
%   to servo rates.  Used for:
%     * Newton FK (fk_3rrs)
%     * servo-speed demand of a planned tilt trajectory
%     * dexterity: cond(J) -> infinity near a singularity
%
%   Computed by central finite differences of the exact IK (accurate to
%   ~1e-9, and automatically consistent with any geometry change).

if nargin < 3
    th0 = ik_3rrs(q, P);   %#ok<NASGU>  (comment silences an editor warning)
end
J = zeros(3,3);
e = [1e-6; 1e-6; 1e-7];            % step sizes: rad, rad, m
for j = 1:3
    d = zeros(3,1); d(j) = e(j);
    J(:,j) = (ik_3rrs(q + d, P) - ik_3rrs(q - d, P)) / (2*e(j));
end
end
