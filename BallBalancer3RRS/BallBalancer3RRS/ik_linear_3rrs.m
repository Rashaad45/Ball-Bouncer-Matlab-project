function [theta, ok] = ik_linear_3rrs(q, P)
% IK_LINEAR_3RRS  BASELINE "hobby" IK used by many ball-balancer builds.
%
%   Assumes the coupler link stays vertical, so the horn tip simply has
%   to rise by the vertical displacement of the ball joint above it:
%       dz_i   = rp*( sin(phi_i) sin(alpha) - cos(phi_i) sin(beta) ) + (h - h0)
%       theta_i = asin( dz_i / L1 )
%   It ignores (a) the inclination of the coupler link, (b) the parasitic
%   x/y/yaw motion.  It is included ONLY as the comparison baseline.

alpha = q(1); beta = q(2); h = q(3);
dz = P.rp * (sin(P.phi(:)) * sin(alpha) - cos(P.phi(:)) * sin(beta)) + (h - P.h0);
s  = dz / P.L1;
ok = all(abs(s) <= 1);             % all() is true if every element is true
s  = max(-1, min(1, s));           % clamp into asin's domain
theta = asin(s);
ok = ok && all(theta >= P.theta_min & theta <= P.theta_max);
end
