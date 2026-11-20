function [theta, side_len, Pj] = ik_fixedcentre_3rrs(q, P)
% IK_FIXEDCENTRE_3RRS  The other common shortcut: "plate normal vector" IK
%   with the plate centre pinned on the z-axis (as in several open-source
%   ball balancers).  Each ball joint is placed at distance rp from the
%   centre, perpendicular to the plate normal n and inside its leg plane.
%
%   Because the true 3-RRS centre DOES move sideways, those three points
%   no longer form a rigid equilateral triangle.  side_len returns the
%   three triangle sides so you can see the violation (rigid = sqrt(3)*rp).

R0 = rot_3rrs(q(1), q(2), 0);
n  = R0(:,3);                       % plate normal = 3rd column of R
Pj = zeros(3,3); theta = zeros(3,1);
for i = 1:3
    u = [cos(P.phi(i)); sin(P.phi(i)); 0];
    w = n(3)*u - (n'*u)*[0;0;1];    % in span{u, z} and perpendicular to n
    w = w / norm(w);
    Pi = [0;0;q(3)] + P.rp * w;
    Pj(:,i) = Pi;
    dr = u'*Pi - P.rb;  dz = Pi(3);  M = hypot(dr,dz);
    K  = (M^2 + P.L1^2 - P.L2^2)/(2*P.L1);
    theta(i) = atan2(dz,dr) - acos(max(-1,min(1,K/M)));
end
side_len = [norm(Pj(:,1)-Pj(:,2)); norm(Pj(:,2)-Pj(:,3)); norm(Pj(:,3)-Pj(:,1))];
end
