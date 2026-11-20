function [theta, ok, info] = ik_3rrs(q, P)
% IK_3RRS  EXACT inverse kinematics of the 3-RRS ball-balancing platform.
%
%   [theta, ok, info] = ik_3rrs(q, P)
%     q     = [alpha; beta; h]  roll, pitch [rad], heave [m]
%     theta = 3x1 servo angles [rad], 0 = horn horizontal, + = horn up
%     ok    = true if all legs are reachable AND inside servo limits
%     info  = struct with parasitic motion, joint positions, diagnostics
%
%   Algorithm (closed form per leg, after the parasitic solve):
%     1. Solve parasitic x, y, gamma  ->  full pose (R, T).
%     2. Ball joint in base frame:  P_i = T + R p_i.
%     3. Express P_i in the leg plane: radial coordinate rho, height z,
%        relative to the servo shaft B_i:   d = [rho - rb ; z].
%     4. Elbow C_i = B_i + L1 [cos(theta); sin(theta)] (in-plane).
%        Closing the link |P_i - C_i| = L2 gives
%           d_r cos(theta) + d_z sin(theta) = K,
%           K = (|d|^2 + L1^2 - L2^2) / (2 L1)
%        i.e.  |d| cos(theta - psi) = K,  psi = atan2(d_z, d_r)
%        =>  theta = psi - acos(K/|d|)   ("elbow-out" branch, theta=0 at home)
%     5. Reachable iff |K| <= |d|  (otherwise the link is too short/long).

alpha = q(1); beta = q(2); h = q(3);

[x, y, gamma, R, res] = parasitic_3rrs(alpha, beta, P);   % step 1
T = [x; y; h];

theta = zeros(3,1);   ok = true;   % preallocate output (faster, clearer)
Pj = zeros(3,3); Cj = zeros(3,3); Bj = zeros(3,3); margin = zeros(3,1);

for i = 1:3
    u  = [cos(P.phi(i)); sin(P.phi(i)); 0];     % radial unit vector of leg i
    pi_plat = P.rp * u;                          % ball joint, platform frame
    Pi = T + R * pi_plat;                        % step 2
    Bi = P.rb * u;                               % servo shaft position

    dr = u' * Pi - P.rb;                         % step 3 (u' * v = dot product)
    dz = Pi(3);
    M  = hypot(dr, dz);                          % |d| without overflow
    K  = (M^2 + P.L1^2 - P.L2^2) / (2*P.L1);     % step 4

    ratio = K / M;
    margin(i) = 1 - abs(ratio);                  % >0 means reachable
    if abs(ratio) > 1                            % step 5
        ok = false;
        ratio = sign(ratio);                     % clamp so acos stays real
    end
    th = atan2(dz, dr) - acos(ratio);
    theta(i) = th;

    if th < P.theta_min || th > P.theta_max      % '||' = short-circuit OR
        ok = false;
    end

    Pj(:,i) = Pi;   Bj(:,i) = Bi;
    Cj(:,i) = Bi + P.L1 * (cos(th) * u + sin(th) * [0;0;1]);
end

if nargout > 2           % nargout = number of outputs the caller asked for
    info.x = x; info.y = y; info.gamma = gamma; info.R = R; info.T = T;
    info.P = Pj; info.C = Cj; info.B = Bj;
    info.reach_margin = margin;     % distance from the link-length singularity
    info.parasitic_res = res;       % should be ~1e-16
    info.link_err = sqrt(sum((Pj - Cj).^2, 1)) - P.L2;  % should be ~0
end
end
