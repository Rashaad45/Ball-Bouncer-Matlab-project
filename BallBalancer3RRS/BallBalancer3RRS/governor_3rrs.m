function [q_cmd, theta_cmd, flags] = governor_3rrs(q_des, q_prev_cmd, theta_prev_cmd, P, W, use_gov)
% GOVERNOR_3RRS  Kinematics-aware command shaping (the "proposed" layer).
%
%   1. WORKSPACE GOVERNOR:  the commanded tilt vector is scaled (direction
%      kept) so its magnitude never exceeds  gov_margin * tilt_max(psi)
%      from the reachable-tilt lookup table W.
%   2. SERVO-RATE GOVERNOR: the per-frame servo step required by the new
%      command is predicted with the IK Jacobian.  If any servo would need
%      more than omega_max*Ts, the WHOLE tilt increment is scaled by the
%      same factor.  Plain per-servo saturation would clip the legs
%      unequally and tilt the plate in the wrong direction (cross-axis
%      coupling); uniform scaling keeps the tilt direction correct.
%   With use_gov = false only the software ceiling P.tilt_max_sw applies.

flags = struct('ws', false, 'rate', false);
q_cmd = q_des(:);

m  = norm(q_cmd(1:2));
lim = P.tilt_max_sw;
if use_gov
    psi = atan2(q_cmd(2), q_cmd(1));
    lim = min(lim, P.gov_margin * interp1(W.psi, W.tilt_max, psi));
end
if m > lim
    q_cmd(1:2) = q_cmd(1:2) * lim / m;
    flags.ws = use_gov;
end

[theta_cmd, ok] = ik_3rrs(q_cmd, P);

if use_gov
    dth   = theta_cmd - theta_prev_cmd(:);
    dmax  = 0.9 * P.omega_max * P.Ts_cam;      % allowed servo step per frame
    worst = max(abs(dth));
    if worst > dmax
        s = dmax / worst;                       % same scale for all legs
        q_cmd = q_prev_cmd(:) + s * (q_cmd - q_prev_cmd(:));
        [theta_cmd, ok] = ik_3rrs(q_cmd, P);
        flags.rate = true;
    end
end
if ~ok                                         % '~' is logical NOT
    theta_cmd = min(max(theta_cmd, P.theta_min), P.theta_max);
end
end
