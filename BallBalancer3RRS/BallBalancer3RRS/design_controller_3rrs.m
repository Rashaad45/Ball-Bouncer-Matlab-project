function C = design_controller_3rrs(P)
% DESIGN_CONTROLLER_3RRS  Per-axis LQR-with-integral + steady-state Kalman.
%
%   Linearised ball dynamics (one axis, plate frame), with the plate tilt
%   modelled as a first-order lag of the commanded tilt u:
%
%       x_dot     = v
%       v_dot     = kb*g*tilt - cv*v          (kb = 5/7 rolling sphere)
%       tilt_dot  = (u - tilt)/tau_s
%
%   The y axis is identical with "virtual tilt" = -alpha, so ONE design
%   serves both axes (the 3-RRS is symmetric).  Discretised exactly with
%   the matrix exponential at the camera rate Ts (zero-order hold).
%
%   Augmented with an integrator  z(k+1) = z(k) - Ts*e(k)  to remove
%   steady-state error from plate offsets, friction and IK model error.

Ts  = P.Ts_cam;
tau = P.tau_s;
Ac = [0 1 0;
      0 -P.cv P.kb*P.g;
      0 0 -1/tau];
Bc = [0; 0; 1/tau];

% Exact ZOH discretisation via one matrix exponential of the block matrix
% [Ac Bc; 0 0]  ->  expm gives [Ad Bd; 0 1].
M  = expm([Ac Bc; zeros(1,4)] * Ts);
Ad = M(1:3,1:3);  Bd = M(1:3,4);
Cd = [1 0 0];

% ---- LQR with integral action (Bryson's rule for weights) ----
Aa = [Ad zeros(3,1); -Ts*Cd 1];
Ba = [Bd; 0];
Q  = diag(1 ./ P.lqr_qmax.^2);           % Bryson: 1/(max acceptable value)^2
R  = 1 / P.lqr_umax^2;
[K, ~] = dlqr_riccati(Aa, Ba, Q, R);

% ---- Steady-state Kalman filter (process noise = unmodelled accel) ----
qa  = 0.25^2;                          % [m/s^2]^2 acceleration noise PSD-ish
Qn  = diag([1e-9, qa*Ts^2, deg2rad(0.5)^2]);
Rn  = P.noise_std^2 + P.quant^2/12;    % quantisation variance = q^2/12
Pk  = eye(3)*1e-4;
for k = 1:5000                         % iterate the filter Riccati equation
    Pp = Ad*Pk*Ad' + Qn;               % predict covariance
    L  = Pp*Cd' / (Cd*Pp*Cd' + Rn);    % Kalman gain  (A/B = A*inv(B))
    Pn = (eye(3) - L*Cd)*Pp;           % update covariance
    if max(abs(Pn(:)-Pk(:))) < 1e-16, Pk = Pn; break; end
    Pk = Pn;
end

C.Ad = Ad; C.Bd = Bd; C.Cd = Cd; C.K = K; C.L = L; C.Ts = Ts;
C.Aa = Aa; C.Ba = Ba;
C.cl_poles = eig(Aa - Ba*K);           % closed-loop poles (inside unit circle)
C.W = workspace_3rrs(P);               % reachable tilt LUT for the governor
end
