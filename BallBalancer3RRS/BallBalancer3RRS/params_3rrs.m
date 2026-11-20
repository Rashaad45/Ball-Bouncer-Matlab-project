function P = params_3rrs()
% PARAMS_3RRS  All physical, sensor and simulation parameters in ONE place.
%
%   P = params_3rrs() returns a struct P.  A struct is MATLAB's "record"
%   type: you create fields with dot-syntax, e.g. P.rb = 0.07, and read
%   them back the same way.  Every other file receives P as an argument,
%   so changing a number here changes the whole project consistently.
%
%   Units: SI everywhere (m, rad, s, kg).  deg2rad() converts degrees.
%
%   Mechanism: 3-RRS parallel manipulator
%     R (active revolute, servo shaft)  ->  R (passive revolute elbow)
%     ->  S (spherical/ball joint under the plate).
%   The two revolute axes of each leg are parallel and tangential, so
%   each leg moves in a vertical plane that contains the z-axis.  That
%   planar constraint is what makes the platform 3-DOF (roll, pitch,
%   heave) and what creates the "parasitic" x, y and yaw motion.

%% ---------------- Geometry (edit to match your CAD) ----------------
P.rb  = 0.070;                 % [m] base radius: centre -> servo shaft axis
P.rp  = 0.090;                 % [m] platform radius: centre -> ball joint
P.L1  = 0.050;                 % [m] servo horn (crank) length
P.L2  = 0.080;                 % [m] coupler link length (elbow -> ball joint)
P.phi = deg2rad([90 210 330]); % [rad] azimuth of legs 1..3 (120 deg apart)

% Home height: plate flat and servo horns horizontal (theta = 0).
% Pythagoras in the leg plane: horizontal gap (rp - rb - L1), link L2.
P.h0  = sqrt(P.L2^2 - (P.rp - P.rb - P.L1)^2);   % [m] ~0.0742 m

% Servo travel limits (typical hobby servo: +/-90 deg, keep margin)
P.theta_min = deg2rad(-60);    % [rad]
P.theta_max = deg2rad( 60);    % [rad]

%% ---------------- Ball ----------------
P.g      = 9.81;               % [m/s^2] gravity
P.mb     = 0.030;              % [kg] ball mass (steel ~20 mm ball ~33 g)
P.rball  = 0.010;              % [m] ball radius
P.kb     = 5/7;                % [-] m/(m + I/r^2) for a SOLID sphere (I=2/5 m r^2)
P.cv     = 0.15;               % [1/s] viscous rolling-resistance coefficient
P.plate_radius = 0.100;        % [m] usable plate radius; ball "falls off" beyond

%% ---------------- Servo actuator model ----------------
% First-order lag + slew-rate limit: a good fit for MG996R/DS3218 class.
P.tau_s     = 0.030;           % [s] position-loop time constant
P.omega_max = deg2rad(400);    % [rad/s] max shaft speed (~0.15 s / 60 deg)

%% ---------------- Vision sensor model ----------------
P.Ts_cam       = 1/60;         % [s] camera frame period (60 FPS)
P.delay_frames = 2;            % [frames] processing + transport latency
P.noise_std    = 0.0015;       % [m] 1-sigma centroid noise
P.quant        = 0.0005;       % [m] pixel quantisation (~0.5 mm / px)

%% ---------------- Controller tuning (Bryson's rule) ----------------
% "Largest acceptable" value of each state / input.  Smaller = tighter.
%            [pos err  vel err  tilt     integral]
P.lqr_qmax = [0.015,   0.30,    deg2rad(30), 0.006];
P.lqr_umax = deg2rad(8);       % [rad] typical tilt command

%% ---------------- Simulation ----------------
P.dt   = P.Ts_cam/10;          % [s] plant integration step (600 Hz)
P.tilt_max_sw = deg2rad(15);   % [rad] software tilt ceiling (all configs)
P.gov_margin  = 0.90;          % [-] governor uses 90% of reachable tilt
end
