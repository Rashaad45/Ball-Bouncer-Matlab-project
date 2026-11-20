function S = scenarios_3rrs()
% SCENARIOS_3RRS  Benchmark scenarios (all positions in metres, plate frame).
%
%   Each reference function returns position, velocity and acceleration,
%   because the feed-forward needs r_ddot.  'deal' copies its inputs to
%   several outputs:  [a,b,c] = deal(1,2,3).

% S1: disturbance rejection -- hold centre, ball is flicked at t = 2 s
S{1}.name = 'Disturbance rejection';
S{1}.T = 5;  S{1}.x0 = [0.05; -0.03; 0; 0];
S{1}.kick_t = 2;  S{1}.kick_v = [0.16; -0.12];
S{1}.ref = @(t) deal([0;0], [0;0], [0;0]);

% S2: circle, radius 50 mm, period 2.5 s (centripetal accel 0.32 m/s^2)
R = 0.05; w = 2*pi/2.5;
S{2}.name = 'Circle tracking';
S{2}.T = 6;  S{2}.x0 = [R; 0; 0; 0];
S{2}.ref = @(t) deal(R*[cos(w*t); sin(w*t)], R*w*[-sin(w*t); cos(w*t)], ...
                     -R*w^2*[cos(w*t); sin(w*t)]);     % '...' continues a line

% S3: figure-8 (Lissajous 1:2), 60 x 40 mm, period 4 s
A = 0.06; B = 0.04; w8 = 2*pi/4;
S{3}.name = 'Figure-8 tracking';
S{3}.T = 8;  S{3}.x0 = [0; 0; A*w8; 0];
S{3}.ref = @(t) deal([A*sin(w8*t); B*sin(2*w8*t)], ...
                     [A*w8*cos(w8*t); 2*B*w8*cos(2*w8*t)], ...
                     [-A*w8^2*sin(w8*t); -4*B*w8^2*sin(2*w8*t)]);

% S4: aggressive point-to-point, smooth quintic from (-60,0) to (+60,0) mm in 0.6 s
S{4}.name = 'Fast point-to-point';
S{4}.T = 4;  S{4}.x0 = [-0.06; 0; 0; 0];
S{4}.ref = @(t) quintic_ref(t, 0.5, 0.6, [-0.06;0], [0.06;0]);
end

function [r, rd, rdd] = quintic_ref(t, t0, Tm, p0, p1)
% Minimum-jerk (quintic) blend: s = 10 tau^3 - 15 tau^4 + 6 tau^5
tau = min(max((t - t0)/Tm, 0), 1);
s   = 10*tau^3 - 15*tau^4 + 6*tau^5;
sd  = (30*tau^2 - 60*tau^3 + 30*tau^4) / Tm;
sdd = (60*tau - 180*tau^2 + 120*tau^3) / Tm^2;
if t < t0 || t > t0 + Tm, sd = 0; sdd = 0; end
r = p0 + s*(p1 - p0);  rd = sd*(p1 - p0);  rdd = sdd*(p1 - p0);
end
