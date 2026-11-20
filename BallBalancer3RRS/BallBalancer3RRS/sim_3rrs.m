function out = sim_3rrs(P, C, cfg, scn, Pplant)
% SIM_3RRS  Closed-loop simulation: camera -> estimator -> LQR -> governor
%           -> IK -> servos -> EXACT forward kinematics -> nonlinear ball.
%
%   out = sim_3rrs(P, C, cfg, scn, Pplant)
%     P      : nominal parameters (what the controller BELIEVES)
%     C      : controller from design_controller_3rrs(P)
%     cfg    : configuration struct
%                .name       label for plots
%                .ik         'exact' | 'linear' | 'fixedcentre'
%                .delay_comp true -> delay-compensated Kalman predictor
%                .governor   true -> workspace + servo-rate governor
%                .ff         true -> flatness feed-forward  tilt = r_ddot/(kb g)
%     scn    : scenario struct
%                .T     duration [s]
%                .x0    initial ball [x; y; vx; vy] (plate frame)
%                .ref   @(t) -> [r(2x1), r_dot(2x1), r_ddot(2x1)]
%                .kick_t, .kick_v  velocity disturbance time and vector
%     Pplant : (optional) TRUE plant parameters -> lets Monte-Carlo runs
%              put geometry/servo/latency errors in the plant only.
%
%   Two time scales:
%     plant integration  dt  = Ts/10  (600 Hz): servos, kinematics, ball
%     controller/camera  Ts  = 1/60 s (60 Hz):  sensing, estimation, control

if nargin < 5 || isempty(Pplant), Pplant = P; end
Ts = P.Ts_cam;  dt = Pplant.dt;
nsub   = round(Ts/dt);                 % plant steps per camera frame
nframe = round(scn.T/Ts);
d      = Pplant.delay_frames;          % TRUE latency (controller assumes P.delay_frames)
dctl   = P.delay_frames;

% ---------------- plant state ----------------
ball  = scn.x0(:);                     % [x; y; vx; vy] in plate frame
theta = zeros(3,1);                    % actual servo angles
q     = [0;0;Pplant.h0];               % actual platform pose
[~,~,inf0] = ik_3rrs(q, Pplant);
T_prev = inf0.T; T_prev2 = inf0.T;     % for finite-difference accelerations
q_prevstep = q;
Jfk   = jacobian_3rrs(q, Pplant);
offs  = zeros(3,1);
if isfield(Pplant,'servo_offset'), offs = Pplant.servo_offset(:); end

% ---------------- controller state ----------------
xh   = [scn.x0(1) 0 0; scn.x0(2) 0 0]';   % 3x2: per-axis estimate [pos;vel;tilt]
xh_del = xh;                           % estimate at the (delayed) measurement time
z    = zeros(1,2);                     % integrators
uhist = zeros(nframe+1, 2);            % virtual-tilt commands per frame
meas  = nan(nframe+1, 2);              % camera samples (true time stamp)
q_cmd_prev = [0;0;P.h0];  th_cmd_prev = zeros(3,1);
theta_cmd  = zeros(3,1);
kicked = false;  failed = false;  t_fail = NaN;

% ---------------- logs (one row per frame) ----------------
L.t = (0:nframe)'*Ts;
Z = @(n) nan(nframe+1, n);             % helper that makes a NaN array
L.ball = Z(2); L.est = Z(2); L.ref = Z(2); L.qcmd = Z(2); L.q = Z(3);
L.th = Z(3); L.thcmd = Z(3); L.kinerr = Z(1); L.par = Z(3); L.gov = Z(2); L.thdot = Z(3);

for k = 0:nframe
    t  = k*Ts;
    kk = k+1;                          % MATLAB arrays are 1-indexed

    % ======== 1. CAMERA: sample now, deliver d frames later ========
    yk = ball(1:2) + Pplant.noise_std*randn(2,1);       % randn = N(0,1)
    meas(kk,:) = (round(yk/Pplant.quant)*Pplant.quant)';  % pixel quantisation
    [R3, Rd3, Rdd3] = scn.ref(t);

    % ======== 2. ESTIMATOR (per axis, same model for x and y) ========
    for a = 1:2
        if cfg.delay_comp
            m = k - dctl;              % frame index the new sample belongs to
            if m >= 0 && k - d >= 0
                if m >= 1                  % predict m-1 -> m with the command of frame m-1
                    xh_del(:,a) = C.Ad*xh_del(:,a) + C.Bd*uhist(m,a);
                end
                ymeas = meas(k-d+1, a);
                xh_del(:,a) = xh_del(:,a) + C.L*(ymeas - C.Cd*xh_del(:,a));
            end
            xp = xh_del(:,a);
            for j = max(m,0):k-1       % roll forward through the latency
                xp = C.Ad*xp + C.Bd*uhist(j+1,a);
            end
            xh(:,a) = xp;
        else                           % naive: treat the delayed sample as current
            if k >= 1, xh(:,a) = C.Ad*xh(:,a) + C.Bd*uhist(k,a); end
            if k - d >= 0
                ymeas = meas(k-d+1, a);
                xh(:,a) = xh(:,a) + C.L*(ymeas - C.Cd*xh(:,a));
            end
        end
    end

    % ======== 3. LQR + integral + feed-forward (virtual tilt) ========
    u = zeros(1,2);
    for a = 1:2
        tff = 0;
        if cfg.ff, tff = (Rdd3(a) + P.cv*Rd3(a)) / (P.kb*P.g); end
        e   = [xh(1,a)-R3(a); xh(2,a)-Rd3(a); xh(3,a)-tff; z(a)];
        u(a) = tff - C.K*e;
    end

    % virtual tilt -> plate angles:  x uses beta,  y uses -alpha
    q_des = [-u(2); u(1); P.h0];
    [q_cmd, th_exact, fl] = governor_3rrs(q_des, q_cmd_prev, th_cmd_prev, P, C.W, cfg.governor);
    switch cfg.ik                      % switch/case selects by string value
        case 'exact',       theta_cmd = th_exact;
        case 'linear',      theta_cmd = ik_linear_3rrs(q_cmd, P);
        case 'fixedcentre', theta_cmd = ik_fixedcentre_3rrs(q_cmd, P);
    end
    u_applied = [q_cmd(2), -q_cmd(1)];
    sat = norm(u_applied - u) > 1e-9;
    for a = 1:2                        % integrate only when not saturated (anti-windup)
        if ~sat, z(a) = z(a) - Ts*(xh(1,a) - R3(a)); end
    end
    uhist(kk,:) = u_applied;
    % static kinematic error of the chosen IK: where would the plate REALLY
    % go if the servos reached theta_cmd exactly?
    qk = fk_3rrs(theta_cmd, P, q_cmd, Jfk);
    L.kinerr(kk,1) = norm(qk(1:2) - q_cmd(1:2));
    q_cmd_prev = q_cmd;  th_cmd_prev = th_exact;

    % log
    L.ball(kk,:) = ball(1:2)'; L.est(kk,:) = xh(1,:); L.ref(kk,:) = R3(:)';
    L.qcmd(kk,:) = q_cmd(1:2)'; L.q(kk,:) = q'; L.th(kk,:) = theta';
    L.thcmd(kk,:) = theta_cmd'; L.gov(kk,:) = [fl.ws fl.rate];
    if k == nframe || failed, break; end

    % ======== 4. PLANT: integrate nsub steps at dt ========
    for s = 1:nsub
        ts = t + (s-1)*dt;
        % --- servo: first-order lag + slew-rate limit + travel limit ---
        thd   = (theta_cmd - theta) / Pplant.tau_s;
        thd   = max(-Pplant.omega_max, min(Pplant.omega_max, thd));
        theta = min(Pplant.theta_max, max(Pplant.theta_min, theta + thd*dt));
        L.thdot(kk,:) = max([abs(thd)'; abs(L.thdot(kk,:))], [], 1);

        % --- exact forward kinematics of the TRUE geometry ---
        [q, okfk, ~, Jfk] = fk_3rrs(theta + offs, Pplant, q, Jfk);
        [~,~,inf_] = ik_3rrs(q, Pplant);
        R = inf_.R;  T = inf_.T;

        % plate angular rates and centre acceleration (finite differences)
        qd   = (q - q_prevstep)/dt;   q_prevstep = q;
        aT   = (T - 2*T_prev + T_prev2)/dt^2;
        T_prev2 = T_prev;  T_prev = T;

        % --- nonlinear ball dynamics in the plate frame ---
        g_eff = R' * ([0;0;-Pplant.g] - aT);  % gravity + parasitic accel, plate axes
        wx = qd(1); wy = qd(2);
        x = ball(1); y = ball(2); vx = ball(3); vy = ball(4);
        ax = Pplant.kb*( g_eff(1) + x*wy^2 + y*wx*wy ) - Pplant.cv*vx;
        ay = Pplant.kb*( g_eff(2) + y*wx^2 + x*wx*wy ) - Pplant.cv*vy;
        vx = vx + ax*dt;  vy = vy + ay*dt;      % semi-implicit Euler
        x  = x + vx*dt;   y  = y + vy*dt;
        ball = [x; y; vx; vy];

        if isfield(scn,'kick_t') && ~kicked && ts >= scn.kick_t
            ball(3:4) = ball(3:4) + scn.kick_v(:);  kicked = true;
        end
        if norm(ball(1:2)) > Pplant.plate_radius || ~okfk
            failed = true; t_fail = ts; break;
        end
    end
    [~,~,ip] = ik_3rrs(q, Pplant);
    L.par(kk,:) = [ip.x ip.y ip.gamma];
end

% ---------------- metrics ----------------
valid = ~isnan(L.ball(:,1));
e  = L.ball - L.ref;  en = sqrt(sum(e.^2,2));
tw = valid & L.t >= min(scn.T*0.25, 2);           % ignore initial transient
out.log = L;  out.cfg = cfg;  out.failed = failed;  out.t_fail = t_fail;
out.rms_mm  = 1e3*sqrt(mean(en(tw).^2));
out.max_mm  = 1e3*max(en(tw));
out.peak_servo_dps = rad2deg(max(L.thdot(:)));
out.rate_sat_pct   = 100*mean(any(L.thdot(valid,:) >= 0.999*Pplant.omega_max, 2));
out.gov_pct        = 100*mean(any(L.gov(valid,:),2));
out.par_max_mm     = 1e3*max(sqrt(sum(L.par(valid,1:2).^2,2)));
out.kin_err_deg    = rad2deg(max(L.kinerr(valid)));      % worst static IK error
end
