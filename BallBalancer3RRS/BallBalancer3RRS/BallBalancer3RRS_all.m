function BallBalancer3RRS_all(varargin)
%% BALLBALANCER3RRS_ALL  Single-file version of the whole 3-RRS ball-balancer project.
%
%  All original .m files (kinematics, controller, simulation, studies,
%  animation, tests) are merged here as local functions.  No toolboxes needed.
%
%  Usage:
%    BallBalancer3RRS_all                      % run everything (= old run_all)
%    BallBalancer3RRS_all('test')              % self-tests only
%    BallBalancer3RRS_all('kinematics','closedloop')   % any subset, in order
%
%  Stages:  'test' 'kinematics' 'closedloop' 'latency' 'montecarlo' 'animate'
%  Outputs go to ./results (figures, CSVs, logs, figure8_proposed.mp4).

all_stages = {'test','kinematics','closedloop','latency','montecarlo','animate'};
if nargin == 0, stages = all_stages; else, stages = lower(varargin); end
bad = setdiff(stages, all_stages);
if ~isempty(bad), error('Unknown stage(s): %s', strjoin(bad, ', ')); end

old_vis = get(0,'defaultfigurevisible');
for k = 1:numel(stages)
    fprintf('\n===== %s =====\n', upper(stages{k}));
    switch stages{k}
        case 'test',       test_3rrs();                 % stops here if the maths is broken
        case 'kinematics', run_kinematics_analysis();   % workspace, parasitic motion, IK-error maps
        case 'closedloop', run_closed_loop();           % 4 scenarios x 4 configurations
        case 'latency',    run_latency_study();         % ablation vs camera latency
        case 'montecarlo', run_monte_carlo();           % robustness to tolerances / servo / camera
        case 'animate'                                  % animation of the proposed controller, figure-8
            if ~exist('results/closed_loop_runs.mat','file'), run_closed_loop(); end
            S = load('results/closed_loop_runs.mat');
            P = params_3rrs();
            animate_3rrs(S.runs{3,4}, P, struct('file', 'results/figure8_proposed.mp4'));
    end
end
set(0,'defaultfigurevisible',old_vis);
end


%% ======================================================================
%% test_3rrs (was script test_3rrs.m)
%% ======================================================================

function test_3rrs()
    %% TEST_3RRS  Self-checks.  Run this first after any change to the maths.
    %  Every assert() stops with an error message if its condition is false.

    P = params_3rrs();
    fprintf('Running 3-RRS self-tests...\n');

    % T1: home pose gives zero servo angles
    [th, ok] = ik_3rrs([0;0;P.h0], P);
    assert(ok && max(abs(th)) < 1e-12, 'T1 home pose failed');

    % T2: parasitic constraints satisfied, closed-form yaw exact
    for m = deg2rad([3 10 20])
        for ps = linspace(0, 2*pi, 7)
            [x,y,g,R,res] = parasitic_3rrs(m*cos(ps), m*sin(ps), P);
            assert(res < 1e-12, 'T2 parasitic residual too large');
        end
    end

    % T3: every coupler has exactly length L2, ball joints are rigid triangle
    q = [deg2rad(12); deg2rad(-7); P.h0 + 0.004];
    [th, ok, I] = ik_3rrs(q, P);
    assert(ok && max(abs(I.link_err)) < 1e-12, 'T3 link length violated');
    sides = [norm(I.P(:,1)-I.P(:,2)), norm(I.P(:,2)-I.P(:,3)), norm(I.P(:,3)-I.P(:,1))];
    assert(max(abs(sides - sqrt(3)*P.rp)) < 1e-12, 'T3 platform not rigid');

    % T4: FK inverts IK
    qf = fk_3rrs(th, P);
    assert(norm(qf - q) < 1e-9, 'T4 FK(IK(q)) ~= q');

    % T5: Jacobian matches a finite-difference velocity check
    J  = jacobian_3rrs(q, P);  dq = [1e-5; -2e-5; 1e-6];
    assert(norm(ik_3rrs(q+dq,P) - th - J*dq) < 1e-8, 'T5 Jacobian mismatch');

    % T6: closed-loop poles stable, KF stable
    C = design_controller_3rrs(P);
    assert(all(abs(C.cl_poles) < 1), 'T6 LQR closed loop unstable');
    assert(all(abs(eig(C.Ad - C.L*C.Cd*C.Ad)) < 1), 'T6 Kalman filter unstable');

    % T7: short closed-loop run keeps the ball on the plate
    S = scenarios_3rrs(); CF = configs_3rrs(); s = S{1}; s.T = 1.5;
    randn('seed', 1); o = sim_3rrs(P, C, CF.P, s);
    assert(~o.failed, 'T7 ball fell off in a 1.5 s run');

    fprintf('All 7 tests passed.\n');
end

%% ======================================================================
%% run_kinematics_analysis (was script run_kinematics_analysis.m)
%% ======================================================================

function run_kinematics_analysis()
    %% RUN_KINEMATICS_ANALYSIS  Kinematic study of the 3-RRS ball balancer.
    %
    %  Produces (in ./results):
    %   K1_workspace_polar.png     reachable tilt vs direction (governor LUT)
    %   K2_parasitic.png           parasitic x/y shift and yaw vs tilt
    %   K3_ik_error_maps.png       static tilt error of the two shortcut IKs
    %   K4_dexterity.png           Jacobian condition number over the tilt disc
    %   kinematics_summary.csv     key numbers for the report
    %
    %  A "script" (no 'function' line) runs top to bottom in the base workspace.

    P = params_3rrs();
    if ~exist('results','dir'), mkdir('results'); end   % make output folder once
    set(0,'defaultfigurevisible','off');                % draw off-screen (faster)

    %% 1. Round-trip accuracy  IK -> FK  over a random set of poses
    rng_seed = 1; rand('seed', rng_seed); randn('seed', rng_seed);
    N = 200; err = zeros(N,1); its = zeros(N,1);
    for n = 1:N
        m   = deg2rad(20)*sqrt(rand);  ps = 2*pi*rand;       % uniform over a disc
        q   = [m*cos(ps); m*sin(ps); P.h0 + 0.01*(rand-0.5)];
        th  = ik_3rrs(q, P);
        [qf, ok, its(n)] = fk_3rrs(th, P);                    % cold start from home
        err(n) = norm(qf - q);
    end
    fprintf('IK->FK round trip: max error %.2e, mean Newton iterations %.1f\n', max(err), mean(its));

    %% 2. Workspace envelope (used online by the governor)
    W = workspace_3rrs(P);
    fprintf('Reachable tilt at h0: min %.1f deg, max %.1f deg\n', rad2deg(W.min_tilt), rad2deg(max(W.tilt_max)));
    figure; polar(W.psi, rad2deg(W.tilt_max), 'b-'); hold on;
    polar(W.psi, rad2deg(P.gov_margin*W.tilt_max), 'r--');
    title('Reachable tilt magnitude [deg] vs tilt direction (blue) and governor limit (red)');
    print('-dpng','-r120','results/K1_workspace_polar.png');

    %% 3. Parasitic motion vs tilt magnitude (worst direction)
    mags = deg2rad(0:1:25); psis = linspace(0, 2*pi, 73);
    par_xy = zeros(size(mags)); par_yaw = zeros(size(mags));
    for i = 1:numel(mags)
        for j = 1:numel(psis)
            [x, y, g] = parasitic_3rrs(mags(i)*cos(psis(j)), mags(i)*sin(psis(j)), P);
            par_xy(i)  = max(par_xy(i), hypot(x,y));
            par_yaw(i) = max(par_yaw(i), abs(g));
        end
    end
    figure;
    subplot(1,2,1); plot(rad2deg(mags), 1e3*par_xy, 'LineWidth', 1.5); grid on;
    xlabel('tilt [deg]'); ylabel('max centre shift [mm]'); title('Parasitic translation');
    subplot(1,2,2); plot(rad2deg(mags), rad2deg(par_yaw), 'LineWidth', 1.5); grid on;
    xlabel('tilt [deg]'); ylabel('max yaw [deg]'); title('Parasitic yaw');
    print('-dpng','-r120','results/K2_parasitic.png');

    %% 4. Static tilt error of shortcut IKs (FK of their servo angles)
    ng = 31; a = linspace(-deg2rad(20), deg2rad(20), ng);
    [AA, BB] = meshgrid(a, a);                % meshgrid: AA(i,j)=a(j) (x), BB(i,j)=a(i) (y)
    % imagesc(x, y, M) paints matrix M as a colour image; NaN = outside the disc
    E_lin = nan(ng); E_fc = nan(ng); CN = nan(ng);
    for i = 1:ng
        for j = 1:ng
            q = [AA(i,j); BB(i,j); P.h0];
            if norm(q(1:2)) > deg2rad(20), continue; end    % only the 20 deg disc
            [th_ex, ok] = ik_3rrs(q, P);
            if ~ok, continue; end
            ql = fk_3rrs(ik_linear_3rrs(q, P), P, q);
            qc = fk_3rrs(ik_fixedcentre_3rrs(q, P), P, q);
            E_lin(i,j) = rad2deg(norm(ql(1:2) - q(1:2)));
            E_fc(i,j)  = rad2deg(norm(qc(1:2) - q(1:2)));
            CN(i,j)    = cond(jacobian_3rrs(q, P));
        end
    end
    figure;
    subplot(1,2,1); imagesc(rad2deg(a), rad2deg(a), E_lin); set(gca,'YDir','normal'); colorbar; axis equal tight;
    xlabel('\alpha roll [deg]'); ylabel('\beta pitch [deg]'); title('Linear IK tilt error [deg]');
    subplot(1,2,2); imagesc(rad2deg(a), rad2deg(a), E_fc); set(gca,'YDir','normal'); colorbar; axis equal tight;
    xlabel('\alpha roll [deg]'); ylabel('\beta pitch [deg]'); title('Fixed-centre IK tilt error [deg]');
    print('-dpng','-r120','results/K3_ik_error_maps.png');

    figure; imagesc(rad2deg(a), rad2deg(a), CN); set(gca,'YDir','normal'); colorbar; axis equal tight;
    xlabel('\alpha [deg]'); ylabel('\beta [deg]'); title('cond(J): dexterity (lower = better)');
    print('-dpng','-r120','results/K4_dexterity.png');

    %% 5. Rigidity violation of the fixed-centre shortcut at 15 deg
    worst_side = 0;
    for j = 1:numel(psis)
        q = [deg2rad(15)*cos(psis(j)); deg2rad(15)*sin(psis(j)); P.h0];
        [~, sl] = ik_fixedcentre_3rrs(q, P);
        worst_side = max(worst_side, max(abs(sl/(sqrt(3)*P.rp) - 1)));
    end

    %% 6. Summary table
    i15 = find(abs(rad2deg(mags) - 15) < 1e-9);
    S = {'roundtrip_max_err', max(err);
         'tilt_reach_min_deg', rad2deg(W.min_tilt);
         'tilt_reach_max_deg', rad2deg(max(W.tilt_max));
         'parasitic_shift_mm_at15', 1e3*par_xy(i15);
         'parasitic_yaw_deg_at15', rad2deg(par_yaw(i15));
         'linearIK_max_err_deg_20disc', max(E_lin(:));
         'fixedcentreIK_max_err_deg_20disc', max(E_fc(:));
         'fixedcentre_side_violation_pct_at15', 100*worst_side;
         'condJ_min', min(CN(:)); 'condJ_max', max(CN(:))};
    fid = fopen('results/kinematics_summary.csv','w');
    for i = 1:size(S,1)
        fprintf(fid, '%s,%.6g\n', S{i,1}, S{i,2});
        fprintf('%-40s %10.4g\n', S{i,1}, S{i,2});
    end
    fclose(fid);
end

%% ======================================================================
%% run_closed_loop (was script run_closed_loop.m)
%% ======================================================================

function run_closed_loop()
    %% RUN_CLOSED_LOOP  Benchmark 4 scenarios x 4 configurations.
    %
    %  Outputs in ./results:
    %   C1_<scenario>.png        XY path + error vs time for every config
    %   C2_servos_proposed.png   servo angles of the proposed controller (circle)
    %   closed_loop_results.csv  metrics table
    %   closed_loop_runs.mat     all logs (load with: load results/closed_loop_runs.mat)

    P  = params_3rrs();
    Cn = design_controller_3rrs(P);
    S  = scenarios_3rrs();
    CF = configs_3rrs();
    cfgs = {CF.B1, CF.B2, CF.A1, CF.P};          % cell array: { } holds mixed items
    if ~exist('results','dir'), mkdir('results'); end
    set(0,'defaultfigurevisible','off');

    fprintf('Closed-loop poles |z| = %s\n', mat2str(abs(Cn.cl_poles'), 3));
    runs = cell(numel(S), numel(cfgs));
    fid = fopen('results/closed_loop_results.csv','w');
    fprintf(fid, 'scenario,config,rms_mm,max_mm,failed,kin_err_deg,peak_servo_dps,rate_sat_pct,gov_pct,parasitic_max_mm\n');

    for s = 1:numel(S)
        for c = 1:numel(cfgs)
            randn('seed', 42);                   % same noise for every config -> fair
            o = sim_3rrs(P, Cn, cfgs{c}, S{s});
            runs{s,c} = o;
            fprintf('%-22s | %-20s | RMS %6.2f mm | max %6.1f mm | fail %d | IKerr %.3f deg\n', ...
                S{s}.name, cfgs{c}.name, o.rms_mm, o.max_mm, o.failed, o.kin_err_deg);
            fprintf(fid, '%s,%s,%.3f,%.3f,%d,%.4f,%.1f,%.2f,%.2f,%.3f\n', S{s}.name, cfgs{c}.name, ...
                o.rms_mm, o.max_mm, o.failed, o.kin_err_deg, o.peak_servo_dps, o.rate_sat_pct, o.gov_pct, o.par_max_mm);
        end

        % ---- figure per scenario ----
        figure('Position',[0 0 1100 420]);
        subplot(1,2,1); hold on; grid on;
        L0 = runs{s,1}.log;
        plot(1e3*L0.ref(:,1), 1e3*L0.ref(:,2), 'k--', 'LineWidth', 1.5);
        for c = 1:numel(cfgs)
            L = runs{s,c}.log;
            plot(1e3*L.ball(:,1), 1e3*L.ball(:,2), 'LineWidth', 1.1);
        end
        th = linspace(0,2*pi,200);
        plot(1e3*P.plate_radius*cos(th), 1e3*P.plate_radius*sin(th), 'Color', [.6 .6 .6]);
        axis([-110 110 -110 110]); axis square;   % same scale on both axes
        xlabel('x [mm]'); ylabel('y [mm]'); title(S{s}.name);
        legend(['reference', cellfun(@(z) z.name, cfgs, 'UniformOutput', false)], 'Location','southoutside');

        subplot(1,2,2); hold on; grid on;
        for c = 1:numel(cfgs)
            L = runs{s,c}.log;
            plot(L.t, 1e3*sqrt(sum((L.ball - L.ref).^2, 2)), 'LineWidth', 1.1);
        end
        xlabel('t [s]'); ylabel('|position error| [mm]'); title('Tracking error');
        fname = sprintf('results/C1_S%d.png', s);
        print('-dpng','-r110', fname);
    end
    fclose(fid);

    % servo angles of the proposed controller in the circle test
    L = runs{2, end}.log;
    figure; plot(L.t, rad2deg(L.th), 'LineWidth', 1.2); grid on;
    xlabel('t [s]'); ylabel('servo angle [deg]'); legend('servo 1','servo 2','servo 3');
    title('Proposed controller: servo angles, circle tracking');
    print('-dpng','-r110','results/C2_servos_proposed.png');

    save('results/closed_loop_runs.mat', 'runs', 'P');

    % plain-text logs (for Python/Excel post-processing or a report)
    if ~exist('results/logs','dir'), mkdir('results/logs'); end
    for s = 1:numel(S), for c = 1:numel(cfgs)
        L = runs{s,c}.log;
        M = [L.t L.ball L.ref L.est rad2deg(L.qcmd) rad2deg(L.q(:,1:2)) 1e3*L.q(:,3) ...
             rad2deg(L.th) rad2deg(L.thcmd) 1e3*L.par(:,1:2) rad2deg(L.par(:,3))];
        fn = sprintf('results/logs/S%d_C%d.csv', s, c);
        fid = fopen(fn, 'w');
        fprintf(fid, 't,x,y,rx,ry,xhat,yhat,alpha_cmd,beta_cmd,alpha,beta,h_mm,th1,th2,th3,th1c,th2c,th3c,par_x_mm,par_y_mm,yaw_deg\n');
        fclose(fid);
        dlmwrite(fn, M, '-append', 'precision', '%.6g');
    end, end
end

%% ======================================================================
%% run_latency_study (was script run_latency_study.m)
%% ======================================================================

function run_latency_study()
    %% RUN_LATENCY_STUDY  Ablation: how much does each layer buy as camera
    %  latency grows?  (The question every vision-based balancer runs into.)
    %
    %  Output: results/L1_latency_ablation.png, results/latency_ablation.csv

    P0 = params_3rrs();  S = scenarios_3rrs();  CF = configs_3rrs();
    cfgs = {CF.A1, CF.A2, CF.A3, CF.P};
    delays = [1 2 3 4 5 6];                          % frames (16.7 ms each)
    scn = S{1};                                      % disturbance rejection
    scn2 = S{2};                                     % circle
    RMS  = nan(numel(delays), numel(cfgs), 2);   % RMS error, NaN when the ball fell off
    FAIL = false(numel(delays), numel(cfgs), 2);
    set(0,'defaultfigurevisible','off');
    if ~exist('results','dir'), mkdir('results'); end

    for i = 1:numel(delays)
        P = P0;  P.delay_frames = delays(i);         % controller knows the latency
        C = design_controller_3rrs(P);
        for c = 1:numel(cfgs)
            randn('seed', 7); o1 = sim_3rrs(P, C, cfgs{c}, scn);
            randn('seed', 7); o2 = sim_3rrs(P, C, cfgs{c}, scn2);
            FAIL(i,c,:) = [o1.failed, o2.failed];
            if ~o1.failed, RMS(i,c,1) = o1.rms_mm; end   % failed run -> stays NaN (gap in plot)
            if ~o2.failed, RMS(i,c,2) = o2.rms_mm; end
            fprintf('delay %d | %-22s | dist %7.2f (fail %d) | circle %7.2f (fail %d) mm\n', ...
                delays(i), cfgs{c}.name, RMS(i,c,1), o1.failed, RMS(i,c,2), o2.failed);
        end
    end

    ms = 1e3*delays*P0.Ts_cam;
    figure('Position',[0 0 1000 400]);
    names = cellfun(@(z) z.name, cfgs, 'UniformOutput', false);
    for k = 1:2
        subplot(1,2,k);
        semilogy(ms, RMS(:,:,k), '-o', 'LineWidth', 1.3); grid on; hold on;
        xlabel('camera latency [ms]'); ylabel('RMS error [mm] (log)');
        if k == 1, title('Disturbance rejection'); else, title('Circle tracking'); end
        legend(names, 'Location', 'northwest');
        % mark ball-fell-off runs with a red x at the top of the axis
        yl = ylim;  [ii, cc] = find(FAIL(:,:,k));
        if ~isempty(ii), plot(ms(ii), yl(2)*ones(size(ii)), 'rx', 'MarkerSize', 10, 'LineWidth', 2); end
    end
    print('-dpng','-r110','results/L1_latency_ablation.png');

    fid = fopen('results/latency_ablation.csv','w');
    fprintf(fid, 'delay_ms,config,dist_rms_mm,dist_failed,circle_rms_mm,circle_failed\n');
    for i = 1:numel(delays), for c = 1:numel(cfgs)
        fprintf(fid, '%.1f,%s,%.3f,%d,%.3f,%d\n', ms(i), cfgs{c}.name, RMS(i,c,1), FAIL(i,c,1), RMS(i,c,2), FAIL(i,c,2));
    end, end
    fclose(fid);
end

%% ======================================================================
%% run_monte_carlo (was script run_monte_carlo.m)
%% ======================================================================

function run_monte_carlo()
    %% RUN_MONTE_CARLO  Robustness: the controller is designed on NOMINAL
    %  parameters, but every trial gets a different TRUE plant:
    %    link lengths / radii  +/- 0.5 mm      (print + assembly tolerance)
    %    servo zero offsets    +/- 1.5 deg     (horn spline calibration)
    %    servo tau             +/- 30 %, speed  -20..+10 %
    %    camera latency        1 .. 4 frames  (controller always assumes 2)
    %    camera noise          1 .. 3 mm
    %  Output: results/M1_monte_carlo.png, results/monte_carlo.csv

    P = params_3rrs();  C = design_controller_3rrs(P);
    S = scenarios_3rrs();  CF = configs_3rrs();
    cfgs = {CF.A1, CF.P};
    Ntr = 40;                                % increase to 200 for a paper-grade study
    rand('seed', 11); randn('seed', 11);
    scn = S{2};  scn.T = 5;                  % circle tracking
    res = nan(Ntr, numel(cfgs));  fail = false(Ntr, numel(cfgs));
    U = @(a,b) a + (b-a)*rand;               % uniform random number in [a,b]
    set(0,'defaultfigurevisible','off');
    if ~exist('results','dir'), mkdir('results'); end

    for n = 1:Ntr
        Pp = P;
        Pp.L1 = P.L1 + U(-5e-4, 5e-4);  Pp.L2 = P.L2 + U(-5e-4, 5e-4);
        Pp.rb = P.rb + U(-5e-4, 5e-4);  Pp.rp = P.rp + U(-5e-4, 5e-4);
        Pp.h0 = sqrt(Pp.L2^2 - (Pp.rp - Pp.rb - Pp.L1)^2);
        Pp.servo_offset = deg2rad(U(-1.5,1.5)*[1 0 0] + U(-1.5,1.5)*[0 1 0] + U(-1.5,1.5)*[0 0 1]);
        Pp.tau_s = P.tau_s * U(0.7, 1.3);
        Pp.omega_max = P.omega_max * U(0.8, 1.1);
        Pp.delay_frames = randi([1 4]);      % randi([a b]) = random integer a..b
        Pp.noise_std = U(0.001, 0.003);
        seed = randi(1e6);
        for c = 1:numel(cfgs)
            randn('seed', seed);
            o = sim_3rrs(P, C, cfgs{c}, scn, Pp);
            res(n,c) = o.rms_mm;  fail(n,c) = o.failed;
        end
        fprintf('trial %3d  delay %d  | A1 %6.2f mm (fail %d) | P %6.2f mm (fail %d)\n', ...
            n, Pp.delay_frames, res(n,1), fail(n,1), res(n,2), fail(n,2));
    end

    ok = ~fail;
    fprintf('\nSuccess rate:  A1 %.0f %%   Proposed %.0f %%\n', 100*mean(ok(:,1)), 100*mean(ok(:,2)));
    fprintf('Median RMS:    A1 %.2f mm  Proposed %.2f mm\n', median(res(ok(:,1),1)), median(res(ok(:,2),2)));
    fprintf('95th pct RMS:  A1 %.2f mm  Proposed %.2f mm\n', prctile_simple(res(ok(:,1),1),95), prctile_simple(res(ok(:,2),2),95));

    figure; hold on; grid on;
    edges = linspace(0, max(res(:))*1.05, 25);
    cnt = @(x) arrayfun(@(i) sum(x >= edges(i) & x < edges(i+1)), 1:numel(edges)-1);
    h1 = cnt(res(ok(:,1),1));  h2 = cnt(res(ok(:,2),2));      % manual histogram
    ctr = 0.5*(edges(1:end-1) + edges(2:end));
    bar(ctr, [h1(:) h2(:)], 'grouped');
    xlabel('circle-tracking RMS error [mm]'); ylabel('trials');
    legend('A1 exact IK, naive KF', 'Proposed');
    title(sprintf('Monte Carlo, %d perturbed plants', Ntr));
    print('-dpng','-r110','results/M1_monte_carlo.png');

    csvwrite('results/monte_carlo.csv', [res fail]);
end

%% ======================================================================
%% params_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% configs_3rrs.m
%% ======================================================================

function C = configs_3rrs()
% CONFIGS_3RRS  Controller/kinematics configurations compared in the study.
%   B1/B2 reproduce what typical open-source builds do; P is the proposed
%   pipeline; A1/A2 are the ablation steps between them.
mk = @(n,ik,dc,gv,ff) struct('name',n,'ik',ik,'delay_comp',dc,'governor',gv,'ff',ff);
C.B1 = mk('B1 Linear IK',        'linear',      false, false, false);
C.B2 = mk('B2 Fixed-centre IK',  'fixedcentre', false, false, false);
C.A1 = mk('A1 Exact IK',         'exact',       false, false, false);
C.A2 = mk('A2 Exact+FF',         'exact',       false, false, true);
C.A3 = mk('A3 Exact+FF+DelayKF', 'exact',       true,  false, true);
C.P  = mk('P  Proposed (full)',  'exact',       true,  true,  true);
end

%% ======================================================================
%% scenarios_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% rot_3rrs.m
%% ======================================================================

function R = rot_3rrs(alpha, beta, gamma)
% ROT_3RRS  Platform orientation  R = Rx(alpha) * Ry(beta) * Rz(gamma).
%
%   alpha : roll  about base x  [rad]
%   beta  : pitch about y       [rad]
%   gamma : yaw   about z       [rad]  (parasitic for a 3-RRS, not commanded)
%
%   The matrix maps a vector written in the PLATFORM frame into the BASE
%   frame:  v_base = R * v_platform.   R' (transpose) does the reverse.

ca = cos(alpha); sa = sin(alpha);   % ';' ends a statement and hides output
cb = cos(beta);  sb = sin(beta);
cg = cos(gamma); sg = sin(gamma);

Rx = [1 0 0; 0 ca -sa; 0 sa ca];    % [a b c; d e f] builds a matrix row by row
Ry = [cb 0 sb; 0 1 0; -sb 0 cb];
Rz = [cg -sg 0; sg cg 0; 0 0 1];

R = Rx * Ry * Rz;                   % '*' is matrix multiplication
end

%% ======================================================================
%% parasitic_3rrs.m
%% ======================================================================

function [x, y, gamma, R, res] = parasitic_3rrs(alpha, beta, P)
% PARASITIC_3RRS  Solve the coupled (parasitic) motions of a 3-RRS platform.
%
%   For a 3-RRS, you can only choose 3 of the 6 pose coordinates
%   (here roll alpha, pitch beta, heave h).  The other three -- lateral
%   shift x, y and yaw gamma -- are FORCED by the constraint that every
%   ball joint must stay inside its own leg plane:
%
%       t_i' * ( T + R * p_i ) = 0 ,   i = 1,2,3          (*)
%
%   t_i = tangential unit vector of leg i (normal to its leg plane)
%   p_i = ball-joint position in the platform frame
%   T   = [x; y; h] platform-centre position (h drops out because t_i(3)=0)
%
%   Method:
%     1. Closed-form yaw for the symmetric 120-deg layout
%           gamma = atan2(-sin(alpha) sin(beta), cos(alpha) + cos(beta))
%        (Tetik, Kalla, Kiper & Bandyopadhyay, eq. 7).
%     2. With gamma fixed, (*) is LINEAR in x, y -> least squares.
%     3. Newton polish on [x y gamma] in case geometry is not symmetric.
%   'res' returns the constraint residual so callers can verify it.

c = [cos(P.phi); sin(P.phi)];          % 2x3, column i = radial unit vector u_i
t = [-sin(P.phi); cos(P.phi)];         % 2x3, column i = tangential unit vector t_i
p = P.rp * [c; zeros(1,3)];            % 3x3, ball joints in platform frame

gamma = atan2(-sin(alpha)*sin(beta), cos(alpha) + cos(beta));   % step 1

% step 2: linear solve for x, y given gamma
R   = rot_3rrs(alpha, beta, gamma);
Rp  = R * p;                           % joints rotated into base frame (3x3)
A   = t';                              % 3x2: row i = t_i'
b   = -sum(t .* Rp(1:2,:), 1)';        % 3x1: -(t_i' * R p_i);  .* is element-wise
xy  = A \ b;                           % backslash = least-squares solve
xyg = [xy; gamma];
r   = A*xy - b;                        % residual of (*) with closed-form gamma

% step 3: Newton polish -- only runs if the closed form is not exact
% (e.g. you changed P.phi to a non-symmetric layout).
if norm(r) > 1e-12
    resfun = @(v) constraint_res(v, alpha, beta, t, p);   % anonymous function
    for k = 1:10
        r = resfun(xyg);
        if norm(r) < 1e-13, break; end % norm() = Euclidean length
        J = zeros(3,3); e = 1e-8;
        for j = 1:3                    % finite-difference Jacobian
            d = zeros(3,1); d(j) = e;
            J(:,j) = (resfun(xyg + d) - r) / e;
        end
        xyg = xyg - J \ r;
    end
    R = rot_3rrs(alpha, beta, xyg(3));
end

x = xyg(1); y = xyg(2); gamma = xyg(3);
res = norm(r);
end

function r = constraint_res(xyg, alpha, beta, t, p)
% Local (sub)function: only visible inside this file.
R  = rot_3rrs(alpha, beta, xyg(3));
Pj = R * p + [xyg(1); xyg(2); 0];      % implicit expansion adds T to every column
r  = sum(t .* Pj(1:2,:), 1)';          % t_i' * P_i for each leg -> 3x1
end

%% ======================================================================
%% ik_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% ik_linear_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% ik_fixedcentre_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% fk_3rrs.m
%% ======================================================================

function [q, ok, iters, J] = fk_3rrs(theta, P, q0, J0)
% FK_3RRS  Forward kinematics by Newton-Raphson on the exact IK.
%
%   [q, ok, iters, J] = fk_3rrs(theta, P, q0, J0)
%     theta : 3x1 measured/actual servo angles [rad]
%     q0    : initial guess [alpha; beta; h] (warm start = previous pose)
%     J0    : (optional) Jacobian to reuse -> "chord" Newton, much faster
%             inside a simulation loop where the pose changes very little
%     q     : platform pose that produces theta;  J = Jacobian used
%
%   FK of a 3-RRS has up to 16 real assembly modes (Tetik et al.), so
%   there is no single closed form.  Newton converges to the mode nearest
%   the warm start, which in a control loop is always the physical one.
%
%      q_{k+1} = q_k - J^{-1} ( IK(q_k) - theta ),   J = d theta / d q

if nargin < 3 || isempty(q0)        % isempty([]) is true -> use default
    q0 = [0; 0; P.h0];
end
q  = q0(:);                         % (:) forces a column vector
ok = false;
if nargin < 4 || isempty(J0)
    J = jacobian_3rrs(q, P);
else
    J = J0;
end

for iters = 1:30
    th = ik_3rrs(q, P);
    r  = th - theta(:);
    if norm(r) < 1e-10
        ok = true; break;
    end
    if iters == 6                   % chord iteration slow? refresh Jacobian
        J = jacobian_3rrs(q, P);
    end
    dq = -J \ r;
    s  = min(1, 0.1 / max(norm(dq(1:2)), eps));   % step limiter (rad)
    q  = q + s * dq;
end
end

%% ======================================================================
%% jacobian_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% workspace_3rrs.m
%% ======================================================================

function W = workspace_3rrs(P, h, npsi, margin_min)
% WORKSPACE_3RRS  Reachable-tilt envelope of the platform at heave h.
%
%   W = workspace_3rrs(P, h, npsi, margin_min)
%     h          : heave at which to evaluate (default P.h0)
%     npsi       : number of tilt directions (default 72 -> every 5 deg)
%     margin_min : required distance from the link singularity (default 0.02)
%
%   For each tilt DIRECTION psi the tilt vector is
%        [alpha; beta] = m * [cos(psi); sin(psi)]
%   and we bisect on the MAGNITUDE m until the exact IK says the pose is
%   no longer feasible (servo limit, link can't close, or too close to a
%   singularity).  Result: W.tilt_max(psi), a lookup table used by the
%   online workspace governor.  It is direction-dependent because the
%   three legs sit at 120 deg: tilting "towards" a leg is not the same as
%   tilting "between" two legs.

if nargin < 2 || isempty(h),          h = P.h0;        end
if nargin < 3 || isempty(npsi),       npsi = 72;       end
if nargin < 4 || isempty(margin_min), margin_min = 0.02; end

psi = linspace(-pi, pi, npsi+1);   % linspace(a,b,n): n evenly spaced points
tilt_max = zeros(size(psi));       % size() returns [rows cols]

for k = 1:numel(psi)               % numel = number of elements
    lo = 0; hi = deg2rad(45);
    for it = 1:30                  % bisection: halves the interval each step
        m = 0.5*(lo + hi);
        q = [m*cos(psi(k)); m*sin(psi(k)); h];
        [~, ok, info] = ik_3rrs(q, P);        % '~' discards an output
        if ok && min(info.reach_margin) > margin_min
            lo = m;
        else
            hi = m;
        end
    end
    tilt_max(k) = lo;
end

W.psi = psi;  W.tilt_max = tilt_max;  W.h = h;
W.min_tilt = min(tilt_max);        % the "guaranteed in every direction" tilt
end

%% ======================================================================
%% dlqr_riccati.m
%% ======================================================================

function [K, S] = dlqr_riccati(A, B, Q, R)
% DLQR_RICCATI  Discrete LQR by iterating the Riccati difference equation.
%
%   Minimises  J = sum( x'Qx + u'Ru ),  x(k+1) = A x(k) + B u(k)
%   Optimal law u = -K x.   Same result as Control Toolbox dlqr(), but
%   written out so (a) no toolbox is needed and (b) you can see the maths:
%
%     S <- Q + A' S A - A' S B (R + B' S B)^{-1} B' S A     (repeat)
%     K  = (R + B' S B)^{-1} B' S A

S = Q;
for k = 1:20000
    K  = (R + B'*S*B) \ (B'*S*A);         % A\B solves A*X = B
    Sn = Q + A'*S*A - A'*S*B*K;
    if max(abs(Sn(:) - S(:))) < 1e-12 * max(1, max(abs(S(:))))
        S = Sn; break;
    end
    S = Sn;
end
K = (R + B'*S*B) \ (B'*S*A);
end

%% ======================================================================
%% design_controller_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% governor_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% sim_3rrs.m
%% ======================================================================

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

%% ======================================================================
%% animate_3rrs.m
%% ======================================================================

function animate_3rrs(out, P, opts)
% ANIMATE_3RRS  3-D animation of the mechanism and ball from a sim log.
%
%   animate_3rrs(out, P)                      play on screen
%   animate_3rrs(out, P, struct('file','results/anim.mp4'))   also save video
%   animate_3rrs(out, P, struct('snap','results/A1_snapshot.png','k',120))
%        render ONE frame k to a PNG (works headless / in Octave)
%
%   Draws: base triangle, servo shafts B_i, horns B_i->C_i, couplers
%   C_i->P_i, the plate disc (tilted, shifted, yawed exactly as the
%   kinematics say), the ball, and the reference point.

if nargin < 3, opts = struct(); end
L = out.log;
n = find(~isnan(L.q(:,1)), 1, 'last');          % last valid frame
ks = 1:2:n;                                      % every 2nd frame (30 FPS)
if isfield(opts,'k'), ks = min(opts.k, n); end

fig = figure('Color','w');
vw = isfield(opts,'file');
if vw && exist('VideoWriter','file')             % MATLAB only
    V = VideoWriter(opts.file, 'MPEG-4'); V.FrameRate = 30; open(V);
else
    vw = false;
end

ang = linspace(0, 2*pi, 60);
[sx, sy, sz] = sphere(16);                       % unit-sphere mesh

for k = ks
    clf; hold on; axis equal; grid on;
    q = L.q(k,:)';
    [~, ~, I] = ik_3rrs(q, P);
    % base
    plot3(P.rb*1.3*cos(P.phi([1:3 1])), P.rb*1.3*sin(P.phi([1:3 1])), zeros(1,4), 'k-', 'LineWidth', 2);
    for i = 1:3
        plot3([I.B(1,i) I.C(1,i)], [I.B(2,i) I.C(2,i)], [I.B(3,i) I.C(3,i)], 'b-', 'LineWidth', 4);
        plot3([I.C(1,i) I.P(1,i)], [I.C(2,i) I.P(2,i)], [I.C(3,i) I.P(3,i)], 'k-', 'LineWidth', 2);
        plot3(I.B(1,i), I.B(2,i), I.B(3,i), 'ks', 'MarkerFaceColor', 'k');
        plot3(I.P(1,i), I.P(2,i), I.P(3,i), 'ro', 'MarkerFaceColor', 'r');
    end
    % plate disc: points in plate frame -> base frame
    % drawn as a 2 x N surface (centre ring + rim) so it renders in MATLAB and Octave
    rim  = I.R * [P.plate_radius*cos(ang); P.plate_radius*sin(ang); zeros(size(ang))] + I.T;
    ctr  = repmat(I.T, 1, numel(ang));                 % repmat tiles the centre N times
    surf([ctr(1,:); rim(1,:)], [ctr(2,:); rim(2,:)], [ctr(3,:); rim(3,:)], ...
         'FaceColor', [0.55 0.75 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot3(rim(1,:), rim(2,:), rim(3,:), 'Color', [0.2 0.4 0.7]);
    % ball (sits on top of the plate along the plate normal)
    bp = I.R * [L.ball(k,1); L.ball(k,2); P.rball] + I.T;
    surf(bp(1) + P.rball*sx, bp(2) + P.rball*sy, bp(3) + P.rball*sz, ...
         'FaceColor', [0.9 0.3 0.1], 'EdgeColor', 'none');
    rp = I.R * [L.ref(k,1); L.ref(k,2); 0.001] + I.T;
    plot3(rp(1), rp(2), rp(3), 'g+', 'MarkerSize', 12, 'LineWidth', 2);
    view(35, 25); xlim([-0.12 0.12]); ylim([-0.12 0.12]); zlim([0 0.13]);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title(sprintf('%s   t = %.2f s   tilt = (%.1f, %.1f) deg', out.cfg.name, L.t(k), ...
          rad2deg(q(1)), rad2deg(q(2))));
    drawnow;
    if vw, writeVideo(V, getframe(fig)); end
end
if vw, close(V); end
if isfield(opts,'snap'), print('-dpng','-r120', opts.snap); end
end

%% ======================================================================
%% prctile_simple.m
%% ======================================================================

function p = prctile_simple(x, pc)
% PRCTILE_SIMPLE  Toolbox-free percentile with linear interpolation.
%   p = prctile_simple(x, 95)  -> 95th percentile of the values in x.
x = sort(x(:));  n = numel(x);       % sort ascending, as a column
if n == 0, p = NaN; return; end      % 'return' leaves the function early
k  = 1 + (n-1)*pc/100;               % fractional rank
lo = floor(k);  hi = ceil(k);
p  = x(lo) + (k-lo)*(x(hi)-x(lo));
end
