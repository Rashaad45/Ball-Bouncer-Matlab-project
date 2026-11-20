%% RUN_CLOSED_LOOP  Benchmark 4 scenarios x 4 configurations.
%
%  Outputs in ./results:
%   C1_<scenario>.png        XY path + error vs time for every config
%   C2_servos_proposed.png   servo angles of the proposed controller (circle)
%   closed_loop_results.csv  metrics table
%   closed_loop_runs.mat     all logs (load with: load results/closed_loop_runs.mat)

clear; clc;
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
