%% RUN_LATENCY_STUDY  Ablation: how much does each layer buy as camera
%  latency grows?  (The question every vision-based balancer runs into.)
%
%  Output: results/L1_latency_ablation.png, results/latency_ablation.csv

clear; clc;
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
