%% RUN_MONTE_CARLO  Robustness: the controller is designed on NOMINAL
%  parameters, but every trial gets a different TRUE plant:
%    link lengths / radii  +/- 0.5 mm      (print + assembly tolerance)
%    servo zero offsets    +/- 1.5 deg     (horn spline calibration)
%    servo tau             +/- 30 %, speed  -20..+10 %
%    camera latency        1 .. 4 frames  (controller always assumes 2)
%    camera noise          1 .. 3 mm
%  Output: results/M1_monte_carlo.png, results/monte_carlo.csv

clear; clc;
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
