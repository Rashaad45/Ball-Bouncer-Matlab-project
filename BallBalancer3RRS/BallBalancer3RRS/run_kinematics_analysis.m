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

clear; clc;                       % clear variables, clear command window
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
