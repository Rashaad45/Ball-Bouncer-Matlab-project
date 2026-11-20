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
