%% TEST_3RRS  Self-checks.  Run this first after any change to the maths.
%  Every assert() stops with an error message if its condition is false.

clear; clc;
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
