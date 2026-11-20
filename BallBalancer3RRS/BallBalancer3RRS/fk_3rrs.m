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
