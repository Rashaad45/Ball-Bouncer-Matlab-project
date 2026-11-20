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
