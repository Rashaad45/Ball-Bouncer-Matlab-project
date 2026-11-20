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
