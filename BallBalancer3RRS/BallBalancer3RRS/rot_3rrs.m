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
