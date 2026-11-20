function W = workspace_3rrs(P, h, npsi, margin_min)
% WORKSPACE_3RRS  Reachable-tilt envelope of the platform at heave h.
%
%   W = workspace_3rrs(P, h, npsi, margin_min)
%     h          : heave at which to evaluate (default P.h0)
%     npsi       : number of tilt directions (default 72 -> every 5 deg)
%     margin_min : required distance from the link singularity (default 0.02)
%
%   For each tilt DIRECTION psi the tilt vector is
%        [alpha; beta] = m * [cos(psi); sin(psi)]
%   and we bisect on the MAGNITUDE m until the exact IK says the pose is
%   no longer feasible (servo limit, link can't close, or too close to a
%   singularity).  Result: W.tilt_max(psi), a lookup table used by the
%   online workspace governor.  It is direction-dependent because the
%   three legs sit at 120 deg: tilting "towards" a leg is not the same as
%   tilting "between" two legs.

if nargin < 2 || isempty(h),          h = P.h0;        end
if nargin < 3 || isempty(npsi),       npsi = 72;       end
if nargin < 4 || isempty(margin_min), margin_min = 0.02; end

psi = linspace(-pi, pi, npsi+1);   % linspace(a,b,n): n evenly spaced points
tilt_max = zeros(size(psi));       % size() returns [rows cols]

for k = 1:numel(psi)               % numel = number of elements
    lo = 0; hi = deg2rad(45);
    for it = 1:30                  % bisection: halves the interval each step
        m = 0.5*(lo + hi);
        q = [m*cos(psi(k)); m*sin(psi(k)); h];
        [~, ok, info] = ik_3rrs(q, P);        % '~' discards an output
        if ok && min(info.reach_margin) > margin_min
            lo = m;
        else
            hi = m;
        end
    end
    tilt_max(k) = lo;
end

W.psi = psi;  W.tilt_max = tilt_max;  W.h = h;
W.min_tilt = min(tilt_max);        % the "guaranteed in every direction" tilt
end
