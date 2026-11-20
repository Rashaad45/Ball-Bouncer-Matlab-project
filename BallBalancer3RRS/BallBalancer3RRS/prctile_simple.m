function p = prctile_simple(x, pc)
% PRCTILE_SIMPLE  Toolbox-free percentile with linear interpolation.
%   p = prctile_simple(x, 95)  -> 95th percentile of the values in x.
x = sort(x(:));  n = numel(x);       % sort ascending, as a column
if n == 0, p = NaN; return; end      % 'return' leaves the function early
k  = 1 + (n-1)*pc/100;               % fractional rank
lo = floor(k);  hi = ceil(k);
p  = x(lo) + (k-lo)*(x(hi)-x(lo));
end
