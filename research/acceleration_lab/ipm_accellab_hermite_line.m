function out = ipm_accellab_hermite_line(values,jets,x,query,cells)
%IPM_ACCELLAB_HERMITE_LINE Value and derivatives of the same fixed C1 cubic.
x = x(:)'; values = values(:)'; jets = jets(:)';
assert(numel(values) == numel(x) && numel(jets) == numel(x) && all(diff(x) > 0));
assert(all(isfinite(query),'all') && all(query >= x(1) & query <= x(end),'all'));
if nargin < 5, cells = discretize(query,x); end
assert(isequal(size(cells),size(query)) && all(cells >= 1 & cells < numel(x),'all'));
left = reshape(x(cells),size(query)); right = reshape(x(cells+1),size(query));
h = right-left; t = (query-left)./h;
assert(all(t >= -128*eps & t <= 1+128*eps,'all'));
v0 = reshape(values(cells),size(query)); v1 = reshape(values(cells+1),size(query));
d0 = reshape(jets(cells),size(query)); d1 = reshape(jets(cells+1),size(query));
c0 = v0; c1 = h.*d0; c2 = -3*v0+3*v1-h.*(2*d0+d1); c3 = 2*v0-2*v1+h.*(d0+d1);
v = ((c3.*t+c2).*t+c1).*t+c0;
dx = (3*c3.*t.^2+2*c2.*t+c1)./h; dxx = (6*c3.*t+2*c2)./h.^2;
v(t == 0) = v0(t == 0); v(t == 1) = v1(t == 1);
dx(t == 0) = d0(t == 0); dx(t == 1) = d1(t == 1);
scale = max(cat(ndims(query)+1,abs(v0),abs(v1),abs(h.*d0),abs(h.*d1), ...
    abs(c0),abs(c1),abs(c2),abs(c3)),[],ndims(query)+1);
scale = max(scale,realmin);
out = struct('value',v,'derivative',dx,'secondDerivative',dxx,'cells',cells, ...
    'fractions',t,'cellWidths',h,'polynomialScales',scale, ...
    'normalizedDerivative',dx.*h./scale,'normalizedCurvature',dxx.*h.^2./scale);
end
