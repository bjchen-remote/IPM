function report = ipm_accellab_hermite_level_roots(values,jets,x,target,window)
%IPM_ACCELLAB_HERMITE_LEVEL_ROOTS All roots of one fixed piecewise cubic.
% Root residuals, local conditioning, duplicate node representations and
% whole threshold cells are retained. The caller chooses connected roots.
x = x(:)'; values = values(:)'; jets = jets(:)';
locations = []; slopes = []; cells = []; fractions = []; residuals = []; normalizedSlopes = [];
tolerances = []; duplicateCounts = []; flatCells = [];
for i = 1:numel(x)-1
    if x(i+1) < window(1) || x(i) > window(2), continue; end
    h = x(i+1)-x(i); v0 = values(i); v1 = values(i+1); d0 = jets(i); d1 = jets(i+1);
    c = [2*v0-2*v1+h*(d0+d1),-3*v0+3*v1-h*(2*d0+d1),h*d0,v0-target];
    scale = max(abs(c));
    if scale == 0, flatCells(end+1) = i; continue; end %#ok<AGROW>
    points = roots(c/scale);
    points = sort(real(points(abs(imag(points)) <= 1e-10 & real(points) >= -128*eps & real(points) <= 1+128*eps)));
    for point = points(:)'
        t = min(1,max(0,point)); z = x(i)+h*t;
        if z < window(1) || z > window(2), continue; end
        q = ipm_accellab_hermite_line(values,jets,x,z,i);
        tolerance = 128*eps*max([abs(z),abs(x(i:i+1)),h])+4*max(eps(abs(z)),eps(h));
        relativeResidual = abs(q.value-target)/max(q.polynomialScales,abs(target));
        if ~isempty(locations) && abs(z-locations(end)) <= max(tolerance,tolerances(end))
            duplicateCounts(end) = duplicateCounts(end)+1;
            if relativeResidual >= residuals(end), continue; end
            index = numel(locations);
        else
            index = numel(locations)+1; duplicateCounts(index) = 1; %#ok<AGROW>
        end
        locations(index) = z; slopes(index) = q.derivative; cells(index) = i; %#ok<AGROW>
        fractions(index) = q.fractions; residuals(index) = relativeResidual; %#ok<AGROW>
        normalizedSlopes(index) = q.normalizedDerivative; tolerances(index) = tolerance; %#ok<AGROW>
    end
end
report = struct('locations',locations,'slopes',slopes,'cells',cells,'fractions',fractions, ...
    'normalizedSlopes',normalizedSlopes,'relativeResiduals',residuals, ...
    'positionTolerances',tolerances,'duplicateCounts',duplicateCounts,'flatCells',flatCells, ...
    'target',target,'window',window,'kind','same_C1_piecewise_cubic_level_roots');
end
