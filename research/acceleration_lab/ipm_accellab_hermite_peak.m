function report = ipm_accellab_hermite_peak(values,forcing,x,Dx,window,user)
%IPM_ACCELLAB_HERMITE_PEAK Continuous maximum of one fixed linear C1 spline.
%   P=max H_h[values] on the fixed closed WINDOW. For an isolated interior
%   maximizer, PPrime=H_h[forcing](a). No time evolution or target P is used.
%   Multiple/near-tied and degenerate maxima are explicit validity failures
%   by default. Optional finite isolated ties use a one-sided Dini derivative.
if nargin < 6, user = struct(); end
opts = struct('activeRelativeTolerance',128*eps,'positionRelativeTolerance',128*eps, ...
    'stationarityRelativeTolerance',1e-10,'minimumRelativeCurvature',1e-10, ...
    'allowIsolatedTies',false,'allowBoundaryMaximum',false);
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:HermitePeakOptions','Unknown peak option.');
for k = 1:numel(names), opts.(names{k}) = user.(names{k}); end
x = x(:)'; values = values(:)'; forcing = forcing(:)';
validateattributes(x,{'numeric'},{'finite','real','increasing'});
validateattributes(values,{'numeric'},{'vector','real','finite','numel',numel(x)});
validateattributes(forcing,{'numeric'},{'vector','real','finite','numel',numel(x)});
validateattributes(window,{'numeric'},{'vector','numel',2,'finite','real','increasing'});
window = window(:)';
assert(window(1) >= x(1) && window(2) <= x(end) && ...
    isequal(size(Dx),[numel(x),numel(x)]),'ipm:HermitePeakDomain','Paired fixed derivative/domain data are required.');
slopes = values*Dx'; forcingSlopes = forcing*Dx';
positions = []; samples = []; sampleForcing = []; derivative = []; curvature = [];
cells = []; fractions = []; boundary = [];
localWidths = []; localScales = []; positionTolerances = [];
for i = 1:numel(x)-1
    a = max(x(i),window(1)); b = min(x(i+1),window(2));
    if b <= a, continue; end
    h = x(i+1)-x(i);
    c = coefficients(values(i),values(i+1),slopes(i),slopes(i+1),h);
    f = coefficients(forcing(i),forcing(i+1),forcingSlopes(i),forcingSlopes(i+1),h);
    low = (a-x(i))/h; high = (b-x(i))/h;
    roots = stationary_roots(c);
    roots = roots(roots > low & roots < high);
    t = [low,high,roots];
    v = polyval(fliplr(c),t); fv = polyval(fliplr(f),t);
    dv = (c(2)+2*c(3)*t+3*c(4)*t.^2)/h;
    cv = (2*c(3)+6*c(4)*t)/h^2;
    % Shared-node representations use the identical nodal value and jet.
    v(t == 0) = values(i); v(t == 1) = values(i+1);
    fv(t == 0) = forcing(i); fv(t == 1) = forcing(i+1);
    dv(t == 0) = slopes(i); dv(t == 1) = slopes(i+1);
    z = x(i)+h*t;
    coefficientScale = max([abs(values(i:i+1)),abs(h*slopes(i:i+1)),abs(c),realmin]);
    localPositionScale = max([abs(z);repmat(max(abs(x(i:i+1))),size(z));repmat(h,size(z))],[],1);
    localPositionTolerance = opts.positionRelativeTolerance*localPositionScale+ ...
        4*max(eps(abs(z)),eps(h));
    positions = [positions,z]; samples = [samples,v]; %#ok<AGROW>
    sampleForcing = [sampleForcing,fv]; derivative = [derivative,dv]; %#ok<AGROW>
    curvature = [curvature,cv]; cells = [cells,repmat(i,size(t))]; %#ok<AGROW>
    fractions = [fractions,t]; %#ok<AGROW>
    boundary = [boundary,z == window(1) | z == window(2)]; %#ok<AGROW>
    localWidths = [localWidths,repmat(h,size(t))]; %#ok<AGROW>
    localScales = [localScales,repmat(coefficientScale,size(t))]; %#ok<AGROW>
    positionTolerances = [positionTolerances,localPositionTolerance]; %#ok<AGROW>
end
assert(~isempty(samples) && all(isfinite(samples)),'ipm:HermitePeakFinite','No finite maximum candidates.');
[P,index] = max(samples);
scale = max(max(abs(values)),realmin);
valueTolerance = opts.activeRelativeTolerance*scale;
normalizedStationarity = abs(derivative).*localWidths./localScales;
normalizedCurvature = curvature.*localWidths.^2./localScales;
eligible = boundary | normalizedStationarity <= opts.stationarityRelativeTolerance;
nearMaximum = samples >= P-valueTolerance;
active = find(nearMaximum & eligible);
[~,order] = sort(positions(active)); active = active(order);
representatives = [];
for i = active
    if isempty(representatives) || abs(positions(i)-positions(representatives(end))) > ...
            max(positionTolerances([i,representatives(end)]))
        representatives(end+1) = i; %#ok<AGROW>
    elseif samples(i) > samples(representatives(end))
        representatives(end) = i;
    end
end
reasons = {};
if P <= 0, reasons{end+1} = 'nonpositive_peak'; end
if isempty(representatives), reasons{end+1} = 'unresolved_stationary_maximum'; end
if any(boundary(active)) && ~opts.allowBoundaryMaximum
    reasons{end+1} = 'maximum_at_search_boundary';
end
interiorActive = active(~boundary(active));
if any(normalizedCurvature(interiorActive) >= -opts.minimumRelativeCurvature)
    reasons{end+1} = 'flat_or_degenerate_maximum';
end
if numel(representatives) > 1 && ~opts.allowIsolatedTies
    reasons{end+1} = 'multiple_or_near_tied_maxima';
end
if ~isempty(representatives)
    [~,j] = max(samples(representatives)); index = representatives(j);
end
PPrime = NaN; derivativeInterval = [NaN,NaN];
if ~isempty(representatives)
    derivativeInterval = [min(sampleForcing(representatives)),max(sampleForcing(representatives))];
    if isempty(reasons), PPrime = derivativeInterval(2); end
end
report = struct('kind','independent_linear_C1_hermite_maximum_functional', ...
    'value',P,'x',positions(index),'valid',isempty(reasons),'invalidReasons',{reasons}, ...
    'PPrime',PPrime,'activeForcingDerivativeInterval',derivativeInterval, ...
    'isolatedPeakCount',numel(representatives),'smoothUniquePeak',isscalar(representatives) && isempty(reasons), ...
    'activePositions',positions(representatives),'activeValues',samples(representatives), ...
    'activeForcingValues',sampleForcing(representatives), ...
    'activeValueSpread',max(samples(nearMaximum))-min(samples(nearMaximum)), ...
    'activeValueTolerance',valueTolerance,'activePositionTolerance',positionTolerances(index), ...
    'selectedCell',cells(index),'selectedFraction',fractions(index), ...
    'selectedStationarityResidual',derivative(index),'selectedSecondDerivative',curvature(index), ...
    'selectedNormalizedStationarityResidual',normalizedStationarity(index), ...
    'selectedNormalizedCurvature',normalizedCurvature(index), ...
    'selectedCellWidth',localWidths(index),'selectedPolynomialScale',localScales(index), ...
    'activeOneSidedSecondDerivatives',curvature(active), ...
    'activeNormalizedCurvatures',normalizedCurvature(active), ...
    'activePositionTolerances',positionTolerances(representatives), ...
    'atNativeNode',fractions(index) == 0 || fractions(index) == 1, ...
    'window',window,'options',opts,'candidatePositions',positions, ...
    'candidateValues',samples,'candidateCells',cells, ...
    'tolerancePolicy','local_cell_polynomial_scale_and_local_coordinate_ulp_v1', ...
    'isProductionGauge',false,'pdeSteps',0, ...
    'interpretation','The maximum value is continuous in the field. PPrime is an envelope derivative on resolved active maxima; pointwise cancellation alone does not certify finite-step conservation.');
end

function c = coefficients(v0,v1,d0,d1,h)
c = [v0,h*d0,-3*v0+3*v1-h*(2*d0+d1),2*v0-2*v1+h*(d0+d1)];
end

function points = stationary_roots(c)
a = 3*c(4); b = 2*c(3); d = c(2);
scale = max(abs([a,b,d])); tolerance = 64*eps*scale;
if scale == 0
    points = [];
elseif abs(a) <= tolerance
    if abs(b) <= tolerance, points = []; else, points = -d/b; end
else
    discriminant = b*b-4*a*d;
    toleranceD = 64*eps*(b*b+abs(4*a*d));
    if discriminant < -toleranceD
        points = [];
    else
        root = sqrt(max(discriminant,0));
        signB = 1; if b < 0, signB = -1; end
        q = -0.5*(b+signB*root);
        if q == 0, points = -b/(2*a); else, points = unique([q/a,d/q]); end
    end
end
end
