function report = green()
%IPMTESTS.FOURTH.GREEN Gate smooth, uncompressed Green quadrature.
%   The oracle independently integrates the analytic image kernels with
%   INTEGRAL2.  The manufactured source vanishes to sixth order on the
%   rectangle boundary and at x=0, where DOUBLE_ODD_OMEGA truncates it.
%   This test deliberately disables source cutoff and aggregation.  It
%   certifies smooth uncompressed Green quadrature only; it does not certify
%   the production default fixed-3000 source compression.

threshold = 3.8;
gridSizes = [17,33,65,129];
gridModes = {'uniform','stretched'};
symmetryModes = {'half_plane','double_odd_omega'};
kappas = [0.5,1,2];
stretch = [0.8,0.6];
greenSourceTolerance = realmin;
greenMaxSources = 1e6;

fprintf(['Running high-order smooth Green verification ' ...
    '(required order %.1f)...\n'],threshold);
report = struct('threshold',threshold,'gridSizes',gridSizes, ...
    'sampleNorm','fixed-boundary-node Linf', ...
    'scope','uncompressed smooth Green quadrature only', ...
    'excludedScope','default fixed-3000 source compression', ...
    'greenSourceTolerance',greenSourceTolerance, ...
    'greenMaxSources',greenMaxSources);
for modeIndex = 1:numel(gridModes)
    mode = gridModes{modeIndex};
    report.(mode) = verify_grid_family(mode,gridSizes,stretch, ...
        symmetryModes,kappas,threshold,greenSourceTolerance, ...
        greenMaxSources);
end
report.passed = true;
fprintf(['High-order smooth, uncompressed Green verification passed. ' ...
    'Fixed-3000 compression was not tested.\n']);
end

function family = verify_grid_family(mode,gridSizes,stretch, ...
        symmetryModes,kappas,threshold,sourceTolerance,maxSources)
h = 1./(gridSizes-1);
metricNames = {'leftLinf','rightLinf','topLinf','sampleLinf'};
numberOfCases = numel(symmetryModes)*numel(kappas);
cases = cell(numberOfCases,1);
caseIndex = 0;
targets = fixed_boundary_targets(mode,stretch);

for symmetryIndex = 1:numel(symmetryModes)
    symmetryMode = symmetryModes{symmetryIndex};
    for kappaIndex = 1:numel(kappas)
        kappa = kappas(kappaIndex);
        caseIndex = caseIndex+1;
        exact = oracle_boundary(targets,symmetryMode,kappa);
        normalizers = side_norms(exact);
        absoluteErrors = zeros(numel(metricNames),numel(gridSizes));
        relativeErrors = zeros(size(absoluteErrors));

        for level = 1:numel(gridSizes)
            numberOfNodes = gridSizes(level);
            ops = green_operators(numberOfNodes,mode,stretch, ...
                symmetryMode,sourceTolerance,maxSources);
            assert_sample_coordinates(ops,numberOfNodes,targets);
            source = smooth_source(ops.X,ops.Y);
            boundary = ipm.field.greenBoundary(source,ops,kappa);
            numerical = sampled_boundary(boundary,numberOfNodes);
            absoluteErrors(:,level) = boundary_errors(numerical,exact);
            relativeErrors(:,level) = absoluteErrors(:,level)./normalizers;
        end

        orders = observed_orders(relativeErrors,h);
        finalOrders = orders(:,end);
        label = sprintf('%s/%s/kappa %.1f smooth Green', ...
            mode,symmetryMode,kappa);
        tailOrders = certify_green_tail( ...
            relativeErrors,orders,threshold,label);

        cases{caseIndex} = struct('symmetryMode',symmetryMode, ...
            'kappa',kappa,'h',h,'metricNames',{metricNames}, ...
            'absoluteErrors',absoluteErrors, ...
            'relativeErrors',relativeErrors,'orders',orders, ...
            'tailOrders',tailOrders,'finalOrders',finalOrders, ...
            'oracle',exact);
        fprintf('  %-9s %-16s kappa %.1f final orders: %s\n', ...
            mode,symmetryMode,kappa, ...
            format_named_values(metricNames,finalOrders));
    end
end

family = struct('gridMode',mode,'stretch',stretch, ...
    'gridSizes',gridSizes,'h',h,'targets',targets, ...
    'metricNames',{metricNames},'cases',{cases});
end

function ops = green_operators(numberOfNodes,mode,stretch, ...
        symmetryMode,sourceTolerance,maxSources)
options = struct('nx',numberOfNodes,'ny',numberOfNodes, ...
    'xlim',[-1,1],'ymax',1,'gridMode',mode, ...
    'gridStretchAutomatic',false,'gridStretch',stretch, ...
    'spatialDiscretization','high_order', ...
    'farBoundaryMode','green','symmetryMode',symmetryMode, ...
    'greenSourceTolerance',sourceTolerance, ...
    'greenMaxSources',maxSources,'adaptiveRemesh',false, ...
    'rescalingMode','physical','saveResults',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'verbose',false);
ops = ipm.mesh.build(ipm.config.resolve(options));
end

function targets = fixed_boundary_targets(mode,stretch)
referenceX = [-0.75,-0.25,0,0.25,0.75];
referenceY = [0.125,0.375,0.625,0.875];
if strcmp(mode,'stretched')
    sampleX = sinh(stretch(1)*referenceX)/sinh(stretch(1));
    sampleY = sinh(stretch(2)*referenceY)/sinh(stretch(2));
else
    sampleX = referenceX;
    sampleY = referenceY;
end
targets = struct( ...
    'leftX',-ones(size(sampleY)),'leftY',sampleY, ...
    'rightX',ones(size(sampleY)),'rightY',sampleY, ...
    'topX',sampleX,'topY',ones(size(sampleX)), ...
    'referenceX',referenceX,'referenceY',referenceY);
end

function assert_sample_coordinates(ops,numberOfNodes,targets)
[xIndices,yIndices] = sample_indices(numberOfNodes);
tolerance = 64*eps;
assert(max(abs(ops.x(xIndices)-targets.topX)) < tolerance && ...
    max(abs(ops.y(yIndices)'-targets.leftY)) < tolerance, ...
    'ipm:HighOrderGreenSampling', ...
    'Green samples must stay at fixed physical points under refinement.');
end

function numerical = sampled_boundary(boundary,numberOfNodes)
[xIndices,yIndices] = sample_indices(numberOfNodes);
numerical = struct('left',boundary.left(yIndices), ...
    'right',boundary.right(yIndices), ...
    'top',boundary.top(xIndices)');
end

function [xIndices,yIndices] = sample_indices(numberOfNodes)
xReference = [-0.75,-0.25,0,0.25,0.75];
yReference = [0.125,0.375,0.625,0.875];
xIndices = round((xReference+1)*(numberOfNodes-1)/2)+1;
yIndices = round(yReference*(numberOfNodes-1))+1;
end

function exact = oracle_boundary(targets,symmetryMode,kappa)
exact = struct( ...
    'left',oracle_side(targets.leftX,targets.leftY, ...
        symmetryMode,kappa), ...
    'right',oracle_side(targets.rightX,targets.rightY, ...
        symmetryMode,kappa), ...
    'top',oracle_side(targets.topX,targets.topY, ...
        symmetryMode,kappa));
end

function values = oracle_side(targetX,targetY,symmetryMode,kappa)
values = zeros(numel(targetX),1);
if strcmp(symmetryMode,'double_odd_omega')
    xLower = 0;
else
    xLower = -1;
end
for index = 1:numel(targetX)
    values(index) = integral2(@(sourceX,sourceY) ...
        oracle_integrand(sourceX,sourceY,targetX(index), ...
        targetY(index),symmetryMode,kappa),xLower,1,0,1, ...
        'Method','iterated','AbsTol',2e-13,'RelTol',2e-12);
end
end

function values = oracle_integrand(sourceX,sourceY,targetX,targetY, ...
        symmetryMode,kappa)
sqrtKappa = sqrt(kappa);
direct = hypot(sqrtKappa*(targetX-sourceX), ...
    (targetY-sourceY)/sqrtKappa);
yImage = hypot(sqrtKappa*(targetX-sourceX), ...
    (targetY+sourceY)/sqrtKappa);
if strcmp(symmetryMode,'double_odd_omega')
    xImage = hypot(sqrtKappa*(targetX+sourceX), ...
        (targetY-sourceY)/sqrtKappa);
    xyImage = hypot(sqrtKappa*(targetX+sourceX), ...
        (targetY+sourceY)/sqrtKappa);
    kernel = (log(xImage)+log(yImage)-log(direct)-log(xyImage)) / ...
        (2*pi);
else
    kernel = (log(yImage)-log(direct))/(2*pi);
end
values = smooth_source(sourceX,sourceY).*kernel;
values(~isfinite(values)) = 0;
end

function values = smooth_source(x,y)
% Sixth-order zeros regularize boundary logarithms and the x=0 cutoff.
values = 1e6*(x.*(1-x.^2).*y.*(1-y)).^6 .* ...
    exp(0.23*x-0.17*y);
end

function errors = boundary_errors(numerical,exact)
left = norm(numerical.left-exact.left,inf);
right = norm(numerical.right-exact.right,inf);
top = norm(numerical.top-exact.top,inf);
allSamples = norm([numerical.left-exact.left; ...
    numerical.right-exact.right;numerical.top-exact.top],inf);
errors = [left;right;top;allSamples];
end

function values = side_norms(exact)
left = max(norm(exact.left,inf),realmin);
right = max(norm(exact.right,inf),realmin);
top = max(norm(exact.top,inf),realmin);
allSamples = max(norm([exact.left;exact.right;exact.top],inf),realmin);
values = [left;right;top;allSamples];
end

function orders = observed_orders(errors,h)
orders = log(errors(:,1:end-1)./errors(:,2:end)) ./ ...
    log(h(1:end-1)./h(2:end));
end

function tailOrders = certify_green_tail(errors,orders,threshold,label)
tailSegmentCount = 2;
tailOrders = orders(:,end-tailSegmentCount+1:end);
assert(all(isfinite(errors),'all') && all(errors > 0,'all') && ...
    all(diff(errors,1,2) < 0,'all') && ...
    all(tailOrders >= threshold,'all'), ...
    'ipm:HighOrderGreenOrder', ...
    ['%s must have finite positive, monotonically decreasing errors and ' ...
    'both tail orders at least %.1f; minimum tail order %.3f. This gate ' ...
    'covers uncompressed smooth quadrature only.'], ...
    label,threshold,min(tailOrders,[],'all'));
end

function text = format_named_values(names,values)
parts = cell(size(names));
for index = 1:numel(names)
    parts{index} = sprintf('%s %.3f',names{index},values(index));
end
text = strjoin(parts,', ');
end
