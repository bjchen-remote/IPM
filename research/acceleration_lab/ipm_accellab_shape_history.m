function report = ipm_accellab_shape_history(states)
%IPM_ACCELLAB_SHAPE_HISTORY Separate inner contraction from shape relaxation.
%   Read-only observations of paired same-grid solver states; no PDE solves.
%   The inner coordinates use instantaneous measured full 90%-level widths.
%   Linear interpolation here is diagnostic and supplies no accuracy order.
assert(iscell(states) && numel(states) >= 2,'ipm:AccelShapeHistory', ...
    'At least two paired states are required.');
count = numel(states);
tau = zeros(1,count);
peak = zeros(1,count);
peakX = zeros(1,count);
curvature = zeros(1,count);
wallWidth = zeros(1,count);
verticalWidth = zeros(1,count);
wallHalfWidth = zeros(1,count);
verticalHalfWidth = zeros(1,count);
wallCoordinates = linspace(-2,2,101)';
verticalCoordinates = linspace(0,3,101)';
wallProfiles = NaN(101,count);
verticalProfiles = NaN(101,count);
for index = 1:count
    state = states{index};
    ops = state.ops;
    assert(isequal(size(state.rho),[numel(ops.y),numel(ops.x)]), ...
        'ipm:AccelShapeGrid','Each state must retain its paired grid.');
    if index > 1
        assert(isequal(ops.x,states{1}.ops.x) && ...
            isequal(ops.y,states{1}.ops.y),'ipm:AccelShapeGrid', ...
            'This diagnostic requires a common fixed grid.');
    end
    source = state.rho*ops.Dx';
    mask = abs(ops.x-ops.rescaling.pinX) <= ops.rescaling.peakTrackingHalfWidth;
    wall = ipm.diagnostics.measureFeatures(source(1,:),ops.x,[0.5,0.9],mask);
    vertical = ipm.diagnostics.measureFeatures( ...
        abs(source(:,wall.peakIndex))',ops.y',[0.5,0.9],[]);
    quadratic = ipm.evolve.quadraticPeakFunctional( ...
        source(1,:),ops.x,wall.peakIndex);
    tau(index) = state.scale.canonicalTime;
    peak(index) = quadratic.value;
    peakX(index) = quadratic.x;
    curvature(index) = quadratic.curvature;
    wallHalfWidth(index) = wall.widths(1);
    verticalHalfWidth(index) = vertical.widths(1);
    wallWidth(index) = wall.widths(2);
    verticalWidth(index) = vertical.widths(2);
    if isfinite(wallWidth(index)) && wallWidth(index) > 0
        wallProfiles(:,index) = interp1(ops.x,source(1,:)/peak(index), ...
            peakX(index)+wallWidth(index)*wallCoordinates,'linear',NaN);
    end
    if isfinite(verticalWidth(index)) && verticalWidth(index) > 0
        verticalTrace = interp1(ops.x,source',peakX(index),'linear')';
        verticalProfiles(:,index) = interp1(ops.y,verticalTrace/peak(index), ...
            verticalWidth(index)*verticalCoordinates,'linear',NaN);
    end
end
assert(all(isfinite(tau)) && all(diff(tau) > 0),'ipm:AccelShapeTime', ...
    'Canonical sample times must be finite and increasing.');
report = struct('kind','inner_scale_shape_diagnostic_only', ...
    'canonicalTime',tau,'quadraticPeak',peak,'peakX',peakX, ...
    'quadraticCurvature',curvature, ...
    'wallCoreWidth',wallWidth,'verticalCoreWidth',verticalWidth, ...
    'wallHalfMaximumWidth',wallHalfWidth, ...
    'verticalHalfMaximumWidth',verticalHalfWidth, ...
    'curvatureRadius',sqrt(peak./abs(curvature)), ...
    'wallCoreLogSlope',log_slope(tau,wallWidth), ...
    'verticalCoreLogSlope',log_slope(tau,verticalWidth), ...
    'wallCoreToHalfWidthRatio',wallWidth./wallHalfWidth, ...
    'verticalCoreToHalfWidthRatio',verticalWidth./verticalHalfWidth, ...
    'wallInnerCoordinates',wallCoordinates,'wallInnerProfiles',wallProfiles, ...
    'verticalInnerCoordinates',verticalCoordinates, ...
    'verticalInnerProfiles',verticalProfiles, ...
    'wallTerminalRelativeShapeChange',terminal_change(wallProfiles), ...
    'verticalTerminalRelativeShapeChange',terminal_change(verticalProfiles), ...
    'interpolation','linear diagnostic only; no discretization order claim', ...
    'interpretation','Contraction or residual rejection may reflect a nonsmooth or multiscale limit.');
end

function slope = log_slope(tau,width)
valid = isfinite(width) & width > 0;
if nnz(valid) < 2
    slope = NaN;
else
    centered = tau(valid)-mean(tau(valid));
    fit = [ones(nnz(valid),1),centered(:)]\log(width(valid))';
    slope = fit(2);
end
end

function relative = terminal_change(profiles)
previous = profiles(:,end-1);
current = profiles(:,end);
valid = isfinite(previous) & isfinite(current);
if nnz(valid) < 10
    relative = NaN;
else
    relative = norm(current(valid)-previous(valid))/max(norm(current(valid)),realmin);
end
end
