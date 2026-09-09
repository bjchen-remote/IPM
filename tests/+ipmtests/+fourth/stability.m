function report = stability(mode)
%IPMTESTS.FOURTH.STABILITY Empirical SSPRK54 stability certificate.
%   REPORT = IPMTESTS.FOURTH.STABILITY() runs the quick, required
%   gates. They include SSPRK54 CFL sweeps for one- and
%   two-dimensional smooth advection, boundedness/TV stress tests for a
%   top hat and a localized high-frequency perturbation, and a 320-step
%   full physical IPM evolution using
%
%       high_order + weno5_fd + ssprk54.
%
%   Both modes require a smoothly mapped-grid probe at the recommended CFL.
%   REPORT = IPMTESTS.FOURTH.STABILITY('heavy') refines the transport
%   grids and advances the IPM case for 480 steps.
%
%   The recommended production CFL is 80 percent of the largest contiguous
%   all-test passing CFL in the empirical sweep, capped by the public solver
%   option limit of 0.5.  Consequently it has at least a 20 percent margin
%   relative to the observed passing edge (and sometimes a larger margin).

if nargin < 1 || isempty(mode)
    mode = 'quick';
end
mode = validatestring(lower(char(string(mode))),{'quick','heavy'}, ...
    mfilename,'mode',1);
settings = stability_settings(mode);

fprintf('Running high-order stability verification (%s mode)...\n',mode);
fprintf('  SSPRK54 CFL candidates: %s\n', ...
    sprintf('%.2f ',settings.cflCandidates));

report = struct('mode',mode,'settings',settings);
report.boundaryContracts = struct( ...
    'advection','outflow_only_linear_expansion', ...
    'oneWayZeroSplit','single_step_extrapolation_regression', ...
    'physicalIpm','closed_impermeable_dirichlet_streamfunction');
report.algebraic.oneWayZeroSplit = verify_one_way_zero_split(settings);
report.advection.oneDimensional = sweep_smooth_1d(settings);
report.advection.twoDimensional = sweep_smooth_2d(settings);
report.nonsmooth.topHat = sweep_top_hat(settings);
report.nonsmooth.highFrequency = sweep_high_frequency(settings);

allPass = report.advection.oneDimensional.pass & ...
    report.advection.twoDimensional.pass & ...
    report.nonsmooth.topHat.pass & ...
    report.nonsmooth.highFrequency.pass;
[certifiedIndex,firstFailedIndex] = contiguous_passing_prefix(allPass);
if certifiedIndex >= 1
    certifiedCfl = settings.cflCandidates(certifiedIndex);
    empiricalMarginCfl = settings.marginFactor*certifiedCfl;
    recommendedCfl = min(settings.publicCflMaximum,empiricalMarginCfl);
else
    certifiedCfl = NaN;
    empiricalMarginCfl = NaN;
    % Continue through the expensive physical-path diagnostic even when a
    % local transport gate fails, so the report identifies whether the
    % defect also blocks the full solver.  The function still fails below.
    recommendedCfl = settings.diagnosticFallbackCfl;
end
report.cfl = struct('candidates',settings.cflCandidates, ...
    'allPass',allPass,'certifiedContiguousMaximum',certifiedCfl, ...
    'firstFailed',candidate_or_nan(settings.cflCandidates, ...
        firstFailedIndex), ...
    'marginFactor',settings.marginFactor, ...
    'empiricalMarginRecommendation',empiricalMarginCfl, ...
    'publicSolverMaximum',settings.publicCflMaximum, ...
    'recommended',recommendedCfl);

if certifiedIndex >= 1
    assert(recommendedCfl <= settings.marginFactor*certifiedCfl+10*eps, ...
        'ipm:StabilityCflMargin', ...
        'The recommended CFL does not retain the required 20%% margin.');
    assert(all(allPass(1:certifiedIndex)), ...
        'ipm:StabilityNonmonotonePrefix', ...
        'The reported certified CFL range contains a failed candidate.');
    fprintf(['  contiguous all-profile CFL certificate: <= %.3f; ' ...
        'recommended production CFL: %.3f.\n'], ...
        certifiedCfl,recommendedCfl);
else
    fprintf(['  no CFL certificate; continuing full-IPM diagnostic at ' ...
        'CFL %.3f.\n'],recommendedCfl);
end

report.physicalIpm = verify_physical_ipm(settings,recommendedCfl);
report.stretched = probe_stretched_advection(settings,recommendedCfl);

assert(certifiedIndex >= 1,'ipm:StabilityNoCertifiedCfl', ...
    ['No CFL candidate passed every smooth-advection and nonsmooth ' ...
    'stability gate.']);
assert(strcmp(report.stretched.status,'passed'), ...
    'ipm:StabilityMappedGrid', ...
    'The required smoothly mapped-grid stability probe did not pass.');
report.passed = true;
fprintf(['High-order stability verification passed: %d-step physical ' ...
    'IPM, max range/mass/C2/C3 drifts %.3e/%.3e/%.3e/%.3e.\n'], ...
    report.physicalIpm.steps,report.physicalIpm.maxRangeViolation, ...
    report.physicalIpm.maxMassDrift, ...
    report.physicalIpm.maxCasimirDrift(1), ...
    report.physicalIpm.maxCasimirDrift(2));
end

function settings = stability_settings(mode)
settings = struct();
settings.cflCandidates = [0.20,0.35,0.50,0.65,0.80,1.00,1.20,1.40, ...
    1.60,1.80,2.00,2.20,2.40,2.80,3.20];
settings.marginFactor = 0.8;
settings.publicCflMaximum = 0.5;
settings.diagnosticFallbackCfl = 0.4;
settings.smoothRangeTolerance = 5e-2;
settings.smoothL2Tolerance = 5e-2;
settings.smoothNormRatioTolerance = 0.10;
settings.nonsmoothRangeTolerance = 2e-2;
settings.nonsmoothTvTolerance = 1.05;
settings.ipmRangeTolerance = 5e-3;
settings.ipmMassTolerance = 2e-11;
settings.ipmCasimirTolerance = [2e-2,3e-2];

if strcmp(mode,'quick')
    settings.oneDimensionalNodes = 129;
    settings.twoDimensionalNodes = [65,49];
    settings.nonsmoothNodes = 161;
    settings.ipmNodes = [41,33];
    settings.ipmSteps = 320;
else
    settings.oneDimensionalNodes = 257;
    settings.twoDimensionalNodes = [97,73];
    settings.nonsmoothNodes = 257;
    settings.ipmNodes = [41,33];
    settings.ipmSteps = 480;
end
settings.smoothFinalTime = 0.25;
settings.nonsmoothFinalTime = 0.30;
settings.ipmMaxDt = 5e-3;
end

function result = verify_one_way_zero_split(settings)
% A strictly positive constant velocity makes the negative LF split zero
% on every physical node.  This catches scale floors that allow roundoff in
% extrapolated ghosts to overflow the normalized WENO smoothness indicators.
x = linspace(-1,1,65);
dx = x(2)-x(1);
velocity = 0.8;
initial = 0.2+0.7*exp(-20*(x+0.25).^2);
rhs = @(q)-ipm.field.weno5FluxDerivative( ...
    q,velocity,ones(size(x)),dx,weno_options());
initialRate = rhs(initial);
dt = min(settings.cflCandidates)*dx/velocity;
[advanced,status] = advance_ssprk54(initial,rhs,dt,dt);
result = struct('finiteRate',all(isfinite(initialRate),'all'), ...
    'finiteStep',status.finite,'cfl',min(settings.cflCandidates), ...
    'steps',status.steps,'rangeViolation', ...
    range_violation(advanced,initial),'message',status.message);
assert(result.finiteRate && result.finiteStep, ...
    'ipm:StabilityOneWayZeroSplit', ...
    'The one-way LF zero-split regression produced non-finite values.');
end

function result = sweep_smooth_1d(settings)
x = linspace(-1,1,settings.oneDimensionalNodes);
dx = x(2)-x(1);
velocity = x;
baseline = 0.15;
initial = baseline+0.75*exp(-35*(x+0.35).^2);
departure = x*exp(-settings.smoothFinalTime);
exact = baseline+0.75*exp(-35*(departure+0.35).^2);
metric = ones(size(x));
rhs = @(q)advective_weno_rhs_1d(q,velocity,metric,dx);

count = numel(settings.cflCandidates);
result = empty_sweep(settings.cflCandidates,count);
result.l2Error = nan(1,count);
result.linfError = nan(1,count);
result.rangeViolation = nan(1,count);
result.normRatioToExact = nan(1,count);
exactNorm = sqrt(mean((exact-baseline).^2));
for index = 1:count
    cfl = settings.cflCandidates(index);
    dtTarget = cfl/(max(abs(velocity))/dx);
    [evolved,status] = advance_ssprk54( ...
        initial,rhs,settings.smoothFinalTime,dtTarget);
    result.steps(index) = status.steps;
    result.message{index} = status.message;
    if ~status.finite
        continue;
    end
    error = evolved-exact;
    result.l2Error(index) = sqrt(mean(error.^2));
    result.linfError(index) = max(abs(error),[],'all');
    result.rangeViolation(index) = range_violation(evolved,initial);
    result.normRatioToExact(index) = ...
        sqrt(mean((evolved-baseline).^2))/max(exactNorm,eps);
    result.pass(index) = result.l2Error(index) <= ...
        settings.smoothL2Tolerance && ...
        result.rangeViolation(index) <= settings.smoothRangeTolerance && ...
        abs(result.normRatioToExact(index)-1) <= ...
        settings.smoothNormRatioTolerance;
end
print_sweep_summary('1-D smooth',result);
end

function result = sweep_smooth_2d(settings)
nx = settings.twoDimensionalNodes(1);
ny = settings.twoDimensionalNodes(2);
x = linspace(-1,1,nx);
y = linspace(-1,1,ny)';
dx = x(2)-x(1);
dy = y(2)-y(1);
[X,Y] = meshgrid(x,y);
growthRate = [0.70,0.45];
velocityX = growthRate(1)*X;
velocityY = growthRate(2)*Y;
baseline = 0.10;
initial = baseline+0.80*exp(-22*((X+0.35).^2+(Y+0.28).^2));
exact = baseline+0.80*exp(-22*( ...
    (X*exp(-growthRate(1)*settings.smoothFinalTime)+0.35).^2 + ...
    (Y*exp(-growthRate(2)*settings.smoothFinalTime)+0.28).^2));
metricX = ones(1,nx);
metricY = ones(1,ny);
rhs = @(q)advective_weno_rhs_2d( ...
    q,velocityX,velocityY,metricX,metricY,dx,dy);
rate = max(abs(velocityX)./dx+abs(velocityY)./dy,[],'all');

count = numel(settings.cflCandidates);
result = empty_sweep(settings.cflCandidates,count);
result.l2Error = nan(1,count);
result.linfError = nan(1,count);
result.rangeViolation = nan(1,count);
result.normRatioToExact = nan(1,count);
exactNorm = sqrt(mean((exact-baseline).^2,'all'));
for index = 1:count
    cfl = settings.cflCandidates(index);
    [evolved,status] = advance_ssprk54( ...
        initial,rhs,settings.smoothFinalTime,cfl/rate);
    result.steps(index) = status.steps;
    result.message{index} = status.message;
    if ~status.finite
        continue;
    end
    error = evolved-exact;
    result.l2Error(index) = sqrt(mean(error.^2,'all'));
    result.linfError(index) = max(abs(error),[],'all');
    result.rangeViolation(index) = range_violation(evolved,initial);
    result.normRatioToExact(index) = ...
        sqrt(mean((evolved-baseline).^2,'all'))/max(exactNorm,eps);
    result.pass(index) = result.l2Error(index) <= ...
        settings.smoothL2Tolerance && ...
        result.rangeViolation(index) <= settings.smoothRangeTolerance && ...
        abs(result.normRatioToExact(index)-1) <= ...
        settings.smoothNormRatioTolerance;
end
print_sweep_summary('2-D smooth',result);
end

function result = sweep_top_hat(settings)
x = linspace(-1,1,settings.nonsmoothNodes);
dx = x(2)-x(1);
velocity = x;
initial = double(x >= -0.62 & x <= -0.18);
rhs = @(q)advective_weno_rhs_1d(q,velocity,ones(size(x)),dx);
result = nonsmooth_sweep(initial,rhs,settings, ...
    settings.cflCandidates/(max(abs(velocity))/dx));
print_sweep_summary('top hat',result);
end

function result = sweep_high_frequency(settings)
x = linspace(-1,1,settings.nonsmoothNodes);
dx = x(2)-x(1);
velocity = x;
envelope = exp(-((x+0.34)/0.42).^8);
initial = 0.5+0.44*envelope.*sin(22*pi*(x+0.34));
rhs = @(q)advective_weno_rhs_1d(q,velocity,ones(size(x)),dx);
result = nonsmooth_sweep(initial,rhs,settings, ...
    settings.cflCandidates/(max(abs(velocity))/dx));
print_sweep_summary('high frequency',result);
end

function result = nonsmooth_sweep(initial,rhs,settings,dtTargets)
count = numel(settings.cflCandidates);
result = empty_sweep(settings.cflCandidates,count);
result.rangeViolation = nan(1,count);
result.totalVariation = nan(1,count);
result.tvRatio = nan(1,count);
initialTv = total_variation(initial);
for index = 1:count
    [evolved,status] = advance_ssprk54(initial,rhs, ...
        settings.nonsmoothFinalTime,dtTargets(index));
    result.steps(index) = status.steps;
    result.message{index} = status.message;
    if ~status.finite
        continue;
    end
    result.rangeViolation(index) = range_violation(evolved,initial);
    result.totalVariation(index) = total_variation(evolved);
    result.tvRatio(index) = result.totalVariation(index)/max(initialTv,eps);
    result.pass(index) = result.rangeViolation(index) <= ...
        settings.nonsmoothRangeTolerance && ...
        result.tvRatio(index) <= settings.nonsmoothTvTolerance;
end
end

function rhs = advective_weno_rhs_1d(q,velocity,metric,spacing)
fluxDerivative = ipm.field.weno5FluxDerivative( ...
    q,velocity,metric,spacing,weno_options());
unit = ones(size(q),'like',q);
velocityDivergence = ipm.field.weno5FluxDerivative( ...
    unit,velocity,metric,spacing,weno_options());
rhs = -(fluxDerivative-q.*velocityDivergence);
end

function rhs = advective_weno_rhs_2d( ...
        q,velocityX,velocityY,metricX,metricY,spacingX,spacingY)
unit = ones(size(q),'like',q);
fluxX = ipm.field.weno5FluxDerivative( ...
    q,velocityX,metricX,spacingX,weno_options());
fluxY = ipm.field.weno5FluxDerivative( ...
    q',velocityY',metricY,spacingY,weno_options())';
divergenceX = ipm.field.weno5FluxDerivative( ...
    unit,velocityX,metricX,spacingX,weno_options());
divergenceY = ipm.field.weno5FluxDerivative( ...
    unit',velocityY',metricY,spacingY,weno_options())';
rhs = -(fluxX+fluxY-q.*(divergenceX+divergenceY));
end

function result = verify_physical_ipm(settings,recommendedCfl)
steps = settings.ipmSteps;
maxDt = settings.ipmMaxDt;
opts = struct( ...
    'nx',settings.ipmNodes(1),'ny',settings.ipmNodes(2), ...
    'xlim',[-4,4],'ymax',4, ...
    'gridMode','uniform','gridStretchAutomatic',false, ...
    'gridStretch',[0,0], ...
    'finalTime',2*steps*maxDt,'physicalFinalTime',Inf, ...
    'cfl',recommendedCfl,'maxDt',maxDt,'minDt',1e-12, ...
    'maxSteps',steps,'outputEvery',1e-12, ...
    'timeIntegrator','ssprk54', ...
    'initialCondition','smooth_blob', ...
    'rescalingMode','physical', ...
    'farBoundaryMode','dirichlet_zero', ...
    'transportBoundaryMode','closed', ...
    'transportScheme','weno5_fd', ...
    'spatialDiscretization','high_order', ...
    'wallTransportMode','conservative_flux', ...
    'symmetryMode','double_odd_omega', ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'gradientStop',1e12, ...
    'rangeStopTolerance',settings.ipmRangeTolerance, ...
    'oscillationTVTolerance',1e6, ...
    'positiveWallNegativeTolerance',Inf, ...
    'storeSnapshots',true,'saveResults',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'verbose',false);
solverResult = ipm.solve(opts);

assert(solverResult.state.steps == steps, ...
    'ipm:StabilityIpmStepCount', ...
    'Physical IPM stopped after %d of %d requested steps (%s).', ...
    solverResult.state.steps,steps,solverResult.state.stopReason);
assert(strcmp(solverResult.state.stopReason,'maximum_steps'), ...
    'ipm:StabilityIpmEarlyStop', ...
    'Physical IPM stopped unexpectedly with reason ''%s''.', ...
    solverResult.state.stopReason);

snapshots = solverResult.snapshots.rho;
assert(numel(snapshots) == steps+1, ...
    'ipm:StabilityIpmSnapshotCount', ...
    'Expected one snapshot per IPM step plus the initial state.');
weights = ipm.mesh.quadrature(solverResult.grid.y) * ...
    ipm.mesh.quadrature(solverResult.grid.x);
milestones = unique([0,round(steps*[1/8,1/4,1/2,3/4]),steps]);
snapshotCount = numel(snapshots);
finite = false(1,snapshotCount);
rhoRange = nan(snapshotCount,2);
rangeViolation = nan(1,snapshotCount);
mass = nan(1,snapshotCount);
casimirPowers = [2,3];
casimir = nan(snapshotCount,numel(casimirPowers));
initial = snapshots{1};
for snapshot = 1:snapshotCount
    rho = snapshots{snapshot};
    finite(snapshot) = all(isfinite(rho),'all');
    if ~finite(snapshot)
        continue;
    end
    rhoRange(snapshot,:) = [min(rho,[],'all'),max(rho,[],'all')];
    rangeViolation(snapshot) = range_violation(rho,initial);
    mass(snapshot) = sum(weights.*rho,'all');
    for powerIndex = 1:numel(casimirPowers)
        casimir(snapshot,powerIndex) = sum( ...
            weights.*rho.^casimirPowers(powerIndex),'all');
    end
end
massDrift = relative_drift(mass);
casimirDrift = relative_drift(casimir);
milestoneIndices = milestones+1;
physicalTime = solverResult.snapshots.physicalTime(:)';

% Also expose the lifecycle's control-volume diagnostic.  The certified
% invariant uses the high-order quadrature paired with this discretization.
historyMassDrift = solverResult.history.common.massDrift;
result = struct('steps',steps,'stopReason',solverResult.state.stopReason, ...
    'cfl',recommendedCfl,'maxDt',maxDt, ...
    'timeStepRange',[min(diff(physicalTime)),max(diff(physicalTime))], ...
    'gridSize',settings.ipmNodes, ...
    'milestones',milestones,'milestoneTimes',physicalTime(milestoneIndices), ...
    'finite',finite(milestoneIndices), ...
    'rhoRange',rhoRange(milestoneIndices,:), ...
    'rangeViolation',rangeViolation(milestoneIndices), ...
    'mass',mass(milestoneIndices), ...
    'massDrift',massDrift(milestoneIndices), ...
    'casimirPowers',casimirPowers, ...
    'casimir',casimir(milestoneIndices,:), ...
    'casimirDrift',casimirDrift(milestoneIndices,:), ...
    'stepSeries',struct('step',0:steps,'physicalTime',physicalTime, ...
        'finite',finite,'rhoRange',rhoRange, ...
        'rangeViolation',rangeViolation,'mass',mass, ...
        'massDrift',massDrift,'casimir',casimir, ...
        'casimirDrift',casimirDrift), ...
    'maxRangeViolation',max(rangeViolation), ...
    'maxMassDrift',max(abs(massDrift)), ...
    'maxCasimirDrift',max(abs(casimirDrift),[],1), ...
    'controlVolumeMassDriftMaximum',max(abs(historyMassDrift)), ...
    'terminalTime',solverResult.state.physicalTime);

assert(all(finite),'ipm:StabilityIpmFinite', ...
    'The full physical IPM run produced a non-finite accepted state.');
assert(result.maxRangeViolation <= settings.ipmRangeTolerance, ...
    'ipm:StabilityIpmRange', ...
    'Physical IPM range violation %.3e exceeds %.3e.', ...
    result.maxRangeViolation,settings.ipmRangeTolerance);
assert(result.maxMassDrift <= settings.ipmMassTolerance, ...
    'ipm:StabilityIpmMass', ...
    'Physical IPM mass drift %.3e exceeds %.3e.', ...
    result.maxMassDrift,settings.ipmMassTolerance);
assert(all(result.maxCasimirDrift <= settings.ipmCasimirTolerance), ...
    'ipm:StabilityIpmCasimir', ...
    'Physical IPM C2/C3 drift exceeds tolerance: %.3e/%.3e.', ...
    result.maxCasimirDrift(1),result.maxCasimirDrift(2));
end

function result = probe_stretched_advection(settings,cfl)
stretch = [1.05,0.85];
nx = settings.twoDimensionalNodes(1);
ny = settings.twoDimensionalNodes(2);
sx = linspace(-1,1,nx);
sy = linspace(-1,1,ny)';
[x,metricX] = sinh_axis(sx,stretch(1));
[y,metricY] = sinh_axis(sy,stretch(2));
[X,Y] = meshgrid(x,y);
growthRate = [0.70,0.45];
velocityX = growthRate(1)*X;
velocityY = growthRate(2)*Y;
baseline = 0.10;
initial = baseline+0.80*exp(-22*((X+0.35).^2+(Y+0.28).^2));
exact = baseline+0.80*exp(-22*( ...
    (X*exp(-growthRate(1)*settings.smoothFinalTime)+0.35).^2 + ...
    (Y*exp(-growthRate(2)*settings.smoothFinalTime)+0.28).^2));
dsx = sx(2)-sx(1);
dsy = sy(2)-sy(1);
rhs = @(q)advective_weno_rhs_2d( ...
    q,velocityX,velocityY,metricX,metricY',dsx,dsy);
rate = max(abs(velocityX)./(metricX*dsx) + ...
    abs(velocityY)./(metricY*dsy),[],'all');
dtTarget = cfl/rate;
try
    [evolved,status] = advance_ssprk54( ...
        initial,rhs,settings.smoothFinalTime,dtTarget);
    l2Error = sqrt(mean((evolved-exact).^2,'all'));
    rangeError = range_violation(evolved,initial);
    passed = status.finite && l2Error <= settings.smoothL2Tolerance && ...
        rangeError <= settings.smoothRangeTolerance;
    if passed
        label = 'passed';
    else
        label = 'experimental_failure';
    end
    result = struct('status',label,'requiredForCertificate',true, ...
        'cfl',cfl,'stretch',stretch,'steps',status.steps, ...
        'finite',status.finite,'l2Error',l2Error, ...
        'rangeViolation',rangeError,'message',status.message);
catch exception
    result = struct('status','experimental_error', ...
        'requiredForCertificate',true,'cfl',cfl,'stretch',stretch, ...
        'identifier',exception.identifier,'message',exception.message);
end
fprintf('  stretched mapped-WENO required probe: %s.\n',result.status);
end

function result = empty_sweep(candidates,count)
result = struct('cfl',candidates,'pass',false(1,count), ...
    'steps',zeros(1,count),'message',{repmat({''},1,count)});
end

function [state,status] = advance_ssprk54(initial,rhs,finalTime,dtTarget)
state = initial;
time = 0;
steps = 0;
finite = true;
message = '';
while time < finalTime
    dt = min(dtTarget,finalTime-time);
    try
        state = explicit_ssprk54_step(state,dt,rhs);
    catch exception
        finite = false;
        message = sprintf('%s: %s',exception.identifier,exception.message);
        break;
    end
    steps = steps+1;
    if any(~isfinite(state),'all') || max(abs(state),[],'all') > 1e8
        finite = false;
        message = 'non-finite or greater-than-1e8 state magnitude';
        break;
    end
    time = time+dt;
end
status = struct('finite',finite,'steps',steps,'time',time, ...
    'message',message);
end

function stateNew = explicit_ssprk54_step(state,dt,rhs)
tableau = ipm.evolve.ssprk54Tableau();
rates = cell(tableau.stages,1);
for stage = 1:tableau.stages
    stageState = state;
    for previous = 1:stage-1
        stageState = stageState + ...
            dt*tableau.A(stage,previous)*rates{previous};
    end
    rates{stage} = rhs(stageState);
end
stateNew = state;
for stage = 1:tableau.stages
    stateNew = stateNew+dt*tableau.b(stage)*rates{stage};
end
end

function options = weno_options()
options = struct('lowerBoundary','extrapolate', ...
    'upperBoundary','extrapolate','extrapolationDegree',5, ...
    'epsilon',1e-12);
end

function value = range_violation(state,reference)
value = max([min(reference,[],'all')-min(state,[],'all'), ...
    max(state,[],'all')-max(reference,[],'all'),0]);
end

function value = total_variation(state)
value = sum(abs(diff(state(:))));
end

function drift = relative_drift(values)
if isvector(values)
    values = values(:);
end
reference = values(1,:);
drift = (values-reference)./max(abs(reference),eps);
if size(drift,2) == 1
    drift = drift.';
end
end

function [lastPassing,firstFailed] = contiguous_passing_prefix(pass)
firstFailed = find(~pass,1,'first');
if isempty(firstFailed)
    lastPassing = numel(pass);
    firstFailed = [];
else
    lastPassing = firstFailed-1;
end
end

function value = candidate_or_nan(candidates,index)
if isempty(index)
    value = NaN;
else
    value = candidates(index);
end
end

function print_sweep_summary(name,result)
passing = result.cfl(result.pass);
if isempty(passing)
    lastPassing = NaN;
else
    lastPassing = max(passing);
end
fprintf('  %-16s largest passing candidate: %.3f.\n',name,lastPassing);
if isnan(lastPassing)
    fprintf('    pass mask: %s\n',sprintf('%d ',result.pass));
    fprintf('    completed steps: %s\n',sprintf('%d ',result.steps));
    if isfield(result,'l2Error')
        fprintf('    L2 errors: %s\n',sprintf('%.3e ',result.l2Error));
    end
    if isfield(result,'rangeViolation')
        fprintf('    range violations: %s\n', ...
            sprintf('%.3e ',result.rangeViolation));
    end
    if isfield(result,'tvRatio')
        fprintf('    TV ratios: %s\n',sprintf('%.3e ',result.tvRatio));
    end
    if ~isempty(result.message{1})
        fprintf('    first-candidate error: %s\n',result.message{1});
    end
end
end

function [mapped,metric] = sinh_axis(reference,stretch)
mapped = sinh(stretch*reference)/sinh(stretch);
metric = stretch*cosh(stretch*reference)/sinh(stretch);
end
