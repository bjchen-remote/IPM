function report = equivalence(originalRoot)
%EQUIVALENCE Compare unchanged physical kernels with the preserved original.
%   Dynamic normalization deliberately differs under the schema-4 exact-gauge
%   contract and is covered by the scaling suite instead.  Physical result
%   projections and remesh kernels remain exact ISEQUALN comparisons.
%   No results, plots, or other output artifacts are written.

copyRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(originalRoot)
    originalRoot = fullfile(fileparts(copyRoot),'sixth_order_integration');
end
savedPath = path;
savedDirectory = pwd;
cleanup = onCleanup(@()restore_environment(savedPath,savedDirectory));
assert(isfolder(originalRoot),'ipm:EquivalenceOriginalMissing', ...
    'The original solver folder does not exist: %s.',originalRoot);
cd(originalRoot);
originalRoot = pwd;
addpath(originalRoot,'-begin');
cd(copyRoot);
assert(strcmp(which('main_ipm'),fullfile(originalRoot,'main_ipm.m')), ...
    'ipm:EquivalenceOriginalPath', ...
    'main_ipm must resolve to the explicitly selected original folder.');
assert(strcmp(which('ipm.solve'),fullfile(copyRoot,'+ipm','solve.m')), ...
    'ipm:EquivalenceCopyPath', ...
    'ipm.solve must resolve to this independent copy.');
assert(~strcmp(originalRoot,copyRoot),'ipm:EquivalenceSameSource', ...
    'The original and copied solver must occupy different folders.');

fprintf('Running legacy physical-kernel equivalence verification...\n');
legacyActive = ipm_active_case();
legacyActive = rmfield(legacyActive,{'lengthScaleGain','maxDynamicRate'});
assert_exact(legacyActive,ipm.config.activeCase(),'activeCase.stableOptions');
assert_exact(ipm_sixth_order_options(),ipm.config.sixthOrder(), ...
    'sixthOrderOptions');
[names,cases] = run_cases();
excludedOperatorDiagnostics = operator_grid_quality_diagnostic_fields();
addedCommonTelemetry = added_common_telemetry_fields();
report = struct('originalRoot',originalRoot,'copyRoot',copyRoot, ...
    'scope','physical_kernel_projection', ...
    'dynamicContract','intentionally_changed_schema4', ...
    'excludedResultFields',{{'metadata','config','fit', ...
        'state.canonicalTime','history.common.canonicalTau', ...
        'snapshots.canonicalTime'}}, ...
    'excludedOperatorGridQualityFields',{excludedOperatorDiagnostics}, ...
    'addedCommonTelemetryFields',{addedCommonTelemetry}, ...
    'addedTelemetryValidation','exact_reconstruction_from_legacy_snapshots', ...
    'projectionNegativeChecks',0, ...
    'runs',struct('name',{},'gridSize',{},'steps',{},'exact',{}), ...
    'telemetry',struct('name',{},'records',{},'fields',{},'exact',{}), ...
    'remesh',struct('order',{},'gridChange',{}, ...
        'relativeMassDefect',{},'exact',{}),'passed',false);
for index = 1:numel(cases)
    options = cases{index};
    original = main_ipm(options);
    copied = ipm.solve(options);
    assert_config_kernel(original.config,copied.config, ...
        [names{index},'.config']);
    [originalKernel,copiedKernel] = paired_result_kernels( ...
        original,copied,names{index});
    assert_exact(originalKernel,copiedKernel,names{index});
    expectedTelemetry = legacy_derived_telemetry(original,copied);
    assert_exact(expectedTelemetry,select_telemetry(copied), ...
        [names{index},'.addedCommonTelemetry']);
    report.telemetry(index) = struct('name',names{index}, ...
        'records',numel(original.history.common.t), ...
        'fields',numel(addedCommonTelemetry),'exact',true);
    if index == 1
        report.projectionNegativeChecks = projection_negative_checks( ...
            original,copied,expectedTelemetry);
    end
    % The last accepted physical step can cross the requested terminal time.
    % This completion bound never relaxes the exact result comparison.
    assert(strcmp(original.state.stopReason,'physical_final_time') && ...
        original.state.steps >= 2 && ...
        original.state.physicalTime >= ...
            options.physicalFinalTime-options.minDt && ...
        original.state.physicalTime <= ...
            options.physicalFinalTime+options.maxDt && ...
        numel(original.snapshots.rho) >= 2, ...
        'ipm:EquivalenceIncompleteRun', ...
        '%s stopped at t=%.17g (%s, %d steps, %d snapshots).', ...
        names{index},original.state.physicalTime, ...
        original.state.stopReason,original.state.steps, ...
        numel(original.snapshots.rho));
    report.runs(index) = struct('name',names{index}, ...
        'gridSize',[options.nx,options.ny], ...
        'steps',original.state.steps,'exact',true);
    fprintf('  %-30s exact (%d steps).\n',names{index},original.state.steps);
end
for order = [2,4,6]
    report.remesh(end+1) = compare_remesh(order);
end
report.passed = true;
fprintf(['Physical-kernel equivalence passed: %d runs and %d remesh ' ...
    'transfers.\n'], ...
    numel(report.runs),numel(report.remesh));
end

function [names,cases] = run_cases()
orders = [2,4,6];
symmetries = {'double_odd_omega','half_plane'};
suffixes = {'DoubleOdd','HalfPlane'};
names = cell(1,numel(orders)*numel(symmetries));
cases = cell(size(names));
cursor = 0;
for order = orders
    if order == 2
        options = short_options();
        prefix = 'baseline';
    else
        options = numerical_options(order);
        prefix = sprintf('order%d',order);
    end
    for symmetryIndex = 1:numel(symmetries)
        cursor = cursor+1;
        candidate = options;
        candidate.rescalingMode = 'physical';
        candidate.symmetryMode = symmetries{symmetryIndex};
        names{cursor} = [prefix,'Physical',suffixes{symmetryIndex}];
        cases{cursor} = candidate;
    end
end
end

function options = short_options()
options = struct('nx',33,'ny',17,'xlim',[-4,4],'ymax',4, ...
    'gridMode','stretched','gridStretch',[1.5,1.5], ...
    'gridStretchAutomatic',false,'targetCenterSpacing',[0.1,0.1], ...
    'initialCondition','degenerate','degeneratePower',4, ...
    'farBoundaryMode','dirichlet_zero','transportBoundaryMode','open', ...
    'transportScheme','muscl_minmod', ...
    'physicalFinalTime',2e-4,'maxDt',1e-4,'minDt',1e-12, ...
    'maxSteps',8,'outputEvery',1e-4,'cfl',0.25, ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'gradientStop',1e9,'rangeStopTolerance',0.1, ...
    'oscillationTVTolerance',10,'positiveWallNegativeTolerance',0.99, ...
    'storeSnapshots',true,'saveResults',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
end

function options = numerical_options(order)
options = short_options();
options.nx = 25;
options.ny = 17;
options.xlim = [-2,2];
options.ymax = 2;
options.gridStretch = [0.7,0.5];
options.initialCondition = 'degenerate_primitive';
options.wallTransportMode = 'conservative_flux';
if order == 4
    options.spatialDiscretization = 'high_order';
    options.transportScheme = 'weno5_fd';
    options.timeIntegrator = 'ssprk54';
    options.remeshTransferScheme = 'high_order';
elseif order == 6
    options.transportScheme = 'weno7_fd';
    options = ipm_sixth_order_options(options);
end
end

function report = compare_remesh(order)
options = numerical_options(order);
options.nx = 25;
options.ny = 25;
options.xlim = [-1,1];
options.ymax = 1;
options.gridStretch = [1,0.75];
options.rescalingMode = 'physical';
options.symmetryMode = 'half_plane';
options.transportBoundaryMode = 'closed';
originalConfig = ipm_resolve_config(options);
copiedConfig = ipm.config.resolve(options);
assert_config_kernel(originalConfig,copiedConfig,'remesh.config');
originalOps = ipm_build_operators_resolved(originalConfig);
copiedOps = ipm.mesh.build(copiedConfig);
assert_exact(operator_kernel(originalOps), ...
    operator_kernel(copiedOps),'remesh.initialOperators');
% A transfer needs the frozen pin, but no dynamic normalization is active.
originalOps.rescaling.pinX = 0;
copiedOps.rescaling.pinX = 0;
rho = 2+exp(0.35*originalOps.X).*(1+0.08*cos(1.1*originalOps.Y));
sx = linspace(-1,1,25);
sy = linspace(0,1,25)';
proposal = struct('x',sinh(0.65*sx)/sinh(0.65), ...
    'y',sinh(0.50*sy)/sinh(0.50));
[originalRho,originalNew,originalMetrics] = ipm_remesh_transfer( ...
    rho,originalOps,originalConfig,proposal);
[copiedRho,copiedNew,copiedMetrics] = ipm.remesh.transfer( ...
    rho,copiedOps,copiedConfig,proposal);
label = sprintf('order%dRemesh',order);
assert_exact(originalRho,copiedRho,[label,'.rho']);
assert_exact(originalMetrics,copiedMetrics,[label,'.metrics']);
assert_exact(operator_kernel(originalNew), ...
    operator_kernel(copiedNew),[label,'.operators']);
% Decomposition objects are checked through their actual solve action.
probe = sin((1:size(originalNew.A,1))');
assert_exact(originalNew.poisson\probe,copiedNew.poisson\probe, ...
    [label,'.factorizedSolve']);
gridChange = [max(abs(originalNew.x-originalOps.x)), ...
    max(abs(originalNew.y-originalOps.y))];
assert(all(gridChange > 1e-3) && ...
    all(isfinite(originalRho),'all') && ~isequaln(rho,originalRho), ...
    'ipm:EquivalenceTrivialRemesh', ...
    '%s must change both grid axes and the transferred samples.',label);
oldMass = sum(rho.*originalOps.integrationWeights,'all');
newMass = sum(originalRho.*originalNew.integrationWeights,'all');
relativeMassDefect = abs(newMass-oldMass)/max(abs(oldMass),1);
if order > 2
    assert(relativeMassDefect <= 5e-13,'ipm:EquivalenceRemeshMass', ...
        '%s lost the high-order conservation contract.',label);
end
report = struct('order',order,'gridChange',gridChange, ...
    'relativeMassDefect',relativeMassDefect,'exact',true);
fprintf('  %-30s exact (field, metrics, operators, solve).\n',label);
end

function assert_config_kernel(original,copied,label)
domains = {'grid','time','physics','elliptic','transport','remesh', ...
    'diagnostics','output'};
for index = 1:numel(domains)
    name = domains{index};
    assert_exact(original.(name),copied.(name),[label,'.',name]);
end
scalingNames = {'rescalingMode','dynamicScaleGeometry'};
for index = 1:numel(scalingNames)
    name = scalingNames{index};
    assert_exact(original.scaling.(name),copied.scaling.(name), ...
        [label,'.scaling.',name]);
end
end

function [originalKernel,copiedKernel] = paired_result_kernels(original,copied,label)
% Never use a generic intersection: every old common observable is retained,
% and every newly unmatched field must belong to this exact known addition.
oldNames = fieldnames(original.history.common);
newNames = fieldnames(copied.history.common);
allowed = added_common_telemetry_fields();
assert(isempty(setdiff(oldNames,newNames)) && ...
    isequal(sort(setdiff(newNames,oldNames)),sort(allowed(:))), ...
    'ipm:EquivalenceTelemetrySchema', ...
    '%s must retain every legacy common field and add exactly the 23 audited telemetry fields.',label);
originalKernel = result_kernel(original);
copiedKernel = result_kernel(copied);
copiedKernel.history.common = rmfield(copiedKernel.history.common,allowed);
end

function projected = result_kernel(result)
projected = rmfield(result,{'metadata','config','fit'});
if isfield(projected.history,'wallCore')
    % Compact wall traces are new output-only telemetry and do not alter
    % the legacy trajectory equivalence kernel.
    projected.history = rmfield(projected.history,'wallCore');
end
% Schema 4 advances canonical tau directly as tau+dt; the legacy code sends
% the constant derivative one through each RK tableau.  In physical mode the
% two forms have the same trajectory but can differ in roundoff.  Remove only
% the three copies of that canonical clock.  Physical time, normalized time,
% every scale, field, velocity, physical observable, and snapshot remain in
% the exact ISEQUALN projection.
projected.state = rmfield(projected.state,'canonicalTime');
projected.history.common = rmfield( ...
    projected.history.common,'canonicalTau');
projected.snapshots = rmfield(projected.snapshots,'canonicalTime');
end

function names = added_common_telemetry_fields()
% record.m adds 1 step index, 11 accepted-step observations, 6 physical-grid
% observations and 5 selected-quadratic-peak observations. These names are
% absent from the preserved original; no shared physical field is removed.
names = {'acceptedStep','acceptedCanonicalDt','transportRate','cflDtLimit', ...
    'maximumDtLimit','canonicalEndpointDtLimit','physicalEndpointDtLimit', ...
    'physicalClockSpeed','realizedCfl','timestepActiveLimiter', ...
    'amplitudeSourceStep','conservativeSourceStep', ...
    'physicalDomainXMinimum','physicalDomainXMaximum','physicalDomainXRadius', ...
    'physicalDomainYMaximum','physicalMinimumDx','physicalMinimumDy', ...
    'physicalQuadraticPeak','inversePhysicalQuadraticPeak', ...
    'quadraticInverseSlopeAlgebraic','quadraticLocalTerminalTime', ...
    'rescaledRhoXInfToQuadraticPeak'};
end

function expected = legacy_derived_telemetry(original,copied)
% Independent observation of already-equal legacy fields. In particular,
% rates come from the original flow implementation, not from new telemetry
% or the new selectTimestep/record helpers. No additional trajectory is run.
h = original.history.common;n = numel(h.t);config = original.config;
assert(strcmp(config.scaling.rescalingMode,'physical') && ...
    ~strcmp(config.scaling.cOmegaGauge,'wall_omega_quadratic_peak') && ...
    original.grid.remeshCount == 0 && copied.state.steps == n-1 && ...
    numel(original.snapshots.rho) == n && ...
    isequaln(original.snapshots.normalizedTime,h.t), ...
    'ipm:EquivalenceTelemetryScope', ...
    'Telemetry reconstruction requires every accepted physical step and its unchanged-grid snapshot.');
assert(all(h.C_l == 1 & h.C_x == 1 & h.C_y == 1 & h.C_omega == 1 & ...
    h.X_shift == 0 & h.c_l == 0 & h.c_omega == 0 & h.timeSpeed == 1), ...
    'ipm:EquivalenceTelemetryScope','Physical-mode identity scales/rates must remain exact.');
ops = ipm_build_operators_resolved(config);
scale = ipm_initial_scale(ops);
names = added_common_telemetry_fields();expected = struct();
for j = 1:numel(names),expected.(names{j}) = NaN(n,1);end
expected.acceptedStep = (0:n-1)';
expected.timestepActiveLimiter = strings(n,1);
canonical = zeros(n,1);
limiterNames = {'max_dt','cfl','canonical_final_time','physical_final_time'};
for k = 1:n
    x = original.snapshots.x{k};y = original.snapshots.y{k};
    assert_exact(x,ops.x,'telemetry.legacySnapshotX');
    assert_exact(y,ops.y,'telemetry.legacySnapshotY');
    expected.physicalDomainXMinimum(k) = (x(1)-h.X_shift(k))/h.C_x(k);
    expected.physicalDomainXMaximum(k) = (x(end)-h.X_shift(k))/h.C_x(k);
    expected.physicalDomainXRadius(k) = max(abs([ ...
        expected.physicalDomainXMinimum(k),expected.physicalDomainXMaximum(k)]));
    expected.physicalDomainYMaximum(k) = y(end)/h.C_y(k);
    expected.physicalMinimumDx(k) = min(diff(x))/h.C_x(k);
    expected.physicalMinimumDy(k) = min(diff(y))/h.C_y(k);
    expected.maximumDtLimit(k) = config.time.maxDt;
    expected.physicalClockSpeed(k) = 1;
    if k == 1
        expected.canonicalEndpointDtLimit(k) = config.time.finalTime;
        expected.physicalEndpointDtLimit(k) = config.time.physicalFinalTime;
        expected.timestepActiveLimiter(k) = "initial";
        continue;
    end
    scale.canonicalTime = h.canonicalTau(k-1);
    scale.physicalTime = h.physicalTime(k-1);
    [~,flow] = ipm_flow_at_state(original.snapshots.rho{k-1},ops,scale);
    assert_exact(flow.timeSpeed,1,'telemetry.legacyPhysicalClock');
    if strcmp(ops.spatialDiscretization,'sixth_order')
        rate = max(abs(flow.transportU1)./ops.transportSpacingX + ...
            abs(flow.transportU2)./ops.transportSpacingY,[],'all');
    else
        rate = max(abs(flow.transportU1)./ops.hx + ...
            abs(flow.transportU2)./ops.hy,[],'all');
    end
    if rate > 0,cflLimit = config.time.cfl/rate;else,cflLimit = config.time.maxDt;end
    limits = [config.time.maxDt,cflLimit, ...
        config.time.finalTime-canonical(k-1), ...
        config.time.physicalFinalTime-h.physicalTime(k-1)];
    [dt,limiter] = min(limits);
    [legacyDt,stopReason] = ipm_select_timestep(flow,scale,ops,config.time);
    assert_exact(dt,legacyDt,'telemetry.legacySelectedDt');
    assert(isempty(stopReason),'ipm:EquivalenceTelemetryScope', ...
        'Every recorded accepted step must have an admissible original dt.');
    canonical(k) = canonical(k-1)+dt;
    expected.acceptedCanonicalDt(k) = dt;
    expected.transportRate(k) = rate;
    expected.cflDtLimit(k) = cflLimit;
    expected.canonicalEndpointDtLimit(k) = limits(3);
    expected.physicalEndpointDtLimit(k) = limits(4);
    expected.realizedCfl(k) = dt*rate;
    expected.timestepActiveLimiter(k) = string(limiterNames{limiter});
    expected.amplitudeSourceStep(k) = dt*abs(h.c_omega(k));
    expected.conservativeSourceStep(k) = dt*abs(2*h.c_l(k)+h.c_omega(k));
end
% The previously excluded new canonical aliases must themselves be the exact
% direct accumulation of the independently reconstructed accepted steps.
assert_exact(canonical,copied.history.common.canonicalTau,'telemetry.directCanonicalHistory');
assert_exact(canonical,copied.snapshots.canonicalTime,'telemetry.directCanonicalSnapshots');
assert_exact(canonical(end),copied.state.canonicalTime,'telemetry.directCanonicalState');
end

function selected = select_telemetry(result)
names = added_common_telemetry_fields();selected = struct();
for k = 1:numel(names),selected.(names{k}) = result.history.common.(names{k});end
end

function count = projection_negative_checks(original,copied,expected)
candidate = copied;candidate.history.common.unknownReadOnlyTelemetry = 0;
must_reject(@()paired_result_kernels(original,candidate,'unknownField'), ...
    'ipm:EquivalenceTelemetrySchema');
candidate = copied;candidate.history.common = rmfield(candidate.history.common,'acceptedStep');
must_reject(@()paired_result_kernels(original,candidate,'missingAllowedField'), ...
    'ipm:EquivalenceTelemetrySchema');
candidate = copied;candidate.history.common = rmfield(candidate.history.common,'physicalRhoXInf');
must_reject(@()paired_result_kernels(original,candidate,'missingLegacyField'), ...
    'ipm:EquivalenceTelemetrySchema');
candidate = copied;candidate.history.common.physicalRhoXInf(end) = ...
    candidate.history.common.physicalRhoXInf(end)+1;
[a,b] = paired_result_kernels(original,candidate,'changedLegacyValue');
must_reject(@()assert_exact(a,b,'changedLegacyValue'),'ipm:EquivalenceMismatch');
candidate = copied;candidate.history.common.transportRate(end) = ...
    candidate.history.common.transportRate(end)+1;
must_reject(@()assert_exact(expected,select_telemetry(candidate),'changedNewValue'), ...
    'ipm:EquivalenceMismatch');
count = 5;
end

function must_reject(action,identifier)
try
    action();
catch exception
    assert(strcmp(exception.identifier,identifier),'ipm:EquivalenceNegativeWrongError', ...
        'Expected %s, received %s.',identifier,exception.identifier);
    return;
end
error('ipm:EquivalenceNegativeAccepted','A corrupted equivalence fixture was accepted.');
end

function projected = operator_kernel(ops)
projected = rmfield(ops,{'poisson','rescaling'});
if isfield(projected,'gridQuality')
    diagnosticFields = operator_grid_quality_diagnostic_fields();
    for axisName = {'x','y'}
        name = axisName{1};
        present = diagnosticFields(isfield( ...
            projected.gridQuality.(name),diagnosticFields));
        if ~isempty(present)
            projected.gridQuality.(name) = rmfield( ...
                projected.gridQuality.(name),present);
        end
    end
end
end

function names = operator_grid_quality_diagnostic_fields()
% Added post-migration diagnostics; none participates in operator assembly.
names = { ...
    'maximumQuadratureWeightToControlWidthAbsoluteDeviationFromUnity', ...
    'maximumQuadratureWeightToControlWidthRatio', ...
    'minimumQuadratureWeightToControlWidthRatio', ...
    'quadratureWeightsStrictlyPositive'};
end

function assert_exact(original,copied,label)
if ~isequaln(original,copied)
    error('ipm:EquivalenceMismatch','First unequal field: %s.', ...
        first_difference(original,copied,label));
end
end

function label = first_difference(original,copied,label)
if ~strcmp(class(original),class(copied)) || ...
        ~isequal(size(original),size(copied))
    label = sprintf('%s (class or size differs)',label);
elseif isstruct(original)
    differentFields = setxor(fieldnames(original),fieldnames(copied));
    if ~isempty(differentFields)
        label = sprintf('%s (missing or added fields: %s)', ...
            label,strjoin(differentFields,', '));
        return;
    end
    names = fieldnames(original);
    for element = 1:numel(original)
        for index = 1:numel(names)
            name = names{index};
            if ~isequaln(original(element).(name),copied(element).(name))
                label = first_difference(original(element).(name), ...
                    copied(element).(name),sprintf('%s(%d).%s', ...
                    label,element,name));
                return;
            end
        end
    end
elseif iscell(original)
    for index = 1:numel(original)
        if ~isequaln(original{index},copied{index})
            label = first_difference(original{index},copied{index}, ...
                sprintf('%s{%d}',label,index));
            return;
        end
    end
elseif isnumeric(original) || islogical(original)
    index = find(original ~= copied & ...
        ~(isnan(original) & isnan(copied)),1);
    if ~isempty(index)
        label = sprintf('%s(%d): %.17g versus %.17g', ...
            label,index,original(index),copied(index));
    end
end
end

function restore_environment(savedPath,savedDirectory)
cd(savedDirectory);
path(savedPath);
end
