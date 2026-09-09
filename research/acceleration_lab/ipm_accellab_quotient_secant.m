function [rhoCandidate,report] = ipm_accellab_quotient_secant(states,user)
%IPM_ACCELLAB_QUOTIENT_SECANT Full-field trial, freshly tested local pullback.
%   Extrapolates complete same-grid rho states, never a pasted local patch.
%   Every trial gets the unchanged full Poisson/RHS solve. This returns an
%   independent profile candidate, not an accepted physical-time trajectory.
if nargin < 2, user = struct(); end
opts = struct('xi',linspace(-2,2,65),'eta',linspace(0,3,49), ...
    'trainingHalfWidthX',1,'trainingHeightY',1.5, ...
    'tailBoxHalfWidthX',2.5,'tailBoxHeightY',2, ...
    'tailDensityTolerance',1e-3,'tailSourceTolerance',1e-3, ...
    'tailVelocityTolerance',1e-3,'secant',struct());
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:QuotientOptions','Unknown quotient option.');
for k = 1:numel(names), opts.(names{k}) = user.(names{k}); end
for name = {'trainingHalfWidthX','trainingHeightY','tailBoxHalfWidthX', ...
        'tailBoxHeightY','tailDensityTolerance','tailSourceTolerance','tailVelocityTolerance'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','finite','positive'});
end
validateattributes(opts.xi,{'numeric'},{'vector','finite','real','increasing'});
validateattributes(opts.eta,{'numeric'},{'vector','finite','real','increasing','nonnegative'});
assert(iscell(states) && numel(states) >= 2 && numel(states) <= 6, ...
    'ipm:QuotientHistory','Supply two to six chronological full solver states.');
base = states{end}; ops = base.ops;
for k = 1:numel(states)
    state = states{k};
    for group = {'grid','physics','elliptic','transport','scaling'}
        assert(isequaln(state.config.(group{1}),base.config.(group{1})), ...
            'ipm:QuotientContract','The full history must share its frozen numerical configuration.');
    end
    assert(isequal(state.ops.x,ops.x) && isequal(state.ops.y,ops.y) && ...
        isequaln(state.ops.rescaling,ops.rescaling) && ...
        state.ops.remeshCount == ops.remeshCount && isequal(size(state.rho),size(base.rho)), ...
        'ipm:QuotientGrid','A remesh or mismatched full-field grid requires a new experiment.');
    if k > 1
        assert(state.scale.canonicalTime > states{k-1}.scale.canonicalTime, ...
            'ipm:QuotientTime','Source canonical times must increase.');
    end
end
lastwarn('');
[baseProbe,baseData] = ipm_accellab_modulated_probe(base,struct('includeExactCoordinates',true));
assert(baseProbe.exactInnerCoordinates.valid,'ipm:QuotientBaseline', ...
    'The baseline exact inner coordinates are invalid.');
targetPeak = baseProbe.exactInnerCoordinates.peak.value;
[XI,ETA] = meshgrid(opts.xi(:)',opts.eta(:));
training = abs(XI) <= opts.trainingHalfWidthX & ETA <= opts.trainingHeightY;
holdout = ~training;
wall = ETA == 0;
assert(nnz(training) >= 12 && nnz(holdout) >= 12 && any(wall(:)), ...
    'ipm:QuotientMasks','Observer training/holdout/wall masks are insufficient.');
observerWeights = trap_weights(opts.eta(:))*trap_weights(opts.xi(:))';
trainingWeights = observerWeights(training)/sum(observerWeights(training));
tail = abs(ops.X) >= opts.tailBoxHalfWidthX | ops.Y >= opts.tailBoxHeightY;
assert(any(tail(:)),'ipm:QuotientTail','The fixed tail matching region is empty.');
tailWeights = ops.integrationWeights(tail);
tailWeights = tailWeights/sum(tailWeights);
tailDensityScale = max(sqrt(sum(tailWeights.*base.rho(tail).^2)),realmin);
tailVelocityScale = max(sqrt(sum(tailWeights.* ...
    (baseData.flow.u1(tail).^2+baseData.flow.u2(tail).^2))),realmin);
rhsCount = 1;
baselineVector = []; baselineAudit = struct();
[baselineVector,baselineAudit] = observed(base.rho,baseProbe,baseData);
history = zeros(numel(base.rho),numel(states));
amplitudeFactors = zeros(1,numel(states));
for k = 1:numel(states)
    amplitudeFactors(k) = targetPeak/peak_value(states{k}.rho);
    history(:,k) = amplitudeFactors(k)*states{k}.rho(:);
end
secant = opts.secant;
assert(~isfield(secant,'project') && ~isfield(secant,'accept'), ...
    'ipm:QuotientOptions','Projection and all independent guards are fixed by the protocol.');
secant.project = @project;
secant.accept = @guards;
if ~isfield(secant,'maxRelativeJump'), secant.maxRelativeJump = 0.05; end
if ~isfield(secant,'svdRelativeFloor'), secant.svdRelativeFloor = 1e-8; end
[vector,report] = ipm_accellab_secant(history,@evaluate,secant);
rhoCandidate = reshape(vector,size(base.rho));
report.kind = 'guarded_full_state_secant_for_local_shape_quotient';
report.resolvedSecantOptions = report.options;
report.options = opts;
for k = 1:numel(report.attempts)
    attempt = report.attempts(k);
    failures = {};
    if isfield(attempt.audit,'valid') && attempt.audit.valid
        [~,failures] = guards(attempt.audit,report.baselineAudit);
    end
    required = (1-report.resolvedSecantOptions.minimumDecrease*attempt.damping)*report.baselineNorm;
    if ~isfinite(attempt.residualNorm) || attempt.residualNorm > required
        failures{end+1} = 'trainingL2_sufficient_decrease'; %#ok<AGROW>
    end
    report.attempts(k).guardFailures = failures;
    report.attempts(k).maximumTrainingResidual = required;
end
report.historyAmplitudeFactors = amplitudeFactors;
report.sourceCanonicalTimes = cellfun(@(s)s.scale.canonicalTime,states);
report.originalRhsEvaluations = rhsCount;
report.pdeAcceptedSteps = 0;
report.isTimeTrajectory = false;
report.artifactDisposition = 'independent_full_field_profile_candidate_never_checkpoint';
report.matchingProtocol = 'Complete-field linear combination; fixed tail rho/source/velocity matching; fresh original Poisson/RHS; common-observer and native holdout gates.';

    function output = project(input)
        rho = reshape(input,size(base.rho));
        output = (targetPeak/peak_value(rho))*input;
    end

    function value = peak_value(rho)
        omega = rho*ops.Dx';
        tracking = abs(ops.x-ops.rescaling.pinX) <= ops.rescaling.peakTrackingHalfWidth;
        indices = find(tracking);
        [~,index] = max(omega(1,tracking));
        peak = ipm.evolve.quadraticPeakFunctional(omega(1,:),ops.x,indices(index));
        value = peak.value;
    end

    function [residual,audit] = evaluate(input)
        rho = reshape(input,size(base.rho));
        if isequal(rho,base.rho)
            residual = baselineVector; audit = baselineAudit; return;
        end
        candidate = base; candidate.rho = rho;
        rhsCount = rhsCount+1;
        lastwarn('');
        [probe,data] = ipm_accellab_modulated_probe(candidate,struct('includeExactCoordinates',true));
        assert(probe.exactInnerCoordinates.valid,'ipm:QuotientCoordinates', ...
            'A trial left the exact inner-coordinate admissibility conditions.');
        [residual,audit] = observed(rho,probe,data);
    end

    function [residual,audit] = observed(rho,probe,data)
        [warningMessage,warningIdentifier] = lastwarn;
        exact = probe.exactInnerCoordinates;
        chart = ipm_accellab_pullback(data.omega,data.forcing,ops.x,ops.y,exact,opts.xi,opts.eta);
        residual = sqrt(trainingWeights).*chart.G(training);
        holdoutWeights = observerWeights(holdout)/sum(observerWeights(holdout));
        native = exact.metrics;
        audit = struct('valid',all(isfinite(chart.G(:))) && isfinite(probe.poissonResidual), ...
            'trainingL2',norm(residual),'trainingInf',max(abs(chart.G(training))), ...
            'holdoutL2',sqrt(sum(holdoutWeights.*chart.G(holdout).^2)), ...
            'holdoutInf',max(abs(chart.G(holdout))), ...
            'wallInf',max(abs(chart.G(wall))), ...
            'nativeCoreL2',native.core.modulatedL2,'nativeCoreInf',native.core.modulatedInf, ...
            'nativeHoldoutL2',native.holdout.modulatedL2, ...
            'nativeHoldoutInf',native.holdout.modulatedInf, ...
            'tailDensityRelativeL2',sqrt(sum(tailWeights.*(rho(tail)-base.rho(tail)).^2))/tailDensityScale, ...
            'tailSourceAbsoluteInfOverPeak',max(abs(data.omega(tail)-baseData.omega(tail)))/targetPeak, ...
            'tailVelocityRelativeL2',sqrt(sum(tailWeights.* ...
                ((data.flow.u1(tail)-baseData.flow.u1(tail)).^2+ ...
                (data.flow.u2(tail)-baseData.flow.u2(tail)).^2)))/tailVelocityScale, ...
            'peak',exact.peak.value,'peakPrime',exact.peakPrime, ...
            'aPrime',exact.translationRate,'beta',exact.logScaleXRate,'gamma',exact.logScaleYRate, ...
            'wallCoreWidth',exact.wallCoreWidth,'verticalCoreWidth',exact.verticalCoreWidth, ...
            'coordinateSignature',exact.signature, ...
            'activePeakMargin',exact.peak.activeNodeMargin, ...
            'crossingMargins',[exact.horizontalLeft.templateInteriorMargin, ...
                exact.horizontalRight.templateInteriorMargin,exact.verticalCrossing.templateInteriorMargin], ...
            'poissonResidual',probe.poissonResidual, ...
            'lastRuntimeWarningIdentifier',warningIdentifier, ...
            'lastRuntimeWarningMessage',warningMessage, ...
            'observerMinimumMovingCellMargin',chart.minimumMovingCellMargin);
    end

    function [accepted,failures] = guards(trial,baseline)
        failures = {};
        tests = [abs(trial.peak/targetPeak-1) <= 1e-8, ...
            trial.tailDensityRelativeL2 <= opts.tailDensityTolerance, ...
            trial.tailSourceAbsoluteInfOverPeak <= opts.tailSourceTolerance, ...
            trial.tailVelocityRelativeL2 <= opts.tailVelocityTolerance];
        labels = {'quadratic_peak_normalization','tailDensityRelativeL2', ...
            'tailSourceAbsoluteInfOverPeak','tailVelocityRelativeL2'};
        failures = [failures,labels(~tests)];
        for field = {'trainingInf','holdoutL2','holdoutInf','wallInf', ...
                'nativeCoreL2','nativeCoreInf','nativeHoldoutL2','nativeHoldoutInf'}
            if trial.(field{1}) > baseline.(field{1})*(1+1e-12)+1e-14
                failures{end+1} = field{1}; %#ok<AGROW>
            end
        end
        accepted = isempty(failures);
    end
end

function weights = trap_weights(axis)
assert(numel(axis) >= 2,'ipm:QuotientObserver','Each observer axis needs at least two nodes.');
weights = [diff(axis(1:2));axis(3:end)-axis(1:end-2);diff(axis(end-1:end))]/2;
end
