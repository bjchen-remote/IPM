function report = ipm_accellab_test_exact_inner()
%IPM_ACCELLAB_TEST_EXACT_INNER Known motion, nonzero PPrime, and time FD.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
x = linspace(-1,1.4,161);
x = x+0.003*sin(1.7*x);
y = linspace(0,1.1,97)'.^1.13;
steps = [1e-3,3e-4,1e-4,3e-5,1e-5,3e-6];
cases = cell(1,2);
labels = {'polynomial','gaussian'};
for index = 1:2
    label = labels{index};
    [omega,forcing,known] = synthetic(label,0,x,y);
    [~,peakIndex] = max(omega(1,:));
    exact = ipm_accellab_exact_inner_rates(omega,forcing,x,y,peakIndex);
    predicted = [exact.peakPrime,exact.translationRate, ...
        exact.wallCoreWidthPrime,exact.verticalCoreWidthPrime];
    derivatives = NaN(numel(steps),4);
    errors = NaN(size(derivatives));
    switched = false(numel(steps),1);
    for k = 1:numel(steps)
        dt = steps(k);
        plus = coordinate_sample(label,dt,x,y,exact.signature);
        minus = coordinate_sample(label,-dt,x,y,exact.signature);
        switched(k) = plus.templateChanged || minus.templateChanged;
        if ~switched(k)
            derivatives(k,:) = (coordinates(plus)-coordinates(minus))/(2*dt);
            errors(k,:) = abs(derivatives(k,:)-predicted);
        end
    end
    assert(~any(switched) && max(errors(end,:)) < 5e-7 && min(max(errors,[],2)) < 5e-8, ...
        'ipm:InnerTestFiniteDifference','Instantaneous derivatives disagree with same-template time differences.');
    amplitude = ipm_accellab_exact_inner_rates(omega,0.07*omega,x,y,peakIndex);
    assert(abs(amplitude.peakPrime-0.07*amplitude.peak.value) < 1e-12 && ...
        max(abs([amplitude.translationRate,amplitude.logScaleXRate,amplitude.logScaleYRate])) < 1e-10, ...
        'ipm:InnerTestAmplitude','Nonzero PPrime did not cancel pure amplitude motion in the widths.');
    if strcmp(label,'polynomial')
        assert(abs(exact.peakPrime-known.amplitudeRate*exact.peak.value) < 1e-11 && ...
            abs(exact.translationRate-known.translationRate) < 1e-11 && ...
            abs(exact.logScaleYRate-known.logScaleYRate) < 1e-10, ...
            'ipm:InnerTestKnownMotion','Exact polynomial peak/vertical motion was not recovered.');
        % The tilted vertical trace makes the aPrime*dXQ term essential.
        assert(abs(exact.verticalCrossing.bracketForcing(1)) > 0, ...
            'ipm:InnerTestVerticalForcing','The chain-rule test became trivial.');
    end
    later = coordinate_sample(label,0.2,x,y,exact.signature);
    assert(later.templateChanged && ~isempty(later.changedTemplateFields), ...
        'ipm:InnerTestSwitch','A real later-time stencil/bracket change was not recorded.');
    nodeLevel = omega(1,peakIndex-3)/exact.peak.value;
    assert_throws(@()ipm_accellab_exact_inner_rates(omega,forcing,x,y,peakIndex, ...
        struct('level',nodeLevel)),'ipm:InnerCoordinateSwitchBoundary');
    cases{index} = struct('profile',label,'knownContinuousMotion',known, ...
        'exactDiscreteMotion',exact,'timeDifferenceSteps',steps, ...
        'predictedCoordinateDerivatives',predicted,'finiteDifferenceDerivatives',derivatives, ...
        'finiteDifferenceAbsoluteErrors',errors,'templateChangesInDifference',switched, ...
        'pureAmplitudeTest',amplitude,'laterTemplateChange',{later.changedTemplateFields});
end

gridSizes = [81,161,321,641];
betaError = zeros(size(gridSizes));
for k = 1:numel(gridSizes)
    axis = linspace(-1,1.4,gridSizes(k));
    axis = axis+0.003*sin(1.7*axis);
    [omega,forcing,known] = synthetic('polynomial',0,axis,y);
    [~,peakIndex] = max(omega(1,:));
    exact = ipm_accellab_exact_inner_rates(omega,forcing,axis,y,peakIndex);
    betaError(k) = abs(exact.logScaleXRate-known.logScaleXRate);
end
assert(betaError(end) < 0.015,'ipm:InnerTestSpatialConsistency', ...
    'The refined discrete horizontal-width rate is inconsistent with known shrinkage.');

% Real original RHS and the identical previously registered LS observation masks.
opts = struct('nx',65,'ny',33,'xlim',[-4,4],'ymax',4, ...
    'initialCondition','degenerate_primitive','degeneratePower',8, ...
    'symmetryMode','double_odd_omega','rescalingMode','dynamic', ...
    'lengthGauge','transport_anchor','transportAnchorX',1, ...
    'cOmegaGauge','wall_omega_quadratic_peak', ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'adaptiveRemesh',false,'saveResults',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
state = ipm.evolve.initialize(opts);
tiny = ipm_accellab_modulated_probe(state,struct('fitHalfWidthX',2, ...
    'fitHeightY',3,'holdoutHalfWidthX',3,'holdoutHeightY',5, ...
    'includeExactCoordinates',true));
assert(tiny.exactInnerCoordinates.valid && tiny.pdeAcceptedSteps == 0 && ...
    tiny.originalRhsEvaluations == 1 && ...
    tiny.exactInnerCoordinates.metrics.holdout.nodeCount == tiny.freeMotion.metrics.holdout.nodeCount, ...
    'ipm:InnerTestTiny','Tiny exact-coordinate residual comparison failed.');
tinyOmega = state.rho*state.ops.Dx';
tracking = abs(state.ops.x-state.ops.rescaling.pinX) <= state.ops.rescaling.peakTrackingHalfWidth;
indices = find(tracking);
[~,local] = max(tinyOmega(1,tracking));
tinyIndex = indices(local);
tinyPeak = ipm.evolve.quadraticPeakFunctional(tinyOmega(1,:),state.ops.x,tinyIndex);
nodeLevel = tinyOmega(1,tinyIndex-3)/tinyPeak.value;
rejection = ipm_accellab_modulated_probe(state,struct('fitHalfWidthX',2, ...
    'fitHeightY',3,'holdoutHalfWidthX',3,'holdoutHeightY',5, ...
    'includeExactCoordinates',true,'exactCoordinateOptions',struct('level',nodeLevel)));
assert(~rejection.exactInnerCoordinates.valid && rejection.freeMotion.valid && ...
    any(strcmp(rejection.exactInnerCoordinates.exceptionIdentifier, ...
        {'ipm:InnerCoordinateSwitchBoundary','ipm:InnerRateBracket'})), ...
    'ipm:InnerTestRejectionRetention', ...
    'Rejection/LS audit: exactValid=%d freeValid=%d exact=%s', ...
    rejection.exactInnerCoordinates.valid,rejection.freeMotion.valid, ...
    jsonencode(rejection.exactInnerCoordinates));
report = struct('status','passed','cases',{cases}, ...
    'horizontalSpatialGridSizes',gridSizes,'horizontalContinuousBetaErrors',betaError, ...
    'tinyOriginalRhsProbe',tiny,'tinyRejectionRetention',rejection,'productionPdeSteps',0, ...
    'switchBoundaryRejected',true,'nonzeroPeakPrimeRetained',true);
fprintf('EXACT INNER PASS: FD errors polynomial %.3e Gaussian %.3e; beta refinement error %.3e.\n', ...
    max(cases{1}.finiteDifferenceAbsoluteErrors(end,:)), ...
    max(cases{2}.finiteDifferenceAbsoluteErrors(end,:)),betaError(end));
end

function sample = coordinate_sample(label,time,x,y,signature)
omega = synthetic(label,time,x,y);
[~,index] = max(omega(1,:));
sample = ipm_accellab_exact_inner_rates(omega,zeros(size(omega)),x,y,index, ...
    struct('previousSignature',signature));
end

function values = coordinates(sample)
values = [sample.peak.value,sample.peak.x,sample.wallCoreWidth,sample.verticalCoreWidth];
end

function [omega,forcing,known] = synthetic(label,time,x,y)
known = struct('translationRate',0.031,'logScaleXRate',-0.65, ...
    'logScaleYRate',-0.67,'amplitudeRate',0.07);
a = 0.137+known.translationRate*time;
lx = 0.45*exp(known.logScaleXRate*time);
ly = 0.50*exp(known.logScaleYRate*time);
amplitude = 1.8*exp(known.amplitudeRate*time);
[X,Y] = meshgrid(x,y);
z = (X-a)/lx; eta = Y/ly; tilt = 0.35;
switch label
    case 'polynomial'
        omega = amplitude*(1-z.^2-eta+tilt*z.*eta);
        omegaX = amplitude*(-2*z+tilt*eta)/lx;
        omegaY = amplitude*(-1+tilt*z)/ly;
    case 'gaussian'
        omega = amplitude*exp(-z.^2-eta.^2+tilt*z.*eta);
        omegaX = omega.*(-2*z+tilt*eta)/lx;
        omegaY = omega.*(-2*eta+tilt*z)/ly;
end
forcing = known.amplitudeRate*omega-known.translationRate*omegaX- ...
    known.logScaleXRate*(X-a).*omegaX-known.logScaleYRate*Y.*omegaY;
end

function assert_throws(callback,identifier)
try
    callback();
catch exception
    assert(strcmp(exception.identifier,identifier),'ipm:InnerTestWrongRejection', ...
        'Unexpected exception: %s.',exception.identifier);
    return;
end
error('ipm:InnerTestMissingRejection','Expected rejection: %s.',identifier);
end
