function [report,artifacts] = ipm_accellab_test_quotient()
%IPM_ACCELLAB_TEST_QUOTIENT Pullback sensitivity, original flow map, and guards.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
x = linspace(-1,1.4,161); x = x+0.003*sin(1.7*x);
y = linspace(0,1.1,97)'.^1.13;
xi = linspace(-2,2,65); eta = linspace(0,3,49);
[omega,forcing] = synthetic(0,x,y);
[~,index] = max(omega(1,:));
exact = ipm_accellab_exact_inner_rates(omega,forcing,x,y,index);
chart = ipm_accellab_pullback(omega,forcing,x,y,exact,xi,eta);
dt = 1e-5;
plus = synthetic_chart(dt,x,y,xi,eta);
minus = synthetic_chart(-dt,x,y,xi,eta);
stable = stable_cells(chart,plus) & stable_cells(chart,minus);
difference = (plus.U-minus.U)/(2*dt);
syntheticError = max(abs(difference(stable)-chart.G(stable)));
assert(nnz(stable) > 0.95*numel(stable) && syntheticError < 1e-7, ...
    'ipm:QuotientTestSensitivity','The specified pullback derivative failed a same-cell time difference.');
amplitude = ipm_accellab_exact_inner_rates(omega,0.07*omega,x,y,index);
amplitudeChart = ipm_accellab_pullback(omega,0.07*omega,x,y,amplitude,xi,eta);
assert(max(abs(amplitudeChart.G),[],'all') < 1e-11, ...
    'ipm:QuotientTestAmplitude','The normalized shape observer retained a spurious amplitude mode.');

opts = struct('nx',129,'ny',65,'xlim',[-4,4],'ymax',4, ...
    'initialCondition','degenerate_primitive','degeneratePower',8, ...
    'symmetryMode','double_odd_omega','rescalingMode','dynamic', ...
    'lengthGauge','transport_anchor','transportAnchorX',1, ...
    'cOmegaGauge','wall_omega_quadratic_peak', ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'adaptiveRemesh',false,'maxDt',0.002,'finalTime',1,'physicalFinalTime',1, ...
    'saveResults',false,'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false);
state = ipm.evolve.initialize(opts);
states = {state};
for step = 1:32
    [state,stop] = ipm.evolve.advance(state);
    assert(isempty(stop),'ipm:QuotientTestHistory','The small original-PDE history stopped early.');
    if mod(step,8) == 0, states{end+1} = state; end %#ok<AGROW>
end
baseRho = state.rho;
[baseProbe,data] = ipm_accellab_modulated_probe(state,struct('includeExactCoordinates',true));
assert(baseProbe.exactInnerCoordinates.valid,'ipm:QuotientTestBaseline','The exact baseline chart is invalid.');
baseChart = ipm_accellab_pullback(data.omega,data.forcing,state.ops.x,state.ops.y, ...
    baseProbe.exactInnerCoordinates,xi,eta);
flowSteps = [1e-3,5e-4,2.5e-4];
flowErrors = zeros(size(flowSteps));
stableFractions = zeros(size(flowSteps));
for k = 1:numel(flowSteps)
    [rho,~,~] = ipm.evolve.stepSsprk54(state.rho,flowSteps(k),state.ops,state.scale);
    omega = rho*state.ops.Dx';
    tracking = abs(state.ops.x-state.ops.rescaling.pinX) <= state.ops.rescaling.peakTrackingHalfWidth;
    indices = find(tracking); [~,index] = max(omega(1,tracking));
    coordinate = ipm_accellab_exact_inner_rates(omega,zeros(size(omega)), ...
        state.ops.x,state.ops.y,indices(index));
    next = ipm_accellab_pullback(omega,zeros(size(omega)),state.ops.x,state.ops.y,coordinate,xi,eta);
    stable = stable_cells(baseChart,next);
    difference = (next.U-baseChart.U)/flowSteps(k);
    flowErrors(k) = max(abs(difference(stable)-baseChart.G(stable)));
    stableFractions(k) = nnz(stable)/numel(stable);
end
assert(all(stableFractions > 0.95) && flowErrors(end) < 0.4*flowErrors(1), ...
    'ipm:QuotientTestFlowMap','The original-PDE pullback difference did not converge to the fresh observer derivative.');
trials = cell(1,3);
candidateFields = cell(1,3);
relaxation = cell(1,3);
for k = 1:3
    memory = [1,2,4];
    [candidate,trials{k}] = ipm_accellab_quotient_secant(states,struct('secant', ...
        struct('memory',memory(k),'damping',[0.25,0.0625,0.015625,0.00390625])));
    candidateFields{k} = candidate;
    assert(isequal(state.rho,baseRho) && all(isfinite(candidate(:))) && ...
        ~trials{k}.isTimeTrajectory && trials{k}.pdeAcceptedSteps == 0, ...
        'ipm:QuotientTestTransaction','An independent quotient trial mutated the source trajectory.');
    fprintf('QUOTIENT trial memory=%d status=%s ratio=%.8g RHS=%d\n', ...
        memory(k),trials{k}.status,trials{k}.residualRatio,trials{k}.originalRhsEvaluations);
    if trials{k}.accepted
        source = state; trial = state; trial.rho = candidate;
        for relaxationStep = 1:4
            [source.rho,source.flow,source.scale] = ipm.evolve.stepSsprk54( ...
                source.rho,0.002,source.ops,source.scale);
            [trial.rho,trial.flow,trial.scale] = ipm.evolve.stepSsprk54( ...
                trial.rho,0.002,trial.ops,trial.scale);
        end
        baselineObservation = observe(source,xi,eta);
        candidateObservation = observe(trial,xi,eta);
        ratios = struct();
        for name = {'trainingL2','trainingInf','holdoutL2','holdoutInf', ...
                'wallInf','nativeCoreL2','nativeCoreInf','nativeHoldoutL2','nativeHoldoutInf'}
            ratios.(name{1}) = candidateObservation.(name{1})/baselineObservation.(name{1});
        end
        relaxation{k} = struct('kind','independent_original_pde_short_relaxation_comparison', ...
            'stepsPerBranch',4,'canonicalIncrement',0.008, ...
            'baseline',baselineObservation,'candidate',candidateObservation, ...
            'candidateToBaselineRatios',ratios,'productionTrajectoryChanged',false);
    end
end
report = struct('status','passed','syntheticPullbackTimeDifferenceError',syntheticError, ...
    'pureAmplitudeRemoved',true,'originalHistorySteps',32,'gridSize',[129,65], ...
    'shortOriginalFlowMapSteps',flowSteps,'originalFlowMapErrors',flowErrors, ...
    'originalFlowMapStableCellFractions',stableFractions,'quotientTrials',{trials}, ...
    'shortRelaxation',{relaxation}, ...
    'productionPdeSteps',0,'productionAccelerationEstablished',false);
artifacts = struct('kind','independent_research_fields_not_native_checkpoint', ...
    'config',state.config,'runtimeRescaling',state.ops.rescaling,'sourceScale',state.scale, ...
    'x',state.ops.x,'y',state.ops.y,'sourceRho',baseRho, ...
    'historyRho',{cellfun(@(s)s.rho,states,'UniformOutput',false)}, ...
    'candidateRho',{candidateFields});
fprintf('QUOTIENT PASS: synthetic derivative error %.3e; original flow-map errors %.3e -> %.3e.\n', ...
    syntheticError,flowErrors(1),flowErrors(end));
end

function result = observe(state,xi,eta)
[probe,data] = ipm_accellab_modulated_probe(state,struct('includeExactCoordinates',true));
assert(probe.exactInnerCoordinates.valid,'ipm:QuotientRelaxationChart','Relaxation left the exact chart.');
chart = ipm_accellab_pullback(data.omega,data.forcing,state.ops.x,state.ops.y, ...
    probe.exactInnerCoordinates,xi,eta);
[XI,ETA] = meshgrid(xi,eta); core = abs(XI) <= 1 & ETA <= 1.5;
weights = ones(size(XI)); weights([1,end],:) = weights([1,end],:)/2;
weights(:,[1,end]) = weights(:,[1,end])/2;
native = probe.exactInnerCoordinates.metrics;
result = struct('trainingL2',sqrt(sum(weights(core).*chart.G(core).^2)/sum(weights(core))), ...
    'trainingInf',max(abs(chart.G(core))), ...
    'holdoutL2',sqrt(sum(weights(~core).*chart.G(~core).^2)/sum(weights(~core))), ...
    'holdoutInf',max(abs(chart.G(~core))),'wallInf',max(abs(chart.G(ETA == 0))), ...
    'nativeCoreL2',native.core.modulatedL2,'nativeCoreInf',native.core.modulatedInf, ...
    'nativeHoldoutL2',native.holdout.modulatedL2,'nativeHoldoutInf',native.holdout.modulatedInf, ...
    'peakPrime',probe.exactInnerCoordinates.peakPrime,'poissonResidual',probe.poissonResidual);
end

function stable = stable_cells(a,b)
stable = a.cellRows == b.cellRows & a.cellColumns == b.cellColumns;
end

function chart = synthetic_chart(time,x,y,xi,eta)
omega = synthetic(time,x,y);
[~,index] = max(omega(1,:));
coordinate = ipm_accellab_exact_inner_rates(omega,zeros(size(omega)),x,y,index);
chart = ipm_accellab_pullback(omega,zeros(size(omega)),x,y,coordinate,xi,eta);
end

function [omega,forcing] = synthetic(time,x,y)
aPrime = 0.031; beta = -0.65; gamma = -0.67; sigma = 0.07;
a = 0.137+aPrime*time; lx = 0.45*exp(beta*time); ly = 0.50*exp(gamma*time);
amplitude = 1.8*exp(sigma*time); [X,Y] = meshgrid(x,y);
z = (X-a)/lx; eta = Y/ly; tilt = 0.35;
omega = amplitude*(1-z.^2-eta+tilt*z.*eta);
omegaX = amplitude*(-2*z+tilt*eta)/lx;
omegaY = amplitude*(-1+tilt*z)/ly;
forcing = sigma*omega-aPrime*omegaX-beta*(X-a).*omegaX-gamma*Y.*omegaY;
end
