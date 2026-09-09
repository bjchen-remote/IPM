function report = ipm_perflab_benchmark(input,options)
%IPM_PERFLAB_BENCHMARK Bounded timings on one explicit current-contract state.
%   INPUT is a flat solver-options structure, checkpoint structure, or native
%   checkpoint filename. No-input calls are deliberately rejected. The
%   benchmark reuses maintained kernels; it implements no time integrator.
%   OPTIONS: repetitions=3, solveSteps=0, profileSolve=false. solveSteps>0
%   runs the sole solver entry ipm.solve with an explicit short step limit.
%   File output, plots and video are disabled. Fresh inputs also disable
%   snapshots; checkpoint inputs retain their frozen snapshot policy.
%   Return REPORT to
%   the caller; persistence is the caller's responsibility under result/.

if nargin < 1 || isempty(input)
    error('ipm:BenchmarkExplicitInput','Pass explicit options or checkpoint.');
end
if nargin < 2
    options = struct();
end
options = settings(options);
controls = struct('saveResults',false,'makePlots',false,'livePlot',false, ...
    'writeVideo',false,'verbose',false, ...
    'checkpoint',struct('enabled',false));
clock = tic;
if ischar(input) || isstring(input) || ...
        (isstruct(input) && isfield(input,'payload'))
    [state,log,cursor] = ipm.output.restoreCheckpoint(input,controls);
    inputKind = 'native_checkpoint';
else
    assert(isstruct(input) && isscalar(input),'Pass scalar solver options.');
    names = fieldnames(controls);
    for index = 1:numel(names)
        input.(names{index}) = controls.(names{index});
    end
    input.storeSnapshots = false;
    state = ipm.evolve.initialize(input);
    log = ipm.output.initializeLog();
    [log,~] = ipm.output.record(log,state);
    cursor = struct('nextOutput',state.config.time.outputEvery, ...
        'nextCheckpoint',Inf,'lastCheckpointStep',-1, ...
        'lastCheckpointCanonicalTime',NaN);
    inputKind = 'explicit_fresh_options';
end
initializationSeconds = toc(clock);
assert(state.config.schemaVersion == 4 && strcmp( ...
    state.config.scaling.scalingContract,'exact_gauge_no_feedback_v1'), ...
    'Benchmark requires the schema-4 exact-gauge contract.');
rho = state.rho;
ops = state.ops;
scale = state.scale;
flow = state.flow;
kappa = 1;
if strcmp(ops.dynamicScaleGeometry,'anisotropic')
    kappa = exp(scale.logC_y-scale.logC_l);
end
source = flow.source;
assemblyBoundary = ops.transportBoundaryMode;
if ops.rescaling.enabled
    assemblyBoundary = 'open';
end
time = state.config.time;
time.finalTime = Inf;
time.physicalFinalTime = Inf;
[dt,stopReason] = ipm.evolve.selectTimestep(flow,scale,ops,time);
assert(isempty(stopReason) && isfinite(dt) && dt > 0, ...
    'The frozen state has no valid positive benchmark timestep.');
method = state.config.time.timeIntegrator;
cachedStep = @()ipm.evolve.stepRk(rho,dt,ops,method,scale,state.rhsCache);
uncachedStep = @()ipm.evolve.stepRk(rho,dt,ops,method,scale);
[rhoCached,flowCached,scaleCached,cacheCached,infoCached] = ...
    ipm.evolve.stepRk(rho,dt,ops,method,scale,state.rhsCache);
[rhoUncached,flowUncached,scaleUncached,cacheUncached,infoUncached] = ...
    ipm.evolve.stepRk(rho,dt,ops,method,scale);
cacheEquivalent = isequaln(rhoCached,rhoUncached) && ...
    isequaln(flowCached,flowUncached) && isequaln(scaleCached,scaleUncached) && ...
    isequaln(cacheCached,cacheUncached);
assert(cacheEquivalent && infoCached.initialStageCacheHit, ...
    'The cached/uncached maintained RK paths did not agree exactly.');
% This in-memory benchmark transaction is installed at the frozen step.
cursor.lastCheckpointStep = state.step;
cursor.lastCheckpointCanonicalTime = scale.canonicalTime;
checkpoint = ipm.output.makeCheckpoint(state,log,cursor);
operations = { ...
    'greenBoundary',@()ipm.field.greenBoundary(source,ops,kappa); ...
    'poisson',@()ipm.field.poisson(source,ops,kappa); ...
    'velocity',@()ipm.field.velocity(rho,ops,kappa); ...
    'rescaledAssembly',@()ipm.evolve.assembleRhs(rho,flow.u1,flow.u2, ...
        flow.c_x,flow.c_y,flow.c_omega,flow.c_r,ops,assemblyBoundary); ...
    'flow',@()ipm.evolve.flow(rho,ops,scale); ...
    'measure',@()ipm.diagnostics.measure(rho,flow,ops,state.mass0); ...
    'record',@()ipm.output.record(log,state); ...
    'makeCheckpoint',@()ipm.output.makeCheckpoint(state,log,cursor); ...
    'stepCached',cachedStep; ...
    'stepUncached',uncachedStep};
timings = struct();
for index = 1:size(operations,1)
    operation = operations{index,2};
    operation();
    samples = zeros(1,options.repetitions);
    for repetition = 1:options.repetitions
        clock = tic;
        operation();
        samples(repetition) = toc(clock);
    end
    timings.(operations{index,1}) = struct( ...
        'medianSeconds',median(samples),'samplesSeconds',samples);
end
report = struct('schemaVersion',1,'kind','ipm_frozen_state_benchmark', ...
    'matlabVersion',version,'computer',computer,'inputKind',inputKind, ...
    'computationalThreads',maxNumCompThreads, ...
    'options',options,'grid',[ops.ny,ops.nx],'timeIntegrator',method, ...
    'spatialDiscretization',ops.spatialDiscretization, ...
    'transportScheme',ops.transportScheme, ...
    'canonicalTime',scale.canonicalTime,'acceptedStep',state.step, ...
    'initializationSeconds',initializationSeconds,'dt',dt, ...
    'cachedRhsEvaluations',infoCached.rhsEvaluations, ...
    'uncachedRhsEvaluations',infoUncached.rhsEvaluations, ...
    'cacheEquivalent',cacheEquivalent,'timings',timings, ...
    'cacheStepSpeedup',timings.stepUncached.medianSeconds / ...
        timings.stepCached.medianSeconds, ...
    'rhoXReuseExact',isequaln(rho*ops.Dx',flow.source));
if options.solveSteps > 0
    controls.maxSteps = state.step+options.solveSteps;
    controls.finalTime = scale.canonicalTime+2*options.solveSteps*dt;
    controls.physicalFinalTime = Inf;
    if options.profileSolve
        profile clear;
        profile on;
        profileCleanup = onCleanup(@()profile('off'));
    end
    clock = tic;
    capturedOutput = evalc('result = ipm.solve(controls,checkpoint);');
    report.solveSeconds = toc(clock);
    if options.profileSolve
        profile off;
        report.profile = profile('info');
        clear profileCleanup;
    end
    report.solveStopReason = result.state.stopReason;
    report.actualSolveSteps = result.history.common.acceptedStep(end)-state.step;
    report.solveCanonicalTime = result.state.canonicalTime;
    report.solveOutput = capturedOutput;
end
fprintf('Frozen benchmark %dx%d, %s, tau=%.8g; cache exact=%d\n', ...
    report.grid,method,scale.canonicalTime,cacheEquivalent);
names = fieldnames(timings);
for index = 1:numel(names)
    fprintf('  %-18s %.6g s\n',names{index},timings.(names{index}).medianSeconds);
end
end

function options = settings(options)
defaults = struct('repetitions',3,'solveSteps',0,'profileSolve',false);
assert(isstruct(options) && isscalar(options),'Pass scalar options.');
names = fieldnames(options);
assert(all(ismember(names,fieldnames(defaults))),'Unknown benchmark option.');
for index = 1:numel(names)
    defaults.(names{index}) = options.(names{index});
end
options = defaults;
validateattributes(options.repetitions,{'numeric'}, ...
    {'scalar','integer','>=',1,'<=',20});
validateattributes(options.solveSteps,{'numeric'}, ...
    {'scalar','integer','>=',0,'<=',20});
validateattributes(options.profileSolve,{'logical'},{'scalar'});
assert(~options.profileSolve || options.solveSteps > 0, ...
    'profileSolve requires solveSteps > 0. It replaces the MATLAB profiler buffer.');
end
