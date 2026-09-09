function report = ipm_accellab_integrate_hermite_peak(frameFile,outputRoot,protocol)
%IPM_ACCELLAB_INTEGRATE_HERMITE_PEAK Registered tiny independent gauge trial.
% Native baselines use unchanged stepSsprk54. Research states never become
% native checkpoints or acquire a production quadratic-gauge identity.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
if nargin < 3, protocol = 'fine'; end
assert(any(strcmp(protocol,{'fine','coarse_order_calibration'})));
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['hermite_peak_integration_',token]);
mkdir(destination);
registration = struct('sourceFrameFile',frameFile,'outputDirectory',destination, ...
    'nx',[49,65],'ny',[25,33],'dt',[0.0005,0.00025,0.000125],'tauEnd',0.08, ...
    'researchEquationId','independent_hermite_C1_peak_amplitude_fixed_window_v1', ...
    'nativeBaselineMethod','maintained ipm.evolve.stepSsprk54', ...
    'researchMethod','same SSPRK54 tableau and arithmetic, with simultaneous five-component scale evolution', ...
    'peakWindowRule','[0,x(end)] on a fixed native grid', ...
    'usesPeakProjection',false,'usesFeedback',false,'nativeCheckpointsWritten',0, ...
    'referencePolicy','Three registered levels provide self differences; finest is a numerical comparison, not exact truth.', ...
    'nativeStageTrace','Maintained step internals are unchanged; native observations are accepted-state only.', ...
    'researchStageTrace','All five RK stage inputs and each terminal accepted state are retained.');
registration.protocol = protocol;
if strcmp(protocol,'coarse_order_calibration')
    registration.dt = [0.004,0.002,0.001];
    registration.referencePolicy = 'Additional coarse self-convergence calibration; the separate original fine levels remain retained at their roundoff floor.';
end
save(fullfile(destination,'registration.mat'),'registration'); write_json(fullfile(destination,'registration.json'),registration);
fprintf('HERMITE_PEAK_INTEGRATION_REGISTERED %s\n',destination);
source = load(frameFile,'config');
schema = ipm.config.schema(); opts = struct();
for k = 1:numel(schema.domainNames)
    domain = source.config.(schema.domainNames{k}); names = fieldnames(domain);
    for j = 1:numel(names), opts.(names{j}) = domain.(names{j}); end
end
opts.finalTime = registration.tauEnd; opts.saveResults = false; opts.storeSnapshots = false;
opts.verbose = false; summaries = cell(2,2,3); runs = cell(2,2,3); parities = cell(2,1);
try
    for g = 1:2
        opts.nx = registration.nx(g); opts.ny = registration.ny(g);
        state = ipm.evolve.initialize(opts);
        assert(~state.config.remesh.adaptiveRemesh && ~state.config.remesh.initialAnalyticRemesh && ...
            strcmp(state.ops.spatialDiscretization,'high_order') && strcmp(state.ops.transportScheme,'weno5_fd'), ...
            'ipm:PeakResearchOperatorScope','The registered spatial scheme or fixed-grid condition changed.');
        parities{g} = parity_check(state,registration.dt(1));
        save(fullfile(destination,sprintf('native_parity_n%d.mat',opts.nx)),'parities');
        fprintf('PEAK_NATIVE_PARITY n%d passed\n',opts.nx);
        for m = 1:2
            modes = {'native_quadratic','research_hermite_C1'}; mode = modes{m};
            for level = 1:3
                dt = registration.dt(level);
                label = sprintf('%s_n%d_level%d',mode,opts.nx,level);
                started = tic;
                run = evolve_case(state,dt,registration.tauEnd,mode,destination,label);
                run.elapsedSeconds = toc(started);
                run.researchEquationId = registration.researchEquationId;
                if m == 1, run.researchEquationId = 'unmodified_native_quadratic_equation'; end
                run.operatorSourceConfig = state.config;
                run.configurationInterpretation = 'operatorSourceConfig names the native operator source; mode and researchEquationId identify the actually integrated RHS.';
                summary = summarize(run);
                save(fullfile(destination,[label,'.mat']),'run','summary','-v7.3');
                write_json(fullfile(destination,[label,'.json']),summary);
                runs{g,m,level} = run; summaries{g,m,level} = summary;
                fprintf('PEAK_RUN n%d %s dt %.8g Hdrift %.4e Qdrift %.4e maxHPrime %.4e seconds %.2f\n', ...
                    opts.nx,mode,dt,summary.maxRelativeHDrift,summary.maxRelativeQuadraticDrift, ...
                    summary.maximumAbsoluteHPrime,run.elapsedSeconds);
            end
        end
        clear state;
    end
    comparisons = cell(2,2);
    for g = 1:2
        for m = 1:2
            comparisons{g,m} = time_comparison(squeeze(runs(g,m,:)));
        end
    end
    if strcmp(protocol,'fine')
        nativeReplay = validate_saved_native(runs{2,1,2},frameFile);
    else
        nativeReplay = struct('performed',false,'reason','Archived replay dt differs from this separately registered coarse protocol; four-step native parity was checked on both grids.');
    end
    report = struct('status','completed_independent_tiny_research','registration',registration, ...
        'nativeParity',{parities},'runSummaries',{summaries},'timeComparisons',{comparisons}, ...
        'nativeSavedReplay',nativeReplay,'allRunsFinite',true,'nativeCheckpointsWritten',0, ...
        'interpretation','Actual finite-time peak drift and field self-convergence are reported separately from instantaneous gauge cancellation. No production or singularity conclusion.');
    save(fullfile(destination,'report.mat'),'report','-v7.3'); write_json(fullfile(destination,'report.json'),report);
catch exception
    failure = struct('status','failed_preserved','registration',registration,'completedSummaries',{summaries}, ...
        'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure'); write_json(fullfile(destination,'failure.json'),failure);
    fprintf('HERMITE_PEAK_INTEGRATION_FAILED %s\n',destination); rethrow(exception);
end
fprintf('HERMITE_PEAK_INTEGRATION_COMPLETE %s\n',destination);
end

function result = parity_check(state,dt)
rho = state.rho; scale = state.scale; z = pack(scale); nativeCache = state.rhsCache;
cache = []; maximumFieldDifference = 0; maximumScaleDifference = 0;
for k = 1:4
    [a,~,s,nativeCache] = ipm.evolve.stepSsprk54(rho,dt,state.ops,scale,nativeCache);
    [b,t,cache] = ipm_accellab_peak_step(rho,z,dt,state.ops,'native_quadratic',cache);
    maximumFieldDifference = max(maximumFieldDifference,max(abs(a-b),[],'all'));
    maximumScaleDifference = max(maximumScaleDifference,max(abs(pack(s)-t)));
    assert(isequaln(a,b) && isequaln(pack(s),t),'ipm:PeakStageParity','Independent stage arithmetic differs from native SSPRK54.');
    rho = a; scale = s; z = t;
end
% The new direct operator path must share the original spatial base RHS.
[original,flow] = ipm.evolve.flow(rho,state.ops);
e = ipm_accellab_peak_rhs(rho,z,state.ops,'research_hermite_C1');
base = ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,flow.c_l,flow.c_l,0,flow.c_r,state.ops,'open');
assert(isequaln(original,base+flow.c_omega*rho) && e.diagnostic.c_l == flow.c_l && flow.c_r == 0);
assert(isequaln(e.rhoRate,base+e.diagnostic.c_omega*rho), ...
    'ipm:PeakBaseOperatorParity','The direct C1 path changed the maintained base spatial operator.');
result = struct('steps',4,'dt',dt,'fieldBitwise',true,'allFiveScaleComponentsBitwise',true, ...
    'maximumFieldDifference',maximumFieldDifference,'maximumScaleDifference',maximumScaleDifference, ...
    'independentBaseRhsBitwise',true);
end

function run = evolve_case(state,dt,tauEnd,mode,destination,label)
count = round(tauEnd/dt); assert(abs(count*dt-tauEnd) < 1e-14);
rho = state.rho; scale = state.scale; z = pack(scale); ops = state.ops;
cache = []; nativeCache = state.rhsCache;
accepted = cell(count+1,1); snapshots = cell(count+1,1); stages = cell(count,1);
snapshots{1} = rho; rhsCount = 0;
if strcmp(mode,'native_quadratic')
    native = struct('rho',rho,'z',z,'rhs',nativeCache.rhoRate,'flow',state.flow);
    e = ipm_accellab_peak_rhs(rho,z,ops,mode,native);
else
    e = ipm_accellab_peak_rhs(rho,z,ops,mode); cache = e; rhsCount = 1;
end
accepted{1} = e.diagnostic;
try
    for step = 1:count
        if strcmp(mode,'native_quadratic')
            [rho,flow,scale,nativeCache,info] = ipm.evolve.stepSsprk54(rho,dt,ops,scale,nativeCache);
            rhsCount = rhsCount+info.rhsEvaluations; z = pack(scale);
            native = struct('rho',rho,'z',z,'rhs',nativeCache.rhoRate,'flow',flow);
            e = ipm_accellab_peak_rhs(rho,z,ops,mode,native);
        else
            [rho,z,cache,trace] = ipm_accellab_peak_step(rho,z,dt,ops,mode,cache);
            stages{step} = trace; e = cache; rhsCount = rhsCount+5;
            assert(all(cellfun(@(s)s.HValid && s.HActiveCount == 1 && ...
                dt*s.transportRate < state.config.time.cfl,trace)), ...
                'ipm:PeakStageGuard','A research stage has a nonunique/invalid maximum or exceeds native CFL.');
        end
        assert(e.diagnostic.HValid && dt*e.diagnostic.transportRate < state.config.time.cfl, ...
            'ipm:PeakAcceptedGuard','An accepted state has an invalid Hermite observer or exceeds native CFL.');
        accepted{step+1} = e.diagnostic; snapshots{step+1} = rho;
        if mod(step,80) == 0
            fprintf('PEAK_PROGRESS %s step %d/%d tau %.6f Hrel %.3e\n', ...
                label,step,count,z(5),(e.diagnostic.H-accepted{1}.H)/accepted{1}.H);
        end
    end
catch exception
    failedState = struct('mode',mode,'dt',dt,'step',step,'rho',rho,'z',z, ...
        'accepted',{accepted},'stageTrace',{stages},'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,[label,'_failed_state.mat']),'failedState','-v7.3'); rethrow(exception);
end
run = struct('mode',mode,'dt',dt,'steps',count,'rho',rho,'z',z, ...
    'scaleComponentNames',{{'logC_l','logC_omega','physicalTime','X_shift','canonicalTime'}}, ...
    'x',ops.x,'y',ops.y,'Dx',ops.Dx,'Dy',ops.Dy,'integrationWeights',ops.integrationWeights, ...
    'accepted',{accepted},'rhoSnapshots',{snapshots},'stageTrace',{stages}, ...
    'newRhsEvaluations',rhsCount,'isNativeCheckpoint',false,'peakProjectionCount',0);
end

function summary = summarize(run)
a = [run.accepted{:}]; H = [a.H]; Q = [a.quadraticP];
allDiagnostics = a;
if strcmp(run.mode,'research_hermite_C1')
    c = [run.stageTrace{:}]; allDiagnostics = [c{:}];
end
qSwitch = find(diff([a.quadraticCenter]) ~= 0)+1;
hSwitch = find(diff([a.HCell]) ~= 0)+1;
summary = struct('mode',run.mode,'researchEquationId',run.researchEquationId,'nx',numel(run.x),'ny',numel(run.y), ...
    'dt',run.dt,'steps',run.steps,'terminalScale',run.z,'elapsedSeconds',run.elapsedSeconds, ...
    'newRhsEvaluations',run.newRhsEvaluations,'initialH',H(1),'terminalH',H(end), ...
    'terminalRelativeHDrift',(H(end)-H(1))/H(1),'maxRelativeHDrift',max(abs(H-H(1)))/H(1), ...
    'initialQuadraticP',Q(1),'terminalQuadraticP',Q(end), ...
    'maxRelativeQuadraticDrift',max(abs(Q-Q(1)))/Q(1), ...
    'maximumAbsoluteHPrime',max(abs([allDiagnostics.HPrime])), ...
    'maximumAbsoluteQuadraticPPrime',max(abs([allDiagnostics.quadraticPPrime])), ...
    'maximumCwAcceptedDifference',max(abs(diff([a.c_omega]))), ...
    'maximumRelativeAcceptedHIncrement',max(abs(diff(H)))/H(1), ...
    'quadraticSwitchRows',qSwitch,'quadraticSwitchTimes',[a(qSwitch).canonicalTime], ...
    'hermiteCellSwitchRows',hSwitch,'hermiteCellSwitchTimes',[a(hSwitch).canonicalTime], ...
    'invalidHCount',sum(~[allDiagnostics.HValid]),'nonuniqueHCount',sum([allDiagnostics.HActiveCount] ~= 1), ...
    'invalidQuadraticObserverCount',sum(~[allDiagnostics.quadraticValid]), ...
    'maximumHSecondDerivative',max([allDiagnostics.HSecondDerivative]), ...
    'minimumHSecondDerivative',min([allDiagnostics.HSecondDerivative]), ...
    'minimumHFraction',min([allDiagnostics.HFraction]),'maximumHFraction',max([allDiagnostics.HFraction]), ...
    'maximumRealizedCfl',run.dt*max([allDiagnostics.transportRate]), ...
    'maximumPoissonResidual',max([allDiagnostics.poissonResidual]), ...
    'stageDiagnosticsIncluded',strcmp(run.mode,'research_hermite_C1'),'peakProjectionCount',0, ...
    'isNativeCheckpoint',false);
end

function comparison = time_comparison(runs)
fine = runs{end}; pairs = cell(2,1); errors = cell(3,1);
for k = 1:3, errors{k} = field_difference(runs{k},fine); end
for k = 1:2, pairs{k} = field_difference(runs{k},runs{k+1}); end
rhoOrder = log2(pairs{1}.rhoAbsolute./pairs{2}.rhoAbsolute);
omegaOrder = log2(pairs{1}.omegaAbsolute./pairs{2}.omegaAbsolute);
comparison = struct('mode',fine.mode,'nx',numel(fine.x),'dt',cellfun(@(r)r.dt,runs), ...
    'finestNumericalComparisonErrors',{errors},'adjacentSelfDifferences',{pairs}, ...
    'rhoObservedSelfOrderL2Infinity',rhoOrder,'omegaObservedSelfOrderL2Infinity',omegaOrder, ...
    'rhoSelfDifferenceAboveRoundoffScale',pairs{2}.rhoAbsolute(2) > 256*eps*max(1,max(abs(fine.rho),[],'all')), ...
    'interpretation','Three-level self-convergence at common canonical time; finest is not an independently certified exact reference. Orders at roundoff are unresolved.');
end

function difference = field_difference(a,b)
assert(isequal(a.x,b.x) && isequal(a.y,b.y));
w = b.integrationWeights; rhoError = a.rho-b.rho;
omega = b.rho*b.Dx'; omegaError = rhoError*b.Dx';
rhoAbsolute = norms(rhoError,w); omegaAbsolute = norms(omegaError,w);
difference = struct('rhoAbsolute',rhoAbsolute,'rhoRelative',rhoAbsolute./norms(b.rho,w), ...
    'omegaAbsolute',omegaAbsolute,'omegaRelative',omegaAbsolute./norms(omega,w), ...
    'scaleAbsolute',abs(a.z-b.z),'scaleSigned',a.z-b.z);
end

function result = validate_saved_native(run,frameFile)
file = fullfile(fileparts(frameFile),'native_replay_with_snapshots.mat');
saved = load(file,'result'); history = saved.result.history;
tau = history.common.canonicalTau; [distance,row] = min(abs(tau-run.z(5)));
assert(distance < 1e-12 && isequaln(saved.result.snapshots.rho{row},run.rho), ...
    'ipm:PeakOriginalReplay','The new native baseline failed bitwise equality against the original verified replay at tau .08.');
result = struct('sourceFile',file,'row',row,'canonicalTime',tau(row), ...
    'terminalRhoBitwise',true,'entireComparedPrefixRhoBitwise',false,'comparedFrames',numel(run.rhoSnapshots));
for k = 1:numel(run.rhoSnapshots)
    assert(isequaln(saved.result.snapshots.rho{k},run.rhoSnapshots{k}), ...
        'ipm:PeakOriginalReplayPrefix','A native baseline prefix field differs from the original verified replay.');
end
result.entireComparedPrefixRhoBitwise = true;
end

function value = norms(field,w)
value = [sqrt(sum(field.^2.*w,'all')/sum(w,'all')),max(abs(field),[],'all')];
end

function z = pack(s)
z = [s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime];
end

function write_json(file,value)
fid = fopen(file,'w'); assert(fid >= 0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
