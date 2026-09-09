function report = ipm_perflab_verify_lawson_remap(outputDirectory,settings)
%IPM_PERFLAB_VERIFY_LAWSON_REMAP Bounded nonlinear time/space research screen.
% Uses the maintained flow and ipm.solve reference. No alternate WENO,
% checkpoint, feedback, or production method is installed by this experiment.
if nargin < 2
    settings = struct();
end
assert(~exist(outputDirectory,'dir'),'Pass a new output directory.');
defaults = struct('gridSize',[33,17],'runGeneratorExpm',true,'runSpatialScreen',true, ...
    'finalCanonicalTime',0.128,'runLargeStepScreen',true);
names = fieldnames(settings);
assert(all(ismember(names,fieldnames(defaults))),'Unknown research setting.');
for k = 1:numel(names)
    defaults.(names{k}) = settings.(names{k});
end
settings = defaults;
assert(numel(settings.gridSize)==2 && all(settings.gridSize<=65) && ...
    all(settings.gridSize>=9),'Only admitted tiny grids up to 65x65.');
validateattributes(settings.finalCanonicalTime,{'numeric'}, ...
    {'scalar','real','finite','positive'});
mkdir(outputDirectory);
thresholds = struct('minimumTailTimeOrder',1.7,'referenceAdequacy',0.02, ...
    'algebraicRhsDefect',1e-12,'maximumGrid',[65,65], ...
    'consistencyDt',[1e-3,1e-4,1e-5,1e-6,1e-7], ...
    'timeSteps',[4,8,16,32],'finalCanonicalTime',0.128, ...
    'fineTailSteps',[64,128,256], ...
    'referenceSteps',512,'checkReferenceSteps',256, ...
    'largeStepMultipliers',[1,2,4,8], ...
    'largeStepRelativeRhoTolerance',0.01, ...
    'interpretation','Research gates only; no production CFL admission.');
thresholds.gridSize = settings.gridSize;
thresholds.runGeneratorExpm = settings.runGeneratorExpm;
thresholds.finalCanonicalTime = settings.finalCanonicalTime;
thresholds.runLargeStepScreen = settings.runLargeStepScreen;
write_json(fullfile(outputDirectory,'frozen_thresholds.json'),thresholds);
report = struct('schemaVersion',1,'kind','nonlinear_lawson_remap_research', ...
    'thresholds',thresholds,'productionModified',false,'largeGridTested',false);
report.nonlinearCF2 = nonlinear_cf2_mms();
fprintf('Nonlinear dilation CF1 orders %s, CF2 orders %s.\n', ...
    mat2str(report.nonlinearCF2.ordersCF1,5),mat2str(report.nonlinearCF2.ordersCF2,5));
if settings.runSpatialScreen
    report.spatialSplit = spatial_split_screen(thresholds);
else
    report.spatialSplit = {};
end
save(fullfile(outputDirectory,'partial_report.mat'),'report');
report.physical = ipm_case(false,outputDirectory,thresholds);
save(fullfile(outputDirectory,'partial_report.mat'),'report');
if report.physical.correctedTimeOrderPassed
    report.exactGauge = ipm_case(true,outputDirectory,thresholds);
else
    report.exactGauge = struct('performed',false, ...
        'reason','Physical semidiscrete time-order gate did not pass.');
end
report.performedExactGauge = report.exactGauge.performed;
report.integrationApproved = false;
report.interpretation = ['A matched-RHS second-order experiment is distinct ' ...
    'from a stable large-step method. Full-speed WENO remains in F-B.'];
save(fullfile(outputDirectory,'lawson_remap_report.mat'),'report','-v7.3');
write_json(fullfile(outputDirectory,'lawson_remap_report.json'),report);
end

function item = nonlinear_cf2_mms()
% Whole-line dilation: rho_t+a*x*rho_x=0, a=int(rho dx)/sqrt(pi).
% For rho=exp(-(beta*x)^2), a=1/beta, a'=a^2, beta'=-a*beta.
% Start beta=2,a=.5; exact beta=2-t, a=1/(2-t). C(y)=diag(a,-a).
counts = [6,12,24,48];
T = 0.6;
exact = [1/(2-T);2-T];
errors = zeros(2,numel(counts));
profiles = cell(2,numel(counts));
x = linspace(-4,4,129);
for variant = 1:2
    for level = 1:numel(counts)
        dt = T/counts(level);
        u = [0.5;2];
        for step = 1:counts(level)
            coefficient = u(1);
            if variant == 2
                middle = exp(0.5*dt*[coefficient;-coefficient]).*u;
                coefficient = middle(1);
            end
            u = exp(dt*[coefficient;-coefficient]).*u;
        end
        errors(variant,level) = norm(u-exact,Inf);
        profiles{variant,level} = exp(-(u(2)*x).^2);
    end
end
orders = log2(errors(:,1:end-1)./errors(:,2:end));
item = struct('equation','rho_t+a*x*rho_x=0; a=integral(rho)/sqrt(pi)', ...
    'coefficientState','a=1/beta; rho=exp(-(beta*x)^2)', ...
    'frozenLinearOperator','diag(a,-a) on [a,beta]', ...
    'exactTerminalState',exact,'dt',T./counts,'errors',errors, ...
    'ordersCF1',orders(1,:),'ordersCF2',orders(2,:), ...
    'profiles',{profiles},'x',x, ...
    'cf2Passed',all(orders(2,end-1:end)>1.9), ...
    'usesExactSubflow',true,'interpolationErrorExcluded',true, ...
    'initialFrozenCoefficientIsFirstOrder',all(orders(1,:)<1.1));
end

function items = spatial_split_screen(thresholds)
items = {};
for stretch = [0,2]
    for nx = [17,33,65]
        opts = tiny_options(nx,(nx+1)/2,false);
        opts.gridMode = 'stretched';
        opts.gridStretch = [stretch,stretch/2];
        try
            state = ipm.evolve.initialize(opts);
        catch exception
            if startsWith(exception.identifier,'ipm:HighOrder') || ...
                    contains(exception.message,'quadrature norm')
                items{end+1} = struct('grid',[(nx+1)/2,nx],'stretch',stretch, ...
                    'admitted',false,'failure',exception.message); %#ok<AGROW>
                fprintf('Spatial %dx%d s=%g: original mesh gate rejected: %s\n', ...
                    (nx+1)/2,nx,stretch,exception.message);
                continue;
            end
            rethrow(exception);
        end
        ops = state.ops;
        rho = state.rho;
        alpha = 0.2*(1-exp(-(ops.y(:)/0.5).^2));
        b = alpha.*ops.X;
        v1 = b+0.1*sin(ops.X).*exp(-ops.Y.^2);
        v2 = 0.07*ops.Y.*exp(-(ops.X.^2+ops.Y.^2));
        fullRhs = ipm.field.transport(rho,v1,v2,ops);
        residualRhs = ipm.field.transport(rho,v1-b,v2,ops);
        backgroundWeno = ipm.field.transport(rho,b,zeros(size(b)),ops);
        g = ipm_perflab_remap_generator(ops.x,6);
        B = @(r)-alpha.*(r*g.D1.').*ops.x;
        bRho = B(rho);
        rSecond = alpha.^2.*((rho*g.D1.').*ops.x+ ...
            (rho*g.D2.').*ops.x.^2);
        continuum = -(v1.*(-2*ops.X.*rho)+v2.*(-2*ops.Y.*rho));
        dtErrors = zeros(numel(thresholds.consistencyDt),3);
        for k = 1:numel(thresholds.consistencyDt)
            h = thresholds.consistencyDt(k);
            R = ipm_perflab_background_remap(rho,ops.x,alpha,h,'poly6');
            derivative = (R-rho)/h;
            dtErrors(k,:) = [norms(derivative-bRho,ops).', ...
                max(abs(derivative+residualRhs-fullRhs),[],'all')];
        end
        item = struct('grid',[ops.ny,ops.nx],'stretch',stretch,'admitted',true, ...
            'wenoVelocityAdditivityDefect',norms(backgroundWeno+residualRhs-fullRhs,ops), ...
            'remapVersusWenoBackground',norms(bRho-backgroundWeno,ops), ...
            'naiveSplitRhsDefect',norms(bRho+residualRhs-fullRhs,ops), ...
            'fullWenoContinuumError',norms(fullRhs-continuum,ops), ...
            'remapSemigroupSecondDerivativeDefect',norms(rSecond-B(bRho),ops), ...
            'correctedAlgebraicRhsDefect',norms(bRho+(fullRhs-bRho)-fullRhs,ops), ...
            'dt',thresholds.consistencyDt,'differenceQuotientErrors',dtErrors, ...
            'differenceQuotientColumns',{{'remapGeneratorL2','remapGeneratorInf','naiveFullRhsInf'}});
        items{end+1} = item; %#ok<AGROW>
        fprintf('Spatial %dx%d s=%g: WENO additivity %.3e, naive split %.3e, R2-B2 %.3e.\n', ...
            ops.ny,ops.nx,stretch,item.wenoVelocityAdditivityDefect(2), ...
            item.naiveSplitRhsDefect(2),item.remapSemigroupSecondDerivativeDefect(2));
    end
end
end

function item = ipm_case(dynamic,outputDirectory,thresholds)
label = 'physical';
alphaMode = 'prescribed_positive';
if dynamic
    label = 'exact_gauge';
    alphaMode = 'anchor';
end
directory = fullfile(outputDirectory,label);
mkdir(directory);
opts = tiny_options(thresholds.gridSize(1),thresholds.gridSize(2),dynamic);
opts.finalTime = thresholds.finalCanonicalTime;
opts.outputEvery = opts.finalTime/thresholds.referenceSteps;
opts.maxDt = opts.outputEvery;
t = tic;
reference = ipm.solve(opts);
referenceSeconds = toc(t);
assert(abs(reference.state.canonicalTime-opts.finalTime)<1e-12, ...
    'Reference did not reach requested canonical time.');
opts.maxDt = opts.finalTime/thresholds.checkReferenceSteps;
checkReference = ipm.solve(opts);
assert(abs(checkReference.state.canonicalTime-opts.finalTime)<1e-12, ...
    'Check reference did not reach the canonical time.');
save(fullfile(directory,'native_references.mat'),'reference','checkReference','-v7.3');
state = ipm.evolve.initialize(opts);
ops = state.ops;
g = ipm_perflab_remap_generator(ops.x,6);
variants = {'raw','corrected','generator_expm'};
if ~thresholds.runGeneratorExpm
    variants = variants(1:2);
end
variantCount = numel(variants);
counts = thresholds.timeSteps;
errors = zeros(variantCount,numel(counts),4);
scales = zeros(variantCount,numel(counts),5);
seconds = zeros(variantCount,numel(counts));
allRuns = cell(variantCount,numel(counts));
for variant = 1:variantCount
    for level = 1:numel(counts)
        dt = opts.finalTime/counts(level);
        timer = tic;
        run = candidate_run(state,dt,counts(level),g,variants{variant},alphaMode);
        seconds(variant,level) = toc(timer);
        eRho = norms(run.rho-reference.state.rho,ops);
        eOmega = norms((run.rho-reference.state.rho)*ops.Dx.',ops);
        eRhoY = norms(ops.Dy*(run.rho-reference.state.rho),ops);
        errors(variant,level,:) = [eRho.',eOmega(1),eRhoY(1)];
        referenceScale = [log(reference.scale.Cx);log(reference.scale.Comega); ...
            reference.state.physicalTime;reference.scale.Xshift;reference.state.canonicalTime];
        scales(variant,level,:) = pack_scale(run.scale)-referenceScale;
        allRuns{variant,level} = run;
        fprintf('%s %s h=%.6g: rho L2 %.3e, inf %.3e, P drift %.3e, time %.2fs.\n', ...
            label,variants{variant},dt,eRho(1),eRho(2),run.P(end)-run.P(1),seconds(variant,level));
    end
end
orders = log2(errors(:,1:end-1,:)./errors(:,2:end,:));
referenceDiscrepancy = norms(reference.state.rho-checkReference.state.rho,ops);
referenceAdequacy = referenceDiscrepancy./reshape(errors(2,end,1:2),2,1);
semidiscreteTimePass = all(orders(2,end-1:end,1:2)>=thresholds.minimumTailTimeOrder,'all') && ...
    all(referenceAdequacy<=thresholds.referenceAdequacy);
referenceTrusted = all(ipm.output.trustedMask(reference.history,reference.config));
checkReferenceTrusted = all(ipm.output.trustedMask(checkReference.history,checkReference.config));
candidateTrusted = all(cellfun(@(r)r.allTrusted,allRuns(2,:)));
candidateEndpointMatched = all(cellfun(@(r) ...
    abs(r.scale.canonicalTime-opts.finalTime)<1e-12,allRuns(2,:)));
timePass = semidiscreteTimePass && referenceTrusted && checkReferenceTrusted && ...
    candidateTrusted && candidateEndpointMatched;
quotient = zeros(numel(thresholds.consistencyDt),variantCount);
[f0,~] = ipm.evolve.flow(state.rho,ops,state.scale);
for k = 1:numel(thresholds.consistencyDt)
    h = thresholds.consistencyDt(k);
    for v = 1:variantCount
        next = ipm_perflab_lawson_remap_step(state.rho,state.scale,ops,h,g,variants{v},alphaMode);
        difference = (next-state.rho)/h-f0;
        quotient(k,v) = max(abs(difference),[],'all');
    end
end
if thresholds.runLargeStepScreen
    stability = large_step_screen(state,g,alphaMode,thresholds,directory);
else
    stability = {};
end
fineRuns = cell(2,numel(thresholds.fineTailSteps));
fineErrors = zeros(2,numel(thresholds.fineTailSteps),2);
for v = 1:2
    for k = 1:numel(thresholds.fineTailSteps)
        count = thresholds.fineTailSteps(k);
        fineRuns{v,k} = candidate_run(state,opts.finalTime/count,count,g,variants{v},alphaMode);
        fineErrors(v,k,:) = norms(fineRuns{v,k}.rho-reference.state.rho,ops);
    end
end
joinedErrors = cat(2,errors(1:2,end,1:2),fineErrors);
fineOrders = log2(joinedErrors(:,1:end-1,:)./joinedErrors(:,2:end,:));
initial = struct('rho',state.rho,'scale',state.scale,'flow',state.flow, ...
    'x',ops.x,'y',ops.y,'config',state.config);
save(fullfile(directory,'candidate_raw_runs.mat'),'allRuns','fineRuns','initial','g','-v7.3');
item = struct('performed',true,'case',label,'grid',[ops.ny,ops.nx], ...
    'canonicalEndpoint',opts.finalTime,'dt',opts.finalTime./counts, ...
    'variants',{variants},'metricNames',{{'rhoL2','rhoInf','omegaL2','rhoYL2'}}, ...
    'errors',errors,'orders',orders,'signedScaleDeviations',scales,'seconds',seconds, ...
    'referenceSeconds',referenceSeconds,'referenceDiscrepancy',referenceDiscrepancy, ...
    'referenceAdequacy',referenceAdequacy,'correctedTimeOrderPassed',timePass, ...
    'semidiscreteTimeOrderPassed',semidiscreteTimePass, ...
    'candidateAllTrusted',candidateTrusted,'checkReferenceTrusted',checkReferenceTrusted, ...
    'candidateCanonicalEndpointMatched',candidateEndpointMatched, ...
    'derivativeDt',thresholds.consistencyDt,'rhsDifferenceQuotientInf',quotient, ...
    'fineTailSteps',thresholds.fineTailSteps,'fineErrors',fineErrors,'fineOrders',fineOrders, ...
    'largeStep',{stability},'referenceTrusted',referenceTrusted, ...
    'referenceHistory',reference.history,'rawCandidateFile',fullfile(directory,'candidate_raw_runs.mat'), ...
    'gaugeProjectionApplied',false,'restoringApplied',false, ...
    'physicalEndpointMatched',~dynamic,'CflGateAppliedToCandidate',false);
end

function run = candidate_run(state,dt,count,g,variant,alphaMode)
rho = state.rho;
scale = state.scale;
P = zeros(count+1,1);
C = P; PPrime = P; alphaMaximum = zeros(count,1);
canonical = (0:count)'*dt;
physical = zeros(count+1,1);
fullRates = zeros(count,1);
residualRates = zeros(count,1);
researchLog = ipm.output.initializeLog();
[researchLog,~] = ipm.output.record(researchLog,state);
for step = 1:count
    [rhoNew,scaleNew,info] = ipm_perflab_lawson_remap_step(rho,scale, ...
        state.ops,dt,g,variant,alphaMode);
    f = info.initialFlow;
    [P(step),C(step),PPrime(step)] = gauge_values(f);
    physical(step) = scale.physicalTime;
    alphaMaximum(step) = info.alphaMaximum;
    fullRates(step) = max(abs(f.transportU1)./state.ops.hx+ ...
        abs(f.transportU2)./state.ops.hy,[],'all');
    if strcmp(alphaMode,'anchor')
        a = f.c_l+interp1(state.ops.x,f.u1.',1,'linear').';
    else
        a = 0.2*(1-exp(-(state.ops.y(:)/0.5).^2));
    end
    residualRates(step) = max(abs(f.transportU1-a.*state.ops.X)./state.ops.hx+ ...
        abs(f.transportU2)./state.ops.hy,[],'all');
    rho = rhoNew; scale = scaleNew;
    assert(all(isfinite(rho),'all'),'Nonfinite candidate.');
    [~,state.flow] = ipm.evolve.flow(rho,state.ops,scale);
    [~,~,details] = ipm.evolve.selectTimestep(f, ...
        state.scale,state.ops,state.config.time);
    details.dt = dt; details.realizedCfl = dt*details.rate;
    details.activeLimiter = 'research_prescribed_step';
    state.rho = rho; state.scale = scale;
    state.step = step; state.normalizedTime = step*dt; state.timeStep = details;
    [researchLog,~] = ipm.output.record(researchLog,state);
end
[~,f] = ipm.evolve.flow(rho,state.ops,scale);
[P(end),C(end),PPrime(end)] = gauge_values(f);
physical(end) = scale.physicalTime;
[trusted,trustChecks] = ipm.output.trustedMask(researchLog.history,state.config);
run = struct('variant',variant,'dt',dt,'steps',count,'rho',rho,'scale',scale, ...
    'P',P,'C',C,'PPrime',PPrime,'canonicalTime',canonical, ...
    'numericalPDerivative',diff(P)/dt,'integratedRawPPrime',trapz(canonical,PPrime), ...
    'PDriftResidual',P(end)-P(1)-trapz(canonical,PPrime), ...
    'physicalTime',physical,'alphaMaximum',alphaMaximum, ...
    'fullVelocityRate',fullRates,'residualVelocityRate',residualRates, ...
    'maxFullCfl',max(dt*fullRates),'maxResidualVelocityCfl',max(dt*residualRates), ...
    'researchHistory',researchLog.history,'trusted',trusted,'allTrusted',all(trusted), ...
    'trustChecks',trustChecks,'registeredProductionTimeMethod',false, ...
    'gaugeTelemetry','Raw maintained flow evaluated on candidate states; no adjustment.');
end

function items = large_step_screen(state,g,alphaMode,thresholds,directory)
[~,~,d] = ipm.evolve.selectTimestep(state.flow,state.scale,state.ops,state.config.time);
baseDt = 0.5/d.rate;
items = cell(1,numel(thresholds.largeStepMultipliers));
for index = 1:numel(items)
    factor = thresholds.largeStepMultipliers(index);
    h = baseDt*factor;
    count = 1;
    opts = tiny_options(state.ops.nx,state.ops.ny,strcmp(alphaMode,'anchor'));
    opts.finalTime = h*count;
    opts.maxDt = min(0.002,h/32);
    opts.maxSteps = 512;
    opts.outputEvery = opts.maxDt;
    ref = ipm.solve(opts);
    terminalReached = abs(ref.state.canonicalTime-opts.finalTime)<1e-11;
    referenceTrusted = all(ipm.output.trustedMask(ref));
    failure = '';
    try
        run = candidate_run(state,h,count,g,'corrected',alphaMode);
        rawRelativeRho = norms(run.rho-ref.state.rho,state.ops)./norms(ref.state.rho,state.ops);
        relativeRho = rawRelativeRho;
        if ~terminalReached || ~referenceTrusted || ~run.allTrusted
            relativeRho(:) = NaN;
        end
        finite = all(isfinite(run.rho),'all');
    catch exception
        run = struct(); relativeRho = [NaN;NaN]; rawRelativeRho = [NaN;NaN]; finite = false;
        failure = sprintf('%s: %s',exception.identifier,exception.message);
    end
    filename = fullfile(directory,sprintf('large_step_%g.mat',factor));
    save(filename,'run','ref','failure','-v7.3');
    items{index} = struct('initialFullCfl',h*d.rate,'dt',h,'steps',count, ...
        'relativeRhoErrors',relativeRho,'finite',finite,'failure',failure, ...
        'rawUnqualifiedRelativeRhoErrors',rawRelativeRho,'referenceTrusted',referenceTrusted, ...
        'referenceReachedEndpoint',terminalReached,'rawFile',filename, ...
        'smallErrorAtThisHorizon',finite && terminalReached && ...
            all(relativeRho<thresholds.largeStepRelativeRhoTolerance));
    fprintf('Large-step initial CFL %.3g: finite %d, rho error %s.\n', ...
        h*d.rate,finite,mat2str(relativeRho.',4));
end
end

function [P,C,PPrime] = gauge_values(flow)
P = flow.omegaGaugeQuadraticPeakValue;
C = flow.omegaGaugeQuadraticCurvature;
PPrime = flow.omegaGaugeQuadraticForcing+flow.c_omega*P;
end

function opts = tiny_options(nx,ny,dynamic)
opts = struct('nx',nx,'ny',ny,'xlim',[-4,4],'ymax',4, ...
    'gridMode','uniform','gridStretchAutomatic',false,'gridStretch',[0,0], ...
    'initialCondition',@(X,Y)-exp(-X.^2-Y.^2),'rescalingMode','physical', ...
    'symmetryMode','double_odd_omega','spatialDiscretization','high_order', ...
    'transportScheme','weno5_fd','timeIntegrator','ssprk54', ...
    'remeshTransferScheme','high_order','farBoundaryMode','dirichlet_zero', ...
    'transportBoundaryMode','open','wallTransportMode','advective_upwind', ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'saveResults',false,'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'verbose',false,'storeSnapshots',false,'cfl',0.5,'maxDt',0.004, ...
    'finalTime',0.128,'physicalFinalTime',Inf,'outputEvery',0.004, ...
    'rangeStopTolerance',0.1,'minimumPeakPoints',2,'minimumLevelPoints',[2,2,2], ...
    'positiveWallNegativeTolerance',0.99,'oscillationTVTolerance',10);
if dynamic
    opts.rescalingMode = 'dynamic';
    opts.lengthGauge = 'transport_anchor'; opts.transportAnchorX = 1;
    opts.cOmegaGauge = 'wall_omega_quadratic_peak';
end
end

function value = norms(error,ops)
w = ops.integrationWeights;
value = [sqrt(sum(abs(error).^2.*w,'all')/sum(w,'all')); ...
    max(abs(error),[],'all')];
end

function z = pack_scale(s)
z = [s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime];
end

function write_json(filename,value)
fid = fopen(filename,'w');
assert(fid>=0,'Cannot create report.');
cleanup = onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(value,PrettyPrint=true),'char');
end
