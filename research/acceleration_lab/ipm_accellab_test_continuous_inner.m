function report = ipm_accellab_test_continuous_inner(frameFile,outputRoot)
%IPM_ACCELLAB_TEST_CONTINUOUS_INNER Exact polynomial and tiny true-PDE tests.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['continuous_inner_tiny_',token]); mkdir(destination);
registration = struct('outputDirectory',destination,'nativeOperatorConfigSource',frameFile, ...
    'timeSteps',[1e-3,5e-4,2.5e-4,1.25e-4], ...
    'coordinateUnitPairs',[1e-3,1e3;1e3,1e-3],'amplitudePerturbation',.37, ...
    'manufacturedRevision','interior width-root phase and larger analytic shape derivative; earlier failed near-node/roundoff test is retained separately', ...
    'peakInterfacePolicy','withhold phase derivative at/within local ULP tolerance of native interfaces');
save(fullfile(destination,'registration.mat'),'registration');
evidence = struct();
try
    z = linspace(-1,1,49); x = 1.4*z+.1*z.^3;
    y = linspace(0,1,41)'.^1.2; Dx = ipm.mesh.fdMatrix(x,1,7); Dy = ipm.mesh.fdMatrix(y,1,7);
    xi = linspace(-2,2,41); eta = linspace(0,3,31); window = [-1.2,1.2];
    parameters = struct('a0',.197,'velocity',.07,'Lx',.55,'Ly',.4, ...
        'beta',-.18,'gamma',-.23,'sigma',.11,'delta',.035,'deltaPrime',.6,'coupling',.35);
    [omega,forcing,truth] = polynomial(0,x,y,parameters);
    rates = ipm_accellab_continuous_inner_rates(omega,forcing,x,y,Dx,Dy,window);
    assert(rates.valid);
    chart = ipm_accellab_pullback_hermite(omega,forcing,x,y,Dx,Dy,rates,xi,eta);
    exactG = polynomial_chart_derivative(xi,eta,truth,parameters);
    rateErrors = abs([rates.peak.value-truth.P,rates.peakPrime-truth.PPrime, ...
        rates.peak.x-truth.a,rates.translationRate-truth.aPrime, ...
        rates.wallCoreWidth-truth.wx,rates.wallCoreWidthPrime-truth.wxPrime, ...
        rates.verticalCoreWidth-truth.wy,rates.verticalCoreWidthPrime-truth.wyPrime]);
    analyticGError = max(abs(chart.G-exactG),[],'all');
    evidence.analyticRateErrors = rateErrors; evidence.analyticGError = analyticGError;
    evidence.parameters = parameters;
    assert(max(rateErrors) < 2e-10 && analyticGError < 2e-10 && abs(rates.peakPrime) > .05, ...
        'ipm:ContinuousInnerPolynomial','Analytic nonzero-amplitude/shape derivatives failed.');
    assert(abs(rates.verticalTranslationForcing) > .001, ...
        'ipm:ContinuousInnerMovingTrace','The test did not exercise the vertical moving-X forcing term.');
    modified = forcing+registration.amplitudePerturbation*omega;
    modifiedRates = ipm_accellab_continuous_inner_rates(omega,modified,x,y,Dx,Dy,window);
    modifiedChart = ipm_accellab_pullback_hermite(omega,modified,x,y,Dx,Dy,modifiedRates,xi,eta);
    amplitudeErrors = abs([modifiedRates.translationRate-rates.translationRate, ...
        modifiedRates.logScaleXRate-rates.logScaleXRate,modifiedRates.logScaleYRate-rates.logScaleYRate, ...
        modifiedRates.peakPrime-rates.peakPrime-registration.amplitudePerturbation*rates.peak.value, ...
        max(abs(modifiedChart.G-chart.G),[],'all')]);
    assert(max(amplitudeErrors) < 2e-10,'ipm:ContinuousInnerAmplitude','A pure amplitude forcing changed geometric/quotient rates.');
    evidence.amplitudeErrors = amplitudeErrors;
    unitErrors = zeros(2,5);
    for k = 1:2
        sx = registration.coordinateUnitPairs(k,1); sy = registration.coordinateUnitPairs(k,2);
        q = ipm_accellab_continuous_inner_rates(omega,forcing,sx*x,sy*y,Dx/sx,Dy/sy,sx*window);
        assert(q.valid); c = ipm_accellab_pullback_hermite(omega,forcing,sx*x,sy*y,Dx/sx,Dy/sy,q,xi,eta);
        unitErrors(k,:) = abs([q.translationRate/sx-rates.translationRate, ...
            q.wallCoreWidthPrime/sx-rates.wallCoreWidthPrime,q.verticalCoreWidthPrime/sy-rates.verticalCoreWidthPrime, ...
            max(abs(c.U-chart.U),[],'all'),max(abs(c.G-chart.G),[],'all')]);
    end
    assert(max(unitErrors,[],'all') < 2e-10,'ipm:ContinuousInnerUnits','Coordinate units changed the normalized shape derivative.');
    evidence.unitErrors = unitErrors;
    timeErrors = zeros(size(registration.timeSteps)); sameTemplates = false(size(timeErrors));
    for k = 1:numel(timeErrors)
        h = registration.timeSteps(k);
        [plus,qp] = manufactured_chart(h,x,y,Dx,Dy,window,xi,eta,parameters);
        [minus,qm] = manufactured_chart(-h,x,y,Dx,Dy,window,xi,eta,parameters);
        sameTemplates(k) = isequal(qp.signature,rates.signature) && isequal(qm.signature,rates.signature);
        timeErrors(k) = max(abs((plus.U-minus.U)/(2*h)-chart.G),[],'all');
    end
    evidence.manufacturedTimeErrors = timeErrors; evidence.manufacturedSameTemplates = sameTemplates;
    evidence.initialGeometry = rates;
    fprintf('CONTINUOUS_INNER_MMS errors %s fixed_templates %s\n',mat2str(timeErrors,12),mat2str(sameTemplates));
    assert(all(sameTemplates) && timeErrors(end) < 1e-7 && timeErrors(end) < timeErrors(1)/10, ...
        'ipm:ContinuousInnerTimeDifference','Manufactured full-profile time differences did not converge.');

    % The phase derivative at a native interface is explicitly withheld,
    % even in this special globally smooth polynomial example.
    eventParameters = parameters; eventParameters.a0 = x(25);
    [o,f] = polynomial(0,x,y,eventParameters);
    phaseEvent = ipm_accellab_continuous_inner_rates(o,f,x,y,Dx,Dy,window);
    assert(phaseEvent.geometryValid && ~phaseEvent.valid && phaseEvent.phase.atOrNumericallyNearNativeInterface, ...
        'ipm:ContinuousInnerPhaseEvent','The native-interface phase limitation was not exposed.');
    sideTimes = [-1e-3,-1e-4,1e-4,1e-3]; sideRates = zeros(size(sideTimes)); sideCells = zeros(size(sideTimes));
    for k = 1:numel(sideTimes)
        [o,f] = polynomial(sideTimes(k),x,y,eventParameters);
        q = ipm_accellab_continuous_inner_rates(o,f,x,y,Dx,Dy,window);
        assert(q.valid); sideRates(k) = q.translationRate; sideCells(k) = q.peak.selectedCell;
    end
    assert(numel(unique(sideCells)) == 2 && max(abs(sideRates-parameters.velocity)) < 1e-10);

    % A transverse width root may cross a native interface: its first
    % derivative exists because the underlying value and first jet are C1.
    rootParameters = parameters;
    [~,~,t0] = polynomial(0,x,y,rootParameters);
    rootParameters.Lx = (x(30)-parameters.a0)/t0.qRight;
    rootParameters.Ly = y(10)/sqrt(.1);
    [o,f,truthRoot] = polynomial(0,x,y,rootParameters);
    rootEvent = ipm_accellab_continuous_inner_rates(o,f,x,y,Dx,Dy,window);
    assert(rootEvent.valid);
    rootEventChart = ipm_accellab_pullback_hermite(o,f,x,y,Dx,Dy,rootEvent,xi,eta);
    rootEventError = max(abs(rootEventChart.G-polynomial_chart_derivative(xi,eta,truthRoot,rootParameters)),[],'all');
    [~,before] = manufactured_chart(-1e-4,x,y,Dx,Dy,window,xi,eta,rootParameters);
    [~,after] = manufactured_chart(1e-4,x,y,Dx,Dy,window,xi,eta,rootParameters);
    assert(rootEventError < 2e-10 && ~isequal(before.signature,after.signature), ...
        'ipm:ContinuousInnerRootEvent','The transverse-root interface test did not exercise a cell crossing.');

    native = native_time_test(frameFile,registration.timeSteps,xi,eta);
    report = struct('status','passed','registration',registration,'analyticRateErrors',rateErrors, ...
        'analyticGError',analyticGError,'polynomialRates',rates,'amplitudeInvarianceErrors',amplitudeErrors, ...
        'coordinateUnitErrors',unitErrors,'manufacturedTimeErrors',timeErrors, ...
        'manufacturedTemplatesFixed',sameTemplates,'phaseInterfaceEvent',phaseEvent, ...
        'phaseSideTimes',sideTimes,'phaseSideRates',sideRates,'phaseSideCells',sideCells, ...
        'transverseRootInterfaceEvent',rootEvent,'transverseRootInterfaceAnalyticGError',rootEventError, ...
        'nativePdeTimeDifference',native, ...
        'interpretation','The same H1/H2 geometry and genuine forcing define every rate and U_tau. Phase-interface derivatives are withheld; polynomial crossing accuracy does not imply general C2 smoothness. No mainline gauge or trajectory changed.');
    save(fullfile(destination,'report.mat'),'report'); write_json(fullfile(destination,'report.json'),report);
catch exception
    failure = struct('registration',registration,'evidence',evidence,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure'); write_json(fullfile(destination,'failure.json'),failure);
    fprintf('CONTINUOUS_INNER_TINY_FAILED %s\n',destination); rethrow(exception);
end
fprintf('CONTINUOUS_INNER_TINY %s\n',destination);
end

function result = native_time_test(frameFile,steps,xi,eta)
data = load(frameFile,'config'); config = data.config; schema = ipm.config.schema(); opts = struct();
for k = 1:numel(schema.domainNames)
    d = config.(schema.domainNames{k}); names = fieldnames(d);
    for j = 1:numel(names), opts.(names{j}) = d.(names{j}); end
end
opts.saveResults = false; opts.storeSnapshots = false; opts.verbose = false;
s = ipm.evolve.initialize(opts); ops = s.ops; assert(ops.nx == 65 && ops.ny == 33);
omega = s.flow.source; forcing = s.rhsCache.rhoRate*ops.Dx'; window = [0,ops.x(end)];
r = ipm_accellab_continuous_inner_rates(omega,forcing,ops.x,ops.y,ops.Dx,ops.Dy,window);
assert(r.valid && abs(r.peakPrime) > 1e-5);
c = ipm_accellab_pullback_hermite(omega,forcing,ops.x,ops.y,ops.Dx,ops.Dy,r,xi,eta);
errors = zeros(size(steps)); geometryErrors = zeros(numel(steps),4); same = false(size(steps));
rhsEvaluations = 1;
for k = 1:numel(steps)
    h = steps(k);
    [rho1,flow1,scale1,cache1,info1] = ipm.evolve.stepSsprk54(s.rho,h,ops,s.scale,s.rhsCache);
    [~,flow2,~,cache2,info2] = ipm.evolve.stepSsprk54(rho1,h,ops,scale1,cache1);
    rhsEvaluations = rhsEvaluations+info1.rhsEvaluations+info2.rhsEvaluations;
    f1 = cache1.rhoRate*ops.Dx'; f2 = cache2.rhoRate*ops.Dx';
    r1 = ipm_accellab_continuous_inner_rates(flow1.source,f1,ops.x,ops.y,ops.Dx,ops.Dy,window);
    r2 = ipm_accellab_continuous_inner_rates(flow2.source,f2,ops.x,ops.y,ops.Dx,ops.Dy,window);
    assert(r1.valid && r2.valid);
    c1 = ipm_accellab_pullback_hermite(flow1.source,f1,ops.x,ops.y,ops.Dx,ops.Dy,r1,xi,eta);
    c2 = ipm_accellab_pullback_hermite(flow2.source,f2,ops.x,ops.y,ops.Dx,ops.Dy,r2,xi,eta);
    errors(k) = max(abs((-3*c.U+4*c1.U-c2.U)/(2*h)-c.G),[],'all');
    derivative = [r.peakPrime,r.translationRate,r.wallCoreWidthPrime,r.verticalCoreWidthPrime];
    geometryErrors(k,:) = abs((-3*geometry(r)+4*geometry(r1)-geometry(r2))/(2*h)-derivative);
    same(k) = isequal(r.signature,r1.signature) && isequal(r.signature,r2.signature);
end
assert(all(same) && errors(end) < 1e-6 && errors(end) < errors(1)/10, ...
    'ipm:ContinuousInnerTruePdeDifference','True native-PDE profile differences did not approach the registered derivative.');
result = struct('gridSize',size(s.rho),'initialRates',r,'timeSteps',steps,'forwardSecondOrderGErrors',errors, ...
    'forwardSecondOrderGeometryErrors',geometryErrors,'sameGeometryTemplates',same, ...
    'newRhsEvaluations',rhsEvaluations,'newTinyPdeSteps',2*numel(steps),'tinyPoissonOperatorBuilds',1, ...
    'method','Two original SSPRK54 forward steps at h and 2h; (-3U0+4Uh-U2h)/(2h). No backward physical-time step or altered gauge.');
end

function [c,r] = manufactured_chart(t,x,y,Dx,Dy,window,xi,eta,p)
[o,f] = polynomial(t,x,y,p); r = ipm_accellab_continuous_inner_rates(o,f,x,y,Dx,Dy,window);
assert(r.valid); c = ipm_accellab_pullback_hermite(o,f,x,y,Dx,Dy,r,xi,eta);
end

function [o,f,truth] = polynomial(t,x,y,p)
[X,Y] = meshgrid(x,y); A = exp(p.sigma*t); a = p.a0+p.velocity*t;
Lx = p.Lx*exp(p.beta*t); Ly = p.Ly*exp(p.gamma*t); delta = p.delta+p.deltaPrime*t;
q = (X-a)/Lx; r = Y/Ly; h = 1-q.^2+delta*q.^3; g = 1-r.^2;
hq = -2*q+3*delta*q.^2; value = h.*g+p.coupling*q.*r.^2;
qx = -p.velocity/Lx-p.beta*q; ry = -p.gamma*r;
o = A*value;
f = p.sigma*o+A*(p.deltaPrime*q.^3.*g+(hq.*g+p.coupling*r.^2).*qx+ ...
    (-2*r.*h+2*p.coupling*q.*r).*ry);
roots0 = roots([delta,-1,0,.1]); roots0 = sort(real(roots0(abs(imag(roots0)) < 1e-12)));
ql = roots0(find(roots0 < 0,1,'last')); qr = roots0(find(roots0 > 0,1,'first'));
dql = -p.deltaPrime*ql^3/(-2*ql+3*delta*ql^2);
dqr = -p.deltaPrime*qr^3/(-2*qr+3*delta*qr^2);
width = Lx*(qr-ql); height = Ly*sqrt(.1);
truth = struct('P',A,'PPrime',p.sigma*A,'a',a,'aPrime',p.velocity, ...
    'wx',width,'wxPrime',p.beta*width+Lx*(dqr-dql),'wy',height,'wyPrime',p.gamma*height, ...
    'qLeft',ql,'qRight',qr,'qDifferencePrime',dqr-dql,'delta',delta);
end

function G = polynomial_chart_derivative(xi,eta,t,p)
[XI,ETA] = meshgrid(xi,eta); q = (t.qRight-t.qLeft)*XI; r = sqrt(.1)*ETA;
G = p.deltaPrime*q.^3.*(1-r.^2)+t.qDifferencePrime*XI.* ...
    ((-2*q+3*t.delta*q.^2).*(1-r.^2)+p.coupling*r.^2);
end

function v = geometry(r)
v = [r.peak.value,r.peak.x,r.wallCoreWidth,r.verticalCoreWidth];
end

function write_json(file,value)
fid = fopen(file,'w'); assert(fid >= 0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
