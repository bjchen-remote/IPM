function report = ipm_accellab_test_hermite_peak()
%IPM_ACCELLAB_TEST_HERMITE_PEAK MMS and guards for a continuous max gauge.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
z = linspace(-1,1,33); x = 1.8*z+0.2*z.^3; Dx = ipm.mesh.fdMatrix(x,1,7);
window = [-1.5,1.5]; speed = 0.37; sigma = 0.13;
times = [-0.3,-0.01,0,0.01,0.3]; errors = zeros(numel(times),3); cells = zeros(size(times));
for k = 1:numel(times)
    [v,f,P,a] = quadratic(times(k),x,speed,sigma);
    p = ipm_accellab_hermite_peak(v,f,x,Dx,window);
    assert(p.valid,'ipm:HermitePeakQuadraticValidity','Manufactured quadratic peak was rejected.');
    errors(k,:) = abs([p.value-P,p.x-a,p.PPrime-sigma*P]); cells(k) = p.selectedCell;
end
assert(max(errors,[],'all') < 2e-11 && numel(unique(cells)) > 1, ...
    'ipm:HermitePeakQuadratic','Moving quadratic peak or its derivative failed across cells.');

% A Gaussian sampled on a fixed grid has its own discrete peak; do not
% replace that peak or its derivative by the continuum amplitude.
[v,f] = gaussian(0,x,speed,sigma);
p = ipm_accellab_hermite_peak(v,f,x,Dx,window); assert(p.valid);
steps = [1e-2,1e-3,1e-4,1e-5]; derivativeErrors = zeros(size(steps));
for k = 1:numel(steps)
    [vp,fp] = gaussian(steps(k),x,speed,sigma);
    [vm,fm] = gaussian(-steps(k),x,speed,sigma);
    pp = ipm_accellab_hermite_peak(vp,fp,x,Dx,window);
    pm = ipm_accellab_hermite_peak(vm,fm,x,Dx,window);
    assert(pp.valid && pm.valid);
    derivativeErrors(k) = abs((pp.value-pm.value)/(2*steps(k))-p.PPrime);
end
assert(derivativeErrors(end) < 1e-7 && derivativeErrors(end) < derivativeErrors(1)/100, ...
    'ipm:HermitePeakTimeDerivative','The discrete envelope derivative did not match finite differences.');
amplitude = ipm_accellab_hermite_amplitude(v,f,x,Dx,window);
assert(amplitude.valid && abs(amplitude.instantaneousPeakDerivative) < 1e-12);
columnAmplitude = ipm_accellab_hermite_amplitude(v',f,x',Dx,window');
assert(isequaln(columnAmplitude,amplitude),'ipm:HermitePeakShape','Row/column input representations changed the result.');
scaled = ipm_accellab_hermite_peak(2.7*v,-0.43*f,x,Dx,window);
homogeneityErrors = abs([scaled.value-2.7*p.value,scaled.x-p.x,scaled.PPrime+0.43*p.PPrime]);
assert(max(homogeneityErrors) < 1e-11);

% Actual maximizer crossing a native node: values remain continuous even
% though the selected cell changes and second derivatives need not match.
crossSteps = [1e-2,1e-4,1e-6]; continuity = zeros(numel(crossSteps),3);
for k = 1:numel(crossSteps)
    h = crossSteps(k); [vp,fp] = quadratic(h,x,speed,sigma); [vm,fm] = quadratic(-h,x,speed,sigma);
    pp = ipm_accellab_hermite_peak(vp,fp,x,Dx,window); pm = ipm_accellab_hermite_peak(vm,fm,x,Dx,window);
    assert(pp.selectedCell ~= pm.selectedCell && pp.valid && pm.valid);
    continuity(k,:) = [abs(pp.value-pm.value),abs(pp.PPrime-pm.PPrime),abs(pp.x-pm.x)];
end
assert(max(continuity(end,:)) < max(continuity(1,:))/1000);

flat = ipm_accellab_hermite_peak(ones(size(x)),x,x,Dx,window);
boundary = ipm_accellab_hermite_peak(3+x,x,x,Dx,window);
% Symmetric separated peaks give equal isolated interior maxima. The
% optional Dini branch is tested separately from the default rejection.
twin = 2-(x.^2-0.5^2).^2;
tied = ipm_accellab_hermite_peak(twin,x,x,Dx,window);
dini = ipm_accellab_hermite_peak(twin,x,x,Dx,window,struct('allowIsolatedTies',true));
assert(~flat.valid && any(strcmp(flat.invalidReasons,'flat_or_degenerate_maximum')));
assert(~boundary.valid && any(strcmp(boundary.invalidReasons,'maximum_at_search_boundary')));
assert(~tied.valid && any(strcmp(tied.invalidReasons,'multiple_or_near_tied_maxima')));
assert(dini.valid && dini.isolatedPeakCount == 2 && dini.PPrime > 0);
diniSteps = [1e-3,1e-4,1e-5]; diniErrors = zeros(size(diniSteps));
for k = 1:numel(diniSteps)
    ph = ipm_accellab_hermite_peak(twin+diniSteps(k)*x,x,x,Dx,window);
    assert(ph.valid); diniErrors(k) = abs((ph.value-dini.value)/diniSteps(k)-dini.PPrime);
end
assert(diniErrors(end) < 1e-5 && diniErrors(end) < diniErrors(1)/50);
report = struct('status','passed','newRhsEvaluations',0,'newPdeSteps',0,'poissonOperatorBuilds',0, ...
    'quadraticTimes',times,'quadraticValuePositionDerivativeErrors',errors,'quadraticSelectedCells',cells, ...
    'gaussianPeak',p,'gaussianTimeSteps',steps,'gaussianDerivativeErrors',derivativeErrors, ...
    'amplitudeCancellation',amplitude,'homogeneityErrors',homogeneityErrors, ...
    'crossCellTimeSteps',crossSteps,'crossCellValueDerivativePositionDifferences',continuity, ...
    'flatGuard',flat,'boundaryGuard',boundary,'tieGuard',tied,'diniExample',dini, ...
    'diniTimeSteps',diniSteps,'diniForwardDerivativeErrors',diniErrors, ...
    'interpretation','Manufactured field paths only. No new research RHS was time integrated, and instantaneous cancellation is not finite-step conservation.');
end

function [v,f,P,a] = quadratic(t,x,speed,sigma)
P = exp(sigma*t); a = speed*t; q = x-a;
v = P*(1-0.4*q.^2); f = sigma*v+0.8*P*speed*q;
end

function [v,f] = gaussian(t,x,speed,sigma)
a = 0.173+speed*t; q = (x-a)/0.43;
v = exp(sigma*t-q.^2); f = sigma*v+2*speed/0.43*q.*v;
end
