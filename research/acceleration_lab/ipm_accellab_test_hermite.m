function report = ipm_accellab_test_hermite()
%IPM_ACCELLAB_TEST_HERMITE Linearity, polynomial, C1 and derivative accuracy.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
[x,y,Dx,Dy] = grid(31,25); [X,Y] = meshgrid(x,y);
[QX,QY] = meshgrid(linspace(-0.83,0.92,27),linspace(0.01,0.97,23));
[field,~,~] = polynomial(X,Y);
[truth,truthX,truthY] = polynomial(QX,QY);
observed = ipm_accellab_tensor_hermite(field,x,y,Dx,Dy,QX,QY);
polynomialErrors = [infinity(observed.value-truth), ...
    infinity(observed.derivativeX-truthX),infinity(observed.derivativeY-truthY)];
assert(max(polynomialErrors) < 2e-10,'ipm:HermitePolynomial', ...
    'Tensor cubic polynomial value or derivative was not reproduced.');
other = sin(1.3*X-0.7*Y);
a = ipm_accellab_tensor_hermite(other,x,y,Dx,Dy,QX,QY);
b = ipm_accellab_tensor_hermite(field+0.43*other,x,y,Dx,Dy,QX,QY);
linearityError = max([infinity(b.value-observed.value-0.43*a.value), ...
    infinity(b.derivativeX-observed.derivativeX-0.43*a.derivativeX), ...
    infinity(b.derivativeY-observed.derivativeY-0.43*a.derivativeY)]);
assert(linearityError < 1e-11,'ipm:HermiteLinearity','The observer is not linear in its field.');

xi = linspace(-2,2,65); eta = linspace(0,3,49);
[x,y,Dx,Dy] = grid(81,65);
[omega,forcing] = moving(0,x,y);
[~,index] = max(omega(1,:));
exact = ipm_accellab_exact_inner_rates(omega,forcing,x,y,index);
chart = ipm_accellab_pullback_hermite(omega,forcing,x,y,Dx,Dy,exact,xi,eta);
timeSteps = [1e-3,3e-4,1e-4,3e-5]; timeErrors = zeros(size(timeSteps));
sameCoordinateTemplate = false(size(timeSteps));
stableCellFractions = zeros(size(timeSteps));
for k = 1:numel(timeSteps)
    [plus,cp] = moving_chart(timeSteps(k),x,y,Dx,Dy,xi,eta);
    [minus,cm] = moving_chart(-timeSteps(k),x,y,Dx,Dy,xi,eta);
    sameCoordinateTemplate(k) = isequal(cp.signature,exact.signature) && isequal(cm.signature,exact.signature);
    timeErrors(k) = infinity((plus.U-minus.U)/(2*timeSteps(k))-chart.G);
    stable = plus.cellRows == minus.cellRows & plus.cellColumns == minus.cellColumns;
    stableCellFractions(k) = nnz(stable)/numel(stable);
end
assert(all(sameCoordinateTemplate) && timeErrors(end) < 1e-7 && ...
    timeErrors(end) < 0.1*timeErrors(1),'ipm:HermiteTimeDerivative', ...
    'The full-observer central time difference did not approach its same-interpolant derivative.');
amplitude = ipm_accellab_exact_inner_rates(omega,0.07*omega,x,y,index);
amplitudeChart = ipm_accellab_pullback_hermite(omega,0.07*omega,x,y,Dx,Dy,amplitude,xi,eta);
amplitudeError = infinity(amplitudeChart.G);
assert(amplitudeError < 1e-11,'ipm:HermiteAmplitude','A pure amplitude mode remains after exact normalization.');

[X,Y] = meshgrid(x,y); field = sin(3*X).*exp(-0.7*Y)+0.2*cos(2*Y+X);
epsilon = min([diff(x),diff(y)'])*[1e-3,1e-5,1e-7];
interfaceJumps = zeros(3,3);
for k = 1:numel(epsilon)
    edgeX = x(41); edgeY = y(33); e = epsilon(k);
    left = ipm_accellab_tensor_hermite(field,x,y,Dx,Dy,[edgeX-e,0.173],[0.413,edgeY-e]);
    right = ipm_accellab_tensor_hermite(field,x,y,Dx,Dy,[edgeX+e,0.173],[0.413,edgeY+e]);
    interfaceJumps(k,:) = [infinity(right.value-left.value), ...
        infinity(right.derivativeX-left.derivativeX),infinity(right.derivativeY-left.derivativeY)];
end
assert(max(interfaceJumps(end,:)) < 1e-6 && ...
    max(interfaceJumps(end,:)) < 0.001*max(interfaceJumps(1,:)), ...
    'ipm:HermiteC1','The two one-sided value/gradient limits did not agree.');
manufactured = struct('valid',true,'peak',struct('x',x(41),'value',1), ...
    'wallCoreWidth',0.1,'verticalCoreWidth',0.1,'translationRate',0.031, ...
    'wallCoreWidthPrime',0,'verticalCoreWidthPrime',0,'peakPrime',0);
center = ipm_accellab_pullback_hermite(field,zeros(size(field)),x,y,Dx,Dy,manufactured,[-1,0,1],[0,1,2]);
crossingSteps = [1e-2,1e-3,1e-4]; crossingErrors = zeros(size(crossingSteps));
for k = 1:numel(crossingSteps)
    plus = manufactured; minus = manufactured;
    plus.peak.x = plus.peak.x+plus.translationRate*crossingSteps(k);
    minus.peak.x = minus.peak.x-minus.translationRate*crossingSteps(k);
    a = ipm_accellab_pullback_hermite(field,zeros(size(field)),x,y,Dx,Dy,plus,[-1,0,1],[0,1,2]);
    b = ipm_accellab_pullback_hermite(field,zeros(size(field)),x,y,Dx,Dy,minus,[-1,0,1],[0,1,2]);
    assert(all(a.cellColumns(:,2) ~= b.cellColumns(:,2)), ...
        'ipm:HermiteCrossingTest','The manufactured moving query did not cross an interface.');
    crossingErrors(k) = infinity((a.U-b.U)/(2*crossingSteps(k))-center.G);
end
assert(crossingErrors(end) < 1e-6 && crossingErrors(end) < 0.1*crossingErrors(1), ...
    'ipm:HermiteCrossingDerivative','A genuinely cross-cell time difference did not converge.');

nx = [21,41,81,161]; ny = [17,33,65,129];
spatialErrors = zeros(4,3); linearValueErrors = zeros(1,4); maximumSpacing = zeros(1,4);
[QX,QY] = meshgrid(linspace(-0.73,0.81,37),linspace(0.025,0.93,31));
[truth,truthX,truthY] = smooth(QX,QY);
for k = 1:4
    [x,y,Dx,Dy] = grid(nx(k),ny(k)); [X,Y] = meshgrid(x,y); field = smooth(X,Y);
    value = ipm_accellab_tensor_hermite(field,x,y,Dx,Dy,QX,QY);
    spatialErrors(k,:) = [infinity(value.value-truth), ...
        infinity(value.derivativeX-truthX),infinity(value.derivativeY-truthY)];
    linear = interp2(x,y,field,QX,QY,'linear');
    linearValueErrors(k) = infinity(linear-truth);
    maximumSpacing(k) = max([diff(x),diff(y)']);
end
assert(spatialErrors(end,1) < spatialErrors(1,1)/100 && ...
    all(spatialErrors(end,2:3) < spatialErrors(1,2:3)/20), ...
    'ipm:HermiteRefinement','Smooth-function value and gradient errors did not decrease with refinement.');
orders = log(spatialErrors(1:end-1,:)./spatialErrors(2:end,:))./ ...
    log(maximumSpacing(1:end-1)'./maximumSpacing(2:end)');
report = struct('status','passed','newRhsEvaluations',0,'newPdeSteps',0,'operatorBuilds',0, ...
    'polynomialValueDxDyErrors',polynomialErrors,'linearityError',linearityError, ...
    'timeSteps',timeSteps,'fullObserverTimeDerivativeErrors',timeErrors, ...
    'sameCoordinateTemplates',sameCoordinateTemplate,'stableObserverCellFractions',stableCellFractions, ...
    'pureAmplitudeResidual',amplitudeError,'interfaceEpsilon',epsilon, ...
    'interfaceValueDxDyJumps',interfaceJumps,'crossCellTimeSteps',crossingSteps, ...
    'crossCellTimeDerivativeErrors',crossingErrors,'gridNx',nx,'gridNy',ny, ...
    'smoothValueDxDyErrors',spatialErrors,'smoothValueDxDyOrders',orders, ...
    'bilinearSmoothValueErrors',linearValueErrors, ...
    'interpretation','Verification of a linear C1 observer only; no original candidate rejection is changed and no q512 interpolation-error explanation is established.');
fprintf('HERMITE PASS polynomial %.3e, time derivative %.3e, final C1 jump %.3e.\n', ...
    max(polynomialErrors),timeErrors(end),max(interfaceJumps(end,:)));
end

function [x,y,Dx,Dy] = grid(nx,ny)
x = linspace(-1.2,1.3,nx); x = x+0.01*sin(1.8*x);
y = linspace(0,1.2,ny)'; y = y+0.006*sin(2*y);
Dx = ipm.mesh.fdMatrix(x,1,7); Dy = ipm.mesh.fdMatrix(y,1,7);
end

function [f,fx,fy] = polynomial(x,y)
f = 1+0.2*x-0.3*y+0.4*x.^2.*y+0.1*x.^3.*y.^2-0.07*x.^2.*y.^3+0.05*x.^3.*y.^3;
fx = 0.2+0.8*x.*y+0.3*x.^2.*y.^2-0.14*x.*y.^3+0.15*x.^2.*y.^3;
fy = -0.3+0.4*x.^2+0.2*x.^3.*y-0.21*x.^2.*y.^2+0.15*x.^3.*y.^2;
end

function [f,fx,fy] = smooth(x,y)
a = exp(-2*x.^2-0.7*y); theta = 0.6*x+0.9*y;
f = a.*cos(theta); fx = a.*(-4*x.*cos(theta)-0.6*sin(theta));
fy = a.*(-0.7*cos(theta)-0.9*sin(theta));
end

function [f,forcing] = moving(time,x,y)
aPrime = 0.031; beta = -0.65; gamma = -0.67; sigma = 0.07;
a = 0.137+aPrime*time; lx = 0.45*exp(beta*time); ly = 0.50*exp(gamma*time);
[X,Y] = meshgrid(x,y); z = (X-a)/lx; eta = Y/ly;
f = 1.8*exp(sigma*time)*exp(-z.^2-eta+0.35*z.*eta);
fx = f.*(-2*z+0.35*eta)/lx; fy = f.*(-1+0.35*z)/ly;
forcing = sigma*f-aPrime*fx-beta*(X-a).*fx-gamma*Y.*fy;
end

function [chart,coordinates] = moving_chart(time,x,y,Dx,Dy,xi,eta)
omega = moving(time,x,y); [~,index] = max(omega(1,:));
coordinates = ipm_accellab_exact_inner_rates(omega,zeros(size(omega)),x,y,index);
chart = ipm_accellab_pullback_hermite(omega,zeros(size(omega)),x,y,Dx,Dy,coordinates,xi,eta);
end

function value = infinity(field)
value = max(abs(field),[],'all');
end
