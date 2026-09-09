function report = ipm_accellab_test_observer_sampling()
%IPM_ACCELLAB_TEST_OBSERVER_SAMPLING Exact polynomial integral and all phases.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
x = linspace(-1.3,1.3,33); x = x+0.013*sin(2*x);
y = linspace(0,1,29)'; y = y+0.006*sin(2*y);
Dx = ipm.mesh.fdMatrix(x,1,7); Dy = ipm.mesh.fdMatrix(y,1,7);
[X,Y] = meshgrid(x,y);
coeff = zeros(4); coeff(1,1) = 1; coeff(2,1) = 0.08;
coeff(1,2) = 0.02; coeff(2,2) = 0.024;
coeff(4,3) = 0.013; coeff(3,4) = -0.009; coeff(4,4) = 0.006;
forcing = polynomial(coeff,X/0.4,Y/0.2);
coordinates = struct('valid',true,'kind','manufactured_fixed_chart_for_quadrature_test', ...
    'peak',struct('x',0,'value',1),'wallCoreWidth',0.4,'verticalCoreWidth',0.2, ...
    'translationRate',0,'wallCoreWidthPrime',0,'verticalCoreWidthPrime',0,'peakPrime',0);
baseline = struct('omega',ones(size(X)),'forcing',forcing,'exactCoordinates',coordinates);
candidate = baseline; candidate.forcing = 1.1*forcing;
[study,~] = ipm_accellab_observer_sampling({baseline,candidate},x,y,Dx,Dy);
integralCore = squared_integral(coeff,[-1,1],[0,1.5]);
integralAll = squared_integral(coeff,[-2,2],[0,3]);
truth = [sqrt(integralCore/3),sqrt((integralAll-integralCore)/9)];
q = study.baseline.splitCellGaussReference.hermite;
error = max(abs([q.trainingL2,q.holdoutL2]-truth));
assert(error < 1e-11,'ipm:SamplingExactIntegral','Split-cell Hermite Gauss norm differs from the analytic polynomial norm.');
for k = 1:numel(study.candidateToBaselineSamplingRatios)
    for method = {'bilinear','hermite'}
        values = struct2array(study.candidateToBaselineSamplingRatios{k}.(method{1}));
        assert(max(abs(values-1.1)) < 1e-10,'ipm:SamplingLinearity','A registered phase/density has inconsistent fixed-field scaling.');
    end
end
for k = 1:numel(study.baseline.nativeRestriction)
    d = study.baseline.nativeRestriction{k}.hermite.derivativeDifferenceFromNative;
    assert(max([d.trainingL2,d.holdoutL2]) < 1e-10, ...
        'ipm:SamplingRestrictionPolynomial','Polynomial Hermite restriction has unexpected error.');
end
assert(~study.accepted && strcmp(study.originalDecision,'rejected_unchanged'), ...
    'ipm:SamplingDecision','Observer diagnostics changed a prior decision.');
report = struct('status','passed','analyticHermiteL2Error',error, ...
    'allPhasesAndDensitiesReported',numel(study.baseline.sampling), ...
    'allRestrictionStridesReported',numel(study.baseline.nativeRestriction), ...
    'study',study,'newRhsEvaluations',0,'newPdeSteps',0,'poissonOperatorBuilds',0);
fprintf('SAMPLING PASS exact Hermite L2 error %.3e; %d registered sampling variants.\n', ...
    error,numel(study.baseline.sampling));
end

function f = polynomial(c,x,y)
f = zeros(size(x));
for i = 1:size(c,1)
    for j = 1:size(c,2)
        f = f+c(i,j)*x.^(i-1).*y.^(j-1);
    end
end
end

function result = squared_integral(c,x,y)
result = 0;
for i = 1:size(c,1)
    for j = 1:size(c,2)
        for k = 1:size(c,1)
            for l = 1:size(c,2)
                px = i+k-1; py = j+l-1;
                result = result+c(i,j)*c(k,l)*(x(2)^px-x(1)^px)/px* ...
                    (y(2)^py-y(1)^py)/py;
            end
        end
    end
end
end
