function report = ipm_perflab_benchmark_weno5(options)
%IPM_PERFLAB_BENCHMARK_WENO5 Compare blocked faces with the maintained kernel.
%   No solver runs, plots, or files are produced. OPTIONS may specify
%   shape=[65,129], repetitions=3, and blockSizes=[1,4,16,64]. Larger shapes
%   are explicit opt-ins; use [513,1025] only after checking available load.

if nargin < 1
    options = struct();
end
options = settings(options);
classes = {'double','single'};
boundaries = {'extrapolate','reflect'};
cases = struct('precision',{},'boundary',{},'family',{},'blockSize',{}, ...
    'exact',{},'maximumAbsoluteDifference',{});
for precision = 1:numel(classes)
    for boundary = 1:numel(boundaries)
        for family = 1:5
            [q,a,J] = data([7,19],classes{precision},family);
            kernelOptions = struct('lowerBoundary',boundaries{boundary}, ...
                'upperBoundary',boundaries{boundary}, ...
                'extrapolationDegree',4+mod(family,2));
            expected = ipm.field.weno5FluxDerivative(q,a,J,1/18,kernelOptions);
            for blockSize = options.blockSizes
                actual = ipm_perflab_weno5_candidate( ...
                    q,a,J,1/18,kernelOptions,blockSize);
                item = struct('precision',classes{precision}, ...
                    'boundary',boundaries{boundary},'family',family, ...
                    'blockSize',blockSize,'exact',isequaln(actual,expected), ...
                    'maximumAbsoluteDifference', ...
                    max(abs(double(actual)-double(expected)),[],'all'));
                cases(end+1) = item; %#ok<AGROW>
            end
        end
    end
end
[q,a,J] = data(options.shape,'double',1);
h = 1/(options.shape(2)-1);
kernelOptions = struct();
baseline = @()ipm.field.weno5FluxDerivative(q,a,J,h,kernelOptions);
baseline();
baselineSeconds = measure(baseline,options.repetitions);
expected = baseline();
timing = struct('blockSize',{},'medianSeconds',{},'speedup',{},'exact',{}, ...
    'samplesSeconds',{});
for blockSize = options.blockSizes
    candidate = @()ipm_perflab_weno5_candidate( ...
        q,a,J,h,kernelOptions,blockSize);
    actual = candidate();
    seconds = measure(candidate,options.repetitions);
    timing(end+1) = struct('blockSize',blockSize, ...
        'medianSeconds',median(seconds), ...
        'speedup',median(baselineSeconds)/median(seconds), ...
        'exact',isequaln(actual,expected),'samplesSeconds',seconds); %#ok<AGROW>
end
report = struct('schemaVersion',1,'kind','ipm_weno5_block_benchmark', ...
    'matlabVersion',version,'computer',computer,'options',options, ...
    'computationalThreads',maxNumCompThreads, ...
    'caseCount',numel(cases),'allCasesExact',all([cases.exact]), ...
    'allTimedShapesExact',all([timing.exact]),'cases',cases, ...
    'baselineSamplesSeconds',baselineSeconds, ...
    'baselineMedianSeconds',median(baselineSeconds),'timing',timing, ...
    'productionSource',which('ipm.field.weno5FluxDerivative'), ...
    'candidateSource',which('ipm_perflab_weno5_candidate'), ...
    'interpretation','kernel timing only; no complete solver speedup claimed');
fprintf('WENO5: %d comparison cases, exact=%d, shape=%dx%d\n', ...
    report.caseCount,report.allCasesExact,options.shape);
for index = 1:numel(timing)
    item = timing(index);
    fprintf('  block=%d median=%.6g s speedup=%.3fx exact=%d\n', ...
        item.blockSize,item.medianSeconds,item.speedup,item.exact);
end
end

function [q,a,J] = data(shape,precision,family)
x = linspace(0,1,shape(2));
y = linspace(0,1,shape(1))';
q = exp(-6*(x-0.4).^2).*(1+0.2*cos(2*pi*y))+0.02*sin(7*x+3*y);
a = 0.3+sin(2*pi*x).*cos(3*y);
J = 1+0.4*x+0.1*y;
switch family
    case 2
        q(:) = 1;
        a = 0.25;
        J = 1;
    case 3
        q(:) = 0;
        a = sin(2*pi*x);
        J = 1+0.3*x;
    case 4
        q = q.*(1e-90+1e60*y);
        a = 0.1+y;
        J = 1+y;
        if strcmp(precision,'single')
            q = 1e-10+q/1e40;
        end
    case 5
        q = double(x > 0.37).*(1+y);
        a = zeros(size(q));
end
q = cast(q,precision);
a = cast(a,precision);
J = cast(J,precision);
end

function seconds = measure(operation,repetitions)
seconds = zeros(1,repetitions);
for repetition = 1:repetitions
    timer = tic;
    operation();
    seconds(repetition) = toc(timer);
end
end

function options = settings(options)
defaults = struct('shape',[65,129],'repetitions',3,'blockSizes',[1,4,16,64]);
assert(isstruct(options) && isscalar(options),'Pass scalar options.');
names = fieldnames(options);
assert(all(ismember(names,fieldnames(defaults))),'Unknown benchmark option.');
for index = 1:numel(names)
    defaults.(names{index}) = options.(names{index});
end
options = defaults;
validateattributes(options.shape,{'numeric'}, ...
    {'vector','numel',2,'integer','positive','finite'});
assert(options.shape(2) >= 6,'At least six nodes are required.');
validateattributes(options.repetitions,{'numeric'}, ...
    {'scalar','integer','>=',1,'<=',20});
validateattributes(options.blockSizes,{'numeric'}, ...
    {'vector','integer','positive','finite'});
options.blockSizes = options.blockSizes(:)';
end
