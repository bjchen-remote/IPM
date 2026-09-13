function report = test_original_zero_inner_forecast(dataFile,outputDirectory)
%TEST_ORIGINAL_ZERO_INNER_FORECAST Fixed 3-train/1-future shape-only holdout.
% No PDE, LU, field reconstruction, candidate selection or solver feedback.
assert(~isfolder(outputDirectory),'Refusing to overwrite forecast evidence.');
data = load(dataFile,'profiles','report');
assert(numel(data.profiles)==4 && numel(data.report.records)==4 && ...
    data.report.sameCaseId && all([data.report.records.fullHistoryTrusted]), ...
    'Four trusted paired source profiles are required.');
times = [data.report.records.tau];
assert(all(diff(times)>0),'The four frames must be chronological.');

wallAxis = linspace(-1,1,401);
verticalAxis = linspace(0,2,401);
wall = zeros(4,numel(wallAxis));
vertical = zeros(4,numel(verticalAxis));
for k = 1:4
    p = data.profiles{k};
    wall(k,:) = interp1(p.innerX,p.wallRX/p.peak,wallAxis,'pchip');
    vertical(k,:) = interp1(p.innerY,p.verticalRX/p.peak,verticalAxis,'pchip');
end
assert(all(isfinite(wall),'all') && all(isfinite(vertical),'all'), ...
    'Every fixed inner observation must be covered.');

report = struct('kind','original_zero_inner_shape_holdout_v1', ...
    'sourceDataFile',dataFile,'trainTimes',times(1:3), ...
    'heldoutTime',times(4),'forecastHorizon',times(4)-times(3), ...
    'sameCaseId',true,'physicalFieldPredicted',false, ...
    'originalTrajectoryModified',false, ...
    'wall',one_shape(wall,times), ...
    'vertical',one_shape(vertical,times));
mkdir(outputDirectory);
fid = fopen(fullfile(outputDirectory,'report.json'),'w');
assert(fid>=0,'Could not create forecast report.');
cleaner = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
save(fullfile(outputDirectory,'report.mat'),'report','-v7.3');
fprintf('INNER_SHAPE_HOLDOUT %s\n',jsonencode(report));
end

function result = one_shape(values,times)
first = values(1,:); second = values(2,:);
source = values(3,:); heldout = values(4,:);
dt = times(3)-times(2);
h = times(4)-times(3);
ratio = h/dt;
constant = source;
linear = source+ratio*(source-second);
increment1 = second-first;
increment2 = source-second;
lambda = dot(increment1,increment2)/dot(increment1,increment1);
cosine = dot(increment1,increment2)/ ...
    (norm(increment1)*norm(increment2));
rankOneValid = isfinite(lambda) && lambda>0;
rankOneError = NaN;
rankOneCoefficient = NaN;
if rankOneValid
    if abs(lambda-1)<1e-8
        rankOneCoefficient=ratio;
    else
        rankOneCoefficient=lambda*(lambda^ratio-1)/(lambda-1);
    end
    rankOne = source+rankOneCoefficient*increment2;
    rankOneError = relative_l2(rankOne,heldout);
end
result = struct('carryForwardRelativeL2',relative_l2(constant,heldout), ...
    'linearRelativeL2',relative_l2(linear,heldout), ...
    'rankOneRelativeL2',rankOneError, ...
    'rankOneValid',rankOneValid,'rankOneLambda',lambda, ...
    'rankOneForecastCoefficient',rankOneCoefficient, ...
    'trainingIncrementCosine',cosine, ...
    'fixedAxisNodeCount',size(values,2));
end

function value = relative_l2(candidate,reference)
value = sqrt(sum((candidate-reference).^2)/sum(reference.^2));
end
