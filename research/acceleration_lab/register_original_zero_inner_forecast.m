function registration = register_original_zero_inner_forecast(dataFile,outputDirectory)
%REGISTER_ORIGINAL_ZERO_INNER_FORECAST Freeze shape predictors before new frames.
% Reads only the first three training profiles. Future field/clock is absent.
assert(~isfolder(outputDirectory),'Refusing to overwrite forecast registration.');
data = load(dataFile,'profiles','report');
assert(numel(data.profiles)>=3 && numel(data.report.records)>=3 && ...
    data.report.sameCaseId && all([data.report.records(1:3).fullHistoryTrusted]));
training = data.report.records(1:3);
times = [training.tau];
assert(all(diff(times)>0));
windows = [7.15,7.25;7.45,7.55];
assert(all(windows(:,1)>times(3)));

axes = {linspace(-1,1,401),linspace(0,2,401)};
shape = cell(1,2);
for kind=1:2
    values=zeros(3,numel(axes{kind}));
    for k=1:3
        p=data.profiles{k};
        if kind==1
            values(k,:)=interp1(p.innerX,p.wallRX/p.peak,axes{kind},'pchip');
        else
            values(k,:)=interp1(p.innerY,p.verticalRX/p.peak,axes{kind},'pchip');
        end
    end
    assert(all(isfinite(values),'all'));
    d1=values(2,:)-values(1,:);
    d2=values(3,:)-values(2,:);
    lambda=dot(d1,d2)/dot(d1,d1);
    shape{kind}=struct('axis',axes{kind},'lastShape',values(3,:), ...
        'lastIncrement',d2,'lambda',lambda, ...
        'realRankOneAvailable',isfinite(lambda)&&lambda>0, ...
        'decayingRankOne',isfinite(lambda)&&lambda>0&&lambda<1);
end
registration=struct('kind','original_zero_frozen_inner_shape_forecast_v1', ...
    'sourceProfileData',dataFile, ...
    'trainingCheckpoints',{ {training.checkpointFile} }, ...
    'trainingTimes',times, ...
    'futureCanonicalWindows',windows, ...
    'futureSelection','first native accepted checkpoint inside each window', ...
    'models',{ {'carry_forward','linear','rank_one_positive_lambda'} }, ...
    'physicalFieldPredicted',false,'originalTrajectoryModified',false, ...
    'wallLambda',shape{1}.lambda,'verticalLambda',shape{2}.lambda, ...
    'createdUtc',char(datetime('now','TimeZone','UTC', ...
        'Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''')));
mkdir(outputDirectory);
save(fullfile(outputDirectory,'registration.mat'),'registration','shape','-v7.3');
fid=fopen(fullfile(outputDirectory,'registration.json'),'w');
assert(fid>=0);
cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(registration,PrettyPrint=true));
fprintf('FROZEN_INNER_FORECAST %s\n',jsonencode(registration));
end
