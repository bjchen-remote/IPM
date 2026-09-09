function summary=collect_frozen_profile_forecasts(campaignDirectories,forecastDirectory)
%COLLECT_FROZEN_PROFILE_FORECASTS Observe earliest qualified future endpoints.
% No rank/rate refit, restore, LU or PDE step. Existing observations are immutable.
assert(iscell(campaignDirectories) && maxNumCompThreads==10);
forecastFile=fullfile(forecastDirectory,'fixed_forecast.mat');
loaded=load(forecastFile,'forecast');registration=loaded.forecast.registration;
items={};
for k=1:numel(campaignDirectories)
    files=dir(fullfile(campaignDirectories{k},'stage_*','endpoint_audit.json'));
    for j=1:numel(files)
        path=fullfile(files(j).folder,files(j).name);a=jsondecode(fileread(path));
        if a.stageAccepted && a.parentEquivalentTau>=registration.futureEquivalentTauInterval(1) && ...
                a.parentEquivalentTau<=registration.futureEquivalentTauInterval(2)
            items{end+1}=struct('auditFile',path,'audit',a); %#ok<AGROW>
        end
    end
end
if isempty(items)
    summary=struct('observationCount',0,'predictorRefitted',false,'pdeAdvanced',false);
    return
end
times=cellfun(@(v)v.audit.parentEquivalentTau,items);[~,order]=sort(times);items=items(order);
selected={};last=-Inf;
for k=1:numel(items)
    a=items{k}.audit;
    if a.parentEquivalentTau-last<registration.minimumSeparationBetweenAcceptedFrames,continue;end
    selected{end+1}=items{k};last=a.parentEquivalentTau; %#ok<AGROW>
    if numel(selected)==registration.maximumFutureFrames,break;end
end
out=fullfile(forecastDirectory,'future_observations');if ~isfolder(out),mkdir(out);end
observations=cell(size(selected));
for k=1:numel(selected)
    a=selected{k}.audit;file=fullfile(out,sprintf('observation_%02d.json',k));
    assert(a.nativeAndResultEntireHistoryTrusted && a.nativeResultStateClocksAndHistoryExactlyPaired && ...
        a.freshLineagePreserved && a.nativeResultScalesPhysicalSamplesSnapshotsAndConfigExactlyPaired);
    if isfile(file)
        r=jsondecode(fileread(file));
        assert(strcmp(r.checkpointFile,a.checkpointFile) && r.parentEquivalentTau==a.parentEquivalentTau && ...
            ~r.predictorRefitted && ~r.pdeAdvanced,'Existing prospective observation identity changed.');
    else
        r=evaluate_frozen_profile_forecast(a.checkpointFile,forecastFile,file);
        assert(r.parentEquivalentTau==a.parentEquivalentTau && ...
            r.absolutePhysicalTime==a.absolutePhysicalTime && r.step==a.step);
    end
    observations{k}=struct('index',k,'stageAuditFile',selected{k}.auditFile,'observationFile',file, ...
        'checkpointFile',a.checkpointFile,'parentEquivalentTau',r.parentEquivalentTau, ...
        'absolutePhysicalTime',r.absolutePhysicalTime,'coreCells',a.coreCells, ...
        'physicalGradientMaximum',a.physicalGradientMaximum,'errors',r.errors, ...
        'distanceToFixedCandidateLimit',r.distanceToFixedCandidateLimit, ...
        'linearToRank1ErrorRatio',r.errors.fixedLinearTime.relativeL2/r.errors.fixedIncrementRank1.relativeL2);
    fprintf('PROSPECTIVE_OBSERVATION %s\n',jsonencode(observations{k}));
end
summary=struct('registration',registration,'observationCount',numel(observations), ...
    'selectionRule','Earliest qualified actual endpoints in time order; fixed separation and count cap.', ...
    'observations',[observations{:}],'predictorRefitted',false,'pdeAdvanced',false, ...
    'asymptoticLimitEstablished',false);
file=fullfile(out,sprintf('summary_count_%02d.json',numel(observations)));
if isfile(file)
    old=jsondecode(fileread(file));assert(isequaln(old,jsondecode(jsonencode(summary))));
else
    fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));
end
end
