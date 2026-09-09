function experiment = ipm_accellab_probe_checkpoints(files,user)
%IPM_ACCELLAB_PROBE_CHECKPOINTS Bounded real-data probe with one operator build.
%   Every checkpoint is validated before a shared-ops research history is
%   formed. Positive history projection is explicit and recorded; projected
%   histories and candidates are not accepted PDE trajectories.
if nargin < 2, user = struct(); end
opts = struct('outputRoot','','allowHistoryProjection',false, ...
    'memories',[1,2,4],'profileOptions',struct());
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:AccelCheckpointOptions', ...
    'Unknown checkpoint-probe option.');
for index = 1:numel(names), opts.(names{index}) = user.(names{index}); end
assert(iscell(files) && numel(files) >= 2 && numel(files) <= 6, ...
    'ipm:AccelCheckpointFiles','Supply two to six chronological checkpoint files.');
assert(islogical(opts.allowHistoryProjection) && isscalar(opts.allowHistoryProjection), ...
    'ipm:AccelCheckpointOptions','allowHistoryProjection must be logical.');
validateattributes(opts.memories,{'numeric'},{'vector','integer','>=',1,'<=',5});
started = tic;
checkpoints = cell(size(files));
for index = 1:numel(files)
    checkpoints{index} = ipm.output.readCheckpoint(files{index});
    fprintf('Q512_ACCEL validated checkpoint %d/%d\n',index,numel(files));
end
lastSaved = checkpoints{end}.payload.state;
lastMetadata = rmfield(lastSaved.runMetadata,'latestCheckpointFile');
for index = 1:numel(files)
    saved = checkpoints{index}.payload.state;
    metadata = rmfield(saved.runMetadata,'latestCheckpointFile');
    assert(isequal(saved.x,lastSaved.x) && isequal(saved.y,lastSaved.y) && ...
        isequaln(saved.rescaling,lastSaved.rescaling) && ...
        isequaln(metadata,lastMetadata), ...
        'ipm:AccelCheckpointLineage', ...
        'The input checkpoints changed grid, gauge references, or run lineage.');
end
fprintf('Q512_ACCEL restoring final checkpoint with one operator build\n');
base = ipm.output.restoreCheckpoint(checkpoints{end},struct( ...
    'finalTime',max(lastSaved.config.time.finalTime,lastSaved.scale.canonicalTime+1), ...
    'physicalFinalTime',max(lastSaved.config.time.physicalFinalTime, ...
        lastSaved.scale.physicalTime+1), ...
    'maxSteps',max(lastSaved.config.time.maxSteps,lastSaved.step+1), ...
    'saveResults',false,'makePlots',false,'livePlot',false, ...
    'writeVideo',false,'verbose',false));
states = cell(size(files));
for index = 1:numel(files)
    saved = checkpoints{index}.payload.state;
    states{index} = base;
    states{index}.rho = saved.rho;
    states{index}.scale = saved.scale;
    states{index}.config = saved.config;
    states{index}.step = saved.step;
    states{index}.normalizedTime = saved.normalizedTime;
    states{index}.flow = [];
    states{index}.rhsCache = [];
end
clear checkpoints;
shape = ipm_accellab_shape_history(states);
relativeDrift = shape.quadraticPeak/shape.quadraticPeak(end)-1;
experiment = struct('kind','q512_independent_stationary_profile_probe', ...
    'status','registered','options',opts,'checkpointFiles',{files}, ...
    'sourceSteps',cellfun(@(s)s.step,states), ...
    'sourceCanonicalTimes',shape.canonicalTime,'operatorBuildCount',1, ...
    'pdeAcceptedSteps',0,'sourceShapeDiagnostics',shape, ...
    'historyRelativePeakDrift',relativeDrift, ...
    'maximumHistoryRelativePeakDrift',max(abs(relativeDrift)), ...
    'lineageNonNumericalFieldExcluded','runMetadata.latestCheckpointFile', ...
    'rhsImplementation',which('ipm.evolve.rhs'), ...
    'historyProjectionApplied',false,'historyAmplitudeFactors',ones(size(files)), ...
    'rawHistoryRejectionIdentifier','','rawHistoryRejectionMessage','', ...
    'trials',{{}},'outputDirectory','','wallSeconds',0);
if ~isempty(opts.outputRoot)
    if ~isfolder(opts.outputRoot), mkdir(opts.outputRoot); end
    [~,unique] = fileparts(tempname);
    experiment.outputDirectory = fullfile(opts.outputRoot,['probe_',unique]);
    mkdir(experiment.outputDirectory);
    save(fullfile(experiment.outputDirectory,'registration.mat'),'experiment','-v7.3');
end
probeStates = states;
try
    [~,raw] = ipm_accellab_profile_probe(states,struct('secant',struct('memory',1)));
    experiment.rawHistoryProbe = raw;
catch exception
    if ~strcmp(exception.identifier,'ipm:AccelProfileGaugeDrift')
        rethrow(exception);
    end
    experiment.rawHistoryRejectionIdentifier = exception.identifier;
    experiment.rawHistoryRejectionMessage = exception.message;
    if ~opts.allowHistoryProjection
        experiment.status = 'history_gauge_drift_rejected';
        experiment.wallSeconds = toc(started);
        save_experiment(experiment);
        return;
    end
    experiment.historyProjectionApplied = true;
    experiment.historyAmplitudeFactors = shape.quadraticPeak(end)./shape.quadraticPeak;
    for index = 1:numel(states)
        probeStates{index}.rho = ...
            experiment.historyAmplitudeFactors(index)*states{index}.rho;
    end
end
for memory = opts.memories
    profileOptions = opts.profileOptions;
    if ~isfield(profileOptions,'secant'), profileOptions.secant = struct(); end
    profileOptions.secant.memory = memory;
    [rhoCandidate,trial] = ipm_accellab_profile_probe(probeStates,profileOptions);
    experiment.trials{end+1} = trial;
    if ~isempty(experiment.outputDirectory)
        file = fullfile(experiment.outputDirectory,sprintf('trial_memory_%d.mat',memory));
        assert(~isfile(file),'ipm:AccelProbeOverwrite','Trial output already exists.');
        save(file,'rhoCandidate','trial','-v7.3');
    end
    fprintf('Q512_ACCEL memory=%d status=%s ratio=%.12g evaluations=%d\n', ...
        memory,trial.status,trial.residualRatio,trial.evaluations);
end
experiment.status = 'completed_independent_profile_probes';
experiment.wallSeconds = toc(started);
save_experiment(experiment);
end

function save_experiment(experiment)
if isempty(experiment.outputDirectory), return; end
file = fullfile(experiment.outputDirectory,'experiment.mat');
assert(~isfile(file),'ipm:AccelProbeOverwrite','Experiment output already exists.');
save(file,'experiment','-v7.3');
jsonFile = fullfile(experiment.outputDirectory,'experiment.json');
assert(~isfile(jsonFile),'ipm:AccelProbeOverwrite','JSON output already exists.');
fid = fopen(jsonFile,'w');
assert(fid >= 0,'ipm:AccelProbeWrite','Cannot write JSON output.');
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(experiment,'PrettyPrint',true));
clear cleanup;
end
