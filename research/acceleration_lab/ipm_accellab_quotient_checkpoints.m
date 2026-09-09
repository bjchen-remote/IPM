function experiment = ipm_accellab_quotient_checkpoints(files,user)
%IPM_ACCELLAB_QUOTIENT_CHECKPOINTS Fixed memory-four, one-damping real probe.
%   Four/five native same-grid frames; one restore; unchanged full RHS for
%   every changed rho. Run only in a parent-approved large-grid memory window.
if nargin < 2, user = struct(); end
opts = struct('outputRoot','','expectedRemeshCount',3);
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:QuotientCheckpointOptions','Unknown option.');
for k = 1:numel(names), opts.(names{k}) = user.(names{k}); end
assert(iscell(files) && any(numel(files) == [4,5]), ...
    'ipm:QuotientCheckpointFiles','Supply four or five chronological native checkpoints.');
validateattributes(opts.expectedRemeshCount,{'numeric'},{'scalar','integer','nonnegative','finite'});
assert(~isempty(opts.outputRoot),'ipm:QuotientCheckpointOutput','An independent output root is required.');
if ~isfolder(opts.outputRoot), mkdir(opts.outputRoot); end
[~,token] = fileparts(tempname);
destination = fullfile(opts.outputRoot,['quotient_checkpoints_',token]);
mkdir(destination);
experiment = struct('kind','fixed_protocol_full_field_inner_shape_secant', ...
    'status','registered','checkpointFiles',{files},'options',opts, ...
    'outputDirectory',destination,'requestedMemory',4,'effectiveMemory',numel(files)-1, ...
    'damping',0.0625,'parameterScan',false,'restoreCalls',0, ...
    'operatorBuildCount','not_instrumented_restore_may_rebuild_native_remeshed_axes', ...
    'pdeAcceptedSteps',0,'signatureValidated',false(size(files)), ...
    'configDifferences',{{}},'metadataDifferences',{{}}, ...
    'exceptionIdentifier','','exceptionMessage','','wallSeconds',0);
save(fullfile(destination,'registration.mat'),'experiment','-v7.3');
started = tic;
try
    saved = cell(size(files));
    for k = 1:numel(files)
        checkpoint = ipm.output.readCheckpoint(files{k});
        saved{k} = checkpoint.payload.state;
        experiment.signatureValidated(k) = true;
        fprintf('QUOTIENT_NATIVE validated %d/%d step %d tau %.12g\n', ...
            k,numel(files),saved{k}.step,saved{k}.scale.canonicalTime);
    end
    reference = saved{end};
    metadataControls = {'latestCheckpointFile','resumeCount','resumedFromStep', ...
        'resumedFromCanonicalTime','resumedFromPhysicalTime','resumedFromCheckpoint', ...
        'resultFile','videoFile','caseMetadata'};
    timeControls = {'finalTime','physicalFinalTime','maxSteps'};
    for k = 1:numel(saved)
        current = saved{k};
        for field = {'x','y','baseX','baseY','rescaling','remeshCount','mass0','rhoRange0'}
            assert(isequaln(current.(field{1}),reference.(field{1})), ...
                'ipm:QuotientCheckpointGrid','Native field %s differs in frame %d.',field{1},k);
        end
        assert(current.remeshCount == opts.expectedRemeshCount, ...
            'ipm:QuotientCheckpointRemesh','The native remesh count differs from registration.');
        configDifference = differences(current.config,reference.config,'config');
        experiment.configDifferences{k} = configDifference;
        numericalCurrent = rmfield(current.config,'output');
        numericalReference = rmfield(reference.config,'output');
        numericalCurrent.time = rmfield(numericalCurrent.time,timeControls);
        numericalReference.time = rmfield(numericalReference.time,timeControls);
        assert(isequaln(numericalCurrent,numericalReference), ...
            'ipm:QuotientCheckpointNumerics', ...
            'A numerical configuration changed; only output and terminal horizons may differ.');
        experiment.metadataDifferences{k} = differences(current.runMetadata, ...
            reference.runMetadata,'runMetadata');
        assert(isequaln(drop_existing(current.runMetadata,metadataControls), ...
            drop_existing(reference.runMetadata,metadataControls)), ...
            'ipm:QuotientCheckpointLineage','Unregistered run lineage metadata differs.');
        if k > 1
            assert(current.scale.canonicalTime > saved{k-1}.scale.canonicalTime && ...
                current.step > saved{k-1}.step, ...
                'ipm:QuotientCheckpointOrder','Native frame times and steps must strictly increase.');
        end
    end
    experiment.allowedNonNumericalDifferences = struct( ...
        'configOutputDomain',true,'configTimeFields',{timeControls}, ...
        'metadataFields',{metadataControls});
    experiment.sourceSteps = cellfun(@(s)s.step,saved);
    experiment.sourceCanonicalTimes = cellfun(@(s)s.scale.canonicalTime,saved);
    fprintf('QUOTIENT_NATIVE same-grid preflight passed; restoring final frame once\n');
    experiment.restoreCalls = 1;
    base = ipm.output.restoreCheckpoint(checkpoint,struct('saveResults',false, ...
        'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false));
    states = cell(size(saved));
    for k = 1:numel(saved)
        states{k} = base;
        states{k}.rho = saved{k}.rho; states{k}.scale = saved{k}.scale;
        states{k}.config = saved{k}.config; states{k}.step = saved{k}.step;
        states{k}.normalizedTime = saved{k}.normalizedTime;
        states{k}.flow = []; states{k}.rhsCache = [];
    end
    clear checkpoint saved;
    experiment.sourceShape = ipm_accellab_shape_history(states);
    experiment.rawRelativePeakDrift = experiment.sourceShape.quadraticPeak/ ...
        experiment.sourceShape.quadraticPeak(end)-1;
    experiment.maximumRawRelativePeakDrift = max(abs(experiment.rawRelativePeakDrift));
    experiment.historyProjectionDeclared = ...
        'Every history and candidate gets positive amplitude normalization to terminal quadratic P before its fresh RHS; no PPrime replacement.';
    protocol = struct('secant',struct('memory',4,'damping',0.0625));
    [rhoCandidate,trial] = ipm_accellab_quotient_secant(states,protocol);
    experiment.trial = trial;
    experiment.status = 'completed_independent_profile_candidate_probe';
    experiment.originalRhsEvaluationsAfterRestore = trial.originalRhsEvaluations;
    artifact = struct('kind','research_full_field_not_native_checkpoint', ...
        'sourceCheckpoint',files{end},'x',base.ops.x,'y',base.ops.y, ...
        'config',base.config,'runtimeRescaling',base.ops.rescaling, ...
        'sourceScale',base.scale,'sourceRho',base.rho,'rhoCandidate',rhoCandidate, ...
        'trial',trial);
    save(fullfile(destination,'candidate.mat'),'artifact','-v7.3');
    fprintf('QUOTIENT_NATIVE status=%s ratio=%.12g originalRHS=%d\n', ...
        trial.status,trial.residualRatio,trial.originalRhsEvaluations);
catch exception
    experiment.status = 'rejected_or_failed';
    experiment.exceptionIdentifier = exception.identifier;
    experiment.exceptionMessage = exception.message;
    experiment.wallSeconds = toc(started);
    persist(experiment);
    rethrow(exception);
end
experiment.wallSeconds = toc(started);
persist(experiment);
end

function output = drop_existing(input,fields)
output = rmfield(input,intersect(fieldnames(input),fields));
end

function output = differences(a,b,prefix)
output = struct('field',{},'sourceExists',{},'referenceExists',{},'source',{},'reference',{});
if isequaln(a,b), return; end
if isstruct(a) && isscalar(a) && isstruct(b) && isscalar(b)
    names = union(fieldnames(a),fieldnames(b));
    for k = 1:numel(names)
        name = names{k}; inA = isfield(a,name); inB = isfield(b,name);
        if inA && inB
            changed = differences(a.(name),b.(name),[prefix,'.',name]);
        else
            av = []; bv = [];
            if inA, av = a.(name); end
            if inB, bv = b.(name); end
            changed = struct('field',[prefix,'.',name], ...
                'sourceExists',inA,'referenceExists',inB,'source',av,'reference',bv);
        end
        output = [output,changed]; %#ok<AGROW>
    end
else
    output = struct('field',prefix,'sourceExists',true,'referenceExists',true, ...
        'source',a,'reference',b);
end
end

function persist(experiment)
destination = experiment.outputDirectory;
assert(~isfile(fullfile(destination,'experiment.mat')),'ipm:QuotientOverwrite','Output already exists.');
save(fullfile(destination,'experiment.mat'),'experiment','-v7.3');
fid = fopen(fullfile(destination,'experiment.json'),'w');
assert(fid >= 0,'ipm:QuotientOutput','Could not open JSON output.');
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(experiment,'PrettyPrint',true));
clear cleanup;
end
