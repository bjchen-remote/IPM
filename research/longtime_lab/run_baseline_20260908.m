% Registered, immutable-source continuation of the schema-4 q512 baseline.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
campaignRoot = fullfile(projectRoot,'result','longtime','20260908_campaign_v1');
sourceRoot = fullfile(campaignRoot,'baseline_source');
addpath(sourceRoot);
addpath(fullfile(sourceRoot,'research','gauge_lab'));
addpath(fullfile(sourceRoot,'research','experiments'));
% Historical numeric signatures include threaded reductions. Use the
% original 10-thread environment; changing it fails exact validation.
assert(maxNumCompThreads == 10,'Historical checkpoint requires its original thread setting.');
parentDirectory = fullfile(projectRoot,'result','verification', ...
    'q512_fresh_exact_gauge_roots_v4', ...
    'root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f');
parentCheckpoint = fullfile(parentDirectory, ...
    'checkpoint_20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7_step0000001818.mat');
assert(isfile(parentCheckpoint));
assert(startsWith(which('ipm.solve'),sourceRoot));
% Each endpoint is a review boundary. Exiting this script is not completion
% of the scientific objective: a mesh/safety exit requires a new experiment.
for segment = 1:100
    cp = ipm.output.readCheckpoint(parentCheckpoint);
    targetTau = cp.payload.state.scale.canonicalTime + 0.2;
    segmentDirectory = fullfile(campaignRoot,sprintf('baseline_segment_%03d',segment));
    assert(~isfolder(segmentDirectory),'Refusing to overwrite a baseline segment.');
    mkdir(segmentDirectory);
    fprintf('LONGTIME_BEGIN targetTau=%.9g source=%s\n',targetTau,which('ipm.solve'));
    opts = struct( ...
        'finalTime',targetTau,'physicalFinalTime',10,'maxSteps',100000, ...
        'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false, ...
        'saveResults',true,'resultFile',fullfile(segmentDirectory,'result.mat'), ...
        'checkpoint',struct('file',fullfile(segmentDirectory,'checkpoint.mat'), ...
        'every',0.1,'atExit',true));
    parentStep = cp.payload.state.step;
    parentHistory = cp.payload.log.history;
    wallClock = tic;
    result = ipm.solve(opts,cp);
    wallSeconds = toc(wallClock);
    cp = ipm.output.readCheckpoint(result.metadata.latestCheckpointFile);
    h = cp.payload.log.history;
    c = h.common; m = h.mesh;
    audit = struct('targetTau',targetTau,'step',cp.payload.state.step, ...
        'tau',c.canonicalTau(end),'physicalTime',c.physicalTime(end), ...
        'physicalRhoXInf',c.physicalRhoXInf(end), ...
        'c_l',c.c_l(end),'c_omega',c.c_omega(end), ...
        'horizontalCoreCells',m.coreGridPoints(end), ...
        'verticalCoreCells',m.verticalCoreGridPoints(end), ...
        'safetyFactor',m.safetyFactor(end), ...
        'stopReason',result.state.stopReason, ...
        'checkpointFile',result.metadata.latestCheckpointFile, ...
        'resultFile',result.metadata.resultFile, ...
        'parentCheckpointFile',parentCheckpoint, ...
        'wallSeconds',wallSeconds, ...
        'segmentSteps',result.state.steps-parentStep);
    groups = fieldnames(parentHistory);
    for gi = 1:numel(groups)
        names = fieldnames(parentHistory.(groups{gi}));
        for ni = 1:numel(names)
            old = parentHistory.(groups{gi}).(names{ni});
            new = h.(groups{gi}).(names{ni});
            assert(isequaln(old,new(1:numel(old))),'Parent history prefix changed.');
        end
    end
    assert(isequal(cp.payload.state.rho,result.state.rho));
    fid = fopen(fullfile(campaignRoot,'baseline_progress.jsonl'),'a');
    assert(fid>=0); fprintf(fid,'%s\n',jsonencode(audit)); fclose(fid);
    fprintf('LONGTIME_ENDPOINT %s\n',jsonencode(audit));
    trusted = ipm.output.trustedMask(h,cp.payload.state.config);
    assert(all(trusted),'Continuous trusted prefix failed.');
    if audit.tau < targetTau-1e-10 || audit.horizontalCoreCells < 17 || ...
            audit.verticalCoreCells < 17 || audit.safetyFactor >= 0.60
        fprintf('LONGTIME_REVIEW_REQUIRED mesh/safety/early endpoint; preserve latest checkpoint.\n');
        break
    end
    parentCheckpoint = result.metadata.latestCheckpointFile;
    clear cp h c m result parentHistory
end
