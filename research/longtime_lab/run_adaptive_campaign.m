function run_adaptive_campaign(projectRoot,bootstrapFile,outputDirectory,maxStages)
%RUN_ADAPTIVE_CAMPAIGN Continue tested schema-4 profiles with audited regrids.
% A batch boundary or admissibility failure does not complete the goal.
assert(maxNumCompThreads==10,'Keep the validated checkpoint thread environment.');
assert(~isfolder(outputDirectory),'Use a new campaign directory.');
validateattributes(maxStages,{'numeric'},{'scalar','integer','positive','<=',200});
variables=whos('-file',bootstrapFile); names={variables.name};
seedResultFile='';
if ismember('checkpoint',names)
    seedCheckpoint=ipm.output.readCheckpoint(bootstrapFile);
    currentFile=char(bootstrapFile); bootstrapKind='native_checkpoint';
elseif ismember('result',names)
    seed=ipm.output.validate(bootstrapFile);
    seedResultFile=char(bootstrapFile); bootstrapKind='native_checkpoint_with_paired_result';
    currentFile=char(seed.metadata.latestCheckpointFile);
    seedCheckpoint=ipm.output.readCheckpoint(currentFile);
else
    loaded=load(bootstrapFile,'report'); boot=loaded.report;
    assert(boot.passed,'A matched-physical-time mesh smoke must pass first.');
    seedResultFile=char(boot.newResultFile); bootstrapKind='passed_mesh_smoke';
    seed=ipm.output.validate(seedResultFile);
    currentFile=char(seed.metadata.latestCheckpointFile);
    seedCheckpoint=ipm.output.readCheckpoint(currentFile);
end
assert(all(ipm.output.trustedMask(seedCheckpoint.payload.log.history, ...
    seedCheckpoint.payload.state.config)));
if ~isempty(seedResultFile)
    assert(all(ipm.output.trustedMask(seed)));
    state=seedCheckpoint.payload.state;
    assert(isequal(state.rho,seed.state.rho) && ...
        isequal(state.x(:),seed.grid.x(:)) && isequal(state.y(:),seed.grid.y(:)) && ...
        state.step==seed.state.steps && ...
        state.scale.canonicalTime==seed.state.canonicalTime && ...
        state.scale.physicalTime==seed.state.physicalTime, ...
        'Bootstrap terminal result and native checkpoint must be paired exactly.');
end
clear loaded seed seedCheckpoint state
mkdir(outputDirectory);
manifest=fullfile(outputDirectory,'campaign_manifest.jsonl');
rootResult=fullfile(projectRoot,'result','verification', ...
    'q512_fresh_exact_gauge_roots_v4', ...
    'root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f', ...
    'result.mat');
families=[24,1;21,1;21,0.75;21,0.5;21,0.35;21,0.25;21,0.18;21,1.5;21,2;21,3];
platform=struct('positiveFineCells',96,'anchorFineCellFraction',.85, ...
    'resolutionPadding',1.10,'verticalResolutionPadding',1.15);
% Independent stage006 frozen-field screens admitted both reserve budgets.
% Running campaigns use frozen copies; this expanded registration applies
% only to a subsequently frozen campaign and is rescored on its actual state.
platform(2)=platform(1);platform(2).positiveFineCells=112;
platform(2).anchorFineCellFraction=.90;
platform(3)=platform(1);platform(3).positiveFineCells=128;
platform(3).anchorFineCellFraction=.90;
emit(manifest,struct('event','registered','bootstrapFile',bootstrapFile, ...
    'bootstrapKind',bootstrapKind,'meshFamiliesTargetSigma',families, ...
    'preferredPlatformControls',platform(1),'registeredPlatformControls',platform, ...
    'preferredPlatformCoreTargets',[21,21], ...
    'seedCheckpoint',currentFile,'solver',which('ipm.solve'), ...
    'regrid',which('ipm_gridlab_regrid_checkpoint'),'maxStages',maxStages, ...
    'deltaTau',0.2,'preferredTargetCoreCells',[21,21], ...
    'fallbackTargetCoreCells',[21,21],'regridBelowCoreCells',20, ...
    'minimumTerminalCoreCells',14,'maximumSafety',0.70, ...
    'maximumSinglePeakJump',2e-3,'maximumCumulativePeakJump',2e-2));
try
    for index=1:maxStages
        cp=ipm.output.readCheckpoint(currentFile);
        h=cp.payload.log.history;
        assert(all(ipm.output.trustedMask(h,cp.payload.state.config)));
        currentCores=[h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)];
        assert(all(currentCores>=14) && h.mesh.safetyFactor(end)<0.70);
        [coreDecay,safetySlope]=local_mesh_trend(h);
        predictedCores=currentCores.*exp(-0.2*coreDecay);
        stageDirectory=fullfile(outputDirectory,sprintf('stage_%03d',index));
        mkdir(stageDirectory);
        cumulativeJump=cumulative_jump(cp);
        assert(cumulativeJump<=2e-2);
        regridFile='';
        if min(currentCores)<20 || min(predictedCores)<16 || ...
                h.mesh.safetyFactor(end)+0.2*safetySlope>=0.65
            emit(manifest,struct('event','design_started','stage',index, ...
                'checkpointFile',currentFile,'tau',cp.payload.state.scale.canonicalTime));
            % A checkpoint-to-result conversion in the existing designer
            % reconstructs operators. A paired terminal result avoids that
            % extra build whenever one is available from this campaign.
            currentInput=currentFile;
            if index==1 && ~isempty(seedResultFile)
                currentInput=seedResultFile;
            elseif index>1 && isfile(lastResultFile)
                currentInput=lastResultFile;
            end
            [candidate,candidateFile,chosenTarget]=design_pair( ...
                currentInput,rootResult,stageDirectory,manifest,index,families,platform);
            assert(candidate.transactionReady,'Frozen regrid candidate rejected.');
            assert(cumulativeJump+candidate.pairScore.aggregate.worst. ...
                rhoXMaximumRelativeChange<=2e-2,'Cumulative peak-jump budget exceeded.');
            limits=struct('maxAdjacentCellRatio',1.08, ...
                'maxLogSpacingCurvature',0.01,'minStencilRcond',1e-9, ...
                'minQuadratureWeightRatio',1e-8);
            options=struct('tag',sprintf('longtime_20260908_stage_%03d',index), ...
                'minimumXCoreCells',chosenTarget-1,'minimumYCoreCells',chosenTarget-1, ...
                'maximumRhoXRelativeChange',2e-3,'maximumMassRelativeDefect',5e-12, ...
                'maximumRelativeRangeViolation',2e-4,'meshLimits',limits, ...
                'outputFile','');
            [cp,regridAudit,~]=ipm_gridlab_regrid_checkpoint( ...
                cp,candidate.candidateX,candidate.candidateY,options);
            assert(regridAudit.passed);
            cumulativeJump=cumulative_jump(cp);
            assert(cumulativeJump<=2e-2);
            [cp,regridFile]=ipm.output.writeCheckpoint(cp, ...
                fullfile(stageDirectory,'regrid_checkpoint.mat'));
            currentFile=char(regridFile);
            save(fullfile(stageDirectory,'regrid_audit.mat'),'regridAudit','candidateFile');
            emit(manifest,struct('event','regrid_accepted','stage',index, ...
                'checkpointFile',currentFile,'cumulativePeakJump',cumulativeJump, ...
                'singlePeakJump',regridAudit.rhoXMaximumRelativeChange));
            clear candidate
        end
        currentMesh=cp.payload.log.history.mesh;
        currentCores=[currentMesh.coreGridPoints(end),currentMesh.verticalCoreGridPoints(end)];
        permitted=0.75*log(currentCores/16)./max(coreDecay,eps);
        deltaTau=min([0.2,permitted]);
        assert(deltaTau>=0.01,'Predicted mesh loss requires a stronger design.');
        targetTau=cp.payload.state.scale.canonicalTime+deltaTau;
        sourceStep=cp.payload.state.step;
        sourceHistory=cp.payload.log.history;
        sourceFile=currentFile;
        options=struct('finalTime',targetTau,'physicalFinalTime',10, ...
            'maxSteps',sourceStep+100000,'saveResults',true, ...
            'resultFile',fullfile(stageDirectory,'result.mat'), ...
            'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false, ...
            'checkpoint',struct('file',fullfile(stageDirectory,'checkpoint.mat'), ...
            'every',0.1,'atExit',true));
        emit(manifest,struct('event','solve_started','stage',index, ...
            'checkpointFile',sourceFile,'targetTau',targetTau));
        timer=tic; result=ipm.solve(options,cp); wallSeconds=toc(timer);
        nextFile=char(result.metadata.latestCheckpointFile);
        terminal=ipm.output.readCheckpoint(nextFile);
        assert(isequal(terminal.payload.state.rho,result.state.rho));
        assert_history_prefix(sourceHistory,result.history);
        trusted=ipm.output.trustedMask(result);
        c=result.history.common; m=result.history.mesh;
        audit=struct('event','endpoint','stage',index,'targetTau',targetTau, ...
            'tau',result.state.canonicalTime,'physicalTime',result.state.physicalTime, ...
            'step',result.state.steps,'segmentSteps',result.state.steps-sourceStep, ...
            'physicalRhoXInf',c.physicalRhoXInf(end), ...
            'growthEFolds',log(c.physicalRhoXInf(end)/c.physicalRhoXInf(1)), ...
            'cL',c.c_l(end),'cOmega',c.c_omega(end), ...
            'coreCells',[m.coreGridPoints(end),m.verticalCoreGridPoints(end)], ...
            'safety',m.safetyFactor(end),'allHistoryTrusted',all(trusted), ...
            'cumulativePeakJump',cumulativeJump,'checkpointFile',nextFile, ...
            'resultFile',char(result.metadata.resultFile),'sourceCheckpoint',sourceFile, ...
            'regridFile',char(regridFile),'wallSeconds',wallSeconds, ...
            'stopReason',result.state.stopReason);
        emit(manifest,audit);
        assert(all(trusted) && all(audit.coreCells>=14) && audit.safety<0.70, ...
            'Terminal resolution/trust gate requires a new design.');
        assert(abs(audit.tau-targetTau)<=1e-10,'Solver stopped before the review endpoint.');
        currentFile=nextFile;
        lastResultFile=char(result.metadata.resultFile);
        clear cp result terminal sourceHistory h c m
    end
    emit(manifest,struct('event','batch_review_required','checkpointFile',currentFile));
catch exception
    emit(manifest,struct('event','review_required','checkpointFile',currentFile, ...
        'identifier',exception.identifier,'message',exception.message));
    rethrow(exception)
end
end

function [decay,safetySlope]=local_mesh_trend(h)
tau=h.common.canonicalTau(:); m=h.mesh;
mask=m.remeshCount(:)==m.remeshCount(end) & tau>=tau(end)-0.35;
decay=[0,0]; safetySlope=0;
if nnz(mask)<4, return; end
time=tau(mask)-tau(end);
names={'coreGridPoints','verticalCoreGridPoints'};
for k=1:2
    values=m.(names{k})(:);
    fit=polyfit(time,log(values(mask)),1);
    decay(k)=max(0,-fit(1));
end
fit=polyfit(time,m.safetyFactor(mask),1);
safetySlope=max(0,fit(1));
end

function [candidate,candidateFile,target]=design_pair( ...
        currentInput,rootResult,stageDirectory,manifest,stage,families,platform)
% Moving the fine platform frees the anchor node index while retaining
% exact x=1, the root provenance, and all independent pair/transfer gates.
target=21;
directory=fullfile(stageDirectory,'platform_design');
controls=struct('targetXCoreCells',target,'targetYCoreCells',target, ...
    'xEqualizationSigma',.5);
mesh_design_stage(currentInput,rootResult,directory,controls);
for variant=1:numel(platform)
    if variant==1
        candidateFile=fullfile(directory,'moving_platform.mat');
    else
        candidateFile=fullfile(directory,sprintf('moving_platform_%02d.mat',variant));
    end
    try
        candidate=mesh_platform_candidate(fullfile(directory,'mesh_design.mat'), ...
            fullfile(directory,'mesh_root_reference.mat'),candidateFile,platform(variant));
        if candidate.transactionReady
            emit(manifest,struct('event','design_accepted','stage',stage, ...
                'family','moving_fine_platform','targetCoreCells',target, ...
                'controls',platform(variant),'candidateFile',candidateFile));
            return
        end
        emit(manifest,struct('event','design_rejected','stage',stage, ...
            'family','moving_fine_platform','controls',platform(variant), ...
            'candidateFile',candidateFile,'status','pair_or_transfer_gate'));
    catch exception
        expected={'ipm:SmoothPeakGridQuadrature','ipm:SmoothPeakGridCellRatio', ...
            'ipm:SmoothPeakGridStencilCondition','ipm:SmoothPeakGridQuadratureMargin', ...
            'ipm:SmoothPeakGridLocalQuadratureMargin','ipm:CoordinatedGridFineInterval', ...
            'ipm:CoordinatedGridCellBudget','ipm:CoordinatedGridInfeasible', ...
            'ipm:CoordinatedGridBranchInfeasible','ipm:CoordinatedGridGrowthOverflow'};
        if ~ismember(exception.identifier,expected), rethrow(exception); end
        emit(manifest,struct('event','design_rejected','stage',stage, ...
            'family','moving_fine_platform','controls',platform(variant), ...
            'identifier',exception.identifier,'message',exception.message));
    end
end
% Every fallback keeps its own artifact and the same hard gates.
for trial=1:size(families,1)
    target=families(trial,1); sigma=families(trial,2);
    directory=fullfile(stageDirectory,sprintf('design_%02d',trial));
    controls=struct('targetXCoreCells',target,'targetYCoreCells',target, ...
        'xEqualizationSigma',sigma);
    report=mesh_design_stage(currentInput,rootResult,directory,controls);
    if report.transactionReady
        candidateFile=fullfile(directory,'anchor1_balanced.mat');
        candidate=mesh_balance_y_candidate(fullfile(directory,'mesh_design.mat'), ...
            fullfile(directory,'mesh_root_reference.mat'),candidateFile,1);
        if candidate.transactionReady
            emit(manifest,struct('event','design_accepted','stage',stage, ...
                'trial',trial,'targetCoreCells',target,'sigma',sigma, ...
                'candidateFile',candidateFile));
            return
        end
    end
    emit(manifest,struct('event','design_rejected','stage',stage, ...
        'trial',trial,'targetCoreCells',target,'sigma',sigma,'status',report.status));
end
error('ipm:LongtimeMeshFamilyExhausted', ...
    'All registered candidates failed unchanged grid/transfer gates; preserve the source checkpoint.');
end

function total=cumulative_jump(cp)
metadata=cp.payload.state.runMetadata;
if isfield(metadata,'gridLabRegrids')
    records=metadata.gridLabRegrids;
    assert(numel(records)==cp.payload.state.remeshCount, ...
        'All prior remesh transactions must have recorded jump provenance.');
    jumps=[records.rhoXMaximumRelativeChange];
    assert(all(isfinite(jumps)) && all(jumps>=0) && all(jumps<=2e-3));
    total=sum(jumps);
else
    assert(cp.payload.state.remeshCount==0,'Missing inherited remesh provenance.');
    total=0;
end
end

function assert_history_prefix(old,new)
groups=fieldnames(old);
for gi=1:numel(groups)
    names=fieldnames(old.(groups{gi}));
    for ni=1:numel(names)
        a=old.(groups{gi}).(names{ni}); b=new.(groups{gi}).(names{ni});
        assert(isequaln(a,b(1:numel(a))),'Same-grid history prefix changed.');
    end
end
end

function emit(file,record)
record.utc=char(datetime('now','TimeZone','UTC', ...
    'Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
fid=fopen(file,'a'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
fprintf('ADAPTIVE %s\n',jsonencode(record));
end
