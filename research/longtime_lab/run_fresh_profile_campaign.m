function campaign=run_fresh_profile_campaign(bootstrapResultFile,initialDataFile,outputDirectory,maxStages,qualificationFile)
%RUN_FRESH_PROFILE_CAMPAIGN Audited independent fresh-box profile continuation.
% Requires a caller-written qualification ledger; preserves legacy false and
% separate prospective auxiliary evidence. Never claims parent-native history.
assert(maxNumCompThreads==10,'Run with the validated original ten threads.');
assert(~isfolder(outputDirectory),'Use a new immutable campaign directory.');
validateattributes(maxStages,{'numeric'},{'scalar','integer','positive','<=',200});
assert(nargin==5 && isfile(qualificationFile),'A caller-provided qualification JSON is mandatory.');
[review,cp,contract]=fresh_profile_campaign_review(bootstrapResultFile,initialDataFile);
qualification=fresh_profile_campaign_admit(qualificationFile,bootstrapResultFile,initialDataFile,cp);
assert(~isfield(qualification,'testFixtureOnly') || ~qualification.testFixtureOnly, ...
    'ipm:FreshCampaignQualification','A test fixture is not a caller campaign qualification.');
assert(all(review.coreCells>=20),'ipm:FreshCampaignBootstrap','Bootstrap must meet the stage endpoint core20 floor.');
firstTargetLocalTau=cp.payload.state.scale.canonicalTime+.2/contract.lineage.canonicalCovarianceFactor;
firstEquivalentTauStep=.2;
if isfield(qualification,'firstTargetLocalCanonicalTime')
    firstTargetLocalTau=qualification.firstTargetLocalCanonicalTime;
    validateattributes(firstTargetLocalTau,{'double'},{'scalar','real','finite','positive'});
    firstEquivalentTauStep=contract.lineage.canonicalCovarianceFactor* ...
        (firstTargetLocalTau-cp.payload.state.scale.canonicalTime);
    assert(isfield(qualification,'budgetLimitedCheckpointFiles') && ...
        any(strcmp(qualification.budgetLimitedCheckpointFiles,review.checkpointFile)) && ...
        firstTargetLocalTau==cp.payload.state.config.time.finalTime && ...
        firstEquivalentTauStep>0 && firstEquivalentTauStep<=.2, ...
        'ipm:FreshCampaignFirstTarget','A shortened first stage must finish the explicitly qualified source scheduling target.');
end
mkdir(outputDirectory);manifest=fullfile(outputDirectory,'campaign_manifest.jsonl');
families=ipm.config.autonomousMeshPolicy(struct());
registration=struct('kind','independent_fresh_profile_campaign_v3','bootstrapResultFile',bootstrapResultFile, ...
    'initialDataFile',initialDataFile,'qualificationFile',qualificationFile,'qualification',qualification, ...
    'initialCheckpoint',review.checkpointFile,'caseId',contract.caseId,'lineage',contract.lineage, ...
    'maxStages',maxStages,'equivalentTauStep',.2,'firstEquivalentTauStep',firstEquivalentTauStep, ...
    'freshCanonicalStep',review.localCanonicalStep, ...
    'targetCoreCells',[32,32],'transactionMinimumCoreCells',[31,31],'minimumLeftFrontCells',20, ...
    'regridBelowCoreCells',26,'predictedEndpointTrigger',22,'trendEquivalentTauWindow',.35, ...
    'minimumStageEndpointCoreCells',20,'maximumSafety',.70,'maximumAdditionalStepsPerStage',4096, ...
    'maximumSinglePeakJump',.002,'maximumCumulativeAbsolutePeakJump',.02,'registeredFamilies',families.search,'meshPolicy',families, ...
    'fixedNodeCount',[895,386],'frozenCFL',contract.CFL,'frozenMaxDt',contract.maxDt, ...
    'parentNativeContinuationClaim',false,'bootstrapLegacyDecisionChanged',false, ...
    'solver',which('ipm.solve'),'designer',which('ipm.remesh.plannedAxisPairs'), ...
    'regrid',which('ipm_gridlab_regrid_checkpoint'),'threads',10);
write_json(fullfile(outputDirectory,'registration.json'),registration);
save(fullfile(outputDirectory,'registration.mat'),'registration','contract','-v7.3');
currentFile=review.checkpointFile;currentResultFile=char(bootstrapResultFile);stageCount=0;
campaign=struct('status','registered','stageCount',0,'latestCheckpointFile',currentFile, ...
    'latestPairedResultFile',currentResultFile,'qualificationFile',qualificationFile, ...
    'parentNativeContinuationClaim',false,'asymptoticGoalCompleted',false);
emit(manifest,struct('event','registered','registrationFile',fullfile(outputDirectory,'registration.json'),'review',review));
try
 for stage=1:maxStages
    if stage>1,[review,cp]=fresh_profile_campaign_review(currentResultFile,initialDataFile,contract);end
    assert(strcmp(currentFile,review.checkpointFile) && all(review.coreCells>=20));
    directory=fullfile(outputDirectory,sprintf('stage_%03d',stage));mkdir(directory);
    save(fullfile(directory,'source_review.mat'),'review','-v7.3');
    write_json(fullfile(directory,'source_review.json'),review);
    if review.regridRequired
        [cp,currentFile,transaction]=try_regrid(cp,currentFile,contract, ...
            directory,families,manifest,stage,review.cumulativeAbsolutePeakJump);
        campaign.latestCheckpointFile=currentFile;campaign.latestPairedResultFile='';
        emit(manifest,struct('event','regrid_accepted','stage',stage,'checkpointFile',currentFile,'transaction',transaction));
    end
    sourceHistory=cp.payload.log.history;sourceReferences=cp.payload.state.rescaling;
    sourceStep=cp.payload.state.step;sourceFile=currentFile;
    targetLocalTau=cp.payload.state.scale.canonicalTime+.2/contract.lineage.canonicalCovarianceFactor;
    if stage==1,targetLocalTau=firstTargetLocalTau;end
    targetEquivalentTau=contract.lineage.parentCanonicalTime+contract.lineage.canonicalCovarianceFactor*targetLocalTau;
    options=struct('finalTime',targetLocalTau,'physicalFinalTime',Inf,'maxSteps',sourceStep+registration.maximumAdditionalStepsPerStage, ...
        'saveResults',true,'resultFile',fullfile(directory,'result.mat'), ...
        'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false, ...
        'checkpoint',struct('enabled',true,'file',fullfile(directory,'checkpoint.mat'), ...
        'every',.1/contract.lineage.canonicalCovarianceFactor,'atExit',true));
    emit(manifest,struct('event','solve_started','stage',stage,'checkpointFile',sourceFile, ...
        'targetLocalCanonicalTime',targetLocalTau,'targetParentEquivalentTau',targetEquivalentTau));
    timer=tic;result=ipm.solve(options,cp);wall=toc(timer);nextResultFile=char(result.metadata.resultFile);
    stopReason=char(result.state.stopReason);clear result cp
    [endpoint,terminal]=fresh_profile_campaign_review(nextResultFile,initialDataFile,contract);
    fresh_profile_history_prefix(sourceHistory,terminal.payload.log.history,false);
    assert(isequaln(sourceReferences,terminal.payload.state.rescaling), ...
        'ipm:FreshCampaignInvariant','Runtime references changed during same-grid evolution.');
    endpoint.stage=stage;endpoint.segmentSteps=endpoint.step-sourceStep;endpoint.wallSeconds=wall;
    endpoint.sourceCheckpoint=sourceFile;endpoint.stopReason=stopReason;
    endpoint.targetParentEquivalentTau=targetEquivalentTau;
    endpoint.parentEquivalentEndpointError=abs(endpoint.parentEquivalentTau-targetEquivalentTau);
    endpoint.stageAccepted=all(endpoint.coreCells>=20) && endpoint.safety<.70 && ...
        endpoint.parentEquivalentEndpointError<=1e-10 && strcmp(stopReason,'final_time') && ...
        endpoint.segmentSteps<=registration.maximumAdditionalStepsPerStage;
    write_json(fullfile(directory,'endpoint_audit.json'),endpoint);save(fullfile(directory,'endpoint_audit.mat'),'endpoint','-v7.3');
    emit(manifest,struct('event','endpoint','audit',endpoint));
    % A normally verified but early/under-buffer endpoint remains recoverable;
    % it is saved before the stricter campaign endpoint decision stops work.
    currentFile=endpoint.checkpointFile;currentResultFile=nextResultFile;
    campaign.latestCheckpointFile=currentFile;campaign.latestPairedResultFile=currentResultFile;
    assert(endpoint.stageAccepted,'ipm:FreshCampaignEndpoint','The stage did not reach its qualified review endpoint.');
    stageCount=stage;campaign.stageCount=stageCount;
    clear terminal sourceHistory sourceReferences endpoint
 end
 campaign.status='batch_review_required';
 emit(manifest,struct('event','batch_review_required','checkpointFile',currentFile,'resultFile',currentResultFile,'stageCount',stageCount));
catch exception
 campaign.status='review_required';campaign.stageCount=stageCount;
 campaign.failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
 write_json(fullfile(outputDirectory,'campaign_status.json'),campaign);
 emit(manifest,struct('event','review_required','campaign',campaign));rethrow(exception)
end
write_json(fullfile(outputDirectory,'campaign_status.json'),campaign);
end

function [accepted,file,audit]=try_regrid(source,sourceFile,contract,directory,families,manifest,stage,priorJump)
limits=struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
    'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8);
o=source.payload.state;Dx=ipm.mesh.fdMatrix(o.x,1,7);
view=struct('rho',o.rho,'x',o.x,'y',o.y,'Dx',Dx,'source',o.rho*Dx','trusted',true);
[axes,axisReport]=ipm.remesh.plannedAxisPairs(view,struct('x',o.baseX,'y',o.baseY), ...
    contract.lineage.physicalAnchorX,families);
save(fullfile(directory,'axis_plan.mat'),'axes','axisReport','-v7.3');
write_json(fullfile(directory,'axis_plan.json'),axisReport);
clear view Dx o
for trial=1:numel(axes)
    trialDirectory=fullfile(directory,sprintf('design_%02d',trial));committed=false;
    try
        [candidate,design]=score_axis_pair(source,sourceFile,axes(trial),families,trialDirectory);
        if ~candidate.transactionReady
            emit(manifest,struct('event','design_rejected','stage',stage,'trial',trial,'candidateFile',design.candidateFile,'reason','frozen_pair_gate'));continue
        end
        assert(candidate.anchorExact && candidate.completePairScoredInOriginalUnits && ...
            candidate.currentTransactionNodeCountSupported && candidate.controls.targetXCoreCells==32 && candidate.controls.targetYCoreCells==32);
        if priorJump+candidate.pairScore.aggregate.worst.rhoXMaximumRelativeChange>.02
            emit(manifest,struct('event','design_rejected','stage',stage,'trial',trial,'reason','predicted_cumulative_jump'));continue
        end
        options=struct('tag',sprintf('fresh_profile_stage_%03d_trial_%02d',stage,trial), ...
            'minimumXCoreCells',31,'minimumYCoreCells',31,'maximumRhoXRelativeChange',.002, ...
            'maximumMassRelativeDefect',5e-12,'maximumRelativeRangeViolation',2e-4,'meshLimits',limits,'outputFile','');
        [proposal,audit]=ipm_gridlab_regrid_checkpoint(source,candidate.candidateX,candidate.candidateY,options);
        assert(audit.passed && audit.xCoreCells>=31 && audit.yCoreCells>=31);
        n=proposal.payload.state;o=source.payload.state;
        assert(n.step==o.step && n.normalizedTime==o.normalizedTime && isequaln(n.scale,o.scale) && ...
            strcmp(n.runMetadata.caseId,contract.caseId) && ...
            isequaln(n.runMetadata.caseMetadata.latePhysicalBoxBranch,contract.lineage) && ...
            isequaln(rmfield(n.rescaling,{'originIndex','pinIndex'}),contract.runtimeReferences), ...
            'ipm:FreshCampaignInvariant','A transaction changed fresh lineage, clocks or references.');
        assert(numel(n.x)==895 && numel(n.y)==386 && isequal(n.baseX,o.baseX) && isequal(n.baseY,o.baseY) && ...
            isequal(n.x([1,end]),contract.initialX([1,end])) && isequal(n.y([1,end]),contract.initialY([1,end])) && ...
            any(n.x==contract.lineage.physicalAnchorX) && any(n.x==-contract.lineage.physicalAnchorX), ...
            'ipm:FreshCampaignInvariant','A transaction changed the registered box/base grid or lost the exact anchor.');
        fresh_profile_history_prefix(source.payload.log.history,proposal.payload.log.history,true);
        assert(all(ipm.output.trustedMask(proposal.payload.log.history,n.config)) && ...
            local_quality(audit.xQuality) && local_quality(audit.yQuality));
        snapshot=struct('datasetIndex',1,'sourceLabel','unwritten actual native transaction', ...
            'trusted',true,'rho',n.rho,'x',n.x,'y',n.y,'scale',struct('Cx',exp(n.scale.logC_l), ...
            'Cy',exp(n.scale.logC_l),'Comega',exp(n.scale.logC_omega)), ...
            'canonicalTime',n.scale.canonicalTime,'physicalTime',n.scale.physicalTime);
        dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','actual_fresh_transaction', ...
            'snapshotCount',1,'snapshots',snapshot,'xLimits',n.x([1,end]),'yLimits',n.y([1,end])');
        actual=ipm_gridlab_score_frozen_pair(dataset,n.x,n.y,struct('trustedOnly',true,'meshLimits',limits));
        assert(actual.admissible && actual.aggregate.minimumXLeftFrontCells>=20 && priorJump+audit.rhoXMaximumRelativeChange<=.02);
        audit.actualLeftFrontCells=actual.aggregate.minimumXLeftFrontCells;audit.cumulativeAbsolutePeakJump=priorJump+audit.rhoXMaximumRelativeChange;
        [accepted,file]=ipm.output.writeCheckpoint(proposal,fullfile(directory,'regrid_checkpoint.mat'));committed=true;
        validated=ipm.output.readCheckpoint(file);assert(fresh_profile_value_equal(validated.payload,accepted.payload), ...
            'ipm:FreshCampaignInvariant','Native transaction round-trip data changed.');
        candidateFile=design.candidateFile;save(fullfile(directory,'regrid_audit.mat'),'audit','candidateFile','file','-v7.3');
        write_json(fullfile(directory,'regrid_audit.json'),audit);file=char(file);return
    catch exception
        failure=struct('event','candidate_failure','stage',stage,'trial',trial,'axisIndices',[axes(trial).xIndex,axes(trial).yIndex], ...
            'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
        if ~isfolder(trialDirectory),mkdir(trialDirectory);end
        write_json(fullfile(trialDirectory,'campaign_failure.json'),failure);emit(manifest,failure);
        if committed || ~any(strcmp(exception.identifier,{'ipm:gridlab:RegridAudit','ipm:gridlab:RegridResolution', ...
                'ipm:gridlab:RegridRhoX','ipm:gridlab:RegridMass','ipm:gridlab:RegridRange', ...
                'ipm:gridlab:RegridMesh', ...
                'ipm:HighOrderQuadratureWeights','ipm:HighOrderGridMetric', ...
                'ipm:FiniteDifferenceScale','ipm:HighOrderRemeshBubble'}))
            rethrow(exception);
        end
    end
end
error('ipm:FreshCampaignMeshFamilyExhausted','All registered same-box candidates failed; preserve the last source native checkpoint.');
end

function [candidate,design]=score_axis_pair(source,sourceFile,axes,policy,directory)
mkdir(directory);n=source.payload.state;
snapshot=struct('datasetIndex',1,'sourceLabel',sourceFile,'trusted',true, ...
    'rho',n.rho,'x',n.x,'y',n.y,'scale',struct('Cx',exp(n.scale.logC_l), ...
    'Cy',exp(n.scale.logC_l),'Comega',exp(n.scale.logC_omega)), ...
    'canonicalTime',n.scale.canonicalTime,'physicalTime',n.scale.physicalTime);
dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','campaign_actual_source', ...
    'snapshotCount',1,'snapshots',snapshot,'xLimits',n.x([1,end]),'yLimits',n.y([1,end])');
limits=rmfield(policy.qualityLimits,{'minWeightToControlWidth','maxWeightToControlWidth'});
score=ipm_gridlab_score_frozen_pair(dataset,axes.x,axes.y,struct('trustedOnly',true, ...
    'meshLimits',limits,'minimumXCoreCells',32,'minimumYCoreCells',32,'minimumXFrontCells',20));
a=score.aggregate;
ready=score.admissible && a.minimumXCoreCells>=32-1e-6 && a.minimumYCoreCells>=32-1e-6 && ...
    a.minimumXLeftFrontCells>=20-1e-6 && a.worst.rhoXMaximumRelativeChange<=.002 && ...
    a.worst.conservationRelativeDefect<=5e-12 && a.worst.relativeRangeViolation<=2e-4;
candidate=struct('kind','campaign_bounded_paired_candidate_v1','transactionReady',ready, ...
    'candidateX',axes.x,'candidateY',axes.y,'anchorExact', ...
    any(axes.x==n.rescaling.transportAnchorX)&&any(axes.x==-n.rescaling.transportAnchorX), ...
    'completePairScoredInOriginalUnits',true,'currentTransactionNodeCountSupported',true, ...
    'controls',struct('targetXCoreCells',32,'targetYCoreCells',32,'minimumXFrontCells',20), ...
    'pairScore',score,'nativeTransactionPerformed',false,'pdeAdvanced',false);
design=struct('candidateFile',fullfile(directory,'candidate.mat'));
save(design.candidateFile,'candidate','axes','policy','sourceFile','-v7.3');
end

function yes=local_quality(q)
yes=q.quadratureWeightsStrictlyPositive && q.minimumQuadratureWeightToControlWidthRatio>=.35 && q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
function emit(file,value)
value.utc=char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
fid=fopen(file,'a');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
fprintf('FRESH_PROFILE %s\n',jsonencode(value));
end
