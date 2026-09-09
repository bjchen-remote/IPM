function report=ipm_perflab_fresh_stage5_native_trial(candidateFile,sourceReviewFile,outDir,executeNative)
%IPM_PERFLAB_FRESH_STAGE5_NATIVE_TRIAL One registered zero-time native trial.
% Default is a no-LU preflight. Only an explicit true executes the maintained
% transaction; it builds one factor at a time and never advances the PDE.
if nargin<4,executeNative=false;end
assert(islogical(executeNative)&&isscalar(executeNative)&&maxNumCompThreads==10);
assert(~isfolder(outDir));mkdir(outDir);
try
    d=load(candidateFile,'candidate','axes','policy','sourceFile','priorJump','gates');
    review=jsondecode(fileread(sourceReviewFile));c=d.candidate;
    assert(strcmp(d.sourceFile,review.checkpointFile)&&review.step==9205&& ...
        d.priorJump==review.cumulativeAbsolutePeakJump&&c.transactionReady&& ...
        c.anchorExact&&c.completePairScoredInOriginalUnits&& ...
        c.currentTransactionNodeCountSupported&&~c.nativeTransactionPerformed&&~c.pdeAdvanced&& ...
        all(structfun(@(v)islogical(v)&&isscalar(v)&&v,d.gates)));
    limits=struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
        'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8);
    assert(isequaln(rmfield(d.policy.qualityLimits, ...
        {'minWeightToControlWidth','maxWeightToControlWidth'}),limits));
    options=struct('tag','fresh_stage005_gap_family_candidate', ...
        'minimumXCoreCells',31,'minimumYCoreCells',31,'maximumRhoXRelativeChange',.002, ...
        'maximumMassRelativeDefect',5e-12,'maximumRelativeRangeViolation',2e-4, ...
        'meshLimits',limits,'outputFile','');
    registration=struct('kind','fresh_stage005_registered_native_trial_v1', ...
        'candidateFile',candidateFile,'sourceReviewFile',sourceReviewFile, ...
        'sourceFile',d.sourceFile,'priorCumulativePeakJump',d.priorJump, ...
        'maximumCumulativePeakJump',.02,'minimumFrontCells',20, ...
        'historyCoreFloor',14,'maximumSafety',.70,'options',options, ...
        'executeNative',executeNative,'advancePDE',false,'sourceFilesReadOnly',true, ...
        'sequentialFactorLifecycle','maintained regrid helper drops old poisson before new transfer/build');
    write_json(fullfile(outDir,'registration.json'),registration);
    issues=checkcode(which(mfilename),'-id');write_json(fullfile(outDir,'checkcode.json'),issues);
    assert(isempty(issues),jsonencode(issues));
    source=ipm.output.readCheckpoint(d.sourceFile);o=source.payload.state;h=source.payload.log.history;
    assert(o.step==9205&&strcmp(o.runMetadata.caseId,review.caseId)&& ...
        o.scale.canonicalTime==review.localCanonicalTime&&o.scale.physicalTime==review.localPhysicalTime&& ...
        all(ipm.output.trustedMask(h,o.config))&&all(h.mesh.coreGridPoints>=14)&& ...
        all(h.mesh.verticalCoreGridPoints>=14)&&all(h.mesh.safetyFactor<.70));
    assert(isequaln([h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)],review.coreCells(:)')&& ...
        all(review.coreCells>=20)&&h.mesh.safetyFactor(end)==review.safety);
    records=o.runMetadata.gridLabRegrids;jumps=[records.rhoXMaximumRelativeChange];
    assert(numel(records)==o.remeshCount&&all(isfinite(jumps))&&all(jumps>=0)&&all(jumps<=.002)&& ...
        sum(abs(jumps))==d.priorJump&&d.priorJump<=.02);
    assert(numel(o.x)==895&&numel(o.y)==386&& ...
        isequal(c.candidateX,d.axes.x)&&isequal(c.candidateY,d.axes.y)&& ...
        isequal(c.candidateX([1,end]),o.x([1,end]))&&isequal(c.candidateY([1,end]),o.y([1,end]))&& ...
        any(c.candidateX==o.rescaling.transportAnchorX)&&any(c.candidateX==-o.rescaling.transportAnchorX));
    assert(~isfield(o.config.remesh,'autonomousMesh'), ...
        'ipm:FreshStage5Scope','This prepared research transaction is not an autonomous-policy continuation.');
    a=c.pairScore.aggregate;
    assert(c.pairScore.admissible&&a.minimumXCoreCells>=32-1e-6&&a.minimumYCoreCells>=32-1e-6&& ...
        a.minimumXLeftFrontCells>=20-1e-6&&a.worst.rhoXMaximumRelativeChange<=.002&& ...
        d.priorJump+a.worst.rhoXMaximumRelativeChange<=.02&& ...
        a.worst.conservationRelativeDefect<=5e-12&&a.worst.relativeRangeViolation<=2e-4&& ...
        local_quality(c.pairScore.quality.x)&&local_quality(c.pairScore.quality.y));
    report=registration;report.sourceNativeSignatureValidated=true;report.preflightPassed=true;
    report.nativeTransactionPerformed=false;report.checkpointWritten=false;
    if ~executeNative
        report.status='prepared_only_no_lu';write_json(fullfile(outDir,'report.json'),report);
        return
    end
    [proposal,audit]=ipm_gridlab_regrid_checkpoint(source,c.candidateX,c.candidateY,options);
    report.nativeTransactionPerformed=true;
    n=proposal.payload.state;nh=proposal.payload.log.history;
    % Preserve all nongeometric state, immutable references and mathematical clocks.
    for name={'step','scale','normalizedTime','mass0','rhoRange0','baseX','baseY'}
        assert(isequaln(n.(name{1}),o.(name{1})));
    end
    assert(n.remeshCount==o.remeshCount+1&&isequal(n.x,c.candidateX)&&isequal(n.y,c.candidateY)&& ...
        isequaln(rmfield(n.rescaling,{'originIndex','pinIndex'}),rmfield(o.rescaling,{'originIndex','pinIndex'})));
    ng=n.config;og=o.config;ng.grid=rmfield(ng.grid,{'customX','customY'});og.grid=rmfield(og.grid,{'customX','customY'});
    assert(fresh_profile_value_equal(ng,og));
    assert(isequaln(n.runMetadata.caseMetadata.latePhysicalBoxBranch,o.runMetadata.caseMetadata.latePhysicalBoxBranch));
    nm=n.runMetadata;om=o.runMetadata;
    for name={'gridLabRegrids','latestGridLabRegridTag'}
        if isfield(nm,name{1}),nm=rmfield(nm,name{1});end
        if isfield(om,name{1}),om=rmfield(om,name{1});end
    end
    assert(fresh_profile_value_equal(nm,om));
    assert(isequaln(n.runMetadata.gridLabRegrids(1:end-1),o.runMetadata.gridLabRegrids));
    fresh_profile_history_prefix(h,nh,true);
    for name={'acceptedStep','t','canonicalTau','physicalTime'}
        assert(isequaln(h.common.(name{1}),nh.common.(name{1})));
    end
    assert(all(ipm.output.trustedMask(nh,n.config))&&all(nh.mesh.coreGridPoints>=14)&& ...
        all(nh.mesh.verticalCoreGridPoints>=14)&&all(nh.mesh.safetyFactor<.70));
    snapshot_prefix(source.payload.log,proposal.payload.log,n.config.output.storeSnapshots,n);
    oldCursor=source.payload.cursor;newCursor=proposal.payload.cursor;
    oldCursor.lastCheckpointStep=n.step;oldCursor.lastCheckpointCanonicalTime=n.scale.canonicalTime;
    assert(isequaln(oldCursor,newCursor));
    assert(audit.passed&&audit.xCoreCells>=31&&audit.yCoreCells>=31&& ...
        local_quality(audit.xQuality)&&local_quality(audit.yQuality));
    snapshot=struct('datasetIndex',1,'sourceLabel','actual unwritten zero-time transaction', ...
        'trusted',true,'rho',n.rho,'x',n.x,'y',n.y, ...
        'scale',struct('Cx',exp(n.scale.logC_l),'Cy',exp(n.scale.logC_l),'Comega',exp(n.scale.logC_omega)), ...
        'canonicalTime',n.scale.canonicalTime,'physicalTime',n.scale.physicalTime);
    dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','actual_fresh_transaction', ...
        'snapshotCount',1,'snapshots',snapshot,'xLimits',n.x([1,end]),'yLimits',n.y([1,end])');
    actual=ipm_gridlab_score_frozen_pair(dataset,n.x,n.y,struct('trustedOnly',true,'meshLimits',limits));
    assert(actual.admissible&&actual.aggregate.minimumXLeftFrontCells>=20&& ...
        d.priorJump+audit.rhoXMaximumRelativeChange<=.02);
    audit.actualLeftFrontCells=actual.aggregate.minimumXLeftFrontCells;
    audit.cumulativeAbsolutePeakJump=d.priorJump+audit.rhoXMaximumRelativeChange;
    save(fullfile(outDir,'unwritten_native_audit.mat'),'proposal','audit','actual','candidateFile','-v7.3');
    [accepted,file]=ipm.output.writeCheckpoint(proposal,fullfile(outDir,'regrid_checkpoint.mat'));
    validated=ipm.output.readCheckpoint(file);
    assert(fresh_profile_value_equal(validated.payload,accepted.payload));
    report.status='native_zero_time_transaction_passed';report.checkpointWritten=true;report.checkpointFile=char(file);
    report.actualAudit=audit;report.allHistoryTrusted=true;report.immutableReferencesAndClocksExact=true;
    report.signedNativeRoundTripPayloadExact=true;report.shortPDEValidationPerformed=false;
    save(fullfile(outDir,'report.mat'),'report','audit','actual','-v7.3');
    write_json(fullfile(outDir,'report.json'),report);
catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(outDir,'failure.mat'),'failure','exception','-v7.3');
    write_json(fullfile(outDir,'failure.json'),failure);rethrow(exception)
end
end
function snapshot_prefix(old,new,stored,state)
assert(isequaln(old.snapshotNormalizedTime,new.snapshotNormalizedTime));
for field={'snapshotRho','snapshotX','snapshotY'}
    a=old.(field{1});b=new.(field{1});assert(isequal(size(a),size(b)));
    if stored,assert(isequaln(a(1:end-1),b(1:end-1)));else,assert(isempty(a)&&isempty(b));end
end
if stored
    assert(isequaln(new.snapshotRho{end},state.rho)&&isequaln(new.snapshotX{end},state.x)&&isequaln(new.snapshotY{end},state.y));
end
exclude={'history','snapshotRho','snapshotX','snapshotY','snapshotNormalizedTime'};
assert(fresh_profile_value_equal(rmfield(old,exclude),rmfield(new,exclude)));
end
function yes=local_quality(q)
yes=q.quadratureWeightsStrictlyPositive&&q.minimumQuadratureWeightToControlWidthRatio>=.35&& ...
    q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
