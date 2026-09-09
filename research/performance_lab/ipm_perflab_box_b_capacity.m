function report=ipm_perflab_box_b_capacity(checkpointFile,resultFile,outDir)
%IPM_PERFLAB_BOX_B_CAPACITY Preserved actual failure; cap-only geometry probe.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
registration=struct('kind','actual_B_H16_T32_capacity_and_registered_cap_counterfactual_v1', ...
    'checkpointFile',checkpointFile,'resultFile',resultFile,'originalCap',110000, ...
    'counterfactualCap',210000,'onlyCounterfactualPolicyChange','nodeFamily.maximumTotalNodes', ...
    'actualRegisteredPhysicalTarget',1.9140859724939365, ...
    'noLU',true,'noPDE',true,'noTransfer',true,'resumeQualification',false, ...
    'historicalFailedFieldReconstructed',false,'newRunFromZeroRequiredForCounterfactual',true);
write_json(fullfile(outDir,'registration.json'),registration);
issues=checkcode(which(mfilename),'-id');write_json(fullfile(outDir,'checkcode.json'),issues);
assert(isempty(issues),jsonencode(issues));
profile clear;profile on;timer=tic;
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;h=cp.payload.log.history;
m=s.runMetadata.autonomousMesh;p=s.config.remesh.autonomousMesh;
loaded=load(resultFile,'result');r=ipm.output.validate(loaded.result);clear loaded;
assert(s.step==1520&&m.currentLevelId==2&&p.version==2&& ...
    p.nodeFamily.maximumTotalNodes==110000&&strcmp(r.state.stopReason,'autonomous_mesh_axis_capacity'));
assert(isequaln(r.state.rho,s.rho)&&isequaln(r.grid.x,s.x)&&isequaln(r.grid.y,s.y)&& ...
    isequaln(r.config,s.config)&&isequaln(r.history,h)&& ...
    r.state.steps==s.step&&r.state.canonicalTime==s.scale.canonicalTime&& ...
    r.state.physicalTime==s.scale.physicalTime&&r.grid.remeshCount==s.remeshCount);
metadataPair=metadata_pair(r.metadata,s.runMetadata,checkpointFile);
write_json(fullfile(outDir,'metadata_pair.json'),metadataPair);assert(metadataPair.passed);
assert(all(ipm.output.trustedMask(r))&&s.scale.physicalTime<registration.actualRegisteredPhysicalTarget);
assert(h.common.acceptedStep(1)==0&&h.common.canonicalTau(1)==0&&h.common.physicalTime(1)==0);
anchor=s.config.scaling.transportAnchorX;
family=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX,m.initialization.selectedBaseY,anchor,p);
assert(isequaln(family,m.referenceFamily));member=family.members(m.currentLevelId);
assert(isequaln(s.baseX,member.baseX)&&isequaln(s.baseY,member.baseY));
eligible=eligible_levels(family,member);assert(isequal(eligible,2));
write_json(fullfile(outDir,'saved_last_failure.json'),m.lastFailure);
write_json(fullfile(outDir,'saved_last_decision.json'),m.lastDecision);
Dx=ipm.mesh.fdMatrix(s.x,1,7);Omega=s.rho*Dx';assert(isequaln(Omega,r.state.omega));
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',Omega,'trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);
assert(isequaln(feature.actualCoreCells,m.window.coreCells(end,:)));
[currentCandidates,currentAxisReport]=ipm.remesh.plannedAxisPairs(view, ...
    struct('x',member.baseX,'y',member.baseY),anchor,p);
counterfactualPolicy=p;counterfactualPolicy.nodeFamily.maximumTotalNodes=registration.counterfactualCap;
assert(isequaln(ipm.config.autonomousMeshPolicy(counterfactualPolicy),counterfactualPolicy));
newFamily=ipm.remesh.referenceAxisFamily(family.rootX,family.rootY,anchor,counterfactualPolicy);
assert(isequaln(rmfield(newFamily,'members'),rmfield(family,'members')));
for k=1:numel(family.members)
    assert(isequaln(rmfield(newFamily.members(k),'resourceAdmitted'), ...
        rmfield(family.members(k),'resourceAdmitted')));
end
newEligible=eligible_levels(newFamily,member);assert(isequal(newEligible,[2;4]));
newMember=newFamily.members(4);
[largerCandidates,largerAxisReport]=ipm.remesh.plannedAxisPairs(view, ...
    struct('x',newMember.baseX,'y',newMember.baseY),anchor,counterfactualPolicy);
report=registration;report.sourceNativeSignatureValidated=true;report.resultNativePayloadPairExact=true;
report.metadataPair=metadataPair;
report.originalStopReasonPreserved=r.state.stopReason;report.registeredEndpointReached=false;
report.allNativeHistoryTrusted=true;report.step=s.step;report.canonicalTime=s.scale.canonicalTime;
report.physicalTime=s.scale.physicalTime;report.remeshCount=s.remeshCount;
report.actualCoreCells=feature.actualCoreCells;report.frontCells=feature.leftFrontCells;
report.safety=m.window.safety(end);report.feature=feature;
report.currentAxisQualityX=ipm.mesh.quality(s.x,7,ipm.mesh.quadrature(s.x));
report.currentAxisQualityY=ipm.mesh.quality(s.y,7,ipm.mesh.quadrature(s.y));
report.current=plan_summary(currentCandidates,currentAxisReport);
report.originalEligibleLevels=eligible;report.originalResourceMask=[family.members.resourceAdmitted];
report.counterfactual=plan_summary(largerCandidates,largerAxisReport);
report.counterfactualEligibleLevels=newEligible;report.counterfactualResourceMask=[newFamily.members.resourceAdmitted];
report.counterfactualNodeCount=newMember.nodeCount;report.counterfactualNodeProduct=prod(newMember.nodeCount);
report.currentLevelId=m.currentLevelId;report.cumulativeAbsolutePeakJump=m.cumulativeAbsolutePeakJump;
report.savedLastFailure=m.lastFailure;report.savedLastDecision=m.lastDecision;
report.failedFieldIsSavedCurrentField=false;
if isfield(m.lastFailure,'sourceStep')
    report.failureStep=m.lastFailure.sourceStep;
    report.failureStepOffsetFromSavedState=m.lastFailure.sourceStep-s.step;
    report.failedFieldIsSavedCurrentField=m.lastFailure.sourceStep==s.step&& ...
        m.lastFailure.sourceCanonicalTime==s.scale.canonicalTime;
end
report.noIndependentRecalculationOfHistoricFailureField=true;
report.remainingGrowthTransitions=p.nodeFamily.maximumAcceptedGrowthTransitions- ...
    nnz([m.transactions.sourceLevelId]~=[m.transactions.targetLevelId]);
report.remainingRemeshCountBudget=s.config.remesh.maxRemeshes-s.remeshCount;
report.remainingPeakJumpBudget=p.maximumCumulativeAbsolutePeakJump-m.cumulativeAbsolutePeakJump;
report.wallSeconds=toc(timer);report.completed=true;
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.output.restoreCheckpoint','ipm.evolve.flow', ...
        'ipm.field.velocity','ipm.field.poisson','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})),'ipm:ResearchNoLU','Forbidden numerical call.');
end
save(fullfile(outDir,'report.mat'),'report','currentCandidates','currentAxisReport', ...
    'largerCandidates','largerAxisReport','profileInfo','-v7.3');
write_json(fullfile(outDir,'report.json'),report);
fprintf('BOX_B_CAPACITY_COMPLETE current=%d/%d cap210=%d/%d proposals=%d noLU=1\n', ...
    report.current.admittedX,report.current.admittedY,report.counterfactual.admittedX, ...
    report.counterfactual.admittedY,report.counterfactual.proposalCount);
end
function result=metadata_pair(a,b,checkpointFile)
% maybeCheckpoint updates this output cursor path only after serializing CP.
names=fieldnames(a);different={};assert(isequal(names,fieldnames(b)));
for k=1:numel(names)
    if ~isequaln(a.(names{k}),b.(names{k})),different{end+1}=names{k};end %#ok<AGROW>
end
result=struct('differentFields',{different},'exception','latestCheckpointFile', ...
    'reason','maybeCheckpoint updates the metadata path after immutable serialization', ...
    'resultLatestCheckpointFile',a.latestCheckpointFile, ...
    'checkpointPriorLatestCheckpointFile',b.latestCheckpointFile, ...
    'passed',isequal(different,{'latestCheckpointFile'})&& ...
    strcmp(a.latestCheckpointFile,checkpointFile)&&isfile(b.latestCheckpointFile)&& ...
    isequaln(rmfield(a,'latestCheckpointFile'),rmfield(b,'latestCheckpointFile')));
end
function index=eligible_levels(family,source)
factor=vertcat(family.members.cellFactors);
index=find(all(factor>=source.cellFactors,2)&[family.members.resourceAdmitted]'&[family.members.qualityPassed]');
end
function row=plan_summary(candidates,a)
row=struct('status',a.status,'admittedX',nnz([a.xTrials.admissible]), ...
    'admittedY',nnz([a.yTrials.admissible]),'totalAxisPairs',a.admittedAxisPairCount, ...
    'proposalCount',numel(candidates),'qualifiedKeepCount',nnz([candidates.unchanged]), ...
    'selectedPairs',a.selectedPairs,'candidateTransferQualified',false, ...
    'xRejections',rejection_counts(a.xTrials),'yRejections',rejection_counts(a.yTrials));
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct('reason',{},'count',{});
for k=1:numel(labels)
    rows(k)=struct('reason',labels{k},'count',nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials)));
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
