function report=ipm_perflab_tau8_capacity(checkpointFile,resultFile,outDir)
%IPM_PERFLAB_TAU8_CAPACITY Actual endpoint and bounded pure geometry only.
% Future factor-three references are research proposals, never native v2
% policy, accepted grids, reconstructed fields or certified LU resources.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
registration=struct('kind','tau8_endpoint_fixed_and_proposed_reference_capacity_v1', ...
    'checkpointFile',checkpointFile,'resultFile',resultFile, ...
    'actualSourceStep',4545,'actualSourceCanonicalTime',8, ...
    'currentRegisteredFactors',[2,2],'proposedAdditionalFactors',[3,2;2,3], ...
    'proposedNodeCap',310000,'proposedMaximumGrowthTransitions',3, ...
    'proposedReferenceRule','direct_initial_selected_base_index_pchip_integer_v1', ...
    'preserveOriginalSelectedKnots',true,'preserveAllIntermediateReferenceKnots',false, ...
    'allOtherSearchAndQualityControlsUnchanged',true,'noFutureFieldDeformation',true, ...
    'noLU',true,'noPDE',true,'noTransfer',true, ...
    'newFamilyNativeContractImplemented',false,'newFamilyMemoryQualified',false);
write_json(fullfile(outDir,'registration.json'),registration);
issues=checkcode(which(mfilename),'-id');write_json(fullfile(outDir,'checkcode.json'),issues);
assert(isempty(issues),jsonencode(issues));
profile clear;profile on;timer=tic;
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;h=cp.payload.log.history;
m=s.runMetadata.autonomousMesh;p=s.config.remesh.autonomousMesh;
loaded=load(resultFile,'result');r=ipm.output.validate(loaded.result);clear loaded;
assert(s.step==4545&&s.scale.canonicalTime==8&&m.currentLevelId==4&& ...
    p.version==2&&p.nodeFamily.maximumTotalNodes==210000&& ...
    strcmp(r.state.stopReason,'final_time'));
assert(isequaln(r.state.rho,s.rho)&&isequaln(r.grid.x,s.x)&&isequaln(r.grid.y,s.y)&& ...
    isequaln(r.config,s.config)&& ...
    isequaln(r.history,h)&&r.state.steps==s.step&&r.state.canonicalTime==s.scale.canonicalTime&& ...
    r.state.physicalTime==s.scale.physicalTime&&r.state.normalizedTime==s.normalizedTime&& ...
    r.grid.remeshCount==s.remeshCount);
metadataPair=metadata_pair(r.metadata,s.runMetadata,checkpointFile);
write_json(fullfile(outDir,'metadata_pair.json'),metadataPair);assert(metadataPair.passed);
assert(all(ipm.output.trustedMask(r))&&~s.config.output.storeSnapshots&& ...
    isempty(r.snapshots.rho)&&isempty(r.snapshots.x)&&isempty(r.snapshots.y));
assert(h.common.acceptedStep(1)==0&&h.common.canonicalTau(1)==0&& ...
    h.common.physicalTime(1)==0&&m.initialization.originalStep==0&& ...
    m.initialization.originalCanonicalTime==0);
anchor=s.config.scaling.transportAnchorX;
family=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX, ...
    m.initialization.selectedBaseY,anchor,p);
assert(isequaln(family,m.referenceFamily));member=family.members(m.currentLevelId);
assert(isequaln(s.baseX,member.baseX)&&isequaln(s.baseY,member.baseY));
counts=vertcat(family.members.nodeCount);factors=vertcat(family.members.cellFactors);
eligible=find(all(factors>=member.cellFactors,2)&[family.members.resourceAdmitted]'& ...
    [family.members.qualityPassed]');assert(isequal(eligible,4));
growth=nnz([m.transactions.sourceLevelId]~=[m.transactions.targetLevelId]);assert(growth==2);
Dx=ipm.mesh.fdMatrix(s.x,1,7);Omega=s.rho*Dx';assert(isequaln(Omega,r.state.omega));
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',Omega,'trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);
assert(isequaln(feature.actualCoreCells,m.window.coreCells(end,:)));
[qx,qxPassed,qxMargin]=quality(s.x,p.qualityLimits);
[qy,qyPassed,qyMargin]=quality(s.y,p.qualityLimits);assert(qxPassed&&qyPassed);
planTimer=tic;
[candidates,axisReport]=ipm.remesh.plannedAxisPairs(view, ...
    struct('x',member.baseX,'y',member.baseY),anchor,p);
currentPlanSeconds=toc(planTimer);
current=plan_summary(candidates,axisReport);current.planningSeconds=currentPlanSeconds;
current.candidateTransferQualified=false;current.currentAxisQualityX=qx;current.currentAxisQualityY=qy;
current.currentQualityMarginX=qxMargin;current.currentQualityMarginY=qyMargin;
current.eligibleNativeLevels=eligible;current.registeredNodeCounts=counts;
current.currentNodeCount=member.nodeCount;current.unusedNodeCap=p.nodeFamily.maximumTotalNodes-prod(member.nodeCount);
current.remainingGrowthTransitions=p.nodeFamily.maximumAcceptedGrowthTransitions-growth;
current.remainingRemeshCountBudget=s.config.remesh.maxRemeshes-s.remeshCount;
current.remainingAbsolutePeakBudget=p.maximumCumulativeAbsolutePeakJump-m.cumulativeAbsolutePeakJump;
current.latestDecision=m.lastDecision;current.futureCapacityHorizonKnown=false;
current.currentTriggerMargins=feature.actualCoreCells-p.regridCoreTrigger;
current.currentEndpointMargins=feature.actualCoreCells-p.endpointCoreFloor;
current.currentHistoryMargins=feature.actualCoreCells-p.historyCoreFloor;
current.safetyMargin=p.maximumSafety-m.window.safety(end);

% Direct initial-root interpolation. Factor two must reproduce the actual
% immutable family bitwise before considering the separate factor-three idea.
rootX=family.rootX;rootY=family.rootY;
x2=refine_x(rootX,2);y2=refine_y(rootY,2);
assert(isequaln(x2,member.baseX)&&isequaln(y2,member.baseY));
x3=refine_x(rootX,3);y3=refine_y(rootY,3);
[qx3,passedX3,marginX3]=quality(x3,p.qualityLimits);
[qy3,passedY3,marginY3]=quality(y3,p.qualityLimits);
references={struct('x',x3,'y',y2),struct('x',x2,'y',y3)};
assert(isequaln(x3(1:3:end),rootX)&&isequaln(y3(1:3:end),rootY)&& ...
    isequal(x3,-fliplr(x3))&&any(x3==anchor)&&any(x3==-anchor));
proposalPolicy=p;proposalPolicy.nodeFamily.maximumTotalNodes=registration.proposedNodeCap;
assert(isequaln(rmfield(proposalPolicy.nodeFamily,'maximumTotalNodes'), ...
    rmfield(p.nodeFamily,'maximumTotalNodes')));
proposed=struct([]);proposedCandidates=cell(1,2);proposedAxisReports=cell(1,2);
for k=1:2
    reference=references{k};targetCount=[numel(reference.x),numel(reference.y)];
    if k==1,referencePassed=passedX3;else,referencePassed=passedY3;end
    row=struct('cellFactors',registration.proposedAdditionalFactors(k,:), ...
        'nodeCount',targetCount,'nodeProduct',prod(targetCount), ...
        'originalSelectedKnotsPreserved',true,'sameBox',true, ...
        'referenceQualityPassed',referencePassed, ...
        'withinProposedNodeCap',prod(targetCount)<=registration.proposedNodeCap, ...
        'nativeV2Eligible',false,'nativeTransferQualified',false, ...
        'LUResourceQualified',false,'plan',struct(),'planningSeconds',0);
    if referencePassed&&row.withinProposedNodeCap
        planTimer=tic;
        [cc,aa]=ipm.remesh.plannedAxisPairs(view,reference,anchor,proposalPolicy);
        row.planningSeconds=toc(planTimer);row.plan=plan_summary(cc,aa);
        proposedCandidates{k}=cc;proposedAxisReports{k}=aa;
    end
    if isempty(proposed),proposed=row;else,proposed(end+1)=row;end %#ok<AGROW>
end
referenceEvidence=struct('rootX',rootX,'rootY',rootY,'factorTwoExact',true, ...
    'x3',x3,'y3',y3,'x3Quality',qx3,'y3Quality',qy3, ...
    'x3Passed',passedX3,'y3Passed',passedY3,'x3Margin',marginX3,'y3Margin',marginY3);
report=registration;report.sourceNativeSignatureValidated=true;
report.resultNativePayloadPairExact=true;report.completeHistoryTrusted=true;
report.metadataPair=metadataPair;
report.originalZeroTimeProvenance=true;report.initialReferencesRegeneratedExactly=true;
report.currentStep=s.step;report.canonicalTime=s.scale.canonicalTime;
report.physicalTime=s.scale.physicalTime;report.remeshCount=s.remeshCount;
report.feature=feature;report.safety=m.window.safety(end);
report.cumulativeAbsolutePeakJump=m.cumulativeAbsolutePeakJump;
report.physicalGradInf=h.common.physicalGradInf(end);report.physicalRhoXInf=h.common.physicalRhoXInf(end);
report.physicalQuadraticPeak=h.common.physicalQuadraticPeak(end);
report.physicalBox=[r.grid.physicalX([1,end]),r.grid.physicalY([1,end])'];
report.current=current;report.proposed=proposed;report.savedLastTransaction=m.transactions(end);
report.currentWindowRecords=numel(m.window.time);
report.currentWindowSpan=m.window.time(end)-m.window.time(1);
report.currentFullWindowSaved=true;report.historicPeakJumpIndependentlyRecomputed=false;
report.inputPayloadUnchanged=isequaln(s,cp.payload.state);assert(report.inputPayloadUnchanged);
report.wallSeconds=toc(timer);report.completed=true;
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.output.restoreCheckpoint','ipm.evolve.flow', ...
        'ipm.field.velocity','ipm.field.poisson','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})),'ipm:ResearchNoLU','Forbidden numerical call.');
end
save(fullfile(outDir,'report.mat'),'report','candidates','axisReport','proposedCandidates', ...
    'proposedAxisReports','referenceEvidence','profileInfo','-v7.3');
write_json(fullfile(outDir,'report.json'),report);
fprintf('TAU8_CAPACITY_COMPLETE step=%d tau=8 x=%d y=%d proposals=%d noLU=1\n', ...
    s.step,current.admittedX,current.admittedY,current.proposalCount);
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
function x=refine_x(root,factor)
positive=root(root>=0);v=pchip(linspace(0,1,numel(positive)),positive, ...
    linspace(0,1,factor*(numel(positive)-1)+1));
v(1:factor:end)=positive;v(1)=0;x=[-fliplr(v(2:end)),v];
end
function y=refine_y(root,factor)
y=pchip(linspace(0,1,numel(root)),root,linspace(0,1,factor*(numel(root)-1)+1))';
y(1:factor:end)=root;y(1)=0;
end
function [q,passed,margin]=quality(axis,limits)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
v=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio, ...
    q.maximumQuadratureWeightToControlWidthRatio];
margin=[limits.maxAdjacentCellRatio-v(1),limits.maxLogSpacingCurvature-v(2), ...
    v(3)-limits.minStencilRcond,v(4)-limits.minQuadratureWeightRatio, ...
    v(5)-limits.minWeightToControlWidth,limits.maxWeightToControlWidth-v(6)];
passed=all(isfinite(v))&&all(margin>=0)&&q.quadratureWeightsStrictlyPositive;
end
function row=plan_summary(candidates,a)
row=struct('status',a.status,'admittedX',nnz([a.xTrials.admissible]), ...
    'admittedY',nnz([a.yTrials.admissible]),'totalAxisPairs',a.admittedAxisPairCount, ...
    'proposalCount',numel(candidates),'qualifiedKeepCount',nnz([candidates.unchanged]), ...
    'migrationProposalCount',nnz(~[candidates.unchanged]), ...
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
