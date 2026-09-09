function report=mesh_audit_second_growth(latestFile,priorFile,growthFile,outDir)
%MESH_AUDIT_SECOND_GROWTH Strict persisted evidence plus current pure axes.
% No native field transfer is replayed and no operator/LU/flow is built.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
issues=checkcode(which(mfilename),'-id');assert(isempty(issues),jsonencode(issues));
registration=struct('kind','second_direction_growth_native_and_axis_capacity_audit_v1', ...
    'latestCheckpoint',latestFile,'preGrowthCheckpoint',priorFile,'postGrowthCheckpoint',growthFile, ...
    'independentHistoricPeakJumpRecomputed',false,'historicTransferReplayed',false, ...
    'remainingGeometryIsFrozenObservationOnly',true,'noLU',true,'noPDE',true);
write_json(fullfile(outDir,'registration.json'),registration);
profile clear;profile on;
cp=ipm.output.readCheckpoint(latestFile);s=cp.payload.state;h=cp.payload.log.history;m=s.runMetadata.autonomousMesh;
pre=ipm.output.readCheckpoint(priorFile);old=pre.payload.state;
post=ipm.output.readCheckpoint(growthFile);z=post.payload.state;
p=s.config.remesh.autonomousMesh;assert(p.version==2&&p.nodeFamily.maximumTotalNodes==210000);
assert(all(ipm.output.trustedMask(h,s.config))&&all(ipm.output.trustedMask(pre.payload.log.history,old.config))&& ...
    all(ipm.output.trustedMask(post.payload.log.history,z.config)));
assert(isequaln(old.config,s.config)&&isequaln(z.config,s.config)&&strcmp(old.runMetadata.caseId,s.runMetadata.caseId)&& ...
    strcmp(z.runMetadata.caseId,s.runMetadata.caseId)&&isequaln(old.runMetadata.caseMetadata,s.runMetadata.caseMetadata)&& ...
    isequaln(z.runMetadata.caseMetadata,s.runMetadata.caseMetadata)&& ...
    isequaln(old.runMetadata.gaugeContract,s.runMetadata.gaugeContract)&& ...
    isequaln(z.runMetadata.gaugeContract,s.runMetadata.gaugeContract));
assert(isequaln(old.runMetadata.autonomousMesh.initialization,m.initialization)&& ...
    isequaln(z.runMetadata.autonomousMesh.initialization,m.initialization)&& ...
    isequaln(old.runMetadata.autonomousMesh.referenceFamily,m.referenceFamily)&& ...
    isequaln(z.runMetadata.autonomousMesh.referenceFamily,m.referenceFamily));
prefixOld=history_prefix(pre.payload.log.history,h);prefixPost=history_prefix(post.payload.log.history,h);
assert(prefixOld&&prefixPost);
assert(isequaln(old.runMetadata.autonomousMesh.transactions,m.transactions(1:old.remeshCount))&& ...
    isequaln(z.runMetadata.autonomousMesh.transactions,m.transactions(1:z.remeshCount)));
referenceNames={'originIndex','pinIndex'};
assert(isequaln(rmfield(old.rescaling,referenceNames),rmfield(s.rescaling,referenceNames))&& ...
    isequaln(rmfield(z.rescaling,referenceNames),rmfield(s.rescaling,referenceNames)));
growthIndices=find([m.transactions.sourceLevelId]~=[m.transactions.targetLevelId]);
assert(numel(growthIndices)==2);gIndex=growthIndices(2);g=m.transactions(gIndex);d=g.controllerDecision;
write_json(fullfile(outDir,'second_growth_raw.json'),g);
assert(g.sourceLevelId==2&&g.targetLevelId==4&& ...
    isequal(g.sourceNodeCount,[641,161])&&isequal(g.targetNodeCount,[641,321]));
assert(old.step<g.sourceStep&&z.step>=g.sourceStep&&s.step>=z.step&& ...
    z.remeshCount>=g.sourceRemeshCount+1&&z.scale.canonicalTime>=g.sourceCanonicalTime&& ...
    z.scale.physicalTime>=g.sourcePhysicalTime&&z.normalizedTime>=d.sourceNormalizedTime);
featurePost=features(z);featureCurrent=features(s);
assert(isequal(featureCurrent.actualCoreCells,m.window.coreCells(end,:)));
sameStepPost=z.step==g.sourceStep;
if sameStepPost
    assert(isequal(featurePost.actualCoreCells,g.coreCells)&&isequal(featurePost.leftFrontCells,g.leftFrontCells));
    assert(isequaln(g.xQuality,ipm.mesh.quality(z.x,7,ipm.mesh.quadrature(z.x)))&& ...
        isequaln(g.yQuality,ipm.mesh.quality(z.y,7,ipm.mesh.quadrature(z.y))));
end
j=find(h.common.acceptedStep==g.sourceStep);historyHasTransactionRow=isscalar(j);
if historyHasTransactionRow
    assert(isequal(h.common.canonicalTau(j),g.sourceCanonicalTime)&& ...
        isequal(h.common.physicalTime(j),g.sourcePhysicalTime)&& ...
        isequal([h.mesh.coreGridPoints(j),h.mesh.verticalCoreGridPoints(j)],g.coreCells)&& ...
        isequal([h.common.nodeCountX(j),h.common.nodeCountY(j)],g.targetNodeCount)&&h.common.meshLevelId(j)==4);
end
before=find(h.common.acceptedStep<g.sourceStep,1,'last');after=find(h.common.acceptedStep>g.sourceStep,1,'first');
neighbors=struct('step',h.common.acceptedStep([before,after]), ...
    'canonicalTime',h.common.canonicalTau([before,after]),'physicalTime',h.common.physicalTime([before,after]), ...
    'levelId',h.common.meshLevelId([before,after]));
% Recompute the trigger algebra from saved source observations. The raw
% pre-remesh per-step window is not persisted after committing the new mesh.
expectedPrediction=d.coreCells.*exp(-p.maximumReviewInterval*d.decayEstimate);
assert(isequal(expectedPrediction,d.predictedCoreCells));
causal=true;
for k=1:numel(d.decayEvidence)
    e=d.decayEvidence(k);
    causal=causal&&e.secantFromStep<=e.secantToStep&&e.secantToStep<=g.sourceStep&& ...
        e.secantFromTime<=e.secantToTime&&e.secantToTime<=g.sourceCanonicalTime&& ...
        e.secantFromTime>=g.sourceCanonicalTime-p.trendWindow;
    if d.trendAvailable
        assert(isequal(d.decayEstimate(k),max([0;e.fitDecay;e.maximumSecant])));
    end
end
assert(causal&&d.requested&&~d.initial&&isempty(d.stopReason));
if any(d.coreCells<p.regridCoreTrigger)
    expectedReason='current_core_trigger';
elseif d.trendAvailable&&any(d.predictedCoreCells<p.predictedCoreBuffer&d.coreCells<p.targetCoreCells)
    expectedReason='forecast_core_trigger';
else
    expectedReason='safety_buffer';
end
assert(strcmp(expectedReason,d.reason));
assert(numel(g.attemptSummary)==g.candidateIndex&&g.attemptSummary(end).passed&& ...
    ~any([g.attemptSummary(1:end-1).passed]));
jumps=[m.transactions.relativePeakJump];assert(isequal(sum(jumps),m.cumulativeAbsolutePeakJump));
assert(isequal(sum(jumps(1:gIndex-1)),g.priorCumulativeAbsolutePeakJump)&& ...
    isequal(g.priorCumulativeAbsolutePeakJump+g.relativePeakJump,g.cumulativeAbsolutePeakJump));
family=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX,m.initialization.selectedBaseY, ...
    s.config.scaling.transportAnchorX,p);assert(isequaln(family,m.referenceFamily));
member=family.members(m.currentLevelId);factors=vertcat(family.members.cellFactors);
eligible=find(all(factors>=member.cellFactors,2)&[family.members.resourceAdmitted]'&[family.members.qualityPassed]');
assert(m.currentLevelId==4&&isequal(eligible,4));
Dx=ipm.mesh.fdMatrix(s.x,1,7);view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
[candidates,axisReport]=ipm.remesh.plannedAxisPairs(view,struct('x',member.baseX,'y',member.baseY), ...
    s.config.scaling.transportAnchorX,p);
qx=ipm.mesh.quality(s.x,7,ipm.mesh.quadrature(s.x));qy=ipm.mesh.quality(s.y,7,ipm.mesh.quadrature(s.y));
remaining=struct('levelId',m.currentLevelId,'actualNodeCount',[numel(s.x),numel(s.y)], ...
    'nodeProduct',numel(s.x)*numel(s.y),'cap',p.nodeFamily.maximumTotalNodes, ...
    'unusedNodeCap',p.nodeFamily.maximumTotalNodes-numel(s.x)*numel(s.y), ...
    'eligibleRegisteredLevels',eligible,'remainingGrowthTransitions',p.nodeFamily.maximumAcceptedGrowthTransitions-numel(growthIndices), ...
    'remainingRemeshCountBudget',s.config.remesh.maxRemeshes-s.remeshCount, ...
    'remainingAbsolutePeakBudget',p.maximumCumulativeAbsolutePeakJump-m.cumulativeAbsolutePeakJump, ...
    'plannerStatus',axisReport.status,'admittedX',nnz([axisReport.xTrials.admissible]), ...
    'admittedY',nnz([axisReport.yTrials.admissible]),'xRejections',rejection_counts(axisReport.xTrials), ...
    'yRejections',rejection_counts(axisReport.yTrials),'proposalCount',numel(candidates), ...
    'qualifiedKeepCount',nnz([candidates.unchanged]),'actualMigrationProposals',nnz(~[candidates.unchanged]), ...
    'predictedCoreFront',vertcat(candidates.predictedCells),'xQuality',qx,'yQuality',qy, ...
    'noLargerRegisteredMember',true,'candidateTransferQualified',false,'futureCapacityHorizonKnown',false);
report=registration;report.nativeSignaturesValidated=true;report.allThreeHistoriesTrusted=true;
report.exactHistoryAndLedgerPrefixes=true;report.originalInitialAxesReferencesAndLineagePreserved=true;
report.step=s.step;report.canonicalTime=s.scale.canonicalTime;report.physicalTime=s.scale.physicalTime;
report.normalizedTime=s.normalizedTime;report.remeshCount=s.remeshCount;report.coreCells=featureCurrent.actualCoreCells;
report.safety=m.window.safety(end);report.currentFeatures=featureCurrent;report.growthTransactions=m.transactions(growthIndices);
report.secondGrowthTransactionIndex=gIndex;report.secondGrowth=g;report.secondGrowthDestinationGeometryRecomputed=sameStepPost;
report.secondGrowthHistoryClocksCoreAndNodeCountExact=historyHasTransactionRow;
report.historyHasExactTransactionStep=historyHasTransactionRow;report.neighboringHistoryObservations=neighbors;
report.postGrowthCheckpointStep=z.step;report.postGrowthCheckpointIsExactTransactionStep=sameStepPost;
report.secondGrowthPostFeature=featurePost;report.expectedReason=expectedReason;report.decisionCausalEndpointsPassed=causal;
report.fullPreGrowthTelemetryWindowAvailable=false;report.completeHistoricalAxisRejectionReportAvailable=false;
report.transactionPeakJumpEvidence='native signed actual transaction audit; no same-step pre-transfer rho was saved';
report.cumulativeAbsolutePeakJump=m.cumulativeAbsolutePeakJump;report.remaining=remaining;
report.currentPhysicalGradient=h.common.physicalGradInf(end);report.currentQuadraticPeak=h.common.physicalQuadraticPeak(end);
report.currentWindowRecords=numel(m.window.time);report.currentWindowTimeSpan=m.window.time(end)-m.window.time(1);
report.latestDecision=m.lastDecision;report.completed=true;
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})));
end
save(fullfile(outDir,'report.mat'),'report','candidates','axisReport','featurePost','featureCurrent','profileInfo','-v7.3');
write_json(fullfile(outDir,'report.json'),report);
fprintf('SECOND_GROWTH_AUDIT_COMPLETE sourceStep=%d currentStep=%d level=4 core=%.8g/%.8g x=%d y=%d proposals=%d newLU=0\n', ...
    g.sourceStep,s.step,report.coreCells,remaining.admittedX,remaining.admittedY,numel(candidates));
end

function f=features(s)
Dx=ipm.mesh.fdMatrix(s.x,1,7);v=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
f=ipm.diagnostics.meshFeatureIntervals(v);
end
function yes=history_prefix(a,b)
yes=isequal(fieldnames(a),fieldnames(b));
for group=fieldnames(a)'
    an=a.(group{1});bn=b.(group{1});yes=yes&&isequal(fieldnames(an),fieldnames(bn));
    for field=fieldnames(an)'
        v=an.(field{1});w=bn.(field{1});n=size(v,1);
        yes=yes&&size(w,1)>=n&&isequaln(v,w(1:n,:));
    end
end
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct('reason',{},'count',{});
for k=1:numel(labels),rows(k)=struct('reason',labels{k},'count',nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials)));end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
