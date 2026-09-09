function report=mesh_test_autonomous_decision(checkpointFile,outputDirectory)
%MESH_TEST_AUTONOMOUS_DECISION Policy fixtures and causal historical replay.
assert(maxNumCompThreads==10 && ~isfolder(outputDirectory));mkdir(outputDirectory);
policy=struct('timeUnit','native_canonical');
base=struct('time',(0:.01:.05)','step',(0:5)','remeshCount',zeros(6,1), ...
    'coreCells',repmat([35,37],6,1),'safety',.25*ones(6,1), ...
    'trusted',true(6,1),'nextStableStep',.01,'cumulativeAbsolutePeakJump',0, ...
    'timeUnit','native_canonical');
records=struct();
d=mesh_autonomous_decision(base,policy);assert(strcmp(d.action,'advance') && d.maximumReviewInterval==.2);records.stable=d;
s=base;s.coreCells=[40*exp(-4*s.time),42*exp(-4*s.time)];
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.reason,'forecast_core_trigger'));records.shrinking=d;
s=base;s.nextStableStep=1e-5;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'advance') && ...
    d.maximumReviewInterval<=512*s.nextStableStep && d.maximumAdvanceSteps<512);records.slowCFL=d;
s=base;s.coreCells=[40*exp(-4*s.time),42*exp(-4*s.time)];s.nextStableStep=1e-5;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'advance') && ...
    all(d.predictedCoreCells>22));records.forecastUsesActualReviewHorizon=d;
s=base;s.time(end)=s.time(end-1)+1e-12;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'advance') && ...
    d.maximumReviewInterval==.2);records.endpointClippingDoesNotLimitFuture=d;
s=base;s.remeshCount(end)=1;
d=mesh_autonomous_decision(s,policy);assert(d.maximumAdvanceSteps==1 && ~d.trendAvailable);records.newEpoch=d;
s=base;s.coreCells(end,1)=25;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.reason,'current_core_trigger'));records.currentTrigger=d;
s=base;s.trusted(2)=false;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'preserve_and_reject'));records.untrusted=d;
s=base;s.cumulativeAbsolutePeakJump=.020001;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'preserve_and_reject'));records.budgetExceeded=d;
s=base;s.coreCells(2,1)=13;
d=mesh_autonomous_decision(s,policy);assert(strcmp(d.action,'preserve_and_reject'));records.historicalFailure=d;
s=base;s.timeUnit='parent_equivalent_canonical';clockRejected=false;
try
    mesh_autonomous_decision(s,policy);
catch ex
    clockRejected=strcmp(ex.identifier,'ipm:AutonomousMeshClock');
end
assert(clockRejected);
% A true zero-time state requires one-step observation, not a late bootstrap.
s=base;names={'time','step','remeshCount','coreCells','safety','trusted'};
for k=1:numel(names),s.(names{k})=s.(names{k})(1,:);end
d=mesh_autonomous_decision(s,policy);assert(d.currentStep==0 && d.maximumAdvanceSteps==1);records.zeroTime=d;
cp=ipm.output.readCheckpoint(checkpointFile);h=cp.payload.log.history;state=cp.payload.state;
assert(all(ipm.output.trustedMask(h,state.config)));
factor=1;timeUnit='native_canonical';
if isfield(state.runMetadata.caseMetadata,'latePhysicalBoxBranch')
    lineage=state.runMetadata.caseMetadata.latePhysicalBoxBranch;
    factor=lineage.canonicalCovarianceFactor;timeUnit='parent_equivalent_canonical';
end
% No absolute physical clock is constructed in this scheduling test. The
% positive time-unit conversion is explicit and does not join IVP histories.
t=factor*h.common.canonicalTau(:);steps=h.common.acceptedStep(:);
epoch=h.mesh.remeshCount(:);cores=[h.mesh.coreGridPoints(:),h.mesh.verticalCoreGridPoints(:)];
safety=h.mesh.safetyFactor(:);trusted=logical(ipm.output.trustedMask(h,state.config));
assert(all(diff(t)>0) && all(diff(steps)>0));
replay=cell(numel(t),1);jump=0;nativeRecords=[];
if isfield(state.runMetadata,'gridLabRegrids'),nativeRecords=state.runMetadata.gridLabRegrids;end
for j=1:numel(t)
    if epoch(j)>0
        assert(numel(nativeRecords)>=epoch(j));
        jump=sum(abs([nativeRecords(1:epoch(j)).rhoXMaximumRelativeChange]));
    end
    if j==1,dt=1e-5;else,dt=(t(j)-t(j-1))/(steps(j)-steps(j-1));end
    sample=struct('time',t(1:j),'step',steps(1:j),'remeshCount',epoch(1:j), ...
        'coreCells',cores(1:j,:),'safety',safety(1:j),'trusted',trusted(1:j), ...
        'nextStableStep',dt,'cumulativeAbsolutePeakJump',jump,'timeUnit',timeUnit);
    r=mesh_autonomous_decision(sample,struct('timeUnit',timeUnit));
    replay{j}=struct('step',steps(j),'time',t(j),'remeshCount',epoch(j), ...
        'core',cores(j,:),'action',r.action,'reason',r.reason, ...
        'maximumAdvanceSteps',r.maximumAdvanceSteps,'maximumReviewInterval',r.maximumReviewInterval);
end
replay=[replay{:}];
assert(~any(strcmp({replay.action},'preserve_and_reject')));
epochs=unique(epoch);epochAudit=cell(numel(epochs),1);
for j=1:numel(epochs)
    indexes=find(epoch==epochs(j));triggers=indexes(strcmp({replay(indexes).action},'regrid'));
    crossing=indexes(any(cores(indexes,:)<22,2));
    detectedBeforeBufferLoss=isempty(crossing) || (~isempty(triggers) && triggers(1)<=crossing(1));
    assert(detectedBeforeBufferLoss);
    epochAudit{j}=struct('remeshCount',epochs(j),'records',numel(indexes), ...
        'firstTriggerStep',first_value(steps,triggers),'firstBelow22Step',first_value(steps,crossing), ...
        'detectedBeforeBufferLoss',detectedBeforeBufferLoss);
end
report=struct('kind','research_autonomous_mesh_decision_test_v1','allPassed',true, ...
    'checkpointFile',checkpointFile,'nativeSignatureValidated',true, ...
    'clockMismatchRejected',clockRejected,'fixtures',records,'epochAudit',[epochAudit{:}], ...
    'replayedRecords',numel(replay),'historicalReplayChangesTrajectory',false, ...
    'nextStableStepInReplay','past accepted average only; causal estimate, not a new CFL solve', ...
    'initialReplayStepEstimate',1e-5,'controllerIntegrated',false,'noLU',true,'noPDE',true);
save(fullfile(outputDirectory,'report.mat'),'report','replay','-v7.3');
fid=fopen(fullfile(outputDirectory,'report.json'),'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('AUTONOMOUS_MESH_POLICY_PASS records=%d epochs=%d noLU=1\n',numel(replay),numel(epochs));
end
function value=first_value(values,indexes)
value=NaN;if ~isempty(indexes),value=values(indexes(1));end
end
