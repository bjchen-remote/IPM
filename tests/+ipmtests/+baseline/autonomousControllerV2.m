function report = autonomousControllerV2(oldCheckpointFile,v1ResultFile)
%IPMTESTS.BASELINE.AUTONOMOUSCONTROLLERV2 No-LU variable-node persistence.
if nargin<1,oldCheckpointFile='';end
if nargin<2,v1ResultFile='';end
[m,c,v,cp]=ipmtests.support.autonomousControllerV2Fixture();
id='ipm:AutonomousMeshController';labels={};
assert(ipm.remesh.validateController(m,c,v).checked);
assert(isequaln(ipm.output.readCheckpoint(cp),cp));
baseLess=rmfield(v,{'baseX','baseY'});assert(ipm.remesh.validateController(m,c,baseLess).checked);
withoutSnapshots=v;withoutSnapshots.snapshots=struct('rho',{{}},'x',{{}},'y',{{}});
off=c;off.output.storeSnapshots=false;assert(ipm.remesh.validateController(m,off,withoutSnapshots).checked);
[sameM,sameV]=same_level_fixture(m,v);
assert(ipm.remesh.validateController(sameM,c,sameV).checked);
forecast=forecast_fixture(m,c,2);assert(ipm.remesh.validateController(forecast,c,v).checked);
underflow=forecast_fixture(m,c,1e6);assert(ipm.remesh.validateController(underflow,c,v).checked);
[postM,postV]=history_state_fixture(m,v,4);assert(ipm.remesh.validateController(postM,c,postV).checked);
[initialM,initialV]=history_state_fixture(m,v,1);assert(ipm.remesh.validateController(initialM,c,initialV).checked);
bad=m;bad.autonomousMesh.version=1;reject(bad,c,v,id);labels{end+1}='memory version';
bad=m;bad.autonomousMesh.currentLevelId=3;reject(bad,c,v,id);labels{end+1}='terminal wrong member';
bad=m;bad.autonomousMesh.referenceFamily.members(4).baseX(2)=bad.autonomousMesh.referenceFamily.members(4).baseX(2)+eps;
reject(bad,c,v,id);labels{end+1}='family changed inserted knot';
bad=m;bad.autonomousMesh.referenceFamily.members(2).resourceAdmitted=1;reject(bad,c,v,id);labels{end+1}='family bool class';
bad=m;bad.autonomousMesh.referenceFamily.rootX(2)=bad.autonomousMesh.referenceFamily.rootX(2)+eps;
reject(bad,c,v,id);labels{end+1}='family root drift';
badV=v;badV.nodeCount=int32(v.nodeCount);reject(m,c,badV,id);labels{end+1}='projected node count class';
badV=v;badV.levelId=int32(4);reject(m,c,badV,id);labels{end+1}='projected level class';
badV=v;badV.baseX=m.autonomousMesh.initialization.selectedBaseX;reject(m,c,badV,id);labels{end+1}='wrong current base';
bad=m;bad.autonomousMesh.transactions(2).sourceLevelId=1;reject(bad,c,v,id);labels{end+1}='broken level prefix';
bad=m;member=m.autonomousMesh.referenceFamily.members(3);
bad.autonomousMesh.transactions(2).targetLevelId=3;
bad.autonomousMesh.transactions(2).targetNodeCount=member.nodeCount;
bad.autonomousMesh.transactions(2).targetCellFactors=member.cellFactors;
reject(bad,c,v,id);labels{end+1}='incomparable level';
bad=m;bad.autonomousMesh.transactions(2).targetNodeCount=[129,65];reject(bad,c,v,id);labels{end+1}='wrong target dimensions';
bad=m;bad.autonomousMesh.transactions(2).sourceCellFactors=[1,1];reject(bad,c,v,id);labels{end+1}='wrong source factors';
names={'sameBox','registeredNodeCounts','matchingReferenceMembers','admittedLevelTransition'};
for k=1:numel(names)
    bad=m;bad.autonomousMesh.transactions(2).(names{k})=false;reject(bad,c,v,id);labels{end+1}=names{k}; %#ok<AGROW>
end
bad=m;bad.autonomousMesh.transactions(2).sameBoxAndNodeCount=true;reject(bad,c,v,id);labels{end+1}='false sameN claim';
bad=m;bad.autonomousMesh.window.nodeCount(1,1)=129;reject(bad,c,v,id);labels{end+1}='window stale dimensions';
bad=m;bad.autonomousMesh.window.levelId(1)=2;reject(bad,c,v,id);labels{end+1}='window stale level';
badV=v;badV.history.common.meshLevelId(2)=1;reject(m,c,badV,id);labels{end+1}='history wrong transaction epoch';
badV=v;badV.history.common.nodeCountX(1)=257;reject(m,c,badV,id);labels{end+1}='history initial dimensions';
badV=v;badV.history.common.nodeCountX=badV.history.common.nodeCountX';reject(m,c,badV,id);labels{end+1}='history dimension shape';
badV=v;badV.snapshots.x{1}=v.x;reject(m,c,badV,id);labels{end+1}='snapshot wrong epoch dimensions';
badV=v;badV.snapshots.rho{2}=zeros(1,1);reject(m,c,badV,id);labels{end+1}='snapshot field shape';
badV=v;badV.snapshots.x{2}(1)=badV.snapshots.x{2}(1)-1e-8;reject(m,c,badV,id);labels{end+1}='snapshot changed box';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.sourceStep=19;reject(bad,c,v,id);labels{end+1}='decision source step';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.requested=false;reject(bad,c,v,id);labels{end+1}='decision not requested';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.reason='safety_buffer';reject(bad,c,v,id);labels{end+1}='unsupported accepted trigger';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.predictedCoreCells=[24,32];reject(bad,c,v,id);labels{end+1}='changed exact forecast';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.reason='forecast_core_trigger';reject(bad,c,v,id);labels{end+1}='current trigger priority';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.coreCells=[32,32];
bad.autonomousMesh.transactions(2).controllerDecision.predictedCoreCells=[32,32];
reject(bad,c,v,id);labels{end+1}='no actual trigger';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.stopReason='autonomous_mesh_resolution_failure';reject(bad,c,v,id);labels{end+1}='stopped decision';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.sourceNormalizedTime=.1;reject(bad,c,v,id);labels{end+1}='decision clock';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.predictedCoreCells=[NaN,32];reject(bad,c,v,id);labels{end+1}='decision nonfinite';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.decayEstimate=[1,0];reject(bad,c,v,id);labels{end+1}='decay unsupported maximum';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.decayEvidence(1).secantToStep=21;reject(bad,c,v,id);labels{end+1}='decay future source';
bad=m;bad.autonomousMesh.transactions(2).controllerDecision.decayEvidence(1).axis=single(1);reject(bad,c,v,id);labels{end+1}='decay wrong class';
bad=m;bad.autonomousMesh.transactions(2).candidateIndex=13;reject(bad,c,v,id);labels{end+1}='attempt budget';
bad=m;bad.autonomousMesh.transactions(2).attemptSummary.passed=false;reject(bad,c,v,id);labels{end+1}='accepted attempt false';
bad=m;bad.autonomousMesh.transactions(2).attemptSummary.errorIdentifier='example:failure';reject(bad,c,v,id);labels{end+1}='accepted attempt error';
bad=m;bad.autonomousMesh.transactions(2).attemptSummary.auditReasons={'rejected'};reject(bad,c,v,id);labels{end+1}='accepted attempt reasons';
bad=m;bad.autonomousMesh.lastDecision.sourceStep=26;reject(bad,c,v,id);labels{end+1}='last decision future source';
bad=m;bad.autonomousMesh.lastDecision.sourceNodeCount=[129,65];reject(bad,c,v,id);labels{end+1}='last decision wrong source member';
bad=m;bad.autonomousMesh.lastDecision=m.autonomousMesh.transactions(end).controllerDecision;
reject(bad,c,v,id);labels{end+1}='stale last decision';
bad=postM;bad.autonomousMesh.lastDecision.sourceRemeshCount=postV.remeshCount;
bad.autonomousMesh.lastDecision.sourceLevelId=postV.levelId;
bad.autonomousMesh.lastDecision.sourceNodeCount=postV.nodeCount;
reject(bad,c,postV,id);labels{end+1}='postcommit decision relabelled target';
bad=postM;bad.autonomousMesh.lastDecision.coreCells=postV.coreCells;
reject(bad,c,postV,id);labels{end+1}='postcommit decision relabelled core';
limited=c;limited.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct('version',2,'nodeFamily',struct('maximumTotalNodes',17000)));
bad=m;bad.autonomousMesh.policy=limited.remesh.autonomousMesh;
bad.autonomousMesh.referenceFamily=ipm.remesh.referenceAxisFamily( ...
    bad.autonomousMesh.initialization.selectedBaseX,bad.autonomousMesh.initialization.selectedBaseY, ...
    limited.scaling.transportAnchorX,limited.remesh.autonomousMesh);
reject(bad,limited,v,id);labels{end+1}='resource inadmissible current';
% A valid signed synthetic variable-N CP exercises read, not numerical solve.
badCP=cp;badCP.payload.state.baseX(2)=badCP.payload.state.baseX(2)+eps;
expect(@()ipm.output.readCheckpoint(badCP),id);
badCP=cp;badCP.payload.state.config.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct());
expect(@()ipm.output.readCheckpoint(badCP),'ipm:CheckpointGrid');
badCP=cp;badCP.payload.state.runMetadata.autonomousMesh.referenceFamily.members(4).baseX(2)=0;
expect(@()ipm.output.readCheckpoint(badCP),id);expect(@()ipm.output.restoreCheckpoint(badCP),id);
actual=struct('v1ResultExact',false,'v2ConfigurationOnlyFixtureExact',false,'resultBridgeRejectedBeforeLU',false,'oldCheckpointExact',false);
if ~isempty(v1ResultFile)
    loaded=load(v1ResultFile,'result');original=loaded.result;
    assert(isequaln(ipm.output.validate(original),original));actual.v1ResultExact=true;
    fixture=initial_only_result_fixture(original);
    assert(isequaln(ipm.output.validate(fixture),fixture));actual.v2ConfigurationOnlyFixtureExact=true;
    expect(@()ipm.output.checkpointFromResult(fixture),'ipm:CheckpointResultVariableNodeFamily');
    actual.resultBridgeRejectedBeforeLU=true;
    fixture.history.common.nodeCountX(1)=fixture.history.common.nodeCountX(1)+2;
    expect(@()ipm.output.validate(fixture),id);
end
if ~isempty(oldCheckpointFile)
    assert(maxNumCompThreads==10);old=ipm.output.readCheckpoint(oldCheckpointFile);
    assert(isequaln(ipm.output.readCheckpoint(old),old));actual.oldCheckpointExact=true;
end
report=struct('allPassed',true,'syntheticFixtureNotPdeEvidence',true,'negativeLabels',{labels}, ...
    'negativeCount',numel(labels),'positiveProjectionCases',9,'variableSizeSignedFixtureRead',true, ...
    'malformedRestoreRejectedBeforeLU',true,'actual',actual,'noLU',true,'noPDE',true, ...
    'positiveRestoreNotExecuted',true,'recordPositivePathNotExecuted',true);
fprintf('AUTONOMOUS_CONTROLLER_V2_PASS negatives=%d noLU=1 noPDE=1\n',numel(labels));
end

function [m,v]=history_state_fixture(m,v,index)
groups=fieldnames(v.history);
for j=1:numel(groups)
    names=fieldnames(v.history.(groups{j}));
    for k=1:numel(names),v.history.(groups{j}).(names{k})=v.history.(groups{j}).(names{k})(1:index);end
end
v.step=v.history.common.acceptedStep(end);v.canonicalTime=v.history.common.canonicalTau(end);
v.normalizedTime=v.canonicalTime;v.physicalTime=v.history.common.physicalTime(end);
v.remeshCount=v.history.mesh.remeshCount(end);v.levelId=v.history.common.meshLevelId(end);
v.nodeCount=[v.history.common.nodeCountX(end),v.history.common.nodeCountY(end)];
v.coreCells=[v.history.mesh.coreGridPoints(end),v.history.mesh.verticalCoreGridPoints(end)];
v.safety=v.history.mesh.safetyFactor(end);
names=fieldnames(v.snapshots);for k=1:numel(names),v.snapshots.(names{k})=v.snapshots.(names{k})(1:index);end
v.x=v.snapshots.x{end};v.y=v.snapshots.y{end};member=m.autonomousMesh.referenceFamily.members(v.levelId);
v.baseX=member.baseX;v.baseY=member.baseY;
memory=m.autonomousMesh;memory.currentLevelId=v.levelId;memory.transactions=memory.transactions(1:v.remeshCount);
memory.cumulativeAbsolutePeakJump=sum([memory.transactions.relativePeakJump]);
memory.window=struct('time',v.canonicalTime,'step',v.step,'remeshCount',v.remeshCount, ...
    'coreCells',v.coreCells,'safety',v.safety,'nodeCount',v.nodeCount,'levelId',v.levelId);
if v.remeshCount>0
    memory.lastDecision=memory.transactions(end).controllerDecision;
else
    d=memory.lastDecision;d.sourceStep=0;d.sourceCanonicalTime=0;d.sourcePhysicalTime=0;d.sourceNormalizedTime=0;
    d.sourceRemeshCount=0;d.sourceLevelId=1;d.sourceNodeCount=v.nodeCount;
    d.initial=true;d.requested=true;d.reason='analytic_zero_time_mesh_selection';
    % Pre-selection source core differs from the accepted initial audit.
    d.coreCells=[16,18];d.predictedCoreCells=d.coreCells;
    for k=1:2
        d.decayEvidence(k).secantFromStep=0;d.decayEvidence(k).secantToStep=0;
        d.decayEvidence(k).secantFromTime=0;d.decayEvidence(k).secantToTime=0;
    end
    memory.lastDecision=d;
end
m.autonomousMesh=memory;
end

function m=forecast_fixture(m,c,decay)
d=m.autonomousMesh.transactions(1).controllerDecision;
d.coreCells=[30,31];d.decayEstimate=[decay,0];
d.predictedCoreCells=d.coreCells.*exp(-c.remesh.autonomousMesh.maximumReviewInterval*d.decayEstimate);
d.trendAvailable=true;d.reason='forecast_core_trigger';
for k=1:2
    d.decayEvidence(k).fitDecay=d.decayEstimate(k);
    d.decayEvidence(k).maximumSecant=0;
    d.decayEvidence(k).secantFromStep=0;d.decayEvidence(k).secantFromTime=0;
end
m.autonomousMesh.transactions(1).controllerDecision=d;
end

function [m,v]=same_level_fixture(m,v)
member=m.autonomousMesh.referenceFamily.members(2);
m.autonomousMesh.currentLevelId=2;
m.autonomousMesh.lastDecision.sourceLevelId=2;m.autonomousMesh.lastDecision.sourceNodeCount=member.nodeCount;
a=m.autonomousMesh.transactions(2);a.targetLevelId=2;a.targetNodeCount=member.nodeCount;
a.targetCellFactors=member.cellFactors;a.sameBoxAndNodeCount=true;a.xQuality=member.xQuality;a.yQuality=member.yQuality;
m.autonomousMesh.transactions(2)=a;
m.autonomousMesh.window.nodeCount=repmat(member.nodeCount,numel(m.autonomousMesh.window.time),1);
m.autonomousMesh.window.levelId=2*ones(numel(m.autonomousMesh.window.time),1);
v.levelId=2;v.nodeCount=member.nodeCount;v.x=member.baseX;v.baseX=member.baseX;v.y=member.baseY;v.baseY=member.baseY;
for k=4:5
    v.history.common.meshLevelId(k)=2;v.history.common.nodeCountX(k)=member.nodeCount(1);v.history.common.nodeCountY(k)=member.nodeCount(2);
    v.snapshots.x{k}=member.baseX;v.snapshots.y{k}=member.baseY;v.snapshots.rho{k}=zeros(member.nodeCount(2),member.nodeCount(1));
end
end

function r=initial_only_result_fixture(r)
% This only relabels a completed no-remesh result for boundary testing; the
% altered fixture is explicitly not evidence of a native version-two run.
assert(r.grid.remeshCount==0);
p=ipm.config.autonomousMeshPolicy(struct('version',2,'nodeFamily',struct('maximumTotalNodes',110000)));
r.config.remesh.autonomousMesh=p;m=r.metadata.autonomousMesh;m.version=2;m.policy=p;
m.referenceFamily=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX,m.initialization.selectedBaseY,r.config.scaling.transportAnchorX,p);
m.currentLevelId=1;nodes=[numel(r.grid.x),numel(r.grid.y)];
m.window.nodeCount=repmat(nodes,numel(m.window.time),1);m.window.levelId=ones(numel(m.window.time),1);
step=r.state.steps;tau=r.state.canonicalTime;
evidence=repmat(struct('axis',0,'fitDecay',0,'maximumSecant',0, ...
    'secantFromStep',step,'secantToStep',step,'secantFromTime',tau,'secantToTime',tau),1,2);
evidence(1).axis=1;evidence(2).axis=2;
m.lastDecision=struct('sourceStep',step,'sourceCanonicalTime',tau,'sourcePhysicalTime',r.state.physicalTime, ...
    'sourceNormalizedTime',tau,'sourceRemeshCount',0,'sourceLevelId',1,'sourceNodeCount',nodes, ...
    'requested',false,'initial',false,'coreCells',m.window.coreCells(end,:), ...
    'predictedCoreCells',m.window.coreCells(end,:),'decayEstimate',[0,0], ...
    'reason','synthetic_configuration_boundary_fixture','stopReason','', ...
    'trendAvailable',false,'decayEvidence',evidence);
r.metadata.autonomousMesh=m;r.metadata.caseMetadata.testFixtureOnly=true;r.metadata.caseMetadata.pdeRunClaim=false;
n=numel(r.history.common.t);r.history.common.nodeCountX=nodes(1)*ones(n,1);
r.history.common.nodeCountY=nodes(2)*ones(n,1);r.history.common.meshLevelId=ones(n,1);
end
function reject(m,c,v,id)
expect(@()ipm.remesh.validateController(m,c,v),id);
end
function expect(action,id)
try
    action();
catch exception
    assert(strcmp(exception.identifier,id),'ipm:ControllerV2UnexpectedError', ...
        'Expected %s, got %s: %s',id,exception.identifier,exception.message);return;
end
error('ipm:ControllerV2MissingError','Expected rejection %s.',id);
end
