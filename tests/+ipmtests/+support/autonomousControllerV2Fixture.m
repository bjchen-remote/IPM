function [metadata,config,current,checkpoint] = autonomousControllerV2Fixture()
%IPMTESTS.SUPPORT.AUTONOMOUSCONTROLLERV2FIXTURE Synthetic persistence only.
% No PDE, physical invariant, or native numerical accuracy claim is made.
[metadata,config,current]=ipmtests.support.autonomousControllerFixture();
policy=ipm.config.autonomousMeshPolicy(struct('version',2, ...
    'nodeFamily',struct('maximumTotalNodes',110000)));
config.remesh.autonomousMesh=policy;config.output.storeSnapshots=true;
memory=metadata.autonomousMesh;memory.version=2;memory.policy=policy;
init=memory.initialization;
family=ipm.remesh.referenceAxisFamily(init.selectedBaseX,init.selectedBaseY, ...
    config.scaling.transportAnchorX,policy);
memory.referenceFamily=family;memory.currentLevelId=4;
levels=[1,2,4];transactions=struct([]);
for k=1:2
    a=memory.transactions(k);old=family.members(levels(k));new=family.members(levels(k+1));
    a.version=2;a.sourceLevelId=old.index;a.targetLevelId=new.index;
    a.sourceNodeCount=old.nodeCount;a.targetNodeCount=new.nodeCount;
    a.sourceCellFactors=old.cellFactors;a.targetCellFactors=new.cellFactors;
    a.sameBox=true;a.registeredNodeCounts=true;a.matchingReferenceMembers=true;
    a.admittedLevelTransition=true;a.sameBoxAndNodeCount=false;
    a.xQuality=new.xQuality;a.yQuality=new.yQuality;
    evidence=repmat(struct('axis',0,'fitDecay',0,'maximumSecant',0, ...
        'secantFromStep',a.sourceStep,'secantToStep',a.sourceStep, ...
        'secantFromTime',a.sourceCanonicalTime,'secantToTime',a.sourceCanonicalTime),1,2);
    evidence(1).axis=1;evidence(2).axis=2;
    a.controllerDecision=struct('sourceStep',a.sourceStep,'sourceCanonicalTime',a.sourceCanonicalTime, ...
        'sourcePhysicalTime',a.sourcePhysicalTime,'sourceNormalizedTime',a.sourceCanonicalTime, ...
        'sourceRemeshCount',a.sourceRemeshCount,'sourceLevelId',old.index,'sourceNodeCount',old.nodeCount, ...
        'requested',true,'initial',false,'coreCells',[25,32],'predictedCoreCells',[25,32], ...
        'decayEstimate',[0,0],'reason','current_core_trigger','stopReason','', ...
        'trendAvailable',false,'decayEvidence',evidence);
    a.candidateIndex=1;
    a.attemptSummary=struct('candidateIndex',1,'passed',true,'errorIdentifier','','auditReasons',{{}});
    if isempty(transactions),transactions=a;else,transactions(end+1)=a;end %#ok<AGROW>
end
memory.transactions=transactions;
member=family.members(4);current.x=member.baseX;current.y=member.baseY;
current.baseX=member.baseX;current.baseY=member.baseY;current.levelId=4;current.nodeCount=member.nodeCount;
memory.window.nodeCount=repmat(member.nodeCount,numel(memory.window.time),1);
memory.window.levelId=4*ones(numel(memory.window.time),1);
memory.lastDecision=transactions(end).controllerDecision;
d=memory.lastDecision;d.sourceStep=current.step;d.sourceCanonicalTime=current.canonicalTime;
d.sourcePhysicalTime=current.physicalTime;d.sourceNormalizedTime=current.normalizedTime;
d.sourceRemeshCount=current.remeshCount;d.sourceLevelId=current.levelId;d.sourceNodeCount=current.nodeCount;
d.requested=false;d.reason='resolved';d.coreCells=current.coreCells;d.predictedCoreCells=current.coreCells;
for k=1:2
    d.decayEvidence(k).secantFromStep=current.step;d.decayEvidence(k).secantToStep=current.step;
    d.decayEvidence(k).secantFromTime=current.canonicalTime;d.decayEvidence(k).secantToTime=current.canonicalTime;
end
memory.lastDecision=d;
historyLevels=[1;2;2;4;4];n=numel(historyLevels);
current.history.common.meshLevelId=historyLevels;
current.history.common.nodeCountX=zeros(n,1);current.history.common.nodeCountY=zeros(n,1);
current.snapshots=struct('rho',{cell(n,1)},'x',{cell(n,1)},'y',{cell(n,1)});
for k=1:n
    member=family.members(historyLevels(k));
    current.history.common.nodeCountX(k)=member.nodeCount(1);
    current.history.common.nodeCountY(k)=member.nodeCount(2);
    current.snapshots.x{k}=member.baseX;current.snapshots.y{k}=member.baseY;
    current.snapshots.rho{k}=zeros(member.nodeCount(2),member.nodeCount(1));
end
metadata.autonomousMesh=memory;
metadata.caseId='synthetic_variable_node_controller_contract_only';
% Add the existing generic checkpoint log contract using explicit synthetic
% identity rates and benign quality observations, not invented solver output.
zero=zeros(n,1);common=current.history.common;
names={'c_l','c_x','c_y','c_omega','c_r','conservativeSource','canonicalCL', ...
    'canonicalCX','canonicalCY','canonicalCOmega','canonicalCR','canonicalConservativeSource', ...
    'wallPeakLocalMaxima','wallPeakTVRatio','positiveWallNegativeRatio','physicalRangeViolation'};
for k=1:numel(names),common.(names{k})=zero;end
common.timeSpeed=ones(n,1);current.history.common=common;
current.history.mesh.maximumCellRatioX=ones(n,1);current.history.mesh.maximumCellRatioY=ones(n,1);
current.history.gauge=struct();current.history.anisotropic=struct();
scale=struct('logC_l',0,'logC_omega',0,'physicalTime',current.physicalTime, ...
    'canonicalTime',current.canonicalTime,'X_shift',0);
state=struct('config',config,'runMetadata',metadata,'rho',current.snapshots.rho{end}, ...
    'x',current.x,'y',current.y,'baseX',current.baseX,'baseY',current.baseY, ...
    'rescaling',struct(),'remeshCount',current.remeshCount,'scale',scale, ...
    'normalizedTime',current.normalizedTime,'step',current.step,'mass0',1,'rhoRange0',[0,1]);
log=struct('history',current.history,'snapshotRho',{current.snapshots.rho}, ...
    'snapshotX',{current.snapshots.x},'snapshotY',{current.snapshots.y}, ...
    'snapshotNormalizedTime',current.history.common.t);
cursor=struct('nextOutput',1,'nextCheckpoint',Inf,'lastCheckpointStep',current.step, ...
    'lastCheckpointCanonicalTime',current.canonicalTime);
payload=struct('state',state,'log',log,'cursor',cursor);
checkpoint=struct('schemaVersion',4,'kind','ipm_accepted_step_checkpoint', ...
    'createdUtc','2026-09-09T00:00:00.000Z','acceptedStep',true,'trustedRecord',true, ...
    'payload',payload,'signature',ipm.output.checkpointSignature(payload));
end
