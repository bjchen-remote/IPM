function report = autonomousController(oldCheckpointFile,actualResultFiles)
%IPMTESTS.BASELINE.AUTONOMOUSCONTROLLER No-LU persistent-ledger boundaries.
if nargin<1,oldCheckpointFile='';end
if nargin<2,actualResultFiles={};end
[m,c,v] = ipmtests.support.autonomousControllerFixture();
expected = ipm.remesh.validateController(m,c,v);assert(expected.checked && expected.transactions==2);
withoutBase = rmfield(v,{'baseX','baseY'});assert(ipm.remesh.validateController(m,c,withoutBase).checked);
bare = rmfield(m,'autonomousMesh');legacy = c;legacy.remesh = rmfield(legacy.remesh,'autonomousMesh');
assert(~ipm.remesh.validateController(bare,legacy,struct()).checked);
off = c;off.remesh.autonomousMesh.enabled = false;assert(~ipm.remesh.validateController(bare,off,struct()).checked);
id = 'ipm:AutonomousMeshController';labels = {};
bad = m;bad.autonomousMesh.version = 2;reject(bad,c,v,id);labels{end+1}='version';
bad = m;bad.autonomousMesh.policy.enabled = 1;reject(bad,c,v,id);labels{end+1}='policy logical class';
bad = m;bad.autonomousMesh.policy.timeUnit = "native_canonical";reject(bad,c,v,id);labels{end+1}='policy text class';
bad = m;bad.autonomousMesh.policy.targetCoreCells = [42,42];reject(bad,c,v,id);labels{end+1}='policy values';
reject(bare,c,v,id);labels{end+1}='missing ledger';reject(m,legacy,struct(),id);labels{end+1}='ledger without policy';
reject(m,off,struct(),id);labels{end+1}='disabled ledger';
fields = {'originalStep','originalCanonicalTime','originalPhysicalTime','originalRemeshCount'};
for k=1:numel(fields)
    bad=m;bad.autonomousMesh.initialization.(fields{k})=1;reject(bad,c,v,id);labels{end+1}=['initial ',fields{k}]; %#ok<AGROW>
end
bad=m;bad.autonomousMesh.initialization.newPhysicalEpochCreated=true;reject(bad,c,v,id);labels{end+1}='new physical epoch';
bad=m;bad.autonomousMesh.initialization=rmfield(bad.autonomousMesh.initialization,'operatorGridRepresentation');reject(bad,c,v,id);labels{end+1}='missing operator representation';
bad=m;bad.autonomousMesh.initialization.operatorGridRepresentation='analytic_axes';reject(bad,c,v,id);labels{end+1}='changed operator representation';
bad=m;bad.autonomousMesh.initialization.audit.passed=false;reject(bad,c,v,id);labels{end+1}='initial audit failed';
bad=m;bad.autonomousMesh.initialization.selectedCandidateIndex=2;reject(bad,c,v,id);labels{end+1}='initial attempt index';
bad=m;bad.autonomousMesh.initialization.attempts.errorMessage='failed';reject(bad,c,v,id);labels{end+1}='initial attempt error';
badV=v;badV.baseX(2)=badV.baseX(2)+eps;reject(m,c,badV,id);labels{end+1}='base axes drift';
bad=m;bad.caseMetadata.latePhysicalBoxBranch=struct();reject(bad,c,v,id);labels{end+1}='late epoch relabel';
fields={'passed','clocksAndReferencesPreserved','sameBoxAndNodeCount'};
for k=1:numel(fields)
    bad=m;bad.autonomousMesh.transactions(2).(fields{k})=false;reject(bad,c,v,id);labels{end+1}=['transaction ',fields{k}]; %#ok<AGROW>
end
bad=m;bad.autonomousMesh.transactions(2).sourceRemeshCount=0;reject(bad,c,v,id);labels{end+1}='source count';
bad=m;bad.autonomousMesh.transactions(2).sourceStep=10;reject(bad,c,v,id);labels{end+1}='source step repeat';
bad=m;bad.autonomousMesh.transactions(2).sourceCanonicalTime=.3;reject(bad,c,v,id);labels{end+1}='future transaction';
bad=m;bad.autonomousMesh.transactions(2).sourcePhysicalTime=.02;reject(bad,c,v,id);labels{end+1}='physical backwards';
bad=m;bad.autonomousMesh.transactions(2).priorCumulativeAbsolutePeakJump=0;reject(bad,c,v,id);labels{end+1}='prior cumulative';
bad=m;bad.autonomousMesh.transactions(2).cumulativeAbsolutePeakJump=.0004;reject(bad,c,v,id);labels{end+1}='row cumulative';
bad=m;bad.autonomousMesh.cumulativeAbsolutePeakJump=.0004;reject(bad,c,v,id);labels{end+1}='terminal cumulative';
bad=m;bad.autonomousMesh.transactions(2).relativePeakJump=.003;reject(bad,c,v,id);labels{end+1}='single peak gate';
bad=m;bad.autonomousMesh.transactions(2).coreCells=[30,40];reject(bad,c,v,id);labels{end+1}='core gate';
bad=m;bad.autonomousMesh.transactions(2).xQuality.minimumStencilRcond=0;reject(bad,c,v,id);labels{end+1}='quality gate';
bad=m;bad.autonomousMesh.transactions(2).massRelativeDefect=1e-8;reject(bad,c,v,id);labels{end+1}='mass gate';
bad=m;bad.autonomousMesh.transactions(2).relativeRangeViolation=.01;reject(bad,c,v,id);labels{end+1}='range gate';
bad=m;bad.autonomousMesh.window.time(2)=.2;reject(bad,c,v,id);labels{end+1}='window time repeated';
bad=m;bad.autonomousMesh.window.step(2)=20;reject(bad,c,v,id);labels{end+1}='window step repeated';
bad=m;bad.autonomousMesh.window.remeshCount(1)=1;reject(bad,c,v,id);labels{end+1}='window old epoch';
bad=m;bad.autonomousMesh.window.coreCells(end,1)=33;reject(bad,c,v,id);labels{end+1}='window terminal core';
bad=m;bad.autonomousMesh.window.safety(end)=.31;reject(bad,c,v,id);labels{end+1}='window terminal safety';
bad=m;bad.autonomousMesh.window.coreCells(2,1)=NaN;reject(bad,c,v,id);labels{end+1}='window nonfinite';
badV=v;badV.history.mesh.remeshCount(3)=0;reject(m,c,badV,id);labels{end+1}='history count mismatch';
badV=v;badV.history.common.physicalTime(1)=1e-6;reject(m,c,badV,id);labels{end+1}='history epoch';
bad=m;bad.autonomousMesh.window.coreCells(1,1)=36;reject(bad,c,v,id);labels{end+1}='shared window history';
bad=m;badV=v;badV.canonicalTime=1;badV.normalizedTime=1;bad.autonomousMesh.window.time(end)=1;
reject(bad,c,badV,id);labels{end+1}='window outside bounded time';
% Exercise the two actual accumulation orders on a longer synthetic ledger.
long=m;longV=v;n=128;jumps=(1+mod(1:n,13))*1e-6;differentReductionRows=0;
transactions=repmat(m.autonomousMesh.transactions(1),1,n);
for k=1:n
    transactions(k).sourceStep=k;transactions(k).sourceCanonicalTime=.001*k;
    transactions(k).sourcePhysicalTime=.0001*k;transactions(k).sourceRemeshCount=k-1;
    transactions(k).relativePeakJump=jumps(k);transactions(k).priorCumulativeAbsolutePeakJump=sum(jumps(1:k-1));
    transactions(k).cumulativeAbsolutePeakJump=transactions(k).priorCumulativeAbsolutePeakJump+jumps(k);
    differentReductionRows=differentReductionRows+(transactions(k).cumulativeAbsolutePeakJump~=sum(jumps(1:k)));
end
long.autonomousMesh.transactions=transactions;long.autonomousMesh.cumulativeAbsolutePeakJump=sum(jumps);
long.autonomousMesh.window=struct('time',[.18;.19;.2],'step',[180;190;200], ...
    'remeshCount',n*ones(3,1),'coreCells',[34,39;33,38;32,37],'safety',.3*ones(3,1));
longV.step=200;longV.canonicalTime=.2;longV.normalizedTime=.2;longV.physicalTime=.04;longV.remeshCount=n;
longV.history=struct('common',struct('acceptedStep',[0;200],'canonicalTau',[0;.2],'physicalTime',[0;.04],'t',[0;.2]), ...
    'mesh',struct('remeshCount',[0;n],'coreGridPoints',[40;32],'verticalCoreGridPoints',[44;37],'safetyFactor',[.3;.3]));
longReport=ipm.remesh.validateController(long,c,longV);assert(longReport.transactions==n);
old = struct('tested',false);
if ~isempty(oldCheckpointFile)
    assert(maxNumCompThreads==10);checkpoint=ipm.output.readCheckpoint(oldCheckpointFile);
    signature=checkpoint.signature;again=ipm.output.readCheckpoint(checkpoint);assert(isequaln(signature,again.signature));
    badCheckpoint=checkpoint;badCheckpoint.payload.state.runMetadata.autonomousMesh=m.autonomousMesh;
    expect(@()ipm.output.readCheckpoint(badCheckpoint),id);expect(@()ipm.output.restoreCheckpoint(badCheckpoint),id);
    old=struct('tested',true,'step',checkpoint.payload.state.step,'signatureExact',true, ...
        'malformedRestoreRejectedBeforeLU',true,'positiveRestoreExecuted',false);
end
actual = struct('file',{},'step',{},'transactions',{},'windowRecords',{});
for k=1:numel(actualResultFiles)
    loaded=load(actualResultFiles{k},'result');r=loaded.result;h=r.history;
    projection=struct('step',r.state.steps,'canonicalTime',r.state.canonicalTime,'physicalTime',r.state.physicalTime, ...
        'normalizedTime',r.state.normalizedTime,'remeshCount',r.grid.remeshCount, ...
        'coreCells',[h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)],'safety',h.mesh.safetyFactor(end), ...
        'x',r.grid.x,'y',r.grid.y,'history',h);
    a=ipm.remesh.validateController(r.metadata,r.config,projection);
    actual(end+1)=struct('file',actualResultFiles{k},'step',r.state.steps,'transactions',a.transactions,'windowRecords',a.windowRecords); %#ok<AGROW>
end
report=struct('allPassed',true,'syntheticFixtureOnly',true,'syntheticLedgerNotPdeEvidence',true, ...
    'negativeLabels',{labels},'negativeCount',numel(labels),'actualResults',actual,'oldCheckpoint',old, ...
    'longLedgerTransactions',n,'rowsWhereAccumulationOrdersDiffer',differentReductionRows, ...
    'noLU',true,'noPDE',true,'postRestorePositiveHookNotExecuted',true);
fprintf('AUTONOMOUS_CONTROLLER_PASS negatives=%d actualResults=%d oldCheckpoint=%d noLU=1\n',numel(labels),numel(actual),old.tested);
end
function reject(m,c,v,id)
expect(@()ipm.remesh.validateController(m,c,v),id);
end
function expect(action,id)
try
    action();
catch exception
    assert(strcmp(exception.identifier,id),'ipm:ControllerUnexpectedError','Expected %s, got %s: %s',id,exception.identifier,exception.message);return;
end
error('ipm:ControllerMissingError','Expected rejection %s.',id);
end
