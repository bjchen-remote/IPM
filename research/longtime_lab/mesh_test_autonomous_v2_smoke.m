function report=mesh_test_autonomous_v2_smoke(probeFile,outDir)
%MESH_TEST_AUTONOMOUS_V2_SMOKE Original analytic t=0 and native exact restart.
% This short test exercises registration/persistence only, not natural growth.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
q=load(probeFile,'o','result');o=q.o;oldResult=q.result;clear q
o.autonomousMesh=struct('version',2,'nodeFamily',struct('maximumTotalNodes',110000));
report=struct('kind','autonomous_v2_zero_time_persistence_smoke', ...
    'naturalGrowthQualification',false,'sourceProbe',probeFile,'tests',struct());
save(fullfile(outDir,'registration.mat'),'o','report');
state=ipm.evolve.initialize(o);log=ipm.output.initializeLog();[log,~]=ipm.output.record(log,state);
m=state.runMetadata.autonomousMesh;
assert(m.version==2 && m.currentLevelId==1 && ...
    isequal([m.referenceFamily.members.resourceAdmitted],[true,true,true,false]));
assert(isequal(state.rho,ipm.field.initialDensity(state.ops,state.config.physics)) && ...
    state.step==0 && state.scale.physicalTime==0 && state.scale.canonicalTime==0);
report.tests.exactAnalyticZeroTimeAndFrozenResourceMask=true;
for k=1:2
    [state,stop,plan]=ipm.evolve.advance(state);assert(isempty(stop)&&isempty(plan));
    [log,~]=ipm.output.record(log,state);
end
cursor=struct('nextOutput',state.normalizedTime+o.outputEvery,'nextCheckpoint',Inf, ...
    'lastCheckpointStep',state.step,'lastCheckpointCanonicalTime',state.scale.canonicalTime);
checkpoint=ipm.output.makeCheckpoint(state,log,cursor);
[checkpoint,checkpointFile]=ipm.output.writeCheckpoint(checkpoint,fullfile(outDir,'split_checkpoint.mat'));
roundtrip=ipm.output.readCheckpoint(checkpointFile);
assert(isequaln(checkpoint.payload,roundtrip.payload));
report.tests.signedNativeCheckpointExact=true;
for k=1:2
    [state,stop,plan]=ipm.evolve.advance(state);assert(isempty(stop)&&isempty(plan));
    [log,~]=ipm.output.record(log,state);
end
baseline=ipm.output.finalize(state,log,'final_time');clear state log m
resumed=ipm.solve(struct(),checkpoint);
groups={'state','physical','history','grid','scale','snapshots'};
for k=1:numel(groups),assert(isequaln(baseline.(groups{k}),resumed.(groups{k})));end
assert(isequaln(baseline.metadata.autonomousMesh,resumed.metadata.autonomousMesh));
report.tests.fourStepSplitTrajectoryAndLedgerBitwise=true;
assert(isequal(oldResult.state.rho,baseline.state.rho) && ...
    isequal(oldResult.grid.x,baseline.grid.x) && isequal(oldResult.grid.y,baseline.grid.y) && ...
    isequaln(oldResult.scale,baseline.scale) && isequaln(oldResult.physical,baseline.physical));
report.tests.firstFourStepsMatchFrozenV1FieldGridScalePhysical=true;
ipm.output.validateV2(oldResult);report.tests.frozenV1ResultStillValid=true;
rejected=false;
try
    ipm.output.checkpointFromResult(baseline);
catch e
    rejected=strcmp(e.identifier,'ipm:CheckpointResultVariableNodeFamily');
    report.bridgeRejectionIdentifier=e.identifier;
end
assert(rejected);report.tests.unsupportedResultCheckpointBridgeExplicitlyRejected=true;
report.checkpointFile=char(checkpointFile);report.allPassed=all(structfun(@(x)isequal(x,true),report.tests));
save(fullfile(outDir,'results.mat'),'baseline','resumed','report','-v7.3');
f=fopen(fullfile(outDir,'report.json'),'w');assert(f>=0);cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(report,PrettyPrint=true));assert(report.allPassed);
fprintf('AUTONOMOUS_V2_SMOKE_PASS originalZeroTime=1 splitBitwise=1 v1PhysicalBitwise=1\n');
end
