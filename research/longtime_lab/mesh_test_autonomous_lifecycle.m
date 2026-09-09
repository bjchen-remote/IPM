function report=mesh_test_autonomous_lifecycle(probeFile,outDir)
%MESH_TEST_AUTONOMOUS_LIFECYCLE Actual small-state transfer/recovery/restart.
% The two migrations are explicitly injected transaction tests. They do not
% claim that the automatic trigger fired along the four-step trajectory.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
fixture=load(probeFile,'o');o=fixture.o;
report=struct('kind','autonomous_actual_lifecycle_regression_v1', ...
    'automaticLongtimeQualification',false,'transactionTriggerManuallyInjected',true, ...
    'sourceProbe',probeFile,'tests',struct());
save(fullfile(outDir,'registration.mat'),'report','o');
state=ipm.evolve.initialize(o);log=ipm.output.initializeLog();[log,~]=ipm.output.record(log,state);
assert(all(ipm.output.trustedMask(log.history,state.config)));
for k=1:2
    [state,stop,automaticPlan]=ipm.evolve.advance(state);
    assert(isempty(stop)&&isempty(automaticPlan));
    q=state.ops;Dx=q.Dx;
    view=struct('rho',state.rho,'x',q.x,'y',q.y,'Dx',Dx,'source',state.rho*Dx','trusted',true);
    reference=struct('x',q.baseX,'y',q.baseY);
    [pairs,~]=ipm.remesh.plannedAxisPairs(view,reference,o.transportAnchorX,state.config.remesh.autonomousMesh);
    pairs=pairs(~[pairs.unchanged]);assert(~isempty(pairs));
    clear q Dx view reference
    if k==1
        % An admissible, monotone axis with broken x symmetry must be
        % rejected after real construction, without committing source data.
        bad=pairs(1);j=(numel(bad.x)+1)/2+3;
        bad.x(j)=bad.x(j)+.001*(bad.x(j+1)-bad.x(j));
        plan=struct('initial',false,'candidates',bad);
        b=sin((1:size(state.ops.A,1))');oldAction=state.ops.poisson\b;
        state.ops=rmfield(state.ops,'poisson');saved=state;
        [state,ok,attempts]=ipm.evolve.applyAutonomousMesh(state,plan);
        assert(~ok && ~any([attempts.passed]) && isequaln(state,saved));
        assert(any(strcmp(attempts(1).audit.reasons,'changed_box_reference_or_node_count')));
        state.ops.poisson=decomposition(state.ops.A,'lu');
        assert(isequal(oldAction,state.ops.poisson\b));
        [recoveredRate,recoveredFlow]=ipm.evolve.flow(state.rho,state.ops,state.scale);
        assert(isequaln(recoveredFlow,saved.flow) && isequal(recoveredRate,saved.rhsCache.rhoRate));
        assert(isequaln(rmfield(state.ops,'poisson'),saved.ops));
        report.tests.recoveredActualRhsFlowAndAllOperatorDataBitwise=true;
        clear recoveredRate recoveredFlow
        report.tests.actualRejectedTransferRollbackExact=true;
        report.tests.originalSparseFactorActionRebuiltBitwise=true;
        clear saved b oldAction bad plan
    end
    source=state;source.ops=rmfield(source.ops,'poisson');
    state.ops=rmfield(state.ops,'poisson');
    plan=struct('initial',false,'candidates',pairs);
    [state,ok,attempts]=ipm.evolve.applyAutonomousMesh(state,plan);
    assert(ok && state.ops.remeshCount==k && numel(state.runMetadata.autonomousMesh.transactions)==k);
    assert(isequal(state.scale,source.scale)&&isequal(state.config,source.config)&& ...
        isequal(state.ops.baseX,source.ops.baseX)&&isequal(state.ops.baseY,source.ops.baseY)&& ...
        state.step==source.step&&state.normalizedTime==source.normalizedTime&& ...
        state.mass0==source.mass0&&isequal(state.rhoRange0,source.rhoRange0));
    assert(numel(attempts)<=3&&sum([attempts.passed])==1);
    [log,~]=ipm.output.record(log,state);clear source pairs plan attempts
end
report.tests.twoActualNativeTransfersPassed=true;
cursor=struct('nextOutput',state.normalizedTime+o.outputEvery,'nextCheckpoint',Inf, ...
    'lastCheckpointStep',state.step,'lastCheckpointCanonicalTime',state.scale.canonicalTime);
checkpoint=ipm.output.makeCheckpoint(state,log,cursor);
[checkpoint,checkpointFile]=ipm.output.writeCheckpoint(checkpoint,fullfile(outDir,'split_checkpoint.mat'));
roundtrip=ipm.output.readCheckpoint(checkpointFile);
assert(isequaln(checkpoint.payload,roundtrip.payload));
report.tests.signedNativeCheckpointRoundtripExact=true;
for k=1:2
    [state,stop,plan]=ipm.evolve.advance(state);assert(isempty(stop)&&isempty(plan));
    [log,~]=ipm.output.record(log,state);
end
baseline=ipm.output.finalize(state,log,'final_time');clear state log
resumed=ipm.solve(struct(),checkpoint);
assert(isequaln(baseline.state,resumed.state)&&isequaln(baseline.physical,resumed.physical)&& ...
    isequaln(baseline.history,resumed.history)&&isequaln(baseline.grid,resumed.grid)&& ...
    isequaln(baseline.scale,resumed.scale)&&isequaln(baseline.snapshots,resumed.snapshots));
assert(isequaln(baseline.metadata.autonomousMesh,resumed.metadata.autonomousMesh));
report.tests.splitTrajectoryFieldsAndLedgerBitwise=true;
report.checkpointFile=char(checkpointFile);report.finalStep=resumed.state.steps;
report.finalCoreCells=resumed.metadata.autonomousMesh.window.coreCells(end,:);
report.allPassed=all(structfun(@(v)isequal(v,true),report.tests));
save(fullfile(outDir,'results.mat'),'baseline','resumed','report','-v7.3');
fid=fopen(fullfile(outDir,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
assert(report.allPassed);fprintf('AUTONOMOUS_LIFECYCLE_PASS twoActualTransfers=1 rollbackExact=1 checkpointSplitBitwise=1\n');
end
