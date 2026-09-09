function report=mesh_test_v2_natural_restart(runDir,oldRunDir,outDir)
%MESH_TEST_V2_NATURAL_RESTART Replay a naturally selected growth transaction.
% Only run after the uninterrupted registered experiment has finished.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
summary=jsondecode(fileread(fullfile(runDir,'report.json')));
assert(summary.registeredHorizonReached && summary.previousCapacityTimeCrossed && ...
    summary.noManualTransactionTriggers && ~isempty(summary.growthTransactions));
report=struct('kind','natural_variable_node_checkpoint_restart_regression', ...
    'sourceRun',runDir,'oldFixedNRun',oldRunDir,'tests',struct());
listing=dir(fullfile(runDir,'checkpoint_*_step*.mat'));rows=struct([]);
for k=1:numel(listing)
    name=fullfile(listing(k).folder,listing(k).name);cp=ipm.output.readCheckpoint(name);s=cp.payload.state;
    row=struct('file',name,'step',s.step,'time',s.scale.canonicalTime, ...
        'level',s.runMetadata.autonomousMesh.currentLevelId,'nodes',[numel(s.x),numel(s.y)]);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    clear cp s
end
[~,order]=sort([rows.step]);rows=rows(order);
first=find([rows.level]>1,1);assert(first>1 && rows(first-1).level==1);
report.nativeCheckpointsRead=numel(rows);report.tests.allNaturalCheckpointsStrictlyReadable=true;
report.sourceCheckpoint=rows(first-1).file;report.targetCheckpoint=rows(first).file;
save(fullfile(outDir,'registration.mat'),'report','rows');

% Independent v1 and v2 runs have different case identities; compare only
% the actual shared numerical prefix before the first natural node growth.
oldListing=dir(fullfile(oldRunDir,'checkpoint_*_step*.mat'));matched=0;
stateFields={'rho','x','y','baseX','baseY','rescaling','remeshCount','scale', ...
    'normalizedTime','step','mass0','rhoRange0'};
for k=1:numel(oldListing)
    old=ipm.output.readCheckpoint(fullfile(oldListing(k).folder,oldListing(k).name));
    match=find([rows.step]==old.payload.state.step & [rows.level]==1,1);
    if isempty(match),clear old;continue;end
    new=ipm.output.readCheckpoint(rows(match).file);
    for j=1:numel(stateFields)
        field=stateFields{j};assert(isequaln(old.payload.state.(field),new.payload.state.(field)));
    end
    groups=fieldnames(old.payload.log.history);
    for g=1:numel(groups)
        names=fieldnames(old.payload.log.history.(groups{g}));
        for j=1:numel(names)
            field=names{j};assert(isequaln(old.payload.log.history.(groups{g}).(field), ...
                new.payload.log.history.(groups{g}).(field)));
        end
    end
    matched=matched+1;clear old new
end
assert(matched>=2);report.v1MatchedNativePrefixCheckpoints=matched;
report.tests.frozenV1SharedFieldsAndEveryExistingHistoryColumnBitwise=true;

source=ipm.output.readCheckpoint(report.sourceCheckpoint);
target=ipm.output.readCheckpoint(report.targetCheckpoint);
overrides=struct('maxSteps',target.payload.state.step,'saveResults',false, ...
    'checkpoint',struct('enabled',true,'file',fullfile(outDir,'replay_checkpoint.mat'), ...
        'every',.2,'atExit',true), ...
    'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false);
resumed=ipm.solve(overrides,source);clear source
s=target.payload.state;log=target.payload.log;
assert(isequal(resumed.state.rho,s.rho) && isequal(resumed.grid.x,s.x) && ...
    isequal(resumed.grid.y,s.y) && resumed.grid.remeshCount==s.remeshCount && ...
    resumed.state.steps==s.step && resumed.state.canonicalTime==s.scale.canonicalTime && ...
    resumed.state.physicalTime==s.scale.physicalTime && ...
    isequaln(resumed.history,log.history) && isequal(resumed.snapshots.rho,log.snapshotRho) && ...
    isequal(resumed.snapshots.x,log.snapshotX) && isequal(resumed.snapshots.y,log.snapshotY) && ...
    isequaln(resumed.metadata.autonomousMesh,s.runMetadata.autonomousMesh));
report.tests.naturalGrowthReplayFieldsClocksHistorySnapshotsControllerBitwise=true;
replayCheckpoint=ipm.output.readCheckpoint(resumed.metadata.latestCheckpointFile);
for j=1:numel(stateFields)
    field=stateFields{j};assert(isequaln(replayCheckpoint.payload.state.(field),s.(field)));
end
report.tests.replayedNativeScalesReferencesAxesAndStateBitwise=true;
report.replayedNativeCheckpoint=char(resumed.metadata.latestCheckpointFile);clear replayCheckpoint
report.actualSnapshotShapeAudit=mesh_snapshot_shape_audit(resumed.snapshots);
assert(report.actualSnapshotShapeAudit.spansMultipleNodeCounts);
report.tests.actualSnapshotsSpanMultipleNodeCounts=report.actualSnapshotShapeAudit.spansMultipleNodeCounts;
save(fullfile(outDir,'replayed_result.mat'),'resumed','-v7.3');clear resumed log s

% Re-sign malformed payloads so semantic rejection, rather than a stale
% signature, is exercised before an operator can be built.
names={'wrong_active_level','changed_initial_budget','wrong_active_base', ...
    'false_resource_mask','wrong_history_node_count','wrong_snapshot_pair'};
for k=1:numel(names)
    bad=target;
    switch names{k}
        case 'wrong_active_level',bad.payload.state.runMetadata.autonomousMesh.currentLevelId=1;
        case 'changed_initial_budget'
            bad.payload.state.config.grid.nx=numel(bad.payload.state.x);
            bad.payload.state.config.grid.ny=numel(bad.payload.state.y);
        case 'wrong_active_base'
            bad.payload.state.baseX=bad.payload.state.runMetadata.autonomousMesh.referenceFamily.rootX;
            bad.payload.state.baseY=bad.payload.state.runMetadata.autonomousMesh.referenceFamily.rootY;
        case 'false_resource_mask',bad.payload.state.runMetadata.autonomousMesh.referenceFamily.members(4).resourceAdmitted=true;
        case 'wrong_history_node_count'
            bad.payload.log.history.common.nodeCountX(end)=numel(bad.payload.state.x)-2;
        case 'wrong_snapshot_pair'
            bad.payload.log.snapshotX{end}=bad.payload.log.snapshotX{1};
            bad.payload.log.snapshotY{end}=bad.payload.log.snapshotY{1};
    end
    bad.signature=ipm.output.checkpointSignature(bad.payload);rejected=false;
    try
        ipm.output.readCheckpoint(bad);
    catch
        rejected=true;
    end
    assert(rejected);clear bad
end
report.semanticNegativeChecks=numel(names);report.tests.signedSemanticCorruptionsRejected=true;
report.allPassed=all(structfun(@(x)isequal(x,true),report.tests));
save(fullfile(outDir,'report.mat'),'report','rows');
f=fopen(fullfile(outDir,'report.json'),'w');assert(f>=0);cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(report,PrettyPrint=true));assert(report.allPassed);
fprintf('AUTONOMOUS_V2_NATURAL_RESTART_PASS prefixCPs=%d acrossGrowthBitwise=1 negatives=%d\n',matched,numel(names));
end
