function report=mesh_test_snapshot_shape_nolu(checkpointFiles,finalResultFile,outDir)
%MESH_TEST_SNAPSHOT_SHAPE_NOLU Independent supplemental actual-frame evidence.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
checkpointFiles=cellstr(string(checkpointFiles));rows=struct([]);
for k=1:numel(checkpointFiles)
    cp=ipm.output.readCheckpoint(checkpointFiles{k});s=cp.payload.state;log=cp.payload.log;
    assert(s.config.output.storeSnapshots && all(ipm.output.trustedMask(log.history,s.config)));
    snapshots=struct('rho',{log.snapshotRho},'x',{log.snapshotX},'y',{log.snapshotY});
    shapes=mesh_snapshot_shape_audit(snapshots);assert(shapes.spansMultipleNodeCounts);
    item=struct('checkpointFile',checkpointFiles{k},'step',s.step,'canonicalTime',s.scale.canonicalTime, ...
        'snapshotAudit',shapes,'nativeSignedReadPassed',true);
    if isempty(rows),rows=item;else,rows(end+1)=item;end %#ok<AGROW>
    if k>1
        assert(isequaln(previous.snapshotRho,log.snapshotRho(1:numel(previous.snapshotRho))) && ...
            isequaln(previous.snapshotX,log.snapshotX(1:numel(previous.snapshotX))) && ...
            isequaln(previous.snapshotY,log.snapshotY(1:numel(previous.snapshotY))) && ...
            isequaln(previous.snapshotNormalizedTime,log.snapshotNormalizedTime(1:numel(previous.snapshotNormalizedTime))));
    end
    previous=log;clear cp s
end
result=ipm.output.validate(finalResultFile);assert(all(ipm.output.trustedMask(result)));
resultAudit=mesh_snapshot_shape_audit(result.snapshots);
assert(resultAudit.spansMultipleNodeCounts && isequaln(result.snapshots.rho,log.snapshotRho) && ...
    isequaln(result.snapshots.x,log.snapshotX) && isequaln(result.snapshots.y,log.snapshotY) && ...
    isequaln(result.snapshots.normalizedTime,log.snapshotNormalizedTime));
bad=result.snapshots;bad.rho{end}=bad.rho{end}(1:end-1,:);caught=false;
try
    mesh_snapshot_shape_audit(bad);
catch exception
    caught=strcmp(exception.identifier,'ipm:SnapshotShapeAudit');
end
assert(caught);emptyAudit=mesh_snapshot_shape_audit(struct('rho',{{}},'x',{{}},'y',{{}}));
assert(emptyAudit.empty && ~emptyAudit.spansMultipleNodeCounts);
report=struct('kind','actual_snapshot_size_supplement_no_lu','allPassed',true, ...
    'checkpoints',rows,'finalResultFile',finalResultFile,'resultAudit',resultAudit, ...
    'allOldSnapshotPrefixesExact',true,'terminalNativeResultSnapshotsExact',true, ...
    'wrongRhoShapeRejected',true,'emptySnapshotsNotClaimedAsMultipleSizes',true, ...
    'oldReportsModified',false,'noLU',true,'noPDE',true);
save(fullfile(outDir,'report.mat'),'report','-v7.3');
f=fopen(fullfile(outDir,'report.json'),'w');assert(f>=0);cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('ACTUAL_SNAPSHOT_SHAPES_PASS checkpoints=%d actualFrames=%d distinctSizes=%d\n', ...
    numel(rows),resultAudit.actualSnapshotCount,size(resultAudit.uniqueActualNodeCounts,1));
end
