function report=test_v5_frozen_ratios(checkpointFile,outputFile)
%TEST_V5_FROZEN_RATIOS Check automatic levels/minimax order on native data.
% It neither builds an elliptic operator nor transfers/advances the field.
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;
assert(all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
initial=[s.config.grid.nx,s.config.grid.ny];
policy=ipm.config.autonomousMeshPolicy(struct('version',5, ...
    'nodeFamily',struct('maximumTotalNodes',1000000, ...
    'initialNodeCount',initial)));
root=s.runMetadata.autonomousMesh.referenceFamily;
family=ipm.remesh.referenceAxisFamily(root.rootX,root.rootY, ...
    s.config.scaling.transportAnchorX,policy);
Dx=ipm.mesh.fdMatrix(s.x,1,7);
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx, ...
    'source',s.rho*Dx','trusted',true);
requested=[2,2;3,2;2,3;4,2;3,3];rows=struct([]);
factors=vertcat(family.members.cellFactors);
for k=1:size(requested,1)
    id=find(all(factors==requested(k,:),2));assert(isscalar(id));
    member=family.members(id);
    [candidates,detail,search]=ipm.remesh.hierarchicalAxisPairs(view, ...
        struct('x',member.baseX,'y',member.baseY), ...
        s.config.scaling.transportAnchorX,policy,policy.nodeFamily.maximumTotalNodes);
    goodX=find([detail.xTrials.admissible]);
    goodY=find([detail.yTrials.admissible]);
    assert(~isempty(goodX)&&~isempty(goodY)&&~isempty(candidates));
    ratios=zeros(numel(goodX)*numel(goodY),1);index=0;
    for i=goodX
        for j=goodY
            index=index+1;
            ratios(index)=max(detail.xTrials(i).quality.maximumAdjacentCellRatio, ...
                detail.yTrials(j).quality.maximumAdjacentCellRatio);
        end
    end
    changed=find(~[candidates.unchanged],1);
    assert(~isempty(changed));c=candidates(changed);
    selected=max(c.quality.x.maximumAdjacentCellRatio, ...
        c.quality.y.maximumAdjacentCellRatio);
    best=min(ratios);
    assert(abs(selected-best)<=1e-12 && selected<=policy.qualityLimits.maxAdjacentCellRatio);
    row=struct('factor',requested(k,:),'levelId',id,'nodeCount',member.nodeCount, ...
        'totalNodes',prod(member.nodeCount),'admittedX',numel(goodX), ...
        'admittedY',numel(goodY),'selectedMaximumAdjacentRatio',selected, ...
        'enumeratedMinimaxRatio',best,'searchStatus',search.status);
    if isempty(rows),rows=row;else,rows(end+1)=row;end
    fprintf('V5_FROZEN factor=%dx%d n=%dx%d ratio=%.9f\n', ...
        requested(k,1),requested(k,2),member.nodeCount(1),member.nodeCount(2),selected);
end
assert(rows(2).selectedMaximumAdjacentRatio<rows(1).selectedMaximumAdjacentRatio && ...
    abs(rows(3).selectedMaximumAdjacentRatio-rows(1).selectedMaximumAdjacentRatio)<1e-12);
assert(rows(5).totalNodes>rows(2).totalNodes && ...
    rows(5).selectedMaximumAdjacentRatio<=rows(2).selectedMaximumAdjacentRatio+1e-12);
report=struct('kind','v5_frozen_minimax_direction_test_v1', ...
    'sourceCheckpoint',checkpointFile,'sourceStep',s.step, ...
    'sourceCanonicalTime',s.scale.canonicalTime,'cap',policy.nodeFamily.maximumTotalNodes, ...
    'registeredLevels',numel(family.members),'rows',rows, ...
    'noLU',true,'noTransfer',true,'noPDE',true,'passed',true);
if nargin>=2&&~isempty(outputFile)
    assert(~isfile(outputFile));fid=fopen(outputFile,'w');assert(fid>=0);
    closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
end
