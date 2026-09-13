function report=audit_frozen_mesh_candidate_levels(checkpointFile,levelIds,outputFile)
%AUDIT_FROZEN_MESH_CANDIDATE_LEVELS Audit registered node allocation only.
% This screens geometry/axis quality on one accepted field, without transfer.
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;
assert(all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
policy=s.config.remesh.autonomousMesh;
family=s.runMetadata.autonomousMesh.referenceFamily;
sourceLevel=s.runMetadata.autonomousMesh.currentLevelId;
sourceFactors=family.members(sourceLevel).cellFactors;
Dx=ipm.mesh.fdMatrix(s.x,1,7);
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx, ...
    'source',s.rho*Dx','trusted',true);
rows=struct([]);
for id=levelIds(:)'
    member=family.members(id);
    assert(member.resourceAdmitted && member.qualityPassed && ...
        all(member.cellFactors>=sourceFactors));
    [candidates,trial]=ipm.remesh.plannedAxisPairs(view, ...
        struct('x',member.baseX,'y',member.baseY), ...
        s.config.scaling.transportAnchorX,policy);
    xGood=find([trial.xTrials.admissible]);
    yGood=find([trial.yTrials.admissible]);
    bestRatio=NaN;bestX=NaN;bestY=NaN;bestCore=[NaN,NaN];
    bestFront=NaN;
    for i=xGood
        for j=yGood
            xr=trial.xTrials(i);yr=trial.yTrials(j);
            ratio=max(xr.quality.maximumAdjacentCellRatio, ...
                yr.quality.maximumAdjacentCellRatio);
            if isnan(bestRatio)||ratio<bestRatio
                bestRatio=ratio;bestX=i;bestY=j;
                bestCore=[xr.coreCells,yr.coreCells];bestFront=xr.frontCells;
            end
        end
    end
    firstChanged=find(~[candidates.unchanged],1);
    selectedRatio=NaN;selectedIndices=[NaN,NaN];
    if ~isempty(firstChanged)
        c=candidates(firstChanged);
        selectedRatio=max(c.quality.x.maximumAdjacentCellRatio, ...
            c.quality.y.maximumAdjacentCellRatio);
        selectedIndices=[c.xIndex,c.yIndex];
    end
    row=struct('levelId',id,'nodeCount',member.nodeCount, ...
        'totalNodes',prod(member.nodeCount), ...
        'xAdmitted',numel(xGood),'yAdmitted',numel(yGood), ...
        'bestMinimaxRatio',bestRatio,'bestIndices',[bestX,bestY], ...
        'bestPredictedCoreCells',bestCore,'bestPredictedFrontCells',bestFront, ...
        'productionFirstChangedRatio',selectedRatio, ...
        'productionFirstChangedIndices',selectedIndices, ...
        'plannedPairCount',numel(candidates),'status',trial.status);
    if isempty(rows),rows=row;else,rows(end+1)=row;end
    fprintf('FROZEN_LEVEL tau=%.6f level=%d n=%dx%d pass=%d/%d minimax=%.6f first=%.6f\n', ...
        s.scale.canonicalTime,id,member.nodeCount(1),member.nodeCount(2), ...
        numel(xGood),numel(yGood),bestRatio,selectedRatio);
end
report=struct('checkpoint',checkpointFile,'step',s.step, ...
    'canonicalTime',s.scale.canonicalTime,'sourceLevel',sourceLevel, ...
    'sourceNodeCount',[numel(s.x),numel(s.y)],'feature',trial.feature, ...
    'levels',rows,'noLU',true,'noTransfer',true,'noPDE',true);
if nargin>=3&&~isempty(outputFile)
    assert(~isfile(outputFile));fid=fopen(outputFile,'w');assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
end
