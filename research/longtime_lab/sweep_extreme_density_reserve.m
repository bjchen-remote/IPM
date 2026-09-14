function report=sweep_extreme_density_reserve(sourceDistribution,policyDistribution,outputFile)
%SWEEP_EXTREME_DENSITY_RESERVE Compare X core reserve at fixed node count.
assert(isfile(sourceDistribution)&&isfile(policyDistribution)&&~isfile(outputFile));
d=jsondecode(fileread(sourceDistribution));p=jsondecode(fileread(policyDistribution));
source=ipm.output.readCheckpoint(d.checkpointFile);
policySource=ipm.output.readCheckpoint(p.checkpointFile);
saved=source.payload.state;
input=policySource.payload.state.config.remesh.autonomousMesh;
input.densityVersion=1;input.search.roundingCells=[8,16,24,32,40];
basePolicy=ipm.config.autonomousMeshPolicy(input);
x=saved.x(:)';y=saved.y(:);Dx=ipm.mesh.fdMatrix(x,1,7);
view=struct('rho',saved.rho,'x',x,'y',y,'Dx',Dx, ...
    'source',saved.rho*Dx','trusted',true);
reference=struct('x',saved.baseX,'y',saved.baseY);
rows=struct([]);parameters=[32,1.10;32,1.20;32,1.30;32,1.40;32,1.50; ...
    32,1.60;32,1.70;34,1.10;36,1.10; ...
    38,1.10;40,1.10;44,1.10;48,1.10];
for i=1:size(parameters,1)
    % Exploratory geometry only: these variants are not registered policies.
    policy=basePolicy;
    policy.targetCoreCells=[parameters(i,1),32];
    policy.search.xPadding=parameters(i,2);
    [pairs,detail,search]=ipm.remesh.hierarchicalAxisPairs(view,reference, ...
        saved.config.scaling.transportAnchorX,policy, ...
        policy.nodeFamily.maximumTotalNodes);
    row=struct('targetX',parameters(i,1),'xPadding',parameters(i,2), ...
        'status',detail.status,'searchStatus',search.status, ...
        'candidateCount',numel(pairs),'predictedCells',[NaN,NaN,NaN], ...
        'maximumRatio',[NaN,NaN],'xTrial',NaN,'yTrial',NaN);
    k=find(~[pairs.unchanged],1);
    if ~isempty(k)
        row.predictedCells=pairs(k).predictedCells;
        row.maximumRatio=[pairs(k).quality.x.maximumAdjacentCellRatio, ...
                          pairs(k).quality.y.maximumAdjacentCellRatio];
        row.xTrial=pairs(k).xIndex;row.yTrial=pairs(k).yIndex;
    end
    if isempty(rows),rows=row;else,rows(end+1)=row;end
end
report=struct('kind','fixed_node_extreme_density_core_reserve_sweep_v1', ...
    'source',sourceDistribution,'tau',d.tau, ...
    'nodeCount',[numel(x),numel(y)],'rows',rows, ...
    'geometryOnly',true,'policyVariantsRegistered',false, ...
    'nativeTransferCertified',false);
fid=fopen(outputFile,'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
