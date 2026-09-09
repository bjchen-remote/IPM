function report=mesh_test_v2_plan_dispatch(checkpointFile,outDir)
%MESH_TEST_V2_PLAN_DISPATCH No-LU controller dispatch fixtures, not PDE runs.
% Only the first case is an exact actual checkpoint observation. Synthetic
% states test family dispatch; their assigned flow counts are not PDE data.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
report=struct('kind','v2_plan_dispatch_readonly_review','noLU',true,'noPDE',true, ...
    'sourceCheckpoint',checkpointFile,'actualV1Exact',false,'syntheticCases',struct([]));
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;h=cp.payload.log.history;
state=s;state.ops=struct('x',s.x,'y',s.y,'nx',numel(s.x),'ny',numel(s.y), ...
    'baseX',s.baseX,'baseY',s.baseY,'Dx',ipm.mesh.fdMatrix(s.x,1,7),'remeshCount',s.remeshCount);
state.flow=struct('coreGridPoints',h.mesh.coreGridPoints(end), ...
    'verticalCoreGridPoints',h.mesh.verticalCoreGridPoints(end),'safetyFactor',h.mesh.safetyFactor(end));
[new,newPlan]=ipm.evolve.planAutonomousMesh(state,false);
[old,oldPlan]=legacy_planAutonomousMesh(state,false);
assert(isequaln(new,old)&&isequaln(newPlan,oldPlan));report.actualV1Exact=true;
policy=ipm.config.autonomousMeshPolicy(struct('version',2,'nodeFamily',struct('maximumTotalNodes',110000)));
x=linspace(-2,2,129);y=linspace(0,2,65)';family=ipm.remesh.referenceAxisFamily(x,y,1,policy);
orders={[1,3,2,4],[2,4],[3,4],4};
for level=1:4
    member=family.members(level);[X,Y]=meshgrid(member.baseX,member.baseY);
    rho=(erf((X-.8)/.6)-erf((X+.8)/.6)).*exp(-4*Y.^2);
    state=struct('rho',rho,'config',s.config,'scale',s.scale,'normalizedTime',.1,'step',1, ...
        'ops',struct('x',member.baseX,'y',member.baseY,'nx',member.nodeCount(1),'ny',member.nodeCount(2), ...
        'baseX',member.baseX,'baseY',member.baseY,'Dx',ipm.mesh.fdMatrix(member.baseX,1,7),'remeshCount',0), ...
        'flow',struct('coreGridPoints',25,'verticalCoreGridPoints',32,'safetyFactor',.3));
    state.config.remesh.autonomousMesh=policy;state.config.grid.nx=numel(x);state.config.grid.ny=numel(y);
    state.config.scaling.transportAnchorX=1;state.scale.canonicalTime=.1;state.scale.physicalTime=.05;
    memory=struct('version',2,'policy',policy,'transactions',struct([]), ...
        'cumulativeAbsolutePeakJump',0,'window',struct('time',[],'step',[], ...
        'remeshCount',[],'coreCells',zeros(0,2),'safety',[],'nodeCount',zeros(0,2),'levelId',[]), ...
        'initialization',struct(),'lastDecision',struct(),'lastFailure',struct(), ...
        'referenceFamily',family,'currentLevelId',level);
    state.runMetadata=struct('caseId','synthetic_geometry_dispatch_only','autonomousMesh',memory);
    before=state;[after,plan]=ipm.evolve.planAutonomousMesh(state,false);
    assert(isequaln(state,before)&&isequaln(after.rho,state.rho)&&isequaln(after.ops,state.ops)&& ...
        isequaln(after.scale,state.scale));
    assert(isequal(plan.axisReport.eligibleLevelOrder,orders{level})&&plan.requested&& ...
        strcmp(plan.reason,'current_core_trigger')&&plan.sourceLevelId==level);
    assert(all([plan.axisReport.levels.proposedPairs]<=3)&&numel(plan.candidates)<=3*numel(orders{level}));
    ids=[];
    for k=1:numel(plan.candidates)
        candidate=plan.candidates(k);target=family.members(candidate.nodeFamilyIndex);
        assert(~candidate.unchanged&&isequal([numel(candidate.x),numel(candidate.y)],target.nodeCount));
        assert(all(target.cellFactors>=member.cellFactors)&&target.resourceAdmitted);
        ids(end+1)=candidate.nodeFamilyIndex; %#ok<AGROW>
    end
    positions=arrayfun(@(id)find(orders{level}==id),ids);assert(all(diff(positions)>=0));
    assert(isequaln(after.runMetadata.autonomousMesh.lastDecision,rmfield(plan,{'candidates','axisReport'})));
    row=struct('sourceLevel',level,'eligibleOrder',orders{level},'candidateLevels',ids, ...
        'proposedPerLevel',[plan.axisReport.levels.proposedPairs],'passed',true, ...
        'nativeStateQualificationClaim',false);
    report.syntheticCases=append_row(report.syntheticCases,row);
    save(fullfile(outDir,sprintf('level_%d.mat',level)),'plan','row','-v7.3');
end
report.allPassed=report.actualV1Exact&&all([report.syntheticCases.passed]);
save(fullfile(outDir,'report.mat'),'report','-v7.3');
fid=fopen(fullfile(outDir,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('PLAN_DISPATCH_COMPLETE actualV1Exact=1 syntheticLevels=4 noLU=1\n');
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
