function report=mesh_prepare_preselected_initial_probes(screenDir,sourceOptionsFile,outDir)
%MESH_PREPARE_PRESELECTED_INITIAL_PROBES Audit analytic observations; no LU/PDE.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
issues=checkcode(which(mfilename),'-id');assert(isempty(issues),jsonencode(issues));
loaded=load(sourceOptionsFile,'o');original=loaded.o;
q=load(fullfile(screenDir,'report.mat'),'report');screen=q.report;
assert(strcmp(screen.observation.rule,'fixed_analytic_probe_v1')&& ...
    screen.observation.maximumObservationNodes==600000&&isequal(screen.nodeCount,[321,161])&& ...
    screen.policy.nodeFamily.maximumTotalNodes==110000&&~screen.nativeInitializationQualification);
registration=struct('kind','externally_preselected_analytic_axis_native_probe_registration_v1', ...
    'analyticScreen',screenDir,'sourceOptionsFile',sourceOptionsFile,'boxes',[32,16;64,32;128,64], ...
    'actualSolverNodes',[321,161],'observationNodeCap',600000,'solverNodeCap',110000, ...
    'maximumSteps',4,'canonicalFinalTime',4e-5,'maxDt',1e-5,'noProbePDEStarted',true, ...
    'productionAnalyticFallbackIntegrated',false,'initialAxesExternallyPreselected',true, ...
    'sourceFieldInterpolated',false,'originalPhysicalDatumTime',0,'nativeInitializationQualified',false);
write_json(fullfile(outDir,'registration.json'),registration);report=registration;report.boxes=struct([]);plans=struct([]);
profile clear;profile on;
for k=1:4
    a=load(fullfile(screenDir,sprintf('box_%d_source_and_design.mat',k)));
    p=a.config.remesh.autonomousMesh;
    assert(isequaln(p,screen.policy)&&isequal([numel(a.reference.x),numel(a.reference.y)],[321,161]));
    assert(isequal(a.reference.x,linspace(a.config.grid.xlim(1),a.config.grid.xlim(2),321))&& ...
        isequal(a.reference.y,linspace(0,a.config.grid.ymax,161)'));
    rho=sample(a.source.x,a.source.y,a.config.physics);Dx=ipm.mesh.fdMatrix(a.source.x,1,7);
    assert(isequal(rho,a.source.rho)&&isequal(Dx,a.source.Dx)&&isequal(rho*Dx',a.source.source));
    observationCount=[numel(a.source.x),numel(a.source.y)];assert(prod(observationCount)<=600000);
    assert(isequal(observationCount,a.row.observation.nodeCount)&& ...
        isequal(a.row.observation.actualSolverReferenceNodes,[321,161])&&prod(observationCount)>110000);
    qx=ipm.mesh.quality(a.source.x,7,ipm.mesh.quadrature(a.source.x));qy=ipm.mesh.quality(a.source.y,7,ipm.mesh.quadrature(a.source.y));
    assert(quality_pass(qx,p.qualityLimits)&&quality_pass(qy,p.qualityLimits)&& ...
        isequaln(qx,a.row.observation.xQuality)&&isequaln(qy,a.row.observation.yQuality));
    positive=linspace(0,4,161);inner=a.source.x(abs(a.source.x)<=4);
    assert(isequal(inner,[-fliplr(positive(2:end)),positive]));
    assert(qx.maximumAdjacentCellRatio<=1.05&&qy.maximumAdjacentCellRatio<=1.05);
    audits=struct([]);selected=a.row.selectedPureCandidate;
    for j=1:numel(a.candidates)
        f=fullfile(screenDir,sprintf('box_%d_candidate_%d.mat',k,j));c=load(f);
        assert(isequaln(c.candidate,a.candidates(j))&&isequal([numel(c.candidate.x),numel(c.candidate.y)],[321,161]));
        rc=sample(c.candidate.x,c.candidate.y,c.config.physics);dc=ipm.mesh.fdMatrix(c.candidate.x,1,7);
        assert(isequal(rc,c.view.rho)&&isequal(dc,c.view.Dx)&&isequal(rc*dc',c.view.source));
        feature=ipm.diagnostics.meshFeatureIntervals(c.view);
        assert(isequal(feature.actualCoreCells,c.measurement.actualCoreCells)&&isequal(feature.leftFrontCells,c.measurement.leftFrontCells));
        family=ipm.remesh.referenceAxisFamily(c.candidate.x,c.candidate.y,1,p);assert(isequaln(family,c.family));
        xq=ipm.mesh.quality(c.candidate.x,7,ipm.mesh.quadrature(c.candidate.x));yq=ipm.mesh.quality(c.candidate.y,7,ipm.mesh.quadrature(c.candidate.y));
        passed=quality_pass(xq,p.qualityLimits)&&quality_pass(yq,p.qualityLimits)&& ...
            all(feature.actualCoreCells>=p.transactionMinimumCoreCells)&&feature.leftFrontCells>=p.minimumFrontCells&& ...
            any(c.candidate.x==1)&&any(c.candidate.x==-1)&&all([family.members.qualityPassed])&&family.members(1).resourceAdmitted;
        assert(passed==c.audit.passed&&isequal([family.members.resourceAdmitted],[true,true,true,false]));
        audits=append_row(audits,struct('candidate',j,'coreCells',feature.actualCoreCells,'frontCells',feature.leftFrontCells, ...
            'xQuality',xq,'yQuality',yq,'allReferenceFamilyGeometryPassed',true,'resourceMask',[family.members.resourceAdmitted], ...
            'actualSourceAndCandidateRhoResampledExactly',true,'passed',passed));
    end
    row=struct('H',a.config.grid.xlim(2),'observationNodeCount',observationCount,'observationNodeProduct',prod(observationCount), ...
        'actualSolverNodeCount',[321,161],'observationQualityPassed',true,'sourceAnalyticRhoAndDxExact',true, ...
        'referenceUniformSolverAxesExact',true,'sourceTailMaximumRatio',[qx.maximumAdjacentCellRatio,qy.maximumAdjacentCellRatio], ...
        'selectedCandidate',selected,'candidateAudits',audits,'nativeProbeRegistered',false);
    if k<=3
        assert(selected>0&&audits(selected).passed);c=load(fullfile(screenDir,sprintf('box_%d_candidate_%d.mat',k,selected)));
        [keep,detail]=ipm.remesh.plannedAxisPairs(c.view,struct('x',c.candidate.x,'y',c.candidate.y),1,p);
        assert(~isempty(keep)&&keep(1).unchanged&&isequal(keep(1).x,c.candidate.x)&&isequal(keep(1).y,c.candidate.y));
        o=original;o.nx=321;o.ny=161;o.xlim=a.config.grid.xlim;o.ymax=a.config.grid.ymax;
        o.customX=c.candidate.x;o.customY=c.candidate.y;o.autonomousMesh=p;
        o.finalTime=4e-5;o.physicalFinalTime=Inf;o.maxDt=1e-5;o.maxSteps=4;o.outputEvery=1e-5;
        o.storeSnapshots=true;o.saveResults=true;o.makePlots=false;o.livePlot=false;o.writeVideo=false;
        label=sprintf('H%d_preselected_initial',row.H);runDir=fullfile(outDir,'runs',label);
        o.resultFile=fullfile(runDir,'result.mat');o.checkpoint=struct('enabled',true,'file',fullfile(runDir,'checkpoint.mat'),'every',1e-5,'atExit',true);
        o.caseMetadata=struct('kind','externally_preselected_analytic_axis_native_probe_v1', ...
            'candidateFile',fullfile(screenDir,sprintf('box_%d_candidate_%d.mat',k,selected)), ...
            'observationRule','fixed_analytic_probe_v1','observationNodes',observationCount, ...
            'solverNodes',[321,161],'productionAnalyticFallbackIntegrated',false, ...
            'initialAxesExternallyPreselected',true,'originalPhysicalTime',0,'originalCanonicalTime',0, ...
            'newPhysicalEpochFromCheckpoint',false,'nativeLongtimeOrBoxQualification',false);
        config=ipm.config.resolve(o);assert(isequal(config.grid.customX,c.candidate.x)&&isequal(config.grid.customY,c.candidate.y));
        plan=struct('label',label,'H',row.H,'options',o,'config',config,'runDirectory',runDir, ...
            'candidateFile',o.caseMetadata.candidateFile,'initialRho',c.view.rho,'family',c.family,'initialKeepPlan',detail);
        plans=append_row(plans,plan);row.nativeProbeRegistered=true;
    else
        assert(selected==0&&isempty(a.candidates)&&a.row.admittedX==0);
    end
    report.boxes=append_row(report.boxes,row);
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for f={'ipm.mesh.build','ipm.evolve.flow','ipm.solve','ipm.remesh.transfer','ipm.evolve.advance'},assert(~any(strcmp(names,f{1})));end
report.noLU=true;report.noPDE=true;report.preparationPassed=true;
save(fullfile(outDir,'registration.mat'),'registration','plans','-v7.3');save(fullfile(outDir,'audit.mat'),'report','profileInfo','-v7.3');
write_json(fullfile(outDir,'audit.json'),report);fprintf('PRESELECTED_INITIAL_PREPARATION_PASS selected=3 H256=capacity_rejected noLU=1\n');
end
function rho=sample(x,y,physics)
[X,Y]=meshgrid(x,y);ops=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y),'symmetryMode',physics.symmetryMode);
rho=ipm.field.initialDensity(ops,physics);
end
function passed=quality_pass(q,p)
passed=q.quadratureWeightsStrictlyPositive&&q.maximumAdjacentCellRatio<=p.maxAdjacentCellRatio&& ...
    q.maximumLogSpacingCurvature<=p.maxLogSpacingCurvature&&q.minimumStencilRcond>=p.minStencilRcond&& ...
    q.minimumQuadratureWeightRatio>=p.minQuadratureWeightRatio&& ...
    q.minimumQuadratureWeightToControlWidthRatio>=p.minWeightToControlWidth&& ...
    q.maximumQuadratureWeightToControlWidthRatio<=p.maxWeightToControlWidth;
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
