function report=mesh_prepare_box_space_factorial(baselineDirectory,sourceDirectory,outDir)
%MESH_PREPARE_BOX_SPACE_FACTORIAL Register three runs; pure t0 screening only.
% SourceDirectory is a frozen complete +ipm tree, identical to baseline code.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
assert(strcmp(which('ipm.solve'),fullfile(sourceDirectory,'+ipm','solve.m')));
loaded=load(fullfile(baselineDirectory,'registration.mat'),'o');original=loaded.o;
br=jsondecode(fileread(fullfile(baselineDirectory,'report.json')));
assert(br.registeredHorizonReached&&br.step==1653&&br.canonicalTime==4);
targetTime=br.physicalTime;
baseline=mesh_audit_box_space_case(fullfile(baselineDirectory,'result.mat'),br.checkpointFile,targetTime);
assert(baseline.passed);
labels={'A_H8_T32','B_H16_T32','C_H8_T42','D_H16_T42'};boxes=[8,4;16,8;8,4;16,8];targets=[32,32,42,42];
limits=struct('cropRhoRelativeL2Inf',1e-3,'cropGradientRelativeL2Inf',2e-3, ...
    'cropPeakRelative',2e-3,'cropPhysicalWidthRelative',.005,'cropRateAbsolute',1e-6,'cropRateRelative',.002);
registration=struct('kind','from_zero_box_space_factorial_v1','sourceDirectory',sourceDirectory, ...
    'baselineDirectory',baselineDirectory,'baselineAudit',baseline,'labels',{labels}, ...
    'boxes',boxes,'targetCore',targets,'physicalFinalTime',targetTime,'initialNodeCount',[321,161], ...
    'maximumTotalNodes',110000,'maximumNewRuns',3,'maximumProposalsPerInitialBox',3, ...
    'maximumAnalyticInitialDesignRounds',1,'limits',limits,'physicalWindows',[0,.8,.25;0,2,1], ...
    'physicalPeakWindow',[.1,.6],'commonQueryShape',[321,161],'physicalRateUnits','d_log_scale_d_physical_time', ...
    'initialDatum','original degenerate_primitive k8 at physical and canonical time zero', ...
    'noSharedInnerGridClaim',true,'pureBoxIsolationClaim',false,'dynamicQualification',false, ...
    'nativeInitializationQualification',false,'noLU',true,'noPDE',true);
write_json(fullfile(outDir,'registration.json'),registration);
report=registration;report.cases=struct([]);plans=struct([]);
profile clear;profile on;
for index=1:4
    o=original;o.xlim=[-boxes(index,1),boxes(index,1)];o.ymax=boxes(index,2);
    o.autonomousMesh=ipm.config.autonomousMeshPolicy(struct('version',2,'targetCoreCells',[targets(index),targets(index)], ...
        'nodeFamily',struct('maximumTotalNodes',110000)));
    policy=o.autonomousMesh;o.finalTime=Inf;o.physicalFinalTime=targetTime;
    runDir=fullfile(outDir,'runs',labels{index});o.resultFile=fullfile(runDir,'result.mat');
    o.checkpoint.file=fullfile(runDir,'checkpoint.mat');
    o.caseMetadata=struct('kind','from_original_zero_box_space_factorial_v1','factorialLabel',labels{index}, ...
        'baselineCaseId',baseline.caseId,'originalPhysicalTime',0,'originalCanonicalTime',0, ...
        'newPhysicalEpochFromCheckpoint',false,'manualMidRunMeshChanges',false,'manualTransactionTriggers',false, ...
        'pureBoxIsolationClaim',false,'maximumTotalNodes',110000, ...
        'memoryBudgetInterpretation','registered node cap only; no LU RAM model');
    config=ipm.config.resolve(o);
    assert(config.grid.nx==321&&config.grid.ny==161&&config.time.maxSteps==30000&&config.time.maxDt==.005&& ...
        strcmp(config.grid.gridMode,'uniform')&&strcmp(config.physics.initialCondition,'degenerate_primitive')&& ...
        config.physics.degeneratePower==8&&config.scaling.transportAnchorX==1);
    x=linspace(config.grid.xlim(1),config.grid.xlim(2),321);y=linspace(0,config.grid.ymax,161)';
    [source,sourceMeasure]=analytic_view(x,y,config);
    [candidates,axisReport]=ipm.remesh.plannedAxisPairs(source,struct('x',x,'y',y),1,policy);
    row=struct('label',labels{index},'reuseExistingBaseline',index==1,'initialSource',sourceMeasure, ...
        'plannerStatus',axisReport.status,'candidateCount',numel(candidates), ...
        'admittedX',nnz([axisReport.xTrials.admissible]),'admittedY',nnz([axisReport.yTrials.admissible]), ...
        'xRejections',rejection_counts(axisReport.xTrials),'yRejections',rejection_counts(axisReport.yTrials), ...
        'selectedPureCandidate',0,'candidateAudits',struct([]),'pureGeometryLaunchable',false, ...
        'nativeInitializationQualified',index==1);
    for k=1:numel(candidates)
        c=candidates(k);[view,measurement]=analytic_view(c.x,c.y,config);
        [qx,rx]=quality(c.x,policy.qualityLimits,'x');[qy,ry]=quality(c.y,policy.qualityLimits,'y');reasons=[rx,ry];
        if any(measurement.actualCoreCells<policy.transactionMinimumCoreCells),reasons{end+1}='actual_core_floor';end %#ok<AGROW>
        if measurement.leftFrontCells<policy.minimumFrontCells,reasons{end+1}='actual_front_floor';end %#ok<AGROW>
        if ~any(c.x==1)||~any(c.x==-1),reasons{end+1}='exact_anchor';end %#ok<AGROW>
        family=struct();failure=struct();familyPassed=false;
        try
            family=ipm.remesh.referenceAxisFamily(c.x,c.y,1,policy);
            familyPassed=all([family.members.qualityPassed])&&family.members(1).resourceAdmitted;
        catch e
            failure=struct('identifier',e.identifier,'message',e.message);reasons{end+1}='reference_family'; %#ok<AGROW>
        end
        audit=struct('index',k,'xIndex',c.xIndex,'yIndex',c.yIndex,'predictedCells',c.predictedCells, ...
            'actual',measurement,'actualMinusPrediction',[measurement.actualCoreCells,measurement.leftFrontCells]-c.predictedCells, ...
            'xQuality',qx,'yQuality',qy,'referenceFamilyPassed',familyPassed,'familyFailure',failure, ...
            'passed',isempty(reasons)&&familyPassed,'reasons',{reasons});
        row.candidateAudits=append_row(row.candidateAudits,audit);
        if audit.passed&&row.selectedPureCandidate==0,row.selectedPureCandidate=k;end
        save(fullfile(outDir,sprintf('%s_candidate_%d.mat',labels{index},k)), ...
            'c','view','measurement','audit','family','config','-v7.3');
    end
    row.pureGeometryLaunchable=row.selectedPureCandidate>0;
    if index==1
        cp=ipm.output.readCheckpoint(br.checkpointFile);initial=cp.payload.state.runMetadata.autonomousMesh.initialization;
        selected=candidates(row.selectedPureCandidate);
        assert(isequal(selected.x,initial.selectedBaseX)&&isequal(selected.y,initial.selectedBaseY));
        clear cp
        row.baselineInitialAxesBitwise=true;
    else
        row.baselineInitialAxesBitwise=false;
    end
    plan=struct('label',labels{index},'reuseExistingBaseline',index==1,'runDirectory',runDir, ...
        'options',o,'config',config,'pureGeometryLaunchable',row.pureGeometryLaunchable, ...
        'pureInitialAudit',row,'baselineResultFile',baseline.resultFile,'baselineCheckpointFile',baseline.checkpointFile);
    plans=append_row(plans,plan);report.cases=append_row(report.cases,row);
    save(fullfile(outDir,[labels{index},'_initial_screen.mat']),'source','sourceMeasure','candidates','axisReport','row','-v7.3');
    write_json(fullfile(outDir,[labels{index},'_initial_screen.json']),row);
    shownCore=[NaN,NaN];if row.selectedPureCandidate>0,shownCore=row.candidateAudits(row.selectedPureCandidate).actual.actualCoreCells;end
    fprintf('BOX_SPACE_PREFLIGHT label=%s target=%d actual=%s feasible=%d newRunStarted=0\n', ...
        labels{index},targets(index),jsonencode(shownCore),row.pureGeometryLaunchable);
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})));
end
report.preparationCompleted=true;report.allPureGeometryLaunchable=all([report.cases.pureGeometryLaunchable]);
save(fullfile(outDir,'registration.mat'),'registration','plans','-v7.3');
save(fullfile(outDir,'preflight_report.mat'),'report','profileInfo','-v7.3');write_json(fullfile(outDir,'preflight_report.json'),report);
end

function [view,measure]=analytic_view(x,y,config)
[X,Y]=meshgrid(x,y);ops=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y), ...
    'symmetryMode',config.physics.symmetryMode);
rho=ipm.field.initialDensity(ops,config.physics);Dx=ipm.mesh.fdMatrix(x,1,7);source=rho*Dx';
assert(all(isfinite(rho),'all')&&isreal(rho));
view=struct('rho',rho,'x',x,'y',y,'Dx',Dx,'source',source,'trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);k=config.physics.degeneratePower;
exact=sign(X).*abs(X).^(k-1)./(1+abs(X).^k+Y.^k);
window=X>=.5&X<=1.8&Y<=1;
measure=struct('actualCoreCells',feature.actualCoreCells,'leftFrontCells',feature.leftFrontCells, ...
    'coreInterval',feature.coreInterval,'frontInterval',feature.frontInterval, ...
    'yCoreWidth',feature.yCoreWidth,'coreCenter',feature.coreCenter, ...
    'wallNodalPeak',max(source(1,:)),'exactContinuousWallPeakX',(k-1)^(1/k), ...
    'exactContinuousWallPeak',(k-1)^((k-1)/k)/k, ...
    'pairedDerivativeInnerRelativeInf',max(abs(source(window)-exact(window)))/max(abs(exact(window))), ...
    'minimumDx',min(diff(x)),'minimumDy',min(diff(y)), ...
    'trustedViewMeansFinitePairedAnalyticDatumOnly',true);
end
function [q,reasons]=quality(axis,limits,label)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
bad=[q.maximumAdjacentCellRatio>limits.maxAdjacentCellRatio,q.maximumLogSpacingCurvature>limits.maxLogSpacingCurvature, ...
    q.minimumStencilRcond<limits.minStencilRcond,q.minimumQuadratureWeightRatio<limits.minQuadratureWeightRatio, ...
    q.minimumQuadratureWeightToControlWidthRatio<limits.minWeightToControlWidth, ...
    q.maximumQuadratureWeightToControlWidthRatio>limits.maxWeightToControlWidth,~q.quadratureWeightsStrictlyPositive];
names={'adjacent_ratio','spacing_curvature','stencil_rcond','global_quadrature','local_weight_min','local_weight_max','positive_weights'};
reasons=cellfun(@(s)[label,':',s],names(bad),'UniformOutput',false);
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct([]);
for k=1:numel(labels)
    rows=append_row(rows,struct('reason',labels{k},'count',nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials))));
end
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
