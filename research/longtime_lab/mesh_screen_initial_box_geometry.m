function report=mesh_screen_initial_box_geometry(originalOptionsFile,outDir,boxes,observationRule)
%MESH_SCREEN_INITIAL_BOX_GEOMETRY Registered analytic t=0 geometry only.
% No mesh.build, flow, LU, transfer, time step, or native checkpoint creation.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
if nargin<3,boxes=[8,4;16,8;32,16];end
if nargin<4,observationRule='configuration_uniform';end
assert(any(strcmp(observationRule,{'configuration_uniform','fixed_analytic_probe_v1'})), ...
    'ipm:InitialBoxObservation','Unknown predeclared analytic observation rule.');
validateattributes(boxes,{'double'},{'real','finite','positive','2d','ncols',2,'nonempty'});
assert(size(boxes,1)<=8 && all(boxes(:,1)>1) && size(unique(boxes,'rows'),1)==size(boxes,1), ...
    'ipm:InitialBoxRegistration','Register at most eight distinct boxes containing the exact anchor.');
issues=checkcode(which(mfilename),'-id');assert(isempty(issues),jsonencode(issues));
loaded=load(originalOptionsFile,'o');original=loaded.o;
policy=ipm.config.autonomousMeshPolicy(struct('version',2, ...
    'nodeFamily',struct('maximumTotalNodes',110000)));
registration=struct('kind','primitive_k8_initial_box_pure_geometry_v1', ...
    'sourceOptionsFile',originalOptionsFile,'initialCondition','degenerate_primitive', ...
    'degeneratePower',8,'nodeCount',[321,161],'halfWidthAndYmax',boxes, ...
    'policy',policy,'maximumProposalsPerBox',3,'actualCoreFloor',policy.transactionMinimumCoreCells, ...
    'actualFrontFloor',policy.minimumFrontCells,'candidateDatum','fresh initialDensity on candidate axes', ...
    'sourceDerivative','paired seven-point fdMatrix','time',0, ...
    'analyticGradientDiagnosticOnly',true,'flowSafetyEvaluated',false, ...
    'nativeInitializationQualification',false,'dynamicBoxErrorQualification',false, ...
    'sourceFieldInterpolation',false,'noMeshBuild',true,'noLU',true,'noFlow',true,'noPDE',true);
if strcmp(observationRule,'fixed_analytic_probe_v1')
    registration.kind='primitive_k8_initial_box_pure_geometry_v2';
    registration.observation=struct('rule',observationRule,'innerHalfWidth',4, ...
        'innerSpacing',.025,'maximumTailRatio',1.05,'tailLogSlopeRampCells',20, ...
        'maximumObservationNodes',600000,'solverNodeBudgetChanged',false);
end
write_json(fullfile(outDir,'registration.json'),registration);
report=registration;report.boxes=struct([]);
profile clear;profile on;
for index=1:size(boxes,1)
    o=original;o.nx=321;o.ny=161;o.xlim=[-boxes(index,1),boxes(index,1)];o.ymax=boxes(index,2);
    o.autonomousMesh=policy;o.saveResults=false;o.makePlots=false;o.livePlot=false;o.writeVideo=false;
    config=ipm.config.resolve(o);
    assert(strcmp(config.physics.initialCondition,'degenerate_primitive')&&config.physics.degeneratePower==8&& ...
        strcmp(config.physics.symmetryMode,'double_odd_omega')&&config.scaling.transportAnchorX==1);
    assert((~isfield(config.grid,'customX')||isempty(config.grid.customX))&& ...
        (~isfield(config.grid,'customY')||isempty(config.grid.customY))&& ...
        (strcmp(config.grid.gridMode,'uniform')||all(config.grid.gridStretch==0)), ...
        'ipm:InitialBoxSourceGrid','This registered protocol must use the original uniform initial-axis branch.');
    % These are precisely mesh.build's uniform branch expressions, before
    % any operator construction. Candidate fields use the same native datum.
    x=linspace(config.grid.xlim(1),config.grid.xlim(2),config.grid.nx);
    y=linspace(0,config.grid.ymax,config.grid.ny)';
    reference=struct('x',x,'y',y);observation=struct();
    if strcmp(observationRule,'fixed_analytic_probe_v1')
        positive=observation_axis(config.grid.xlim(2));x=[-fliplr(positive(2:end)),positive];
        y=observation_axis(config.grid.ymax)';
        assert(numel(x)*numel(y)<=registration.observation.maximumObservationNodes);
        [qx,rx]=quality(x,policy.qualityLimits,'x');[qy,ry]=quality(y,policy.qualityLimits,'y');
        assert(isempty(rx)&&isempty(ry),'ipm:InitialBoxObservation','Analytic probe axes must themselves pass original axis quality.');
        observation=struct('nodeCount',[numel(x),numel(y)],'xQuality',qx,'yQuality',qy, ...
            'actualSolverReferenceNodes',[numel(reference.x),numel(reference.y)],'noInterpolation',true);
    end
    [source,sourceMeasure]=analytic_view(x,y,config);
    [candidates,axisReport]=ipm.remesh.plannedAxisPairs(source,reference,1,policy);
    row=struct('index',index,'xlim',config.grid.xlim,'ymax',config.grid.ymax, ...
        'sourceMeasurement',sourceMeasure,'plannerStatus',axisReport.status, ...
        'admittedX',nnz([axisReport.xTrials.admissible]),'admittedY',nnz([axisReport.yTrials.admissible]), ...
        'xRejections',rejection_counts(axisReport.xTrials),'yRejections',rejection_counts(axisReport.yTrials), ...
        'proposalCount',numel(candidates),'candidateAudits',struct([]), ...
        'selectedPureCandidate',0,'pureInitialSelectionFeasible',false, ...
        'nativeInitializationQualification',false,'dynamicBoxErrorQualification',false);
    if strcmp(observationRule,'fixed_analytic_probe_v1'),row.observation=observation;end
    for k=1:numel(candidates)
        candidate=candidates(k);[view,measurement]=analytic_view(candidate.x,candidate.y,config);
        limits=policy.qualityLimits;
        [qx,rx]=quality(candidate.x,limits,'x');[qy,ry]=quality(candidate.y,limits,'y');
        reasons=[rx,ry];
        if any(measurement.actualCoreCells<policy.transactionMinimumCoreCells),reasons{end+1}='actual_analytic_core_floor';end %#ok<AGROW>
        if measurement.leftFrontCells<policy.minimumFrontCells,reasons{end+1}='actual_analytic_front_floor';end %#ok<AGROW>
        if ~any(candidate.x==1)||~any(candidate.x==-1),reasons{end+1}='exact_anchor';end %#ok<AGROW>
        family=struct();familyPassed=false;failure=struct();
        try
            family=ipm.remesh.referenceAxisFamily(candidate.x,candidate.y,1,policy);
            familyPassed=all([family.members.qualityPassed])&&family.members(1).resourceAdmitted;
        catch e
            failure=struct('identifier',e.identifier,'message',e.message);
            reasons{end+1}='reference_family_rejected'; %#ok<AGROW>
        end
        passed=isempty(reasons)&&familyPassed;
        difference=struct('actualMinusSourcePredictedCells', ...
            [measurement.actualCoreCells,measurement.leftFrontCells]-candidate.predictedCells, ...
            'xCoreWidthRelativeChange',diff(measurement.coreInterval)/diff(sourceMeasure.coreInterval)-1, ...
            'yCoreWidthRelativeChange',measurement.yCoreWidth/sourceMeasure.yCoreWidth-1, ...
            'frontWidthRelativeChange',diff(measurement.frontInterval)/diff(sourceMeasure.frontInterval)-1, ...
            'coreCenterChange',measurement.coreCenter-sourceMeasure.coreCenter, ...
            'wallNodalPeakRelativeChange',measurement.wallNodalPeak/sourceMeasure.wallNodalPeak-1);
        audit=struct('candidateIndex',k,'xIndex',candidate.xIndex,'yIndex',candidate.yIndex, ...
            'predictedCellsFromCoarseSource',candidate.predictedCells,'analyticMeasurement',measurement, ...
            'sourceToCandidateDifference',difference,'xQuality',qx,'yQuality',qy, ...
            'referenceFamilyPassed',familyPassed,'referenceFamilyFailure',failure, ...
            'passed',passed,'reasons',{reasons});
        row.candidateAudits=append_row(row.candidateAudits,audit);
        if passed&&row.selectedPureCandidate==0,row.selectedPureCandidate=k;end
        save(fullfile(outDir,sprintf('box_%d_candidate_%d.mat',index,k)), ...
            'candidate','view','measurement','audit','family','config','-v7.3');
        fprintf('INITIAL_BOX_CANDIDATE H=%g k=%d predicted=%.8g/%.8g actual=%.8g/%.8g front=%.8g family=%d purePass=%d\n', ...
            boxes(index,1),k,candidate.predictedCells(1:2),measurement.actualCoreCells,measurement.leftFrontCells,familyPassed,passed);
    end
    row.pureInitialSelectionFeasible=row.selectedPureCandidate>0;
    report.boxes=append_row(report.boxes,row);
    save(fullfile(outDir,sprintf('box_%d_source_and_design.mat',index)), ...
        'config','source','sourceMeasure','reference','candidates','axisReport','row','-v7.3');
    write_json(fullfile(outDir,sprintf('box_%d.json',index)),row);
    fprintf('INITIAL_BOX_COMPLETE H=%g x=%d y=%d candidates=%d pureFeasible=%d\n', ...
        boxes(index,1),row.admittedX,row.admittedY,numel(candidates),row.pureInitialSelectionFeasible);
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})),'ipm:InitialBoxUnexpectedSolver','The registered screen is analytic geometry only.');
end
report.protocolCompleted=true;report.allBoxesPureFeasible=all([report.boxes.pureInitialSelectionFeasible]);
save(fullfile(outDir,'report.mat'),'report','profileInfo','-v7.3');write_json(fullfile(outDir,'report.json'),report);
fprintf('INITIAL_BOX_GEOMETRY_COMPLETE boxes=%d allPureFeasible=%d dynamicQualification=0 noLU=1\n', ...
    size(boxes,1),report.allBoxesPureFeasible);
end

function axis=observation_axis(endpoint)
% Fixed inner observations plus a smooth bounded-ratio tail. This is never
% a solver grid or an interpolated evolved field; it samples the t=0 datum.
inner=4;spacing=.025;ramp=20;maximumLogRatio=log(1.05);
assert(endpoint>inner,'ipm:InitialBoxObservation','The registered analytic probe requires endpoint>4.');
span=endpoint-inner;n=ceil(log1p(.05*span/spacing)/maximumLogRatio)+ramp;
indices=(0:n-1)';exponent=indices-ramp/2;
mask=indices<ramp;exponent(mask)=indices(mask).^2/(2*ramp);
assert(spacing*n<span && spacing*sum(exp(maximumLogRatio*exponent))>=span);
rate=fzero(@(value)spacing*sum(exp(value*exponent))-span,[0,maximumLogRatio]);
tail=inner+cumsum(spacing*exp(rate*exponent));tail(end)=endpoint;
axis=[linspace(0,inner,161),tail'];
assert(all(diff(axis)>0)&&any(axis==1));
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
