function report=mesh_test_anchor_grid_tiny(outputDirectory,tinyFreshReport)
%MESH_TEST_ANCHOR_GRID_TINY Unit covariance plus native transaction contract.
% Synthetic paired arrays exercise only derivatives/quadrature/transfers;
% native identity transactions below are restricted to <=20000 field entries.
assert(~isfolder(outputDirectory));mkdir(outputDirectory);
originalThreads=maxNumCompThreads(10);cleanup=onCleanup(@()maxNumCompThreads(originalThreads));
controls=struct('positiveFineCells',16,'anchorFineCellFraction',.5, ...
    'roundingCells',8,'targetXCoreCells',24,'targetYCoreCells',12, ...
    'minimumXFrontCells',20,'verticalResolutionPadding',1.15);
scales=[1,.25,4,.37];cases=cell(size(scales));
for k=1:numel(scales)
    scale=scales(k);x=scale*linspace(-4,4,257);y=scale*linspace(0,4,129)';
    [X,Y]=meshgrid(x/scale,y/scale);
    rho=(1-exp(-X.^2/2)).*exp(-Y.^2/(2*.6^2));
    snapshot=struct('datasetIndex',1,'sourceLabel','analytic tiny covariance fixture', ...
        'trusted',true,'rho',rho,'x',x,'y',y,'scale',struct('Cx',1,'Cy',1,'Comega',1), ...
        'canonicalTime',0,'physicalTime',0);
    dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset', ...
        'tag','analytic_unit_covariance_not_a_solver_history','snapshotCount',1, ...
        'snapshots',snapshot,'xLimits',x([1,end]),'yLimits',y([1,end])');
    reference=struct('kind','analytic_fixture_initial_axes','x',x,'y',y, ...
        'provenance',struct('syntheticFixture',true,'scale',scale));
    cases{k}=mesh_anchor_grid_candidate(dataset,reference,scale, ...
        fullfile(outputDirectory,sprintf('scale_%d.mat',k)),controls);
    assert(cases{k}.pairScore.admissible,'The covariance fixture must receive a full pair score.');
end
axisErrors=zeros(numel(scales)-1,2);qualityErrors=zeros(numel(scales)-1,1);metricErrors=zeros(numel(scales)-1,1);
base=cases{1};
for k=2:numel(scales)
    axisErrors(k-1,:)=[max(abs(cases{k}.candidateX/scales(k)-base.candidateX)), ...
        max(abs(cases{k}.candidateY/scales(k)-base.candidateY))];
    qualityNames={'maximumAdjacentCellRatio','maximumLogSpacingCurvature', ...
        'minimumStencilRcond','minimumQuadratureWeightRatio', ...
        'minimumQuadratureWeightToControlWidthRatio','maximumQuadratureWeightToControlWidthRatio'};
    for axisName={'x','y'}
        a=base.pairScore.quality.(axisName{1});b=cases{k}.pairScore.quality.(axisName{1});
        for n=qualityNames
            qualityErrors(k-1)=max(qualityErrors(k-1),abs(a.(n{1})-b.(n{1})));
        end
    end
    metricNames={'fieldRelativeL2','rhoXRelativeL2','rhoYRelativeL2', ...
        'rhoXMaximumRelativeChange','conservationRelativeDefect','relativeRangeViolation'};
    for n=metricNames
        metricErrors(k-1)=max(metricErrors(k-1),abs( ...
            base.pairScore.aggregate.worst.(n{1})-cases{k}.pairScore.aggregate.worst.(n{1})));
    end
    assert(cases{k}.anchorExact && max(axisErrors(k-1,:))<1e-12 && ...
        qualityErrors(k-1)<1e-10 && metricErrors(k-1)<1e-10);
end
fresh=jsondecode(fileread(tinyFreshReport));cp=ipm.output.readCheckpoint(fresh.checkpointFile);
s=cp.payload.state;assert(numel(s.rho)<=20000 && s.config.scaling.transportAnchorX~=1);
limits=struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
    'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8);
options=struct('minimumXCoreCells',2,'minimumYCoreCells',2, ...
    'maximumRhoXRelativeChange',2e-3,'maximumMassRelativeDefect',5e-12, ...
    'maximumRelativeRangeViolation',2e-4,'meshLimits',limits,'tag','arbitrary_anchor_tiny_identity');
[identity,audit]=ipm_gridlab_regrid_checkpoint(cp,s.x,s.y,options);
assert(audit.passed && identity.payload.state.config.scaling.transportAnchorX== ...
    s.config.scaling.transportAnchorX && isequal(identity.payload.state.rescaling,s.rescaling));
clear identity
newCountRejected=false;
try
    ipm_gridlab_regrid_checkpoint(cp,linspace(s.x(1),s.x(end),numel(s.x)+2),s.y,options);
catch exception
    newCountRejected=strcmp(exception.identifier,'ipm:gridlab:RegridNodeCount');
    assert(newCountRejected,exception.message);
end
assert(newCountRejected);
report=struct('scales',scales,'maximumNormalizedAxisErrors',axisErrors, ...
    'dimensionlessQualityErrors',qualityErrors,'relativePairMetricErrors',metricErrors, ...
    'nativeIdentityAnchor',s.config.scaling.transportAnchorX, ...
    'nativeIdentityPassed',audit.passed,'nativeReferencesPreserved',true, ...
    'nativeChangedCountRejected',newCountRejected,'productionHardMeshLimits',limits, ...
    'tinyOnlyCoreFloor',2,'noPDEAdvanced',true,'q512OperatorsBuilt',false,'passed',true);
save(fullfile(outputDirectory,'tiny_report.mat'),'report','-v7.3');
fid=fopen(fullfile(outputDirectory,'tiny_report.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
fprintf('ANCHOR_GRID_TINY_PASS covariance=1 nonUnitAnchor=1 countContract=1 noPDE=1\n');
end
