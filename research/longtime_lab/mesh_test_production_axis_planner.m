function report=mesh_test_production_axis_planner(genericRunDirectory,outDir)
%MESH_TEST_PRODUCTION_AXIS_PLANNER No-LU migration/closure and actual replay.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
names={'ipm.remesh.plannedAxisPairs','ipm.remesh.roundedAxis','ipm.remesh.corePatchAxis', ...
    'ipm.remesh.equalizedAxis','ipm.diagnostics.meshFeatureIntervals'};
for name=names
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
files=matlab.codetools.requiredFilesAndProducts(which('ipm.remesh.plannedAxisPairs'));
assert(all(contains(files,[filesep '+ipm' filesep])) && ~any(contains(files,[filesep 'research' filesep])), ...
    'ipm:PlannerResearchDependency','Production planner closure must contain only +ipm files.');
policy=registered_policy();
report=struct('kind','production_axis_planner_migration_nolu_v1','policy',policy, ...
    'dependencyFiles',{files},'noResearchDependency',true,'noLU',true,'noPDE',true, ...
    'actualCases',struct([]),'tinyCases',struct([]),'negativeChecks',struct([]));
write_json(fullfile(outDir,'registration.json'),report);
for label={'original_t0','full_box_step3326','fresh_step3668'}
    old=load(fullfile(genericRunDirectory,[label{1} '.mat']));
    view=accepted_view(old.snapshot);
    [candidates,axisReport]=ipm.remesh.plannedAxisPairs(view,old.reference,old.anchor,policy);
    geometry=axisReport.feature;expected=old.design.feature;
    sameFeatures=isequal(geometry.coreInterval,expected.coreInterval)&& ...
        isequal(geometry.frontInterval,expected.frontInterval)&&geometry.yCoreWidth==expected.yCoreWidth&& ...
        geometry.coreCenter==expected.coreCenter;
    sameX=isequaln(strip_rows(axisReport.xTrials),strip_rows(old.design.xTrials));
    sameY=isequaln(strip_rows(axisReport.yTrials),strip_rows(old.design.yTrials));
    if old.candidate.transactionReady
        selected=find(arrayfun(@(v)isequal(v.x,old.candidate.candidateX)&&isequal(v.y,old.candidate.candidateY),candidates),1);
        selectedFound=~isempty(selected);
    else
        selected=[];selectedFound=isempty(candidates)&&strcmp(axisReport.status,old.design.status);
    end
    row=struct('label',label{1},'sameFeatureIntervals',sameFeatures,'sameAllAxisScores',sameX&&sameY, ...
        'genericSelectedPairPresentBitwise',selectedFound,'selectedProductionIndex',selected, ...
        'productionCandidateCount',numel(candidates),'status',axisReport.status, ...
        'actualCoreCells',geometry.actualCoreCells,'leftFrontCells',geometry.leftFrontCells, ...
        'passed',sameFeatures&&sameX&&sameY&&selectedFound);
    report.actualCases=append_row(report.actualCases,row);
    save(fullfile(outDir,[label{1} '.mat']),'candidates','axisReport','row','-v7.3');
    write_json(fullfile(outDir,[label{1} '.json']),row);
    assert(row.passed,'ipm:PlannerMigrationMismatch','Actual geometry/axes must preserve the research arithmetic.');
    fprintf('PRODUCTION_AXIS_CASE %s candidates=%d exact=1\n',label{1},numel(candidates));
end
% Tiny analytic input: field derivatives are recomputed by the same pure Dx.
x=linspace(-4,4,321);y=linspace(0,1,129)';[X,Y]=meshgrid(x,y);
rho=(erf((X-2)/1.8)-erf((X+2)/1.8)).*exp(-Y.^2);
snapshot=struct('rho',rho,'x',x,'y',y,'scale',struct('Cx',1,'Cy',1,'Comega',1), ...
    'canonicalTime',0,'physicalTime',0,'trusted',true,'provenance',struct('kind','analytic_tiny_noPDE_fixture'));
reference=struct('x',x,'y',y,'provenance',struct('kind','same_analytic_tiny_initial_axes'));
view=accepted_view(snapshot);feature=ipm.diagnostics.meshFeatureIntervals(view);
dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','tiny_planner_reference', ...
    'snapshotCount',1,'snapshots',snapshot,'xLimits',x([1,end]),'yLimits',y([1,end])');
dataset.snapshots.datasetIndex=1;dataset.snapshots.sourceLabel='tiny';
profiles=ipm_gridlab_feature_profiles(dataset);p=profiles.records;
core=[max(0,p.safetyCoreCenter-p.safetyCoreWidthLeft),min(x(end),p.safetyCoreCenter+p.safetyCoreWidthRight)];
front=[max(0,p.leftFrontCenter-p.leftFrontHalfWidth),min(x(end),p.leftFrontCenter+p.leftFrontHalfWidth)];
assert(isequal(core,feature.coreInterval)&&isequal(front,feature.frontInterval)&&p.yCoreWidth==feature.yCoreWidth);
[base,baseReport]=ipm.remesh.plannedAxisPairs(view,reference,1,policy);
report.tinyCases=append_row(report.tinyCases,struct('kind','analytic_tiny','featureExact',true,'candidateCount',numel(base),'status',baseReport.status,'passed',true));
for factor=[.37,4]
    u=snapshot;u.x=factor*x;u.y=factor*y;uv=accepted_view(u);
    r=reference;r.x=factor*x;r.y=factor*y;
    [unit,unitReport]=ipm.remesh.plannedAxisPairs(uv,r,factor,policy);
    assert(numel(unit)==numel(base));error=0;
    for k=1:numel(base)
        error=max(error,max([max(abs(unit(k).x/factor-base(k).x))/max(abs(x)), ...
            max(abs(unit(k).y/factor-base(k).y))/max(abs(y))]));
    end
    row=struct('kind',sprintf('unit_%g',factor),'featureExact',false,'candidateCount',numel(unit), ...
        'status',unitReport.status,'passed',error<=1e-8);
    report.tinyCases=append_row(report.tinyCases,row);assert(row.passed);
end
bad=view;bad.source(2,2)=bad.source(2,2)+1;
report.negativeChecks=append_row(report.negativeChecks,reject('source_pair',@()ipm.diagnostics.meshFeatureIntervals(bad),'ipm:AutonomousMeshSourcePair'));
bad=view;bad.trusted=false;
report.negativeChecks=append_row(report.negativeChecks,reject('untrusted',@()ipm.diagnostics.meshFeatureIntervals(bad),'ipm:AutonomousMeshView'));
badReference=reference;badReference.y(end)=2;
report.negativeChecks=append_row(report.negativeChecks,reject('domain_change',@()ipm.remesh.plannedAxisPairs(view,badReference,1,policy),'ipm:AutonomousMeshReference'));
report.allPassed=all([report.actualCases.passed])&&all([report.tinyCases.passed])&&all([report.negativeChecks.passed]);
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report','-v7.3');assert(report.allPassed);
fprintf('PRODUCTION_AXIS_PLANNER_COMPLETE noLU=1 closure=%d actual=3 tiny=3 negative=3\n',numel(files));
end

function p=registered_policy()
p=struct('targetCoreCells',[32,32],'minimumFrontCells',20, ...
    'qualityLimits',struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
    'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8,'minWeightToControlWidth',.35,'maxWeightToControlWidth',1.65), ...
    'search',struct('fineCells',[32,48,64,80,96,112,128],'roundingCells',[16,24,32,40,48], ...
    'coreFineCellFractions',[.5,.65],'ySigma',[.18,.25,.35,.5],'maximumAxisCandidates',70, ...
    'maximumPairCandidates',3,'xPadding',1.10,'yPadding',1.15,'maximumWarpFraction',.25));
end
function view=accepted_view(s)
Dx=ipm.mesh.fdMatrix(s.x,1,7);
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
end
function rows=strip_rows(rows)
for name={'construction','failure','equalizerInfo'},if isfield(rows,name{1}),rows=rmfield(rows,name{1});end,end
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function row=reject(label,f,expected)
actual='';
try
    f();
catch e
    actual=e.identifier;
end
row=struct('label',label,'expected',expected,'actual',actual,'passed',strcmp(actual,expected));assert(row.passed);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
