function report=ipm_perflab_preflight_fresh_campaign(plan,reportDirectory)
%IPM_PERFLAB_PREFLIGHT_FRESH_CAMPAIGN No-LU launch prerequisites, not acceptance.
% Required: sourceRoot, outputDirectory, startKind, entryPoint. time_zero
% requires explicit customX/customY in options. fresh_checkpoint requires
% bootstrapResultFile, initialDataFile and qualificationFile. Optional
% axisProbes are {label,functionName,arguments,geometry}, calling only the
% three explicit pure-axis factories below. No solver/restore is executed.
assert(maxNumCompThreads==10,'ipm:PreflightThreads','Native signatures require the registered ten threads.');
assert(java.io.File(reportDirectory).isAbsolute(),'ipm:PreflightReportPath','Use an absolute new report directory.');
reportDirectory=canonical(reportDirectory);
assert(isstruct(plan)&&isscalar(plan)&&~isfolder(reportDirectory));mkdir(reportDirectory);
save(fullfile(reportDirectory,'registered_plan.mat'),'plan','-v7.3');
oldPath=path;oldDirectory=pwd;cleanup=onCleanup(@()restore_environment(oldPath,oldDirectory));
report=struct('kind','no_LU_campaign_launch_preflight','passed',false,'poissonBuilds',0, ...
    'poissonSolves',0,'pdeSteps',0,'checkpointWrites',0,'nativeTransferTests',0, ...
    'longTimeStabilityEstablished',false,'originalTimeZeroInputRecognized',false, ...
    'requiresActualInitialFlowAndLongTimeValidation',true,'failures',{{}});
try
    required={'sourceRoot','outputDirectory','startKind','entryPoint'};assert(all(isfield(plan,required)));
    sourceRoot=canonical(plan.sourceRoot);assert(isfolder(sourceRoot));cd(sourceRoot);addpath(sourceRoot);
    for suffix={'research/longtime_lab','research/grid_lab','research/experiments'}
        directory=fullfile(sourceRoot,suffix{1});if isfolder(directory),addpath(directory);end
    end
    report.sourceRoot=sourceRoot;report.entryPoint=plan.entryPoint;report.startKind=plan.startKind;
    entries={plan.entryPoint,'ipm.config.resolve','ipm.field.initialDensity','ipm.mesh.quality','ipm.mesh.quadrature'};
    if isfield(plan,'additionalEntries'),entries=[entries,plan.additionalEntries(:)'];end
    if isfield(plan,'callbackDependencies'),entries=[entries,plan.callbackDependencies(:)'];end
    if isfield(plan,'axisProbes')&&~isempty(plan.axisProbes)
        entries=[entries,{plan.axisProbes.functionName}];
    end
    report.dependencies=ipm_perflab_preflight_dependencies(sourceRoot,unique(entries,'stable'));
    assert(report.dependencies.passed,'ipm:PreflightDependencies','Registered frozen dependency closure is incomplete.');
    assert(java.io.File(plan.outputDirectory).isAbsolute(),'ipm:PreflightOutputPath','Use an absolute campaign output directory.');
    output=canonical(plan.outputDirectory);
    assert(~isfolder(output)&&~isfile(output),'ipm:PreflightOutputExists','Campaign output must be new.');
    parent=fileparts(output);while ~isfolder(parent),parent=fileparts(parent);end
    probe=[tempname(parent),'.preflight'];fid=fopen(probe,'w');assert(fid>=0,'ipm:PreflightOutputWrite','Cannot write the campaign output ancestor.');
    closeProbe=onCleanup(@()delete_owned_probe(probe));fwrite(fid,'preflight','char');fclose(fid);delete(probe);clear closeProbe
    report.output=struct('requested',output,'testedExistingAncestor',parent,'newDirectoryAvailable',true, ...
        'writeProbePassed',true,'campaignDirectoryCreated',false);
    if strcmp(plan.startKind,'time_zero')
        assert(strcmp(plan.entryPoint,'ipm.solve'),'ipm:PreflightTimeZeroEntry','The late-fresh driver cannot establish original t=0 lineage.');
        assert(isfield(plan,'options')&&isstruct(plan.options));
        config=ipm.config.resolve(plan.options);
        assert(~isempty(config.grid.customX)&&~isempty(config.grid.customY), ...
            'ipm:PreflightExplicitAxes','Supply the actual precomputed initial axes; preflight never calls mesh.build to obtain them.');
        x=config.grid.customX(:)';y=config.grid.customY(:);
        initial=config.physics.initialCondition;
        assert(~isa(initial,'function_handle'),'ipm:PreflightCallbackScope', ...
            'Time-zero sample checks execute predefined native initial data only; arbitrary callbacks need a separate no-LU review.');
        assert(numel(x)*numel(y)<=2e6,'ipm:PreflightResourceBound','Initial sample check exceeds this bounded no-LU preflight.');
        [X,Y]=meshgrid(x,y);sampleOps=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y),'symmetryMode',config.physics.symmetryMode);
        rho=ipm.field.initialDensity(sampleOps,config.physics);
        assert(isequal(size(rho),[numel(y),numel(x)])&&isreal(rho)&&all(isfinite(rho),'all'), ...
            'ipm:PreflightInitialData','Initial density samples are invalid.');
        report.initialData=struct('source','original configured time-zero datum','shape',size(rho), ...
            'finite',true,'range',[min(rho,[],'all'),max(rho,[],'all')],'flowNotEvaluated',true);
        report.originalTimeZeroInputRecognized=true;clear X Y rho
    elseif strcmp(plan.startKind,'fresh_checkpoint')
        assert(strcmp(plan.entryPoint,'run_fresh_profile_campaign'));
        [review,cp,contract]=fresh_profile_campaign_review(plan.bootstrapResultFile,plan.initialDataFile);
        qualificationFile=plan.qualificationFile;registeredQualification=jsondecode(fileread(qualificationFile));
        if strcmp(registeredQualification.kind,'fresh_profile_campaign_bootstrap_chain_qualification_v1')
            % Admission rechecks the native chain and writes evidence. Keep
            % this audit inside the new preflight directory, never in the
            % caller's existing qualification/output directory.
            auditQualification=registeredQualification;
            auditQualification.newChainEvidenceDirectory=fullfile(reportDirectory,'native_chain_admission');
            qualificationFile=fullfile(reportDirectory,'qualification_for_readonly_recheck.json');
            write_json(qualificationFile,auditQualification);
            report.qualificationRecheck=struct('originalFile',plan.qualificationFile, ...
                'derivedFile',qualificationFile,'onlyChangedField','newChainEvidenceDirectory', ...
                'originalFileRewritten',false,'nativeDataWritten',false);
        end
        qualification=fresh_profile_campaign_admit(qualificationFile,plan.bootstrapResultFile,plan.initialDataFile,cp);
        assert(~isfield(qualification,'testFixtureOnly')||~qualification.testFixtureOnly);
        assert(all(review.coreCells>=20),'ipm:PreflightBootstrap','Current bootstrap fails the original core20 floor.');
        config=cp.payload.state.config;x=cp.payload.state.x;y=cp.payload.state.y;
        report.bootstrap=review;save(fullfile(reportDirectory,'native_contract.mat'),'contract','-v7.3');
        report.contract=struct('file',fullfile(reportDirectory,'native_contract.mat'),'caseId',contract.caseId, ...
            'initialNodeCount',[numel(contract.initialX),numel(contract.initialY)], ...
            'fullContractPreservedInMat',true,'pairedAndReferencesChecked',true);
        report.initialData=struct('source','late fresh IVP; never relabelled original time zero','nativeReadOnly',true,'operatorRestored',false);
        clear cp
    else
        error('ipm:PreflightStartKind','Unknown startup kind.');
    end
    assert(config.schemaVersion==4&&strcmp(config.scaling.scalingContract,'exact_gauge_no_feedback_v1'));
    assert(strcmp(config.transport.spatialDiscretization,'high_order'), ...
        'ipm:PreflightScope','This campaign preflight is scoped to the maintained high_order policy.');
    report.config=struct('resolved',true,'schemaVersion',config.schemaVersion,'spatialDiscretization',config.transport.spatialDiscretization, ...
        'transportScheme',config.transport.transportScheme,'timeIntegrator',config.time.timeIntegrator, ...
        'adaptiveRemesh',config.remesh.adaptiveRemesh,'initialAnalyticRemesh',config.remesh.initialAnalyticRemesh);
    report.grid=pair_quality(x,y,config);
    assert(report.grid.passed,'ipm:PreflightAxisQuality','Initial/current axes fail the original campaign geometric gates.');
    factories=struct('label',{},'functionName',{},'geometry',{},'passed',{},'quality',{},'failure',{});
    if isfield(plan,'axisProbes')
        allowed={'mesh_core_patch_axis','mesh_general_rounded_axis','ipm_gridlab_equalize_axis'};
        for k=1:numel(plan.axisProbes)
            p=plan.axisProbes(k);q=struct('label',p.label,'functionName',p.functionName,'geometry',p.geometry,'passed',false,'quality',[],'failure',[]);
            try
                assert(ismember(p.functionName,allowed),'ipm:PreflightFactoryScope','Only registered pure-axis factories may be executed.');
                axis=feval(p.functionName,p.arguments{:});
                if strcmp(p.geometry,'x')
                    q.quality=pair_quality(axis,y,config);
                elseif strcmp(p.geometry,'y')
                    q.quality=pair_quality(x,axis,config);
                else
                    error('ipm:PreflightFactoryGeometry','Factory geometry must be x or y.');
                end
                q.passed=q.quality.passed;
            catch exception
                q.failure=error_view(exception);
            end
            factories(end+1)=q; %#ok<AGROW>
        end
    end
    report.axisFactories=factories;
    report.axisFactoryAssurance='Actual registered constructors and full axis-quality gates only; no field transfer, peak jump, front/core or native transaction acceptance is inferred.';
    if ~isempty(factories)
        geometries=unique({factories.geometry});
        for k=1:numel(geometries)
            selected=strcmp({factories.geometry},geometries{k});
            assert(any([factories(selected).passed]),'ipm:PreflightAxisFamiliesExhausted', ...
                'All registered %s-axis factory probes failed; do not begin an unattended batch.',geometries{k});
        end
    end
    report.passed=true;
catch exception
    report.failure=error_view(exception);report.failures={report.failure};
end
report.interpretation='Launch prerequisites only. Complete registered dependency closure, resolved configuration, actual axes and write access do not certify future candidate feasibility, trustworthy initial RHS, or continuous long-time stability.';
save(fullfile(reportDirectory,'report.mat'),'report','-v7.3');write_json(fullfile(reportDirectory,'report.json'),report);
startKind='unspecified';if isfield(plan,'startKind'),startKind=plan.startKind;end
fprintf('CAMPAIGN_PREFLIGHT passed=%d start=%s noLU=1\n',report.passed,startKind);
end
function p=pair_quality(x,y,config)
x=x(:)';y=y(:);anchor=config.scaling.transportAnchorX;
geometry=all(isfinite(x))&&all(isfinite(y))&&isreal(x)&&isreal(y)&&all(diff(x)>0)&&all(diff(y)>0)&& ...
    numel(x)==config.grid.nx&&numel(y)==config.grid.ny&&isequal(x([1,end]),config.grid.xlim)&& ...
    y(1)==0&&y(end)==config.grid.ymax&&isequal(x,-fliplr(x))&&any(x==anchor)&&any(x==-anchor);
qx=axis_quality(x);qy=axis_quality(y);
p=struct('nodeCount',[numel(x),numel(y)],'geometryAndExactAnchors',geometry,'xQuality',qx,'yQuality',qy, ...
    'passed',geometry&&qx.passed&&qy.passed);
end
function q=axis_quality(axis)
weights=ipm.mesh.quadrature(axis);v=ipm.mesh.quality(axis,7,weights);
q=struct('diagnostic',v,'limits',struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
    'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8,'minWeightToControlWidth',.35,'maxWeightToControlWidth',1.65));
q.passed=v.maximumAdjacentCellRatio<=1.08&&v.maximumLogSpacingCurvature<=.01&&v.minimumStencilRcond>=1e-9&& ...
    v.minimumQuadratureWeightRatio>=1e-8&&v.quadratureWeightsStrictlyPositive&& ...
    v.minimumQuadratureWeightToControlWidthRatio>=.35&&v.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function e=error_view(exception)
e=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
end
function restore_environment(oldPath,oldDirectory)
path(oldPath);cd(oldDirectory);
end
function delete_owned_probe(file)
if isfile(file),delete(file);end
end
function value=canonical(value)
value=char(java.io.File(value).getCanonicalPath());
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
