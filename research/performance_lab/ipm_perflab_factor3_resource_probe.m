function report=ipm_perflab_factor3_resource_probe(bundleFile,outDir,allowLarge)
% One operator-only factor + frozen RHS measurement; never a native state.
% Default refuses LU. Root must grant the externally watched resource slot.
if nargin<3,allowLarge=false;end
assert(islogical(allowLarge)&&isscalar(allowLarge)&&maxNumCompThreads==10);
assert(~isfolder(outDir));mkdir(outDir);
try
    d=load(bundleFile,'bundle');b=d.bundle;clear d;
    assert(strcmp(b.kind,'operator_only_factor3_frozen_input_v1')&&b.readyForOperatorResourceProbe&& ...
        ~b.nativeTransferQualified&&~b.LUResourceQualified&&~b.PDEQualified&& ...
        isequal(size(b.rhoFrozen),[321,961])&&all(isfinite(b.rhoFrozen),'all'));
    report=struct('kind','factor3_operator_resource_probe_v1','bundleFile',bundleFile, ...
        'sourceFile',b.sourceFile,'sourceStep',b.sourceStep,'sourceCanonicalTime',b.sourceScale.canonicalTime, ...
        'sourcePhysicalTime',b.sourceScale.physicalTime,'nodeCount',[961,321], ...
        'candidateParameters',b.candidateParameters,'allowLarge',allowLarge, ...
        'nativePolicyVersionThreeImplemented',false,'nativeTransactionPerformed',false, ...
        'checkpointWritten',false,'PDETimeAdvanced',false,'operatorBuilt',false,'flowEvaluated',false, ...
        'interpretation','Frozen operator-only input: constant-conservative six-point X/Y interpolation; not the native wall-trace-preserving transfer.', ...
        'externalMemoryMeasurementRequired',true);
    write_json(fullfile(outDir,'registration.json'),report);
    if ~allowLarge
        report.status='prepared_only_no_lu';write_json(fullfile(outDir,'report.json'),report);return
    end
    assert(isequal(getenv('IPM_FACTOR3_RESOURCE_WINDOW'),'root_authorized_watched_single_factor'), ...
        'ipm:ResearchFactor3Window','Run only through the root-authorized resource watchdog.');
    cp=ipm.output.readCheckpoint(b.sourceFile);s=cp.payload.state;
    assert(s.step==4545&&s.scale.canonicalTime==8&&isequaln(s.scale,b.sourceScale)&& ...
        isequaln(s.config,b.sourceConfig)&&isequaln(s.rescaling,b.sourceReferences)&& ...
        all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
    [draft,family]=ipm_perflab_factor3_family(s.config.remesh.autonomousMesh,s.runMetadata.autonomousMesh.referenceFamily);
    assert(isequaln(draft,b.futureDraftPolicy)&&isequaln(family,b.futureReferenceFamily));
    assert(strcmp(s.config.scaling.rescalingMode,'dynamic')&&strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic')&& ...
        strcmp(s.config.transport.spatialDiscretization,'high_order')&&strcmp(s.config.transport.transportScheme,'weno5_fd')&& ...
        strcmp(s.config.scaling.scalingContract,'exact_gauge_no_feedback_v1'));
    % This lower-level operator override is explicit experimental data. It
    % does not change the immutable CP config or pass any native commit path.
    grid=s.config.grid;grid.nx=961;grid.ny=321;grid.customX=b.candidate.x;grid.customY=b.candidate.y;
    rho=b.rhoFrozen;scale=s.scale;sourceReferences=s.rescaling;sourceConfig=s.config;
    fprintf('FACTOR3_PHASE begin_single_factor\n');timer=tic;
    ops=ipm.mesh.build(s.config,grid);report.factorAndGridBuildSeconds=toc(timer);report.operatorBuilt=true;
    assert(isa(ops.poisson,'decomposition')&&isequal(ops.x,b.candidate.x)&&isequal(ops.y,b.candidate.y));
    ops.baseX=family.members(5).baseX;ops.baseY=family.members(5).baseY;
    refs=s.rescaling;[~,refs.originIndex]=min(abs(ops.x));[~,refs.pinIndex]=min(abs(ops.x-refs.pinX));
    ops.rescaling=refs;ops.remeshCount=s.remeshCount;
    report.factorClass=class(ops.poisson);report.matrixSize=size(ops.A);report.matrixNnz=nnz(ops.A);
    report.actualRuntimeReferences=refs;report.sourceConfigKeptExact=isequaln(s.config,b.sourceConfig);
    assert(isequaln(rmfield(refs,{'originIndex','pinIndex'}),rmfield(sourceReferences,{'originIndex','pinIndex'})));
    fprintf('FACTOR3_PHASE single_original_flow\n');timer=tic;
    [rhoRate,flow]=ipm.evolve.flow(rho,ops,scale);report.flowSeconds=toc(timer);report.flowEvaluated=true;
    assert(all(isfinite(rhoRate),'all')&&ipm.evolve.isFinite(struct('rho',rho,'flow',flow,'scale',scale)));
    report.coreCells=[flow.coreGridPoints,flow.verticalCoreGridPoints];report.safety=flow.safetyFactor;
    report.cL=flow.c_l;report.cOmega=flow.c_omega;report.cR=flow.c_r;
    report.poissonResidual=flow.poissonResidual;
    report.poissonResidualIsForwardAccuracyCertificate=false;
    report.operatorFlowFinite=true;report.singleFactorConstructed=true;report.oldFactorEverConstructed=false;
    clear ops cp s
    fprintf('FACTOR3_PHASE factor_released_before_output\n');
    report.factorReleasedBeforeOutput=true;report.status='operator_measurement_completed_resource_verdict_external';
    candidate=b.candidate;futureDraftPolicy=draft;futureReferenceFamily=family;
    save(fullfile(outDir,'frozen_operator_output.mat'),'rho','rhoRate','flow','scale', ...
        'candidate','futureDraftPolicy','futureReferenceFamily','sourceConfig','sourceReferences','-v7.3');
    save(fullfile(outDir,'report.mat'),'report','-v7.3');write_json(fullfile(outDir,'report.json'),report);
    fprintf('FACTOR3_OPERATOR_COMPLETE nodes=961x321 oneFactor=1 flow=1 nativeCommit=0 pdeSteps=0\n');
catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(outDir,'failure.mat'),'failure','exception','-v7.3');write_json(fullfile(outDir,'failure.json'),failure);rethrow(exception)
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
