function report=ipm_accellab_H64_directional_resource_probe(bundleFile,outDir,executeLU)
% Explicitly dormant by default. One new-grid factor/flow, no native commit.
if nargin<3,executeLU=false;end
assert(islogical(executeLU)&&isscalar(executeLU)&&maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
try
 d=load(bundleFile,'bundle');b=d.bundle;clear d;
 assert(strcmp(b.kind,'actual_H64_directional_operator_only_input_v1')&& ...
  any(b.targetMemberId==[7,8])&&b.readyForOperatorResourceProbe&& ...
  ~b.nativeTransferQualified&&~b.LUResourceQualified&&~b.PDEQualified);
 member=b.futureReferenceFamily.members(b.targetMemberId);N=member.nodeCount;
 assert(isequal(size(b.rhoFrozen),fliplr(N))&&all(isfinite(b.rhoFrozen),'all'));
 report=struct('kind','H64_directional_single_factor_resource_probe_v1','bundleFile',bundleFile, ...
  'sourceFile',b.sourceFile,'sourceStep',b.sourceStep,'sourceCanonicalTime',b.sourceScale.canonicalTime, ...
  'sourcePhysicalTime',b.sourceScale.physicalTime,'nodeCount',N,'targetMemberId',b.targetMemberId, ...
  'cellFactors',member.cellFactors,'box',[b.candidate.x([1,end]),b.candidate.y([1,end])'], ...
  'runtimeSource',fileparts(fileparts(which('ipm.solve'))),'executeLU',executeLU, ...
  'actualLUExecuted',false,'flowEvaluated',false,'operatorFlowFinite',false,'operatorCleared',false, ...
  'nativeTransferQualified',false,'nativeTransactionPerformed',false,'checkpointWritten',false,'PDETimeAdvanced',false, ...
  'resourceQualified',false,'externalWatchdogVerdictRequired',true, ...
  'interpretation','Six-point constant-conservative frozen input, not native Y-bubble transfer. Operator-only measurement is not a trajectory or a future LU-memory bound.');
 write_json(fullfile(outDir,'registration.json'),report);
 if ~executeLU
  report.status='prepared_only_no_LU';write_json(fullfile(outDir,'report.json'),report);return
 end
 assert(strcmp(getenv('IPM_H64_DIRECTIONAL_RESOURCE_WINDOW'),'root_authorized_watched_single_factor'), ...
  'ipm:ResearchH64ResourceWindow','Explicit root resource window and external watchdog required.');
 cp=ipm.output.readCheckpoint(b.sourceFile);s=cp.payload.state;
 assert(s.step==b.sourceStep&&isequaln(s.scale,b.sourceScale)&&isequaln(s.config,b.sourceConfig)&& ...
  isequaln(s.rescaling,b.sourceReferences)&&all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
 m=s.runMetadata.autonomousMesh;draft=s.config.remesh.autonomousMesh;draft.version=3;draft.nodeFamily=struct('maximumTotalNodes',310000);
 draft=ipm.config.autonomousMeshPolicy(draft);
 family=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX,m.initialization.selectedBaseY,m.referenceFamily.anchor,draft);
 assert(isequaln(draft,b.futureDraftPolicy)&&isequaln(family,b.futureReferenceFamily)&& ...
  isequaln(family.members(1:4),m.referenceFamily.members));
 assert(strcmp(s.config.scaling.rescalingMode,'dynamic')&&strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic')&& ...
  strcmp(s.config.transport.spatialDiscretization,'high_order')&&strcmp(s.config.transport.transportScheme,'weno5_fd')&& ...
  strcmp(s.config.scaling.scalingContract,'exact_gauge_no_feedback_v1'));
 grid=s.config.grid;grid.nx=N(1);grid.ny=N(2);grid.customX=b.candidate.x;grid.customY=b.candidate.y;
 rho=b.rhoFrozen;scale=s.scale;sourceConfig=s.config;sourceReferences=s.rescaling;
 fprintf('H64_DIRECTIONAL_PHASE begin_single_factor member=%d\n',b.targetMemberId);timer=tic;
 ops=ipm.mesh.build(s.config,grid);report.buildSeconds=toc(timer);report.actualLUExecuted=true;
 assert(isa(ops.poisson,'decomposition')&&isequal(ops.x,b.candidate.x)&&isequal(ops.y,b.candidate.y));
 ops.baseX=member.baseX;ops.baseY=member.baseY;refs=s.rescaling;
 [~,refs.originIndex]=min(abs(ops.x));[~,refs.pinIndex]=min(abs(ops.x-refs.pinX));
 ops.rescaling=refs;ops.remeshCount=s.remeshCount;
 assert(isequaln(rmfield(refs,{'originIndex','pinIndex'}),rmfield(sourceReferences,{'originIndex','pinIndex'})));
 report.matrixSize=size(ops.A);report.matrixNnz=nnz(ops.A);report.actualRuntimeReferences=refs;
 fprintf('H64_DIRECTIONAL_PHASE single_original_flow\n');timer=tic;
 [rhoRate,flow]=ipm.evolve.flow(rho,ops,scale);report.flowSeconds=toc(timer);report.flowEvaluated=true;
 report.operatorFlowFinite=all(isfinite(rhoRate),'all')&&ipm.evolve.isFinite(struct('rho',rho,'flow',flow,'scale',scale));
 assert(report.operatorFlowFinite);report.coreCells=[flow.coreGridPoints,flow.verticalCoreGridPoints];
 report.safety=flow.safetyFactor;report.cL=flow.c_l;report.cOmega=flow.c_omega;report.cR=flow.c_r;
 report.poissonResidual=flow.poissonResidual;report.poissonResidualIsForwardAccuracyCertificate=false;
 assert(isequaln(s.config,sourceConfig)&&isequaln(s.scale,scale));
 clear ops cp s
 report.operatorCleared=true;report.status='operator_complete_external_resource_verdict_pending';
 fprintf('H64_DIRECTIONAL_PHASE factor_released_before_output\n');
 candidate=b.candidate;save(fullfile(outDir,'operator_output.mat'),'rho','rhoRate','flow','scale','candidate', ...
  'sourceConfig','sourceReferences','family','draft','-v7.3');
 save(fullfile(outDir,'report.mat'),'report','-v7.3');write_json(fullfile(outDir,'report.json'),report);
catch e
 failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
 write_json(fullfile(outDir,'failure.json'),failure);rethrow(e)
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
