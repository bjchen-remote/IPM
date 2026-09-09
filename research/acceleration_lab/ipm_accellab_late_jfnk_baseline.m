function report=ipm_accellab_late_jfnk_baseline(registrationFile,outDir,executeBaseline)
% Research-only, default dormant: strict native read, one saved-grid build,
% one original flow, factor release, immutable array cache and readback.
% No state projection, candidate, timestep or native checkpoint is produced.
if nargin<3,executeBaseline=false;end
assert(islogical(executeBaseline)&&isscalar(executeBaseline)&&maxNumCompThreads==10);
assert(~isfolder(outDir));mkdir(outDir);
reg=jsondecode(fileread(registrationFile));
report=struct('kind','late9205_JFNK_baseline_only_v1','executeBaseline',executeBaseline, ...
 'status','started','nativeSignatureValidated',false,'actualLUExecuted',false, ...
 'freshRHSCalls',0,'operatorCleared',false,'candidateConstructed',false, ...
 'PDETimeAdvanced',false,'nativeCheckpointWritten',false,'accelerationQualified',false);
try
 assert(strcmp(reg.kind,'late9205_JFNK_baseline_registration_v1')&&reg.step==9205&& ...
  ~reg.candidateAllowed&&~reg.PDEAllowed&&reg.baselineFreshRHSBudget==1);
 verify_sources(reg);
 profile clear;profile on;
 cp=ipm.output.readCheckpoint(reg.checkpointFile);s=cp.payload.state;h=cp.payload.log.history;
 assert(s.step==reg.step&&strcmp(s.runMetadata.caseId,reg.caseId)&& ...
  isequal([numel(s.x),numel(s.y)],[895,386])&&s.remeshCount==10&& ...
  all(ipm.output.trustedMask(h,s.config))&&s.scale.X_shift==0);
 assert(strcmp(s.config.scaling.rescalingMode,'dynamic')&& ...
  strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic')&& ...
  strcmp(s.config.scaling.scalingContract,'exact_gauge_no_feedback_v1')&& ...
  strcmp(s.config.scaling.lengthGauge,'transport_anchor')&& ...
  strcmp(s.config.scaling.cOmegaGauge,'wall_omega_quadratic_peak')&& ...
  strcmp(s.config.transport.spatialDiscretization,'high_order')&& ...
  strcmp(s.config.transport.transportScheme,'weno5_fd')&& ...
  strcmp(s.config.transport.wallTransportMode,'conservative_flux'));
 lineage=s.runMetadata.caseMetadata.latePhysicalBoxBranch;
 parent=ipm.output.readCheckpoint(lineage.parentCheckpoint);ps=parent.payload.state;
 assert(ps.step==lineage.parentStep&&strcmp(ps.runMetadata.caseId,lineage.parentCaseId)&& ...
  ps.scale.physicalTime==lineage.absolutePhysicalEpoch&& ...
  ps.scale.canonicalTime==lineage.parentCanonicalTime&& ...
  exp(ps.scale.logC_l)==lineage.parentCx&&exp(ps.scale.logC_omega)==lineage.parentComega&& ...
  lineage.canonicalCovarianceFactor==lineage.parentCx/lineage.parentComega&& ...
  ~lineage.parentHistoryInherited&&~lineage.parentCheckpointModified);
 clear parent ps
 sourceSignature=cp.signature;sourceLineage=lineage;
 report.nativeSignatureValidated=true;report.sourceSignatureRetainedInMAT=true;
 report.sourceCheckpointSHA256=sha256_file(reg.checkpointFile);
 report.parentCheckpointSHA256=sha256_file(lineage.parentCheckpoint);
 report.checkpointFile=reg.checkpointFile;report.runtimeSource=reg.runtimeSource;
 report.step=s.step;report.nodeCount=[numel(s.x),numel(s.y)];report.remeshCount=s.remeshCount;
 report.caseId=s.runMetadata.caseId;report.allHistoryTrusted=true;report.historyRecords=numel(h.common.t);
 report.nativeCanonicalTime=s.scale.canonicalTime;
 report.parentEquivalentTau=lineage.parentCanonicalTime+lineage.canonicalCovarianceFactor*s.scale.canonicalTime;
 report.localPhysicalTime=s.scale.physicalTime;report.physicalEpoch=lineage.absolutePhysicalEpoch;
 report.absolutePhysicalTime=lineage.absolutePhysicalEpoch+s.scale.physicalTime;
 report.scale=s.scale;
 report.lineage=struct('parentCheckpoint',lineage.parentCheckpoint,'parentCaseId',lineage.parentCaseId, ...
  'parentStep',lineage.parentStep,'parentCanonicalTime',lineage.parentCanonicalTime, ...
  'absolutePhysicalEpoch',lineage.absolutePhysicalEpoch,'parentCx',lineage.parentCx, ...
  'parentComega',lineage.parentComega,'canonicalCovarianceFactor',lineage.canonicalCovarianceFactor, ...
  'parentHistoryInherited',lineage.parentHistoryInherited,'parentCheckpointModified',lineage.parentCheckpointModified);
 report.nativeCoreCells=[h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)];
 report.nativeInstantaneousRates=[h.common.c_l(end),h.common.c_omega(end),h.common.c_r(end)];
 report.symmetryMode=s.config.physics.symmetryMode;
 report.futureCandidateFreshRHSBudget=reg.totalResearchFreshRHSBudget;
 report.futureBudgetIncludesThisBaseline=true;
 if ~executeBaseline
  profile off;profileInfo=profile('info');report.callAudit=call_audit(profileInfo,reg,0);
  report.status='strict_native_prepared_only_no_LU';verify_sources(reg);
  save(fullfile(outDir,'preparation.mat'),'report','profileInfo','sourceSignature','sourceLineage','-v7.3');
  write_json(fullfile(outDir,'report.json'),report);return
 end
 assert(strcmp(getenv('IPM_LATE_JFNK_BASELINE_WINDOW'),'root_authorized_single_factor_baseline_only'), ...
  'ipm:ResearchBaselineWindow','Root-authorized, watched, one-factor baseline window required.');
 % Avoid restoreCheckpoint's config-grid build followed by a saved-grid build.
 % The frozen numerical config itself is neither resolved again nor modified.
 grid=s.config.grid;grid.nx=numel(s.x);grid.ny=numel(s.y);grid.customX=s.x;grid.customY=s.y;
 fprintf('LATE_JFNK_PHASE single_saved_grid_build\n');timer=tic;
 ops=ipm.mesh.build(s.config,grid);report.buildSeconds=toc(timer);report.actualLUExecuted=true;
 assert(isa(ops.poisson,'decomposition')&&isequal(ops.x,s.x)&&isequal(ops.y,s.y));
 ops.baseX=s.baseX;ops.baseY=s.baseY;ops.remeshCount=s.remeshCount;
 ops.rescaling=ipm.output.restoreRuntimeReferences(ops.rescaling,s.rescaling,ops.x,true);
 assert(isequaln(ops.rescaling,s.rescaling)&&isequal(ops.baseX,s.baseX)&&isequal(ops.baseY,s.baseY));
 assert(isequal(ops.Dx,ipm.mesh.fdMatrix(s.x,1,7))&&isequal(ops.Dy,ipm.mesh.fdMatrix(s.y,1,7)));
 exactWeights=ipm.mesh.quadrature(s.y)*ipm.mesh.quadrature(s.x);
 assert(isequal(ops.integrationWeights,exactWeights)&&all(exactWeights>0,'all'));
 report.savedAxesBaseReferencesExact=true;report.pairedDerivativesWeightsExact=true;
 report.matrixSize=size(ops.A);report.matrixNnz=nnz(ops.A);
 fprintf('LATE_JFNK_PHASE single_original_flow\n');timer=tic;
 [rhoRate,flow]=ipm.evolve.flow(s.rho,ops,s.scale);
 report.flowSeconds=toc(timer);report.freshRHSCalls=1;
 assert(all(isfinite(rhoRate),'all')&&ipm.evolve.isFinite(struct('rho',s.rho,'flow',flow,'scale',s.scale)));
 Omega=s.rho*ops.Dx';FX=rhoRate*ops.Dx';assert(isequal(Omega,flow.source));
 % Reconstruct only the recorded timestep telemetry for a read-only one-row
 % record replay. This does not advance state or invent an accepted record.
 state=s;state.ops=ops;state.flow=flow;state.timeStep=recorded_timestep(h.common);
 blank=struct('history',struct('common',struct(),'gauge',struct(),'mesh',struct(),'anisotropic',struct()), ...
  'snapshotRho',{{}},'snapshotNormalizedTime',[],'snapshotX',{{}},'snapshotY',{{}});
 [oneRow,~]=ipm.output.record(blank,state);
 report.nativeTerminalRecordParity=terminal_pair(oneRow.history,h);
 clear state oneRow blank
 k=find(s.x==flow.omegaGaugeQuadraticStencilCenterX);assert(isscalar(k));
 peak=ipm.evolve.quadraticPeakFunctional(Omega(1,:),s.x,k);
 assert(isequaln(peak.value,flow.omegaGaugeQuadraticPeakValue)&&isequaln(peak.x,flow.omegaGaugeQuadraticPeakX));
 DPF=sum(peak.weights.*FX(1,peak.indices));
 stencilDPF=sum(peak.weights.*(rhoRate(1,:)*ops.Dx(peak.indices,:)'));
 report.quadratic=struct('functional',peak,'selectorIndex',k, ...
  'actualFullProductDPF',DPF,'actualStencilProductDPF',stencilDPF, ...
  'nativeAlgebraicDPF',flow.omegaGaugeQuadraticRateResidual, ...
  'fullMinusStencilDPF',DPF-stencilDPF,'fullMinusNativeAlgebraicDPF',DPF-flow.omegaGaugeQuadraticRateResidual, ...
  'actualPeakLogRate',DPF/peak.value,'nativeRateScale',flow.omegaGaugeQuadraticRateScale, ...
  'activeSelectionMargin',flow.omegaGaugeQuadraticActiveSelectionMargin, ...
  'interpretation','All three floating-point paths retained. Zero algebraic rate is not cross-selector peak conservation.');
 Cx=exp(s.scale.logC_l);Cw=exp(s.scale.logC_omega);clock=Cx/Cw;
 report.units=struct('Cx',Cx,'Comega',Cw,'dPhysicalTime_dNativeTau',1/clock, ...
  'nativeRates',report.nativeInstantaneousRates,'physicalRates',clock*report.nativeInstantaneousRates, ...
  'parentCanonicalRates',report.nativeInstantaneousRates/lineage.canonicalCovarianceFactor);
 report.anchor=struct('x',ops.rescaling.transportAnchorX, ...
  'nativeVelocityAtAnchor',interp1(s.x,flow.u1(1,:),ops.rescaling.transportAnchorX,'linear'), ...
  'nativeConstraintResidual',flow.lengthGaugeResidual,'reportedTransportVelocity',flow.transportAnchorVelocity, ...
  'meaning','Velocity in native coordinates; physical units require division by Comega. The anchor is an algebraic rate constraint, not a fixed peak-position constraint.');
 W=ops.integrationWeights;report.fullRhs=struct('nativeWeightedRms',sqrt(sum(W.*rhoRate.^2,'all')/sum(W,'all')), ...
  'nativeInfinity',max(abs(rhoRate),[],'all'),'rhoWeightedRms',sqrt(sum(W.*s.rho.^2,'all')/sum(W,'all')), ...
  'nativeOmegaRhsInfinity',max(abs(FX),[],'all'),'physicalTimeCoordinateRhsMultiplier',clock, ...
  'meaning','F=d rho/d native tau at fixed computational X,Y. clock*F remains a computational-coordinate derivative, not the physical Eulerian density RHS.');
 report.instantaneousBaselinePoissonResidual=flow.poissonResidual;
 report.poissonResidualIsForwardErrorCertificate=false;
 % Save a pure-array operator view; explicitly omit global matrix and factor.
 pureOps=rmfield(ops,{'A','poisson'});assert(~contains_object(pureOps));
 assert(isequaln(ipm.output.checkpointSignature(cp.payload),cp.signature));
 cache=struct('kind','late9205_original_true_RHS_baseline_arrays_v1','sourceState',s, ...
  'sourceHistory',h,'sourceSignature',cp.signature,'rhoRate',rhoRate,'flow',flow, ...
  'Omega',Omega,'FX',FX,'operators',pureOps,'sourceRegistration',reg, ...
  'noNativeCheckpoint',true,'noCandidate',true,'noPDETimeAdvanced',true);
 clear ops cp pureOps rhoRate flow Omega FX W exactWeights
 report.operatorCleared=true;
 fprintf('LATE_JFNK_PHASE factor_released_before_cache_and_observers\n');
 tmp=fullfile(outDir,'true_rhs_cache.writing.mat');target=fullfile(outDir,'true_rhs_cache.mat');
 assert(~isfile(target));save(tmp,'cache','-v7.3');readback=load(tmp,'cache');
 assert(same_data(cache,readback.cache),'ipm:ResearchCacheReadback','Every saved array/config/reference must read back exactly.');
 clear readback
 [ok,msg]=movefile(tmp,target);assert(ok,msg);report.cacheReadbackExact=true;
 report.cacheFile=target;report.cacheSHA256=sha256_file(target);
 assert(report.nativeTerminalRecordParity.passed,'ipm:ResearchNativePair','Original terminal record replay mismatch; unqualified cache retained.');
 window=s.config.scaling.transportAnchorX*reg.peakWindowAnchorFactors(:)';
 timer=tic;coordinates=ipm_accellab_continuous_inner_rates(cache.Omega,cache.FX,s.x,s.y, ...
  cache.operators.Dx,cache.operators.Dy,window);
 report.continuousCoordinates=coordinates;report.continuousObservationQualified=coordinates.valid;
 if coordinates.valid
  [residual,details]=ipm_accellab_continuous_residual(cache.Omega,cache.FX,s.x,s.y, ...
   cache.operators.Dx,cache.operators.Dy,coordinates);
  report.continuousResidual=residual;
  report.residualUnitMultipliers=struct('nativeTauToParentTau',1/lineage.canonicalCovarianceFactor, ...
   'nativeTauToPhysicalTime',clock,'nativeGToIntrinsicPhysicalPeakTime',1/coordinates.peak.value, ...
   'meaning','G=U_nativeTau. G/P=U_physicalTime/(actual physical C1 Omega peak). PPrime is measured, never set to zero.');
  save(fullfile(outDir,'continuous_observation.mat'),'coordinates','residual','details','-v7.3');
  report.sampleInfinityRefinement=struct('coarse',residual.sampledMetrics{1}.G, ...
   'fine',residual.sampledMetrics{2}.G,'independentSpatialConvergenceValidated',false);
 end
 report.observerSeconds=toc(timer);
 profile off;profileInfo=profile('info');report.callAudit=call_audit(profileInfo,reg,1);
 verify_sources(reg);report.status='baseline_cache_complete_candidate_not_authorized';
 report.spatialResidualQualification='Single native grid only; sampling refinement is not PDE/grid convergence. Previous nonmonotone multi-grid U_tau evidence remains unresolved.';
 save(fullfile(outDir,'audit.mat'),'report','profileInfo','-v7.3');write_json(fullfile(outDir,'report.json'),report);
 fprintf('LATE_JFNK_BASELINE complete calls=%d candidate=0 PDE=0\n',report.freshRHSCalls);
catch e
 profile off;clear ops state
 report.status='failed_preserved';report.failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
 write_json(fullfile(outDir,'failure.json'),report);rethrow(e)
end
end

function t=recorded_timestep(h)
map={'dt','acceptedCanonicalDt';'rate','transportRate';'dtCfl','cflDtLimit'; ...
 'maxDtLimit','maximumDtLimit';'canonicalLimit','canonicalEndpointDtLimit'; ...
 'physicalLimit','physicalEndpointDtLimit';'physicalClockSpeed','physicalClockSpeed'; ...
 'realizedCfl','realizedCfl';'activeLimiter','timestepActiveLimiter'};
t=struct();for j=1:size(map,1),t.(map{j,1})=h.(map{j,2})(end);end
end
function p=terminal_pair(actual,expected)
p=struct('passed',true,'comparedFields',0,'mismatches',{{}});
for g={'common','gauge','mesh','anisotropic'}
 key=g{1};a=actual.(key);b=expected.(key);names=union(fieldnames(a),fieldnames(b));
 for j=1:numel(names)
  n=names{j};ok=isfield(a,n)&&isfield(b,n)&&isequaln(a.(n)(end),b.(n)(end));
  p.comparedFields=p.comparedFields+1;
  if ~ok,p.passed=false;p.mismatches{end+1}=[key,'.',n];end
 end
end
p.telemetryMeaning='One-row record replay uses the real saved timestep scalar fields only; no accepted history is emitted or altered.';
end
function audit=call_audit(p,reg,expected)
audit=struct();
for item={'mesh/build','evolve/flow','evolve/rhs','field/velocity','field/poisson'}
 parts=strsplit(item{1},'/');file=fullfile(reg.runtimeSource,'+ipm',['+',parts{1}],[parts{2},'.m']);
 match=strcmp({p.FunctionTable.FileName},file)&~contains({p.FunctionTable.FunctionName},'>');
 n=sum([p.FunctionTable(match).NumCalls]);audit.(strjoin(parts,'_'))=n;assert(n==expected);
end
for item={'output/restoreCheckpoint','evolve/initialize','evolve/initializeScaling','evolve/advance','evolve/stepRk','solve'}
 parts=strsplit(item{1},'/');
 if numel(parts)==1,file=fullfile(reg.runtimeSource,'+ipm',[parts{1},'.m']);
 else,file=fullfile(reg.runtimeSource,'+ipm',['+',parts{1}],[parts{2},'.m']);end
 assert(~any(strcmp({p.FunctionTable.FileName},file)),['Forbidden call: ',item{1}]);
end
audit.forbiddenLifecycleCallsAbsent=true;
end
function verify_sources(reg)
assert(~isfolder(fullfile(pwd,'+ipm')),'ipm:ResearchCwd','Run from an artifact directory without +ipm.');
for item={'ipm_accellab_late_jfnk_baseline','ipm_accellab_continuous_inner_rates', ...
 'ipm_accellab_continuous_residual','ipm_accellab_tensor_hermite', ...
 'ipm_accellab_hermite_peak','ipm_accellab_hermite_line','ipm_accellab_hermite_level_roots'}
 assert(strcmp(which(item{1}),fullfile(reg.helperSource,[item{1},'.m'])));
end
for item={'ipm.solve','ipm.output.readCheckpoint','ipm.mesh.build','ipm.output.restoreRuntimeReferences', ...
 'ipm.evolve.flow','ipm.evolve.rhs','ipm.evolve.isotropicGauge','ipm.field.velocity','ipm.field.poisson'}
 parts=strsplit(item{1},'.');path=reg.runtimeSource;
 for j=1:numel(parts)-1,path=fullfile(path,['+',parts{j}]);end
 path=fullfile(path,[parts{end},'.m']);assert(strcmp(which(item{1}),path));
end
manifest=jsondecode(fileread(reg.sourceManifest));
for j=1:numel(manifest.entries)
 a=manifest.entries(j);assert(strcmp(sha256_file(a.path),a.sha256),['Changed source/input: ',a.path]);
end
end
function yes=same_data(a,b)
if isa(a,'function_handle')&&isa(b,'function_handle'),yes=same_data(functions(a),functions(b));return;end
if isstruct(a)&&isstruct(b)
 yes=isequal(size(a),size(b))&&isequal(fieldnames(a),fieldnames(b));if ~yes,return;end
 names=fieldnames(a);for k=1:numel(a),for j=1:numel(names),if ~same_data(a(k).(names{j}),b(k).(names{j})),yes=false;return;end,end,end
elseif iscell(a)&&iscell(b)
 yes=isequal(size(a),size(b));if ~yes,return;end
 for k=1:numel(a),if ~same_data(a{k},b{k}),yes=false;return;end,end
else,yes=isequaln(a,b);
end
end
function yes=contains_object(a)
if isstruct(a),yes=false;names=fieldnames(a);for k=1:numel(a),for j=1:numel(names),yes=yes||contains_object(a(k).(names{j}));end,end
elseif iscell(a),yes=any(cellfun(@contains_object,a));
else,yes=isobject(a)&&~isstring(a);
end
end
function value=sha256_file(file)
md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));
while true,b=fread(fid,1048576,'*uint8');if isempty(b),break;end;md.update(typecast(b,'int8'));end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
