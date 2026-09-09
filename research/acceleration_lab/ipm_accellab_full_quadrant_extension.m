function report=ipm_accellab_full_quadrant_extension(registrationFile,outDir,executeObservation)
% Dormant, research-only full-angle extension; never native evolution.
% The original model and original wall-sector numbers remain immutable.
if nargin<3,executeObservation=false;end
assert(islogical(executeObservation)&&isscalar(executeObservation)&&~isfolder(outDir));mkdir(outDir);
reg=jsondecode(fileread(registrationFile));assert(maxNumCompThreads==10);
assert(strcmp(reg.kind,'fixed_full_quadrant_H64_extension_v1')&& ...
 isequal(reg.steps(:)',[2213,4291,5874])&&isequal(reg.radiusBand(:)',[4,5])&& ...
 reg.radialQueryCount==17&&reg.angularQueryCount==129&&reg.modelGaussOrder==32&& ...
 isequal(reg.modelRKSteps(:)',[32,64,128])&&~reg.prospectiveExperiment);
verify_sources(reg,registrationFile);
report=struct('kind','full_quadrant_H64_frozen_model_extension','prospectiveExperiment',false, ...
 'executeObservation',executeObservation,'nativeCPRead',false,'LUCalls',0,'PDECalls',0, ...
 'originalModelModified',false,'originalWallSectorGateModified',false,'fullExteriorBoundProven',false, ...
 'status','prepared_only_not_executed');
if ~executeObservation,write_json(fullfile(outDir,'report.json'),report);return;end
assert(strcmp(getenv('IPM_FULL_QUADRANT_OBSERVATION_WINDOW'),'root_authorized_noLU_observation'), ...
 'ipm:ResearchNoLUWindow','An allocated no-LU analysis slot is required.');
profile clear;profile on;
try
 angular=angular_audit(reg);save(fullfile(outDir,'angular_audit.mat'),'angular');
 report.original32FullAngleQuadratureAudit=rmfield(angular,{'theta','gauss','reference','absoluteError'});
 cases=cell(3,1);raw=cell(3,1);
 for k=1:3
  cp=ipm.output.readCheckpoint(reg.inputs(k).path);s=cp.payload.state;h=cp.payload.log.history;
  assert(s.step==reg.steps(k)&&strcmp(s.runMetadata.caseId,reg.caseId)&& ...
   all(ipm.output.trustedMask(h,s.config))&&h.common.acceptedStep(1)==0&&h.common.physicalTime(1)==0&& ...
   s.scale.X_shift==0&&isequal(s.config.grid.xlim,[-64,64])&&s.config.grid.ymax==32&& ...
   strcmp(s.config.physics.symmetryMode,'double_odd_omega')&& ...
   strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic')&& ...
   strcmp(s.config.transport.spatialDiscretization,'high_order'));
  C=exp(s.scale.logC_l);W=exp(s.scale.logC_omega);x=s.x/C;y=s.y/C;t=s.scale.physicalTime;
  Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);om=s.rho*Dx';
  assert(C*(1/W)*max(abs(om),[],'all')==h.common.physicalRhoXInf(end));
  omega=exp(s.scale.logC_l-s.scale.logC_omega)*om;
  [X,Y]=meshgrid(x,y);R=hypot(X,Y);lead=X.^7./(X.^8+Y.^8);
  support=false(size(X));support(1:end-3,4:end-3)=true;
  weights=(ipm.mesh.quadrature(s.y)*ipm.mesh.quadrature(s.x))/C^2;
  assert(isequal(size(weights),size(omega))&&all(isfinite(weights),'all')&&all(weights>0,'all'));
  positiveBand=X>0 & R>=4 & R<=5 & support;
  axisBand=X==0 & Y>=4 & Y<=5 & support;
  [qr,qt]=meshgrid(linspace(4,5,17),linspace(0,pi/2,129));
  qx=qr.*cos(qt);qy=qr.*sin(qt);qx(end,:)=0;qy(1,:)=0;
  covered=qx>=x(4)&qx<=x(end-3)&qy>=0&qy<=y(end-3);
  % Omega already differentiates rho once in X. Hermite X/mixed jets add
  % another maintained X derivative; retain a distinct composed support.
  hermiteCovered=qx>=x(7)&qx<=x(end-6)&qy>=0&qy<=y(end-3);
  [iy,ix]=find(positiveBand);theta=atan2(Y(positiveBand),X(positiveBand));
  angularPopulation=histcounts(theta,linspace(0,pi/2,9));
  coverage=struct('radiusBand',[4,5],'angularRange',[0,pi/2], ...
   'registeredQueryCount',numel(qx),'nativeSourceCoverageFraction',mean(covered,'all'), ...
   'hermiteComposedCoverageFraction',mean(hermiteCovered,'all'), ...
   'completeNativeSourceQuerySupport',all(covered,'all'),'completeHermiteQuerySupport',all(hermiteCovered,'all'), ...
   'nativePositivePointCount',nnz(positiveBand),'distinctX',numel(unique(ix)),'distinctY',numel(unique(iy)), ...
   'eightAngularBinPopulation',angularPopulation,'axisPointCount',nnz(axisBand), ...
   'nativePositiveWeightsSum',sum(weights(positiveBand)),'exactQuarterAnnulusArea',pi*(25-16)/4, ...
   'nativeSupportRule','Original rows1:Ny-3, columns4:Nx-3; exact axis nodes recorded separately.', ...
   'hermiteSupportRule','Queries between x7 and x(Nx-6), and y0 to y(Ny-3), account for composed rho-to-Omega-to-Hermite derivatives.');
  coverage.eligible=coverage.completeNativeSourceQuerySupport&&coverage.distinctX>=5&& ...
   coverage.distinctY>=3&&all(angularPopulation>0)&&nnz(positiveBand)<=reg.maximumNativeModelPoints;
  q=struct('step',s.step,'caseId',s.runMetadata.caseId,'canonicalTime',s.scale.canonicalTime, ...
   'physicalTime',t,'physicalBox',[x(1),x(end),y(end)],'coverage',coverage, ...
   'nativeSignatureValidated',true,'allNativeHistoryTrusted',true, ...
   'scale',s.scale,'angularQuadratureQualified',angular.passed);
  [sector,sectorData]=original_sector(reg,k,X,Y,R,omega,lead,support,t);
  q.originalWallSector=sector;
  native=struct('x',X(positiveBand),'y',Y(positiveBand),'actualOmega',omega(positiveBand), ...
   'leadingOmega',lead(positiveBand),'physicalPositiveWeights',weights(positiveBand), ...
   'model',[],'metrics',[],'failure',[]);
  if nnz(positiveBand)>0&&nnz(positiveBand)<=reg.maximumNativeModelPoints
   [native.model,native.failure]=model_values(native.x,native.y,t,reg);
   if isempty(native.failure)
    native.metrics=band_metrics(native.actualOmega,native.model.omega128,native.leadingOmega,native.physicalPositiveWeights);
   end
  else
   native.failure=struct('identifier','ipm:ObservationSupport','message','No positive native band nodes or registered point budget exceeded.');
  end
  axis=struct('x',X(axisBand),'y',Y(axisBand),'actualOmega',omega(axisBand), ...
   'modelAnalyticOddSymmetryTrace',zeros(nnz(axisBand),1), ...
   'interpretation','Exact x=0 trace of the continuous odd-source model, separately labeled; the locked X>0 numerical helper is not evaluated or modified on this axis.');
  if any(axisBand,'all'),axis.actualAbsoluteInfinity=max(abs(axis.actualOmega));else,axis.actualAbsoluteInfinity=NaN;end
  polar=struct('radius',qr,'theta',qt,'x',qx,'y',qy,'nativeCovered',covered,'hermiteCovered',hermiteCovered, ...
   'actualHermite',NaN(size(qx)),'actualLinear',NaN(size(qx)),'model',[],'failure',[], ...
   'leadingOmega',cos(qt).^7./(cos(qt).^8+sin(qt).^8)./qr);
  polar.leadingOmega(end,:)=0;
  if any(hermiteCovered,'all')
   v=ipm_accellab_tensor_hermite(omega,x,y,C*Dx,C*Dy,qx(hermiteCovered),qy(hermiteCovered));
   polar.actualHermite(hermiteCovered)=v.value;
   polar.actualLinear(hermiteCovered)=interp2(x,y,omega,qx(hermiteCovered),qy(hermiteCovered),'linear');
  end
  positiveQuery=qx>0;
  [polar.model,polar.failure]=model_values(qx(positiveQuery),qy(positiveQuery),t,reg);
  polar.modelOmega=NaN(size(qx));polar.modelOmega(~positiveQuery)=0;
  if isempty(polar.failure),polar.modelOmega(positiveQuery)=polar.model.omega128;end
  polarWeights=trap_weights(qt(:,1))*trap_weights(qr(1,:))'.*qr;
  assert(all(polarWeights>0,'all'));
  polar.weights=polarWeights;
  polar.weightsMeaning='Independent fixed-polar trapezoid weights r dr dtheta, not native solver quadrature and not an exact integral certificate.';
  mask=hermiteCovered&isfinite(polar.modelOmega);
  if any(mask,'all')
   polar.supportedSubsetHermiteMetrics=band_metrics(polar.actualHermite(mask),polar.modelOmega(mask),polar.leadingOmega(mask),polarWeights(mask));
   polar.supportedSubsetLinearMetrics=band_metrics(polar.actualLinear(mask),polar.modelOmega(mask),polar.leadingOmega(mask),polarWeights(mask));
   polar.observerDifference=band_metrics(polar.actualHermite(mask),polar.actualLinear(mask),polar.leadingOmega(mask),polarWeights(mask));
  end
  q.nativeBandMetrics=native.metrics;q.nativeModelFailure=native.failure;
  q.polarModelFailure=polar.failure;q.exactAxis=rmfield(axis,{'x','y','actualOmega','modelAnalyticOddSymmetryTrace'});
  q.axisAndPositiveNativeRelativeInfinity=NaN;
  if ~isempty(native.metrics)&&isfinite(axis.actualAbsoluteInfinity)
   q.axisAndPositiveNativeRelativeInfinity=max(native.metrics.absoluteInfinity,axis.actualAbsoluteInfinity)/native.metrics.globalLeadingInfinity;
  end
  q.extendedObservationQualified=coverage.eligible&&angular.passed&&isempty(native.failure)&& ...
   ~isempty(native.metrics)&&native.model.integrationPassed&& ...
   q.axisAndPositiveNativeRelativeInfinity<=reg.relativeErrorThreshold&&native.metrics.relativeWeightedL2<=reg.relativeErrorThreshold;
  q.fullHermiteObservationQualified=q.extendedObservationQualified&&coverage.completeHermiteQuerySupport&& ...
   isempty(polar.failure)&&polar.model.integrationPassed&&isfield(polar,'supportedSubsetHermiteMetrics')&& ...
   polar.supportedSubsetHermiteMetrics.relativeInfinity<=reg.relativeErrorThreshold&& ...
   polar.supportedSubsetHermiteMetrics.relativeWeightedL2<=reg.relativeErrorThreshold;
  q.interpretation='Post-development extension only. Unsupported portions remain missing; subset errors never confer full-angle support. No full exterior source or boundary closure follows.';
  cases{k}=q;raw{k}=struct('native',native,'axis',axis,'polar',polar,'originalSector',sectorData);
  save(fullfile(outDir,sprintf('case_%d.mat',s.step)),'q','native','axis','polar','sectorData','-v7.3');
  assert(isequaln(ipm.output.checkpointSignature(cp.payload),cp.signature));
  clear cp s h Dx Dy X Y R om omega lead weights native polar axis sectorData
 end
 profile off;profileInfo=profile('info');files={profileInfo.FunctionTable.FileName};
 for suffix={'+mesh/build.m','+field/poisson.m','+field/velocity.m','+evolve/flow.m','+evolve/advance.m','+output/restoreCheckpoint.m'}
  assert(~any(endsWith(files,suffix{1})),['Forbidden numerical call: ',suffix{1}]);
 end
 report.status='extension_observed_no_PDE_or_LU';report.nativeCPRead=true;report.cases=cases;
 report.allExtendedObservationsQualified=all(cellfun(@(c)c.extendedObservationQualified,cases));
 report.allOriginalWallSectorNumbersExact=all(cellfun(@(c)c.originalWallSector.originalActualAndCoordinatesBitwise,cases));
 report.noLUCallGraphVerified=true;verify_sources(reg,registrationFile);
 save(fullfile(outDir,'observation.mat'),'report','raw','angular','profileInfo','-v7.3');
 write_json(fullfile(outDir,'report.json'),report);
catch e
 profile off;failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
 save(fullfile(outDir,'failure.mat'),'failure','report');write_json(fullfile(outDir,'failure.json'),failure);rethrow(e)
end
end

function [r,raw]=original_sector(reg,k,X,Y,R,omega,lead,support,t)
entry=reg.inputs(k);old=load(entry.originalModelFile,'report','raw');
oldCase=old.report.cases{entry.originalModelCase};oldRaw=old.raw{entry.originalModelCase};
assert(oldCase.step==entry.step&&oldCase.physicalTime==t);
for j=1:2
 b=oldCase.bands{j}.radiusBand;mask=support&X>0&R>=b(1)&R<=b(2)&Y<=.25*R;
 assert(isequaln(X(mask),oldRaw{j}.x)&&isequaln(Y(mask),oldRaw{j}.y)&& ...
  isequaln(omega(mask),oldRaw{j}.actualOmega)&&isequaln(lead(mask),oldRaw{j}.leadingOmega));
end
r=struct('sourceFile',entry.originalModelFile,'originalActualAndCoordinatesBitwise',true, ...
 'originalNumbersRetainedVerbatim',true,'modelRerun',false,'bands',{oldCase.bands}, ...
 'meaning','Prior registered sector predictions and gates are copied exactly under an identical model SHA; they are not reclassified using new band-global denominators.');
raw=oldRaw;
end
function [v,failure]=model_values(x,y,t,reg)
failure=[];v=[];
try
 [~,o32]=frozen_homogeneous_characteristics(x,y,t,32);
 [~,o64]=frozen_homogeneous_characteristics(x,y,t,64);
 [rho,o128,state]=frozen_homogeneous_characteristics(x,y,t,128);
 ref=max(abs(x.^7./(x.^8+y.^8)));
 assert(isfinite(ref)&&ref>0&&all(isfinite(o128))&&all(isfinite(rho)));
 changes=[max(abs(o32-o64))/ref,max(abs(o64-o128))/ref];
 v=struct('rho128',rho,'omega32',o32,'omega64',o64,'omega128',o128,'backwardState128',state, ...
  'globalLeadingInfinityDenominator',ref,'integrationDifferences',changes, ...
  'integrationPassed',all(changes<=reg.integrationLimits(:)'));
catch e
 failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack, ...
  'pointCount',numel(x),'xRange',[min(x),max(x)],'yRange',[min(y),max(y)],'t',t);
end
end
function m=band_metrics(actual,model,leading,w)
assert(numel(actual)==numel(model)&&numel(model)==numel(w)&&all(w>0)&&all(isfinite(w)));
actual=actual(:);model=model(:);leading=leading(:);w=w(:);d=actual-model;
denInf=max(abs(leading));denL2=sqrt(sum(w.*leading.^2));total=sum(w);
assert(denInf>0&&denL2>0&&all(isfinite([actual;model;leading])));
m=struct('absoluteInfinity',max(abs(d)),'absoluteWeightedRms',sqrt(sum(w.*d.^2)/total), ...
 'globalLeadingInfinity',denInf,'globalLeadingWeightedL2',denL2, ...
 'relativeInfinity',max(abs(d))/denInf,'relativeWeightedL2',sqrt(sum(w.*d.^2))/denL2, ...
 'pointCount',numel(d),'positiveWeightSum',total,'pointwiseDivisionByLeadingUsed',false, ...
 'meaning','Finite weighted samples on the stated support, not a rigorous full-band norm or independent grid convergence.');
end
function a=angular_audit(reg)
theta=linspace(0,pi/2,reg.angularQueryCount)';k=(1:31)';b=k./sqrt(4*k.^2-1);
[q,l]=eig(diag(b,1)+diag(b,-1));[nodes,i]=sort(diag(l));weights=2*q(1,i)'.^2;
ang=.5*(nodes+1)*theta';co=cos(ang);si=sin(ang);ff=co.^7./(co.^8+si.^8);
ic=(.5*theta').*(weights'*(co.*ff));is=(.5*theta').*(weights'*(si.*ff));
c=cos(theta);z=sin(theta);p=z.*(pi/4-ic')+c.*is';dp=c.*(pi/4-ic')-z.*is';
u1=-z.*p-c.*dp;u2=c.*p-z.*dp;icr=zeros(size(theta));isr=icr;
for j=2:numel(theta)
 icr(j)=integral(@(v)cos(v).^8./(cos(v).^8+sin(v).^8),0,theta(j),'AbsTol',1e-13,'RelTol',1e-13);
 isr(j)=integral(@(v)sin(v).*cos(v).^7./(cos(v).^8+sin(v).^8),0,theta(j),'AbsTol',1e-13,'RelTol',1e-13);
end
pr=z.*(pi/4-icr)+c.*isr;dpr=c.*(pi/4-icr)-z.*isr;
g=[ic',is',p,dp,u1,u2];ref=[icr,isr,pr,dpr,-z.*pr-c.*dpr,c.*pr-z.*dpr];err=abs(g-ref);
a=struct('theta',theta,'gauss',g,'reference',ref,'absoluteError',err, ...
 'fields',{{'Ic','Is','p','pPrime','U1','U2'}},'maximumAbsoluteDifference',max(err,[],1), ...
 'registeredGate',reg.angularQuadratureGate,'passed',all(max(err,[],1)<=reg.angularQuadratureGate), ...
 'originalOrderUnchanged',true,'axisComputedPAndU1',[p(end),u1(end)], ...
 'fullAnglePropagationErrorBoundClaim',false);
end
function w=trap_weights(x)
x=x(:);d=diff(x);w=[d(1)/2;(d(1:end-1)+d(2:end))/2;d(end)/2];
end
function verify_sources(reg,registrationFile)
assert(~isfolder(fullfile(pwd,'+ipm')));
for item={'output/readCheckpoint','output/checkpointSignature','output/trustedMask','mesh/fdMatrix','mesh/quadrature'}
 parts=strsplit(item{1},'/');name=['ipm.',parts{1},'.',parts{2}];
 assert(strcmp(which(name),fullfile(reg.runtimeSource,'+ipm',['+',parts{1}],[parts{2},'.m'])));
end
assert(strcmp(which('frozen_homogeneous_characteristics'),fullfile(reg.helperSource,'frozen_homogeneous_characteristics.m')));
assert(strcmp(which('ipm_accellab_full_quadrant_extension'),fullfile(reg.helperSource,'ipm_accellab_full_quadrant_extension.m')));
assert(strcmp(which('ipm_accellab_tensor_hermite'),fullfile(reg.helperSource,'ipm_accellab_tensor_hermite.m')));
manifest=jsondecode(fileread(reg.sourceManifest));paths={manifest.entries.path};
assert(any(strcmp(paths,registrationFile))&&all(ismember({reg.inputs.path},paths)));
for j=1:numel(manifest.entries),e=manifest.entries(j);assert(strcmp(sha256_file(e.path),e.sha256),['Changed file: ',e.path]);end
end
function value=sha256_file(file)
md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));
while true,b=fread(fid,1048576,'*uint8');if isempty(b),break;end;md.update(typecast(b,'int8'));end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
