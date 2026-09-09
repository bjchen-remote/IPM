function report=ipm_accellab_analyze_t0_tail_pair(registrationFile)
% Saved fixed-operator responses plus independent smooth exterior integrals.
% No Poisson operator/factor, restore, flow, RHS or PDE is evaluated here.
reg=jsondecode(fileread(registrationFile));assert(maxNumCompThreads==10&&~isfolder(reg.outputDirectory));mkdir(reg.outputDirectory);
profile clear;profile on;series=jsondecode(fileread(reg.initialIntegralSeriesFile));
report=struct('kind','actual_t0_full_tail_operator_response_analysis','registration',reg,'cases',{{}}, ...
 'noLU',true,'noPDE',true,'nativeGaugeModified',false,'evolvedClosureValidated',false,'directCLCorrectionUsed',false);
[X,Y]=meshgrid(linspace(0,4,21),linspace(0,2,11));query=[X(:),Y(:)];wallX=linspace(0,4,161)';
peakExact=7^(1/8);allQuery=[query;wallX,zeros(size(wallX));1,0;peakExact,0];
fieldData=cell(1,2);
for k=1:2
 original=jsondecode(fileread(reg.reportFiles{k}));d=load(reg.fieldFiles{k});H=reg.H(k);
 assert(original.LUCount==1&&original.PoissonSolveCount==2&&original.PDECount==0&& ...
  original.givenFlowAssemblyBitwiseNativeParity&&original.initialRatesBitwiseNativeParity&& ...
  original.referencesRhoOmegaAndPeakUnchanged&&isequal(d.x([1,end]),[-H,H])&&d.y(end)==H/2&& ...
  d.scale.logC_l==0&&d.scale.logC_omega==0&&d.scale.X_shift==0);
 assert(isequaln(d.Omega,d.rho*d.Dx')&&isequaln(d.baseline.peak,d.trial.peak));
 old=d.baseline.physical;new=d.trial.physical;deltaPsi=new.psi-old.psi;du1=new.u1-old.u1;du2=new.u2-old.u2;
 Dxx=ipm.mesh.fdMatrix(d.x,2,7);Dyy=ipm.mesh.fdMatrix(d.y,2,7);
 hx=deltaPsi*Dxx';hy=Dyy*deltaPsi;interior=hx(2:end-1,2:end-1)+hy(2:end-1,2:end-1);
 boundaryError=max([max(abs(deltaPsi(:,1)-d.tailBoundary.left(:))),max(abs(deltaPsi(:,end)-d.tailBoundary.right(:))), ...
  max(abs(deltaPsi(1,:)-d.tailBoundary.bottom(:)')),max(abs(deltaPsi(end,:)-d.tailBoundary.top(:)'))]);
 divergence=du1*d.Dx'+d.Dy*du2;
 harmonic=struct('fullBoundaryTraceAbsoluteInfError',boundaryError,'interiorDiscreteLaplacianAbsoluteInf',max(abs(interior),[],'all'), ...
  'relativeToSeparatedLaplacianTerms',max(abs(interior),[],'all')/max(max(abs(hx(2:end-1,2:end-1)),[],'all')+max(abs(hy(2:end-1,2:end-1)),[],'all'),realmin), ...
  'velocityFromDeltaPsiAbsoluteInf',max([max(abs(du1+d.Dy*deltaPsi),[],'all'),max(abs(du2-deltaPsi*d.Dx'),[],'all')]), ...
  'discreteDivergenceAbsoluteInf',max(abs(divergence),[],'all'),'notForwardAccuracyCertificate',true);
 integ=cell(1,2);for n=1:2,integ{n}=exterior_integral(allQuery,H,reg.quadratureOrders(n));end
 fine=integ{2};coarse=integ{1};quad=struct('orders',reg.quadratureOrders,'absoluteChanges',struct());
 for name={'psi','u1','u2','wallU1x'}
  field=name{1};mask=isfinite(fine.(field))&isfinite(coarse.(field));quad.absoluteChanges.(field)=max(abs(fine.(field)(mask)-coarse.(field)(mask)));
 end
 quad.passed=max(structfun(@(v)v,quad.absoluteChanges))<=reg.absoluteQuadratureGate;
 assert(quad.passed,'ipm:TailResponseQuadrature','Predeclared exterior-integral convergence gate failed.');
 nQuery=size(query,1);nWall=numel(wallX);idxAnchor=nQuery+nWall+1;idxPeak=idxAnchor+1;
 exactTailCL=-fine.u1(idxAnchor);exactTailCW=exactTailCL+fine.wallU1x(idxPeak);
 if iscell(series.cases),allH=cellfun(@(v)v.H,series.cases);else,allH=[series.cases.H];end
 selected=find(allH==H);assert(isscalar(selected));
 oldExactTail=series.tailCL(selected,end);exactInside=series.analyticInsideCL(selected);exactInfinite=series.infiniteCLComparator;
 rates0=d.baseline.rates;rates1=d.trial.rates;diffRates=rates1-rates0;
 row=struct('H',H,'actualNodeCount',size(d.rho([1,end],:)), ...
  'originalReportFile',reg.reportFiles{k},'quadrature',quad,'harmonic',harmonic, ...
  'baselineRates',rates0,'withTailRates',rates1,'actualRateDifference',diffRates,'exactExteriorCL',exactTailCL, ...
  'previousIndependentExactExteriorCL',oldExactTail,'newVsPreviousExactCL',exactTailCL-oldExactTail, ...
  'infiniteExactInitialCL',exactInfinite,'exactInsideCL',exactInside, ...
  'nativeBaselineMinusExactInsideCL',rates0(1)-exactInside, ...
  'discreteHarmonicAnchorResponseMinusExactTailCL',diffRates(1)-exactTailCL, ...
  'correctedNativeMinusExactInfiniteCL',rates1(1)-exactInfinite, ...
  'continuousExactTailCWChange',exactTailCW,'actualCWResponseMinusContinuousTailCW',diffRates(2)-exactTailCW, ...
  'CWComparatorMeaning','Continuous exact k8 wall-peak gauge response only; not a full exact Cw value or a replacement for native quadratic/WENO gauge.', ...
  'CWLeadingAHContributionCancelsInContinuousIsotropicGauge',true);
 row.actualNodeCount=[numel(d.x),numel(d.y)];
 row.CLDecompositionClosure=(rates0(1)-exactInside)+(diffRates(1)-exactTailCL)-(rates1(1)-exactInfinite);
 [NX,NY]=meshgrid(d.x,d.y);exactSource=NX.^7./(1+NX.^8+NY.^8);sourceError=d.Omega-exactSource;
 row.nativeSourceDerivativeError=struct('wholeBoxInf',max(abs(sourceError),[],'all'),'coreInf',max(abs(sourceError(NX>=0&NX<=2&NY<=1))));
 a=ipm_accellab_tensor_hermite(deltaPsi,d.x,d.y,d.Dx,d.Dy,query(:,1),query(:,2));
 b=ipm_accellab_tensor_hermite(du1,d.x,d.y,d.Dx,d.Dy,query(:,1),query(:,2));
 c=ipm_accellab_tensor_hermite(du2,d.x,d.y,d.Dx,d.Dy,query(:,1),query(:,2));
 observation=struct('psi',a.value,'u1',b.value,'u2',c.value,'psiDerivedU1',-a.derivativeY,'psiDerivedU2',a.derivativeX);
 row.responseComparedWithExactExterior=struct();
 for name={'psi','u1','u2'}
  f=name{1};v=observation.(f);e=fine.(f)(1:nQuery);row.responseComparedWithExactExterior.(f)=struct('absoluteInf',max(abs(v-e)), ...
   'relativeInf',max(abs(v-e))/max(abs(e)),'meaning','Common finite inner grid, includes paired Hermite observation error.');
 end
 row.responseComparedWithExactExterior.velocityRepresentationDiscrepancy=max([max(abs(observation.u1-observation.psiDerivedU1)),max(abs(observation.u2-observation.psiDerivedU2))]);
 oldQ=struct();newQ=struct();
 for name={'u1','u2'}
  f=name{1};a=ipm_accellab_tensor_hermite(old.(f),d.x,d.y,d.Dx,d.Dy,query(:,1),query(:,2));
  b=ipm_accellab_tensor_hermite(new.(f),d.x,d.y,d.Dx,d.Dy,query(:,1),query(:,2));oldQ.(f)=a.value;newQ.(f)=b.value;
 end
 wold=ipm_accellab_tensor_hermite(old.u1,d.x,d.y,d.Dx,d.Dy,wallX,zeros(size(wallX)));
 wnew=ipm_accellab_tensor_hermite(new.u1,d.x,d.y,d.Dx,d.Dy,wallX,zeros(size(wallX)));
 fieldData{k}=struct('query',query,'exactExterior',fine,'operatorResponse',observation,'baseline',oldQ,'withTail',newQ, ...
  'wallX',wallX,'wallU1Baseline',wold.value,'wallU1WithTail',wnew.value,'exactTailWallU1',fine.u1(nQuery+(1:nWall)), ...
  'exactTailWallU1x',fine.wallU1x(nQuery+(1:nWall)));
 report.cases{k}=row;
 fprintf('TAIL_RESPONSE H=%g dCL=%.15g exact=%.15g harmonicGap=%.5g dCW=%.7g exactCW=%.7g\n',H,diffRates(1),exactTailCL,row.discreteHarmonicAnchorResponseMinusExactTailCL,diffRates(2),exactTailCW);
end
r1=report.cases{1};r2=report.cases{2};cross=struct('H',reg.H,'baselineRateDifference',r2.baselineRates-r1.baselineRates, ...
 'withTailRateDifference',r2.withTailRates-r1.withTailRates,'exactExteriorCLDifference',r2.exactExteriorCL-r1.exactExteriorCL);
cross.CLGapReductionFactor=abs(cross.baselineRateDifference(1)/cross.withTailRateDifference(1));
cross.CLRelativeCrossGapBefore=abs(cross.baselineRateDifference(1))/r1.baselineRates(1);
cross.CLRelativeCrossGapAfter=abs(cross.withTailRateDifference(1))/r1.withTailRates(1);
cross.remainingGapHasJointDiscretizationAndInitialAxesEffects=true;
for branch={'baseline','withTail'}
 name=branch{1};a=fieldData{1}.(name);b=fieldData{2}.(name);e=hypot(b.u1-a.u1,b.u2-a.u2);
 cross.(name)=struct('commonQueryVelocityAbsoluteInf',max(e),'commonQueryVelocityRelativeInfToH8',max(e)/max(hypot(a.u1,a.u2)));
end
report.crossCase=cross;profile off;profileInfo=profile('info');files={profileInfo.FunctionTable.FileName};
for suffix={'+mesh/build.m','+field/poisson.m','+field/velocity.m','+evolve/flow.m','+output/restoreCheckpoint.m','+evolve/advance.m'}
 assert(~any(endsWith(files,suffix{1})),['Forbidden call: ',suffix{1}]);
end
report.noLUCallGraphVerified=true;save(fullfile(reg.outputDirectory,'analysis.mat'),'report','fieldData','profileInfo','-v7.3');write_json(fullfile(reg.outputDirectory,'report.json'),report);
end
function values=exterior_integral(query,H,n)
% Independent origin-polar integration with r=r_box(theta)/u, u in(0,1).
% All targets are strictly inside the box: no near-source singular quadrature.
assert(all(abs(query(:,1))<H)&all(query(:,2)<H/2));
[u,w]=gauss01(n);angles=[0,atan(.5),pi/2];xi=[];eta=[];ws=[];
for j=1:2
 theta=angles(j)+(angles(j+1)-angles(j))*u';wt=(angles(j+1)-angles(j))*w';
 c=cos(theta);s=sin(theta);r0=H./max(c,2*s);R=r0./u;
 x=R.*c;y=R.*s;omega=x.^7./(1+x.^8+y.^8);weight=(w*wt).*r0.^2./u.^3;
 xi=[xi;x(:)];eta=[eta;y(:)];ws=[ws;omega(:).*weight(:)]; %#ok<AGROW>
end
nq=size(query,1);values=struct('psi',zeros(nq,1),'u1',zeros(nq,1),'u2',zeros(nq,1),'wallU1x',NaN(nq,1));
for k=1:nq
 x=query(k,1);y=query(k,2);D=(x-xi).^2+(y-eta).^2;O=(x+xi).^2+(y+eta).^2;
 Z=16*x*y*xi.*eta./(D.*O);G=log1p(Z)/(4*pi);
 gx=(16*y*xi.*eta./(D.*O)-Z.*(2*(x-xi)./D+2*(x+xi)./O))./(4*pi*(1+Z));
 gy=(16*x*xi.*eta./(D.*O)-Z.*(2*(y-eta)./D+2*(y+eta)./O))./(4*pi*(1+Z));
 values.psi(k)=sum(ws.*G);values.u1(k)=-sum(ws.*gy);values.u2(k)=sum(ws.*gx);
 if y==0
  kernel=-4*xi.*eta./(pi*D.*O).*(1-x*(2*(x-xi)./D+2*(x+xi)./O));values.wallU1x(k)=sum(ws.*kernel);
 end
end
end
function [x,w]=gauss01(n)
k=(1:n-1)';b=k./sqrt(4*k.^2-1);[Q,T]=eig(diag(b,1)+diag(b,-1),'vector');[t,idx]=sort(T);x=(t+1)/2;w=Q(1,idx)'.^2;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
