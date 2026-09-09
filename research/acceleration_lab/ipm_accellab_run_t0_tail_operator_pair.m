function report=ipm_accellab_run_t0_tail_operator_pair(contextFile,out,executeNative)
%IPM_ACCELLAB_RUN_T0_TAIL_OPERATOR_PAIR Research instantaneous BC comparison.
% Default is no LU. Explicit true builds ONE fixed native operator; baseline
% and exact-initial-tail boundary use that same factor, rho and original refs.
% Never advances, projects, writes native checkpoints, or edits C formulas.
if nargin<3,executeNative=false;end
assert(islogical(executeNative)&&isscalar(executeNative)&&maxNumCompThreads==10&&~isfolder(out));mkdir(out);
report=struct('kind','original_t0_fixed_operator_full_tail_pair_v1','contextFile',contextFile, ...
    'executeNative',executeNative,'LUCount',0,'PoissonSolveCount',0,'PDECount',0, ...
    'nativeCheckpointWritten',false,'realEvolvedClosureValidated',false, ...
    'coreObservationBox',[0,2,0,1],'holdoutObservationBox',[2,4,0,2], ...
    'changedOnlyExternalGreenTrace',true,'directCLCorrectionUsed',false,'acceptedTrajectoryClaim',false);
write_json(fullfile(out,'registration.json'),report);
try
    loaded=load(contextFile,'context');c=loaded.context;
    assert(c.audit.quadraturePassed&&c.mathematicalTime==0&&~c.isNativeCheckpoint&&c.LUCount==0&&c.PDECount==0);
    source=ipm.output.readCheckpoint(c.sourceCheckpoint);s=source.payload.state;l=source.payload.log;
    assert(s.remeshCount==0&&isequaln(c.config,s.config)&&isequaln(c.runtimeReferences,s.rescaling));
    assert(isequaln(c.x,s.x)&&isequaln(c.y,s.y)&&isequaln(c.rho,l.snapshotRho{1})&& ...
        isequaln(c.x,l.snapshotX{1})&&isequaln(c.y,l.snapshotY{1})&&l.snapshotNormalizedTime(1)==0);
    assert(all(ipm.output.trustedMask(l.history,s.config)));
    assert(numel(c.rho)<=60000&&isequaln(c.x,s.baseX)&&isequaln(c.y,s.baseY));
    assert(strcmp(c.config.transport.spatialDiscretization,'high_order')&& ...
        strcmp(c.config.transport.transportScheme,'weno5_fd')&& ...
        strcmp(c.config.physics.symmetryMode,'double_odd_omega')&& ...
        strcmp(c.config.elliptic.farBoundaryMode,'green')&& ...
        strcmp(c.config.scaling.dynamicScaleGeometry,'isotropic'));
    assert(isequaln(c.Dx,ipm.mesh.fdMatrix(c.x,1,7))&& ...
        isequaln(c.Dy,ipm.mesh.fdMatrix(c.y,1,7))&&isequaln(c.Omega,c.rho*c.Dx'));
    assert(isequaln(c.integrationWeights,ipm.mesh.quadrature(c.y)*ipm.mesh.quadrature(c.x))&& ...
        isequaln(c.mass0,s.mass0)&&isequaln(c.rhoRange0,s.rhoRange0));
    assert(isequaln([c.initialNativeCL,c.initialNativeCOmega,c.initialNativeCR], ...
        [l.history.common.c_l(1),l.history.common.c_omega(1),l.history.common.c_r(1)])&& ...
        c.initialNativeQuadraticPeak==l.history.gauge.omegaGaugeQuadraticPeakValue(1));
    assert(all(c.tailBoundary.bottom==0)&&c.tailBoundary.left(1)==0&&c.tailBoundary.right(1)==0&& ...
        c.tailBoundary.left(end)==c.tailBoundary.top(1)&&c.tailBoundary.right(end)==c.tailBoundary.top(end));
    report.label=c.label;report.sourceStrictlyRevalidated=true;report.contextStrictlyPaired=true;
    report.boundaryQuadratureAudit=c.audit;report.actualNodeCount=[numel(c.x),numel(c.y)];
    report.forwardVelocityErrorCertified=false;
    if ~executeNative
        report.status='prepared_only_no_lu';write_json(fullfile(out,'report.json'),report);return
    end
    clear source s l loaded
    timer=tic;report.rssBeforeBuildBytes=rss_bytes();grid=c.config.grid;grid.customX=c.x;grid.customY=c.y;
    assert(isequaln([grid.nx,grid.ny],report.actualNodeCount));
    buildTimer=tic;ops=ipm.mesh.build(c.config,grid);report.buildSeconds=toc(buildTimer);report.LUCount=1;
    ops.rescaling=ipm.output.restoreRuntimeReferences(ops.rescaling,c.runtimeReferences,ops.x,true);
    assert(isequaln(ops.rescaling,c.runtimeReferences)&&isequaln(ops.x,c.x)&&isequaln(ops.y,c.y)&& ...
        isequaln(ops.Dx,c.Dx)&&isequaln(ops.Dy,c.Dy)&&isequaln(ops.integrationWeights,c.integrationWeights));
    report.rssAfterBuildBytes=rss_bytes();rho=c.rho;Omega=c.Omega;scale=ipm.evolve.initialScale(ops);
    assert(isequaln(rho*ops.Dx',Omega));
    t=tic;[rhsBaseline,nativeFlow]=ipm.evolve.rhs(rho,ops,scale);report.baselineRHSSeconds=toc(t);report.PoissonSolveCount=1;
    assert(isequaln([nativeFlow.c_l,nativeFlow.c_omega,nativeFlow.c_r], ...
        [c.initialNativeCL,c.initialNativeCOmega,c.initialNativeCR]));
    basePhysical=physical_fields(nativeFlow);[checkRhs,baseline]=given_flow(rho,basePhysical,ops);
    assert(isequaln(checkRhs,rhsBaseline)&&isequaln(baseline.rates,[nativeFlow.c_l,nativeFlow.c_omega,nativeFlow.c_r]));
    report.givenFlowAssemblyBitwiseNativeParity=true;report.initialRatesBitwiseNativeParity=true;
    baseBoundary=ipm.field.greenBoundary(Omega,ops,1);tailBoundary=c.tailBoundary;combinedBoundary=baseBoundary;
    for key={'left','right','bottom','top'}
        name=key{1};assert(isequal(size(baseBoundary.(name)),size(tailBoundary.(name)))&&all(isfinite(tailBoundary.(name)),'all'));
        combinedBoundary.(name)=baseBoundary.(name)+tailBoundary.(name);
    end
    t=tic;[psi,~,actualBoundary,solveInfo]=ipm.field.poisson(Omega,ops,1,combinedBoundary);
    report.tailPoissonSeconds=toc(t);report.PoissonSolveCount=2;
    assert(isequaln(actualBoundary,combinedBoundary));
    u1=-(ops.Dy*psi);u2=psi*ops.Dx';u2(1,:)=0;
    solveInfo.algebraicResidual=solveInfo.relativeResidual;
    trialPhysical=struct('u1',u1,'u2',u2,'psi',psi,'source',Omega, ...
        'poissonResidual',solveInfo.relativeResidual,'poissonSolveInfo',solveInfo);
    t=tic;[rhsTrial,trial]=given_flow(rho,trialPhysical,ops);report.tailGaugeAssemblySeconds=toc(t);
    assert(isequaln(ops.rescaling,c.runtimeReferences)&&isequaln(rho,c.rho)&&isequaln(Omega,c.Omega));
    assert(isequaln(baseline.peak,trial.peak)&&baseline.peak.value==c.initialNativeQuadraticPeak);
    report.referencesRhoOmegaAndPeakUnchanged=true;report.sameCachedFactorUsed=true;
    report.boundaryTraceInf=struct();
    for key={'left','right','bottom','top'}
        name=key{1};report.boundaryTraceInf.(name)=struct('native',max(abs(baseBoundary.(name))), ...
            'exactInitialTail',max(abs(tailBoundary.(name))),'combined',max(abs(combinedBoundary.(name))));
    end
    baseline.physical=basePhysical;trial.physical=trialPhysical;
    report.baseline=compact_observation(baseline,rhsBaseline,ops,report);
    report.withExactInitialExteriorTail=compact_observation(trial,rhsTrial,ops,report);
    report.rateDifference=trial.rates-baseline.rates;
    report.coreDifference=field_difference(basePhysical,trialPhysical,rhsBaseline,rhsTrial,ops,report.coreObservationBox);
    report.holdoutDifference=field_difference(basePhysical,trialPhysical,rhsBaseline,rhsTrial,ops,report.holdoutObservationBox);
    report.wholeBoxDifference=field_difference(basePhysical,trialPhysical,rhsBaseline,rhsTrial,ops, ...
        [ops.x(1),ops.x(end),ops.y(1),ops.y(end)]);
    report.rssAfterBothBranchesBytes=rss_bytes();report.wallSecondsBeforeSave=toc(timer);
    x=ops.x;y=ops.y;Dx=ops.Dx;Dy=ops.Dy; %#ok<NASGU>
    save(fullfile(out,'fields.mat'),'x','y','rho','Omega','scale','Dx','Dy','baseline','trial', ...
        'rhsBaseline','rhsTrial','baseBoundary','tailBoundary','combinedBoundary','solveInfo','-v7.3');
    clear ops
    report.rssAfterFactorReleaseBytes=rss_bytes();report.status='instantaneous_comparison_completed';
    report.resultMeaning='At t0 only, the same original discretization/rho/Crefs with the exact analytic initial exterior Green trace added. Differences are not an evolved-tail closure or overall PDE accuracy certification.';
    save(fullfile(out,'report.mat'),'report');write_json(fullfile(out,'report.json'),report);disp(jsonencode(report));
catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack,'partialReport',report);
    save(fullfile(out,'failure.mat'),'failure','exception','-v7.3');write_json(fullfile(out,'failure.json'),failure);rethrow(exception)
end
end

function p=physical_fields(f)
p=struct('u1',f.u1,'u2',f.u2,'psi',f.psi,'source',f.source, ...
    'poissonResidual',f.poissonResidual,'poissonSolveInfo',f.poissonSolveInfo);
end
function [rhs,o]=given_flow(rho,physical,ops)
assert(ops.rescaling.enabled&&strcmp(ops.dynamicScaleGeometry,'isotropic')&& ...
    strcmp(ops.rescaling.lengthGauge,'transport_anchor')&&strcmp(ops.rescaling.cOmegaGauge,'wall_omega_quadratic_peak'));
details=ipm.diagnostics.trackFeatures(rho,[],physical,ops);
[cl,cw,cr,baseRhs,details,tx,ty]=ipm.evolve.isotropicGauge(rho,[],physical,ops,details);
rhs=baseRhs+cw*rho;rhsX=rhs*ops.Dx';
index=find(ops.x==details.omegaGaugeQuadraticStencilCenterX);assert(isscalar(index));
peak=ipm.evolve.quadraticPeakFunctional(physical.source(1,:),ops.x,index);
actualPPrime=sum(peak.weights.*rhsX(1,peak.indices));
o=struct('rates',[cl,cw,cr],'peak',peak,'actualPPrimeFromCompleteRHS',actualPPrime, ...
    'algebraicGaugePPrime',details.omegaGaugeQuadraticRateResidual, ...
    'details',details,'transportU1',tx,'transportU2',ty);
assert(all(isfinite([rhs(:);o.rates(:);actualPPrime])));
end
function o=compact_observation(v,F,ops,registration)
p=v.physical;grad1=p.u1*ops.Dx';grad2=ops.Dy*p.u1;grad3=p.u2*ops.Dx';grad4=ops.Dy*p.u2;
o=struct('canonicalAndPhysicalRatesAtIdentityScale',v.rates, ...
    'actualPPrimeFromCompleteRHS',v.actualPPrimeFromCompleteRHS,'algebraicGaugePPrime',v.algebraicGaugePPrime, ...
    'quadraticPeak',v.peak,'poissonResidual',p.poissonResidual, ...
    'velocityInf',[max(abs(p.u1),[],'all'),max(abs(p.u2),[],'all')], ...
    'velocityDerivativeInf',[max(abs(grad1),[],'all'),max(abs(grad2),[],'all'),max(abs(grad3),[],'all'),max(abs(grad4),[],'all')], ...
    'completeRHSInf',max(abs(F),[],'all'));
for name={'coreObservationBox','holdoutObservationBox'}
    box=registration.(name{1});mask=ops.X>=box(1)&ops.X<=box(2)&ops.Y>=box(3)&ops.Y<=box(4);
    o.(name{1})=struct('box',box,'pointCount',nnz(mask),'u1Inf',max(abs(p.u1(mask))), ...
        'u2Inf',max(abs(p.u2(mask))),'u1xInf',max(abs(grad1(mask))),'u1yInf',max(abs(grad2(mask))), ...
        'u2xInf',max(abs(grad3(mask))),'u2yInf',max(abs(grad4(mask))),'rhsInf',max(abs(F(mask))));
end
end
function o=field_difference(a,b,Fa,Fb,ops,box)
mask=ops.X>=box(1)&ops.X<=box(2)&ops.Y>=box(3)&ops.Y<=box(4);
dx=b.u1-a.u1;dy=b.u2-a.u2;dxx=dx*ops.Dx';dxy=ops.Dy*dx;dyx=dy*ops.Dx';dyy=ops.Dy*dy;df=Fb-Fa;
o=struct('box',box,'pointCount',nnz(mask),'u1Inf',max(abs(dx(mask))),'u2Inf',max(abs(dy(mask))), ...
    'u1xInf',max(abs(dxx(mask))),'u1yInf',max(abs(dxy(mask))),'u2xInf',max(abs(dyx(mask))), ...
    'u2yInf',max(abs(dyy(mask))),'rhsInf',max(abs(df(mask))), ...
    'velocityRelativeInf',max(hypot(dx(mask),dy(mask)))/max(hypot(a.u1(mask),a.u2(mask))));
end
function n=rss_bytes()
[status,output]=system(sprintf('/bin/ps -o rss= -p %d',feature('getpid')));
assert(status==0);n=1024*str2double(strtrim(output));assert(isfinite(n)&&n>0);
end
function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
