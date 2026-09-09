function report=ipm_accellab_three_rates_from_cache(cacheFile,observationReportFile,outputRoot)
%IPM_ACCELLAB_THREE_RATES_FROM_CACHE Actual combined transport, no Poisson.
% This is an instantaneous fixed-physical-flow diagnostic, not a new IVP.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['three_rates_cached_',token]);mkdir(destination);
p=struct('cacheFile',cacheFile,'observationReportFile',observationReportFile,'outputDirectory',destination, ...
    'constraintTolerance',5e-11,'minimumRcond',1e-8,'rateDifferenceStep',1e-5, ...
    'maximumNewtonIterations',12,'maximumBacktracks',10,'newPoissonEvaluations',0,'pdeSteps',0, ...
    'constraints','Actual PPrime/P, aPrime/wx, beta of F(rates); gamma remains measured.', ...
    'rateParameterization','q=[cl,cw,(cr+cl*a)/wx]; actual cX=cY=cl, cOmega=cw, cR=q3*wx-cl*a.', ...
    'sourceSymmetryChanged',false,'candidateWrittenAsCheckpoint',false,'peakProjectionCount',0, ...
    'supportedWallModes',{{'advective_upwind','conservative_flux'}},'sourceWallModeChanged',false);
save(fullfile(destination,'registration.mat'),'p');write_json(fullfile(destination,'registration.json'),p);
records={};q=[];r=[];calls=0;
try
    d=load(cacheFile);prior=load(observationReportFile,'report');prior=prior.report;
    assert(strcmp(prior.status,'completed_audited_cache_no_LU') && isequaln(prior.cacheAudit,d.cacheAudit));
    assert(d.cacheAudit.fullFlowAndRhoRateSaved && d.cacheAudit.noEllipticFactorSaved);
    ops=d.transportOps;rho=d.rho;flow=d.flow;Omega=d.Omega;window=prior.localWindow;
    assert(strcmp(ops.spatialDiscretization,'high_order') && strcmp(ops.transportScheme,'weno5_fd') && ...
        strcmp(ops.transportBoundaryMode,'open') && ismember(ops.wallTransportMode,{'advective_upwind','conservative_flux'}));
    assert(isequal(ops.x,d.x) && isequal(ops.y,d.y) && isequal(ops.Dx,d.Dx) && isequal(ops.Dy,d.Dy));
    assert(isequal(Omega,flow.source) && isequal(Omega,rho*ops.Dx') && isequal(d.FX,d.rhoRate*ops.Dx'));
    native=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,flow.c_l,flow.c_l,flow.c_omega,flow.c_r,ops,ops.transportBoundaryMode);
    nativeReproduction=struct('bitwise',isequal(native,d.rhoRate),'maximumAbsoluteDifference',max(abs(native-d.rhoRate),[],'all'));
    assert(nativeReproduction.bitwise,'ipm:CachedTransportParity','Cached native complete RHS was not reproduced bitwise.');
    r0=ipm_accellab_continuous_inner_rates(Omega,d.FX,d.x,d.y,d.Dx,d.Dy,window);assert(r0.valid);
    assert(isequaln(r0,prior.coordinates),'ipm:CachedGeometryParity','The original complete geometry changed.');
    a=r0.peak.x;w=r0.wallCoreWidth;
    coordinate=@(F) functional(Omega,F,ops,window,w);
    F0=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,0,0,0,0,ops,ops.transportBoundaryMode);
    b=coordinate(F0);rx=rho*ops.Dx';ry=ops.Dy*rho;
    generators=cat(3,-(ops.X-a).*rx-ops.Y.*ry,rho,-w*rx);A=zeros(3);
    for k=1:3,rr=coordinate(generators(:,:,k));A(:,k)=rr(1:3);end
    assert(rcond(A)>=p.minimumRcond,'ipm:CachedLinearCondition','Centered generator Jacobian is ill-conditioned.');
    seed=-A\b(1:3);q=seed;evaluate=@(v) actual(v,rho,flow,ops,coordinate,a,w);
    [r,seedF,~]=evaluate(q);calls=calls+1;seedResidual=r;
    for iteration=1:p.maximumNewtonIterations
        if norm(r,inf)<=p.constraintTolerance,break;end
        [J,n]=jacobian(evaluate,q,p.rateDifferenceStep);calls=calls+n;
        assert(rcond(J)>=p.minimumRcond,'ipm:CachedActualCondition','Actual centered rate Jacobian is ill-conditioned.');
        direction=-J\r;attempts={};accepted=false;
        for backtrack=0:p.maximumBacktracks
            damping=2^-backtrack;candidate=q+damping*direction;[next,~,detail]=evaluate(candidate);calls=calls+1;
            attempts{end+1}=struct('damping',damping,'q',candidate,'residual',next, ...
                'actualPPrime',detail.rates.peakPrime,'gamma',detail.rates.logScaleYRate); %#ok<AGROW>
            if norm(next,2)<norm(r,2)*(1-1e-4*damping) || norm(next,inf)<=p.constraintTolerance
                accepted=true;break
            end
        end
        records{end+1}=struct('iteration',iteration,'oldQ',q,'oldResidual',r,'jacobian',J, ...
            'rcond',rcond(J),'direction',direction,'attempts',{attempts},'accepted',accepted); %#ok<AGROW>
        assert(accepted,'ipm:CachedRateLineSearch','All registered actual transport Newton backtracks failed.');
        q=candidate;r=next;
    end
    [r,F,detail]=evaluate(q);calls=calls+1;
    [J1,n1]=jacobian(evaluate,q,p.rateDifferenceStep);[J2,n2]=jacobian(evaluate,q,p.rateDifferenceStep/2);calls=calls+n1+n2;
    jDifference=norm(J1-J2,'fro')/max(norm(J2,'fro'),realmin);
    cr=q(3)*w-q(1)*a;
    [candidateResidual,candidateDetails]=ipm_accellab_continuous_residual(Omega,F*ops.Dx',d.x,d.y,d.Dx,d.Dy,detail.rates);
    lambda=d.cacheAudit.lineage.canonicalCovarianceFactor;cx0=d.cacheAudit.lineage.parentCx;
    delta=F-native;reconstructed=F+q(1)*ops.X.*rx+q(1)*ops.Y.*ry-q(2)*rho+cr*rx;
    defect=reconstructed-F0;weights=ops.integrationWeights;
    rate=struct('cl',q(1),'cw',q(2),'cr',cr);
    rateParent=struct('cl',q(1)/lambda,'cw',q(2)/lambda,'cr',cr*cx0/lambda, ...
        'remainingGamma',detail.rates.logScaleYRate/lambda,'actualPPrime',detail.rates.peakPrime/lambda^2);
    comparison=struct();
    for key={'coreL2','holdoutL2','fullLocalL2'}
        name=key{1};comparison.(name)=struct('nativeParent',prior.residual.gauss.G.(name)/lambda, ...
            'threeGaugeParent',candidateResidual.gauss.G.(name)/lambda, ...
            'ratio',candidateResidual.gauss.G.(name)/prior.residual.gauss.G.(name));
    end
    report=struct('status','completed_fixed_flow_actual_three_rate_probe','registration',p, ...
        'cacheAudit',d.cacheAudit,'nativeReproduction',nativeReproduction,'sourceSymmetryMode',ops.symmetryMode, ...
        'sourceWallTransportMode',ops.wallTransportMode, ...
        'sourceGeometry',r0,'sourceNativeRates',d.nativeRates,'sourcePhysicalNormalWallVelocityInfinity',max(abs(flow.u2(1,:))), ...
        'sourceDensitySymmetryRelativeInf',max(abs(rho-fliplr(rho)),[],'all')/max(abs(rho),[],'all'), ...
        'linearGeneratorMatrix',A,'linearGeneratorRcond',rcond(A),'linearSeed',seed,'seedActualResidual',seedResidual, ...
        'newtonRecords',{records},'actualRates',rate,'parentRates',rateParent,'actualConstraintResidual',r, ...
        'actualContinuousGeometry',detail.rates,'finalJacobian',J2,'finalJacobianRcond',rcond(J2), ...
        'jacobianRelativeDifference',jDifference,'rateAssemblyCalls',calls,'fixedFlowAssemblyOutsideNewton',2, ...
        'newPoissonEvaluations',0,'pdeSteps',0,'candidateResidual',candidateResidual,'residualComparison',comparison, ...
        'coordinateChainDefect',struct('wholeBoxRelativeInf',max(abs(defect),[],'all')/max(abs(F0),[],'all'), ...
            'wholeBoxRelativeL2',sqrt(sum(weights.*defect.^2,'all')/sum(weights.*F0.^2,'all'))), ...
        'completeRhsRelativeChange',max(abs(delta),[],'all')/max(abs(native),[],'all'), ...
        'threeGaugeDensityRhsEvenSymmetryDefect',max(abs(F-fliplr(F)),[],'all')/max(abs(F),[],'all'), ...
        'transportNormalWallVelocityInfinity',max(abs(detail.u2(1,:))), ...
        'passedInstantaneousConstraintProtocol',norm(r,inf)<=p.constraintTolerance && rcond(J2)>=p.minimumRcond && jDifference<1e-5, ...
        'productionPromotionAllowed',false,'futureIsotropicLuReuseAspect',1, ...
        'interpretation','The native physical velocity is fixed, with the actual nonlinear WENO assembly recomputed for every rate. This tests instantaneous coordinate feasibility only. Translation breaks the original double-odd trajectory restriction; a new half-plane physical IVP and matched physical-time/boundary/velocity-tail validation would be required. No new Poisson means a changed-state tail velocity has not been evaluated. Nonzero coordinate consistency and shape residual differences are retained, not called acceleration.');
    save(fullfile(destination,'fields.mat'),'F','F0','seedF','native','defect','delta','candidateDetails','-v7.3');
    save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
        'newtonRecords',{records},'lastQ',q,'lastResidual',r,'rateAssemblyCalls',calls,'newPoissonEvaluations',0,'pdeSteps',0);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('THREE_CACHED %s pass=%d residual=%.3e gammaParent=%.9g calls=%d coreRatio=%.9g holdRatio=%.9g\n', ...
    destination,report.passedInstantaneousConstraintProtocol,norm(r,inf),rateParent.remainingGamma,calls,comparison.coreL2.ratio,comparison.holdoutL2.ratio);
end
function [r,F,d]=actual(q,rho,flow,ops,coordinate,a,w)
[F,~,u1,u2]=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,q(1),q(1),q(2),q(3)*w-q(1)*a,ops,ops.transportBoundaryMode);
[v,rates]=coordinate(F);r=v(1:3);d=struct('rates',rates,'u1',u1,'u2',u2);
end
function [v,r]=functional(Omega,F,ops,window,w)
r=ipm_accellab_continuous_inner_rates(Omega,F*ops.Dx',ops.x,ops.y,ops.Dx,ops.Dy,window);
assert(r.valid,'ipm:CachedGeometry','Continuous rates cannot be evaluated.');
v=[r.peakPrime/r.peak.value;r.translationRate/w;r.logScaleXRate;r.logScaleYRate];
end
function [J,n]=jacobian(evaluate,q,h0)
J=zeros(3);n=0;
for k=1:3
    h=h0*max(1,abs(q(k)));unit=zeros(3,1);unit(k)=h;
    plus=evaluate(q+unit);minus=evaluate(q-unit);n=n+2;J(:,k)=(plus-minus)/(2*h);
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
