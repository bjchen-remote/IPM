function report=ipm_accellab_audit_local_lf_cache(cacheFile,smoothReportFile,costReportFile,outputRoot)
%IPM_ACCELLAB_AUDIT_LOCAL_LF_CACHE Original-rate full-array operator audit.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['local_lf_2d_native_rate_',token]);mkdir(destination);
p=struct('cacheFile',cacheFile,'smoothReportFile',smoothReportFile,'costReportFile',costReportFile, ...
    'outputDirectory',destination,'primaryTheta',.5,'holdoutTheta',[.3,.7], ...
    'telescopingRelativeTolerance',5e-11,'symmetryNativeMultiplier',5,'symmetryRelativeAllowance',1e-8, ...
    'localAssemblySecondsCeiling',180,'normalization','conditional_global','rateChanges',0, ...
    'newPoissonEvaluations',0,'pdeSteps',0,'productionModified',false,'oldDecisionsChanged',false);
save(fullfile(destination,'registration.mat'),'p');write_json(fullfile(destination,'registration.json'),p);
Fglobal=[];Flocal=[];
try
    q=load(costReportFile,'report');cost=q.report;assert(cost.withinRegisteredCeiling && strcmp(cost.registration.cacheFile,cacheFile));
    d=load(cacheFile);s=load(smoothReportFile,'report');s=s.report;
    assert(d.cacheAudit.rhsCachePairedExactly && d.cacheAudit.fullFlowAndRhoRateSaved && ...
        strcmp(d.cacheAudit.sourceCheckpoint,s.cases{1}.sourceCheckpoint));
    ops=d.transportOps;rho=d.rho;rates=[d.flow.c_l,d.flow.c_omega,d.flow.c_r];
    assert(strcmp(ops.wallTransportMode,'conservative_flux') && strcmp(ops.transportBoundaryMode,'open'));
    assert(isequal(ops.x(:)',-flip(ops.x(:)')), ...
        'ipm:LocalLfNativeReflection','Source x axes are not exactly reflection paired.');
    assert(isequal(d.Omega,rho*ops.Dx') && isequal(d.FX,d.rhoRate*ops.Dx'));
    timer=tic;[Fglobal,globalDetails]=ipm_accellab_local_lf_rhs(rho,d.flow,rates,ops,'line');globalSeconds=toc(timer);
    globalBitwise=isequal(Fglobal,d.rhoRate);
    preflight=struct('globalBitwiseNative',globalBitwise,'maximumDifference',max(abs(Fglobal-d.rhoRate),[],'all'),'seconds',globalSeconds);
    save(fullfile(destination,'global_preflight.mat'),'preflight');write_json(fullfile(destination,'global_preflight.json'),preflight);
    assert(globalBitwise,'ipm:LocalLfFullNativeParity','The full two-dimensional line-mode assembly is not bitwise native.');
    fprintf('LOCAL_LF_GLOBAL_PARITY passed=1 seconds=%.3f\n',globalSeconds);
    timer=tic;[Flocal,localDetails]=ipm_accellab_local_lf_rhs(rho,d.flow,rates,ops,'face');localSeconds=toc(timer);
    FXglobal=Fglobal*ops.Dx';FXlocal=Flocal*ops.Dx';deltaF=Flocal-Fglobal;deltaFX=FXlocal-FXglobal;
    globalFlux=flux_balance(globalDetails,ops);localFlux=flux_balance(localDetails,ops);
    fluxPassed=max([globalFlux.maximumRelativeError,localFlux.maximumRelativeError])<=p.telescopingRelativeTolerance;
    nativeSymmetry=max(abs(Fglobal-fliplr(Fglobal)),[],'all');localSymmetry=max(abs(Flocal-fliplr(Flocal)),[],'all');
    rhsScale=max(abs(Fglobal),[],'all');symmetryBound=p.symmetryNativeMultiplier*nativeSymmetry+p.symmetryRelativeAllowance*rhsScale;
    symmetry=struct('nativeRhsAbsoluteDefect',nativeSymmetry,'localRhsAbsoluteDefect',localSymmetry,'registeredAbsoluteBound',symmetryBound, ...
        'passed',localSymmetry<=symmetryBound,'sourceRhoExactlyEven',isequal(rho,fliplr(rho)),'xAxesExactlyPaired',true, ...
        'sourceRhoEvenAbsoluteDefect',max(abs(rho-fliplr(rho)),[],'all'), ...
        'physicalU1OddDefect',max(abs(d.flow.u1+fliplr(d.flow.u1)),[],'all'), ...
        'physicalU2EvenDefect',max(abs(d.flow.u2-fliplr(d.flow.u2)),[],'all'), ...
        'sourceOmegaOddDefect',max(abs(d.Omega+fliplr(d.Omega)),[],'all'));
    old=load(s.cases{1}.sourceReport,'report');window=old.report.localWindow;
    observers=cell(2,3);residualFields=cell(2,3);
    for mode=1:2
        FX=FXglobal;if mode==2,FX=FXlocal;end
        for k=1:3
            calibration=s.unitCalibration(k);
            rr=ipm_accellab_smooth_inner_rates(d.Omega,FX,d.x,d.y,d.Dx,d.Dy,window,calibration.theta);assert(rr.valid);
            if mode==1,assert(isequaln(rr,s.cases{1}.observers{k}.rawRates));end
            coordinates=chart_coordinates(rr,calibration);
            [residual,residualFields{mode,k}]=ipm_accellab_continuous_residual(d.Omega,FX,d.x,d.y,d.Dx,d.Dy,coordinates);
            lambda=d.cacheAudit.lineage.canonicalCovarianceFactor;
            observers{mode,k}=struct('theta',calibration.theta,'rawRates',rr,'calibration',calibration, ...
                'residual',residual,'parentRates',[rr.peakPrime/rr.peak.value,rr.translationRate/coordinates.wallCoreWidth,rr.logScaleXRate,rr.logScaleYRate]/lambda, ...
                'parentCoreResidualL2',residual.gauss.G.coreL2/lambda,'parentHoldoutResidualL2',residual.gauss.G.holdoutL2/lambda);
        end
    end
    differences=cell(1,3);
    for k=1:3
        b=residualFields{1,k};l=residualFields{2,k};
        assert(isequal(b.quadratureX,l.quadratureX) && isequal(b.quadratureY,l.quadratureY) && isequal(b.quadratureWeights,l.quadratureWeights));
        [XI,ETA]=meshgrid(b.quadratureX,b.quadratureY);core=abs(XI)<1 & ETA<1.5;
        differences{k}=struct('theta',s.unitCalibration(k).theta,'exactSameGaussQueries',true, ...
            'metrics',region_difference(b.quadratureFields.G,l.quadratureFields.G,b.quadratureWeights,core));
    end
    cr=s.cases{1}.observers{1}.coordinatesUsedForResidual;XI=(ops.X-cr.peak.x)/cr.wallCoreWidth;ETA=ops.Y/cr.verticalCoreWidth;
    core=abs(XI)<=1 & ETA<=1.5;hold=abs(XI)<=2 & ETA<=3 & ~core;
    spatial=struct('densityRhs',native_difference(Fglobal,Flocal,ops.integrationWeights,core,hold), ...
        'omegaRhs',native_difference(FXglobal,FXlocal,ops.integrationWeights,core,hold), ...
        'largestDensityChange',largest(deltaF,ops,d.scale),'largestOmegaChange',largest(deltaFX,ops,d.scale), ...
        'differentiatedDifferenceRoundoff',max(abs(deltaF*ops.Dx'-deltaFX),[],'all'));
    boundary=struct('physicalNormalWallVelocityInfinity',max(abs(d.flow.u2(1,:))), ...
        'transportNormalWallVelocityInfinity',max(abs(localDetails.transportU2(1,:))), ...
        'lowerYDensityFaceDifferenceInfinity',max(abs(localDetails.yDensityFlux.faceFlux(:,1)-globalDetails.yDensityFlux.faceFlux(:,1))), ...
        'upperYDensityFaceDifferenceInfinity',max(abs(localDetails.yDensityFlux.faceFlux(:,end)-globalDetails.yDensityFlux.faceFlux(:,end))), ...
        'leftXDensityFaceDifferenceInfinity',max(abs(localDetails.xDensityFlux.faceFlux(:,1)-globalDetails.xDensityFlux.faceFlux(:,1))), ...
        'rightXDensityFaceDifferenceInfinity',max(abs(localDetails.xDensityFlux.faceFlux(:,end)-globalDetails.xDensityFlux.faceFlux(:,end))), ...
        'nativeQuadratureMassRhs',sum(ops.integrationWeights.*Fglobal,'all'),'localQuadratureMassRhs',sum(ops.integrationWeights.*Flocal,'all'), ...
        'meaning','Numerical lower-y face is a ghost-cell face, not the physical wall plane. The original conservative-wall y term and open-boundary/free-stream-corrected assembly remain intact; no zero-total-mass claim is made.');
    report=struct('status','completed_fixed_native_rate_local_alpha_operator_audit','registration',p,'cacheAudit',d.cacheAudit, ...
        'nativeRates',rates,'globalModeBitwiseNative',globalBitwise,'globalSeconds',globalSeconds,'localSeconds',localSeconds, ...
        'timingCeilingPassed',localSeconds<=p.localAssemblySecondsCeiling,'globalFluxBalance',globalFlux,'localFluxBalance',localFlux, ...
        'fluxContractsPassed',fluxPassed,'symmetry',symmetry,'boundary',boundary,'spatialRhsDifference',spatial, ...
        'observers',{observers},'smoothResidualDifferences',{differences}, ...
        'auditContractsPassed',globalBitwise&&fluxPassed&&symmetry.passed&&localSeconds<=p.localAssemblySecondsCeiling, ...
        'productionQualification',false,'interpretation','One native field/physicalflow/rate state, two spatial operators. Actual PPrime and smooth geometry rates are recomputed from each full RHS. This quantifies discretization sensitivity, not the continuum truth, a new trajectory, a gauge root or an acceleration. Native C and all earlier failed decisions remain unchanged.');
    save(fullfile(destination,'fields.mat'),'Fglobal','Flocal','FXglobal','FXlocal','deltaF','deltaFX','globalDetails','localDetails','residualFields','-v7.3');
    save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure','Fglobal','Flocal','-v7.3');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('LOCAL_LF_2D_AUDIT %s contracts=%d localSeconds=%.3f primaryCoreDifference=%.8g primaryHoldDifference=%.8g\n', ...
    destination,report.auditContractsPassed,localSeconds,differences{1}.metrics.core.relativeL2Difference,differences{1}.metrics.holdout.relativeL2Difference);
end
function report=flux_balance(d,ops)
keys={'xDensityFlux','yDensityFlux','xUnitFlux','yUnitFlux'};
derivatives={d.xDensityFluxDerivative,d.yDensityFluxDerivative',d.xVelocityDivergence,d.yVelocityDivergence'};
metric={ops.metricX,ops.metricY',ops.metricX,ops.metricY'};
spacing=[ops.computationalSpacingX,ops.computationalSpacingY,ops.computationalSpacingX,ops.computationalSpacingY];
rows=cell(1,4);maximum=0;
for k=1:4
    f=d.(keys{k}).faceFlux;lhs=sum(derivatives{k}.*metric{k},2)*spacing(k);rhs=f(:,end)-f(:,1);
    error=abs(lhs-rhs)./max(1,max(abs(f),[],2));maximum=max(maximum,max(error));
    rows{k}=struct('kernel',keys{k},'maximumRelativeError',max(error),'maximumAbsoluteError',max(abs(lhs-rhs)), ...
        'leftFlux',f(:,1),'rightFlux',f(:,end),'mappedDerivativeSum',lhs);
end
report=struct('rows',{rows},'maximumRelativeError',maximum,'meaning','Mapped one-dimensional flux telescoping before free-stream/source conversion; not conservation of open-boundary physical density mass.');
end
function c=chart_coordinates(r,calibration)
c=struct('valid',true,'peak',r.peak,'peakPrime',r.peakPrime,'translationRate',r.translationRate, ...
    'wallCoreWidth',calibration.horizontalFactor*r.rawWallWidth,'verticalCoreWidth',calibration.verticalFactor*r.rawVerticalWidth, ...
    'wallCoreWidthPrime',calibration.horizontalFactor*r.rawWallWidthPrime, ...
    'verticalCoreWidthPrime',calibration.verticalFactor*r.rawVerticalWidthPrime, ...
    'logScaleXRate',r.logScaleXRate,'logScaleYRate',r.logScaleYRate,'signature',r.signature);
end
function report=region_difference(b,l,w,core)
report=struct();
for name={'core','holdout','fullLocal'}
    key=name{1};mask=true(size(core));if strcmp(key,'core'),mask=core;elseif strcmp(key,'holdout'),mask=~core;end
    ww=w(mask);bb=b(mask);ll=l(mask);
    report.(key)=struct('baselineL2',l2(bb,ww),'localL2',l2(ll,ww),'normRatio',l2(ll,ww)/l2(bb,ww), ...
        'relativeL2Difference',l2(ll-bb,ww)/l2(bb,ww),'differenceSampledInfinity',max(abs(ll-bb)));
end
end
function report=native_difference(b,l,w,core,hold)
report=struct();
for name={'wholeBox','core','holdout'}
    key=name{1};mask=true(size(core));if strcmp(key,'core'),mask=core;elseif strcmp(key,'holdout'),mask=hold;end
    ww=w(mask);bb=b(mask);ll=l(mask);
    report.(key)=struct('baselineL2',l2(bb,ww),'differenceRelativeL2',l2(ll-bb,ww)/max(l2(bb,ww),realmin), ...
        'differenceRelativeInf',max(abs(ll-bb))/max(max(abs(bb)),realmin),'differenceAbsoluteInf',max(abs(ll-bb)));
end
end
function r=largest(delta,ops,scale)
[value,i]=max(abs(delta),[],'all','linear');[row,column]=ind2sub(size(delta),i);
r=struct('absoluteChange',value,'signedChange',delta(i),'row',row,'column',column, ...
    'computationalX',ops.x(column),'computationalY',ops.y(row), ...
    'physicalX',(ops.x(column)-scale.X_shift)/exp(scale.logC_l),'physicalY',ops.y(row)/exp(scale.logC_l));
end
function v=l2(f,w)
v=sqrt(sum(w.*f.^2)/sum(w));
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
