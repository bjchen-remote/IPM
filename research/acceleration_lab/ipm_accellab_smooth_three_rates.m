function report=ipm_accellab_smooth_three_rates(cacheFile,smoothReportFile,oldFailureReportFile,outputRoot)
%IPM_ACCELLAB_SMOOTH_THREE_RATES Fixed-flow actual global-LF research rates.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['smooth_three_rates_',token]);mkdir(destination);
p=struct('cacheFile',cacheFile,'smoothReportFile',smoothReportFile,'oldFailureReportFile',oldFailureReportFile, ...
    'outputDirectory',destination,'primaryTheta',.5,'holdoutTheta',[.3,.7],'power',4, ...
    'constraintTolerance',5e-11,'minimumRcond',1e-8,'rateDifferenceStep',1e-5, ...
    'maximumIterations',10,'maximumBacktracks',10,'directionModelRelativeTolerance',.01, ...
    'nonincreaseRelativeAllowance',1e-10,'coordinateRhsRelativeTolerance',1e-3, ...
    'fixedPhysicalWindow',[-3,3,0,3],'tailThreshold',2.5, ...
    'rateParameterization','increments q=[deltaCL,deltaCW,deltaVPeak/w]; cr=cr0+q3*w-q1*a', ...
    'newPoissonEvaluations',0,'pdeSteps',0,'productionModified',false,'peakProjectionCount',0, ...
    'jacobianContract','One-sided secant models and actual combined directional checks; no differentiability assumed.', ...
    'protocolFile',fullfile(directory,'SMOOTH_THREE_RATE_PROTOCOL.md'));
save(fullfile(destination,'registration.mat'),'p');write_json(fullfile(destination,'registration.json'),p);
records={};q=[];r=[];calls=0;F=[];
try
    d=load(cacheFile);s=load(smoothReportFile,'report');smooth=s.report;
    o=load(oldFailureReportFile,'report');old=o.report;
    assert(strcmp(old.status,'original_rate_solve_rejected_unchanged') && ~old.originalDecisionReplaced);
    assert(strcmp(d.cacheAudit.sourceCheckpoint,smooth.cases{1}.sourceCheckpoint) && strcmp(d.cacheAudit.sourceCheckpoint,old.sourceCheckpoint));
    assert(d.cacheAudit.rhsCachePairedExactly && d.cacheAudit.fullFlowAndRhoRateSaved && d.cacheAudit.noEllipticFactorSaved);
    assert(smooth.allDirectionalChecksPassed && smooth.registration.primaryTheta==.5);
    ops=d.transportOps;rho=d.rho;Omega=d.Omega;flow=d.flow;
    assert(strcmp(ops.spatialDiscretization,'high_order') && strcmp(ops.transportScheme,'weno5_fd') && ...
        strcmp(ops.wallTransportMode,'conservative_flux') && strcmp(ops.transportBoundaryMode,'open'));
    assert(isequal(Omega,rho*ops.Dx') && isequal(d.FX,d.rhoRate*ops.Dx') && isequal(Omega,flow.source));
    oldSource=load(smooth.cases{1}.sourceReport,'report');oldSource=oldSource.report;
    window=oldSource.localWindow;calibration=smooth.unitCalibration;
    assert(isequal([calibration.theta],[.5,.3,.7]));
    r0=ipm_accellab_smooth_inner_rates(Omega,d.FX,d.x,d.y,d.Dx,d.Dy,window,.5);assert(r0.valid);
    assert(isequaln(r0,smooth.cases{1}.observers{1}.rawRates));
    a=r0.peak.x;w=calibration(1).horizontalFactor*r0.rawWallWidth;
    coordinate=@(rhs) functional(Omega,rhs,ops,window,w);
    evaluate=@(v) actual(v,rho,flow,ops,coordinate,a,w);
    [nativeR,nativeF,nativeDetail]=evaluate(zeros(3,1));calls=calls+1;
    assert(isequal(nativeF,d.rhoRate),'ipm:SmoothRateNativeParity','Zero increments do not reproduce native RHS bitwise.');
    [nativeSides,n]=side_matrices(evaluate,zeros(3,1),nativeR,p.rateDifferenceStep);calls=calls+n;
    [nativeSidesHalf,n]=side_matrices(evaluate,zeros(3,1),nativeR,p.rateDifferenceStep/2);calls=calls+n;
    rx=rho*ops.Dx';ry=ops.Dy*rho;
    generators=cat(3,-(ops.X-a).*rx-ops.Y.*ry,rho,-w*rx);A=zeros(3);
    for k=1:3,A(:,k)=coordinate(generators(:,:,k));end
    assert(rcond(A)>=p.minimumRcond,'ipm:SmoothRateSeedCondition','Continuum-generator initial matrix is ill-conditioned.');
    seed=-A\nativeR;q=seed;[r,F,detail]=evaluate(q);calls=calls+1;seedR=r;seedF=F;
    stopReason='maximum_iterations';
    for iteration=1:p.maximumIterations
        if norm(r,inf)<=p.constraintTolerance,stopReason='root_tolerance';break;end
        [sides,n]=side_matrices(evaluate,q,r,p.rateDifferenceStep);calls=calls+n;
        candidates={};accepted=false;
        for side=1:2
            J=sides.forward;if side==2,J=sides.backward;end
            attempt=struct('side',side,'rcond',rcond(J),'modelMatrix',J,'accepted',false);
            if rcond(J)<p.minimumRcond
                attempt.reason='one_sided_model_condition';candidates{end+1}=attempt;continue %#ok<AGROW>
            end
            direction=-J\r;
            h=p.rateDifferenceStep*max(1,norm(q,inf))/max(1,norm(direction,inf));
            [dr1,~,dt1]=evaluate(q+h*direction);[dr2,~,dt2]=evaluate(q+(h/2)*direction);calls=calls+2;
            response1=(dr1-r)/h;response2=(dr2-r)/(h/2);predicted=J*direction;
            denominator=max([norm(predicted,2),norm(response2,2),1e-12]);
            discrepancy=max(norm(response1-predicted,2),norm(response2-predicted,2))/denominator;
            attempt.direction=direction;attempt.directionStep=h;attempt.predictedResponse=predicted;
            attempt.actualDirectionalResponses=[response1,response2];attempt.directionModelRelativeError=discrepancy;
            attempt.directionProbeSelectors={dt1.selector,dt2.selector};
            if discrepancy>p.directionModelRelativeTolerance
                attempt.reason='combined_direction_disagrees_with_one_sided_model';candidates{end+1}=attempt;continue %#ok<AGROW>
            end
            backtracks={};
            for backtrack=0:p.maximumBacktracks
                damping=2^-backtrack;candidate=q+damping*direction;
                [next,nextF,nextDetail]=evaluate(candidate);calls=calls+1;
                backtracks{end+1}=struct('damping',damping,'q',candidate,'actualResidual',next, ...
                    'selector',nextDetail.selector,'actualPPrime',nextDetail.rates.peakPrime); %#ok<AGROW>
                if norm(next,2)<norm(r,2)*(1-1e-4*damping) || norm(next,inf)<=p.constraintTolerance
                    accepted=true;break
                end
            end
            attempt.backtracks=backtracks;attempt.accepted=accepted;
            if accepted,attempt.reason='actual_residual_descent';else,attempt.reason='all_backtracks_failed';end
            candidates{end+1}=attempt; %#ok<AGROW>
            if accepted,break;end
        end
        records{end+1}=struct('iteration',iteration,'oldQ',q,'oldResidual',r,'oldSelector',detail.selector, ...
            'oneSidedModels',sides,'attempts',{candidates},'accepted',accepted); %#ok<AGROW>
        if ~accepted,stopReason='registered_one_sided_direction_or_backtracking_rejection';break;end
        q=candidate;r=next;F=nextF;detail=nextDetail;
    end
    [r,F,detail]=evaluate(q);calls=calls+1;
    if norm(r,inf)<=p.constraintTolerance,stopReason='root_tolerance';end
    [finalSides,n]=side_matrices(evaluate,q,r,p.rateDifferenceStep);calls=calls+n;
    [finalSidesHalf,n]=side_matrices(evaluate,q,r,p.rateDifferenceStep/2);calls=calls+n;
    oldFields=load(fullfile(fileparts(oldFailureReportFile),'fixed_fields.mat'),'F','FX');
    assert(isequal(oldFields.F*ops.Dx',oldFields.FX));
    oldQ=old.fixedQ;oldW=oldSource.coordinates.wallCoreWidth;
    oldRates=[oldQ(1),oldQ(2),oldQ(3)*oldW-oldQ(1)*a];
    newRates=[flow.c_l+q(1),flow.c_omega+q(2),flow.c_r+q(3)*w-q(1)*a];
    allF={nativeF,oldFields.F,F};rateRows=[flow.c_l,flow.c_omega,flow.c_r;oldRates;newRates];
    labels={'original_native_rates','old_90percent_rejected_rates_unchanged','new_smooth_terminal_rates'};
    audits=cell(1,3);residualFields=cell(3,3);nativeFields=cell(1,3);rhsDefects=cell(1,3);
    for j=1:3
        rows=cell(1,3);FX=allF{j}*ops.Dx';
        for k=1:3
            rr=ipm_accellab_smooth_inner_rates(Omega,FX,d.x,d.y,d.Dx,d.Dy,window,calibration(k).theta);assert(rr.valid);
            coordinates=chart_coordinates(rr,calibration(k));
            [metric,residualFields{j,k}]=ipm_accellab_continuous_residual(Omega,FX,d.x,d.y,d.Dx,d.Dy,coordinates);
            rows{k}=struct('theta',calibration(k).theta,'rawRates',rr,'fixedCalibration',calibration(k), ...
                'coordinates',coordinates,'residual',metric);
        end
        [nativeMetric,nativeFields{j}]=native_metric(Omega,FX,ops,rows{1}.coordinates);
        delta=rateRows(j,:)-rateRows(1,:);
        rhsDefects{j}=allF{j}-nativeF+delta(1)*(ops.X.*rx+ops.Y.*ry)-delta(2)*rho+delta(3)*rx;
        physicalX=(ops.X-d.scale.X_shift)/exp(d.scale.logC_l);physicalY=ops.Y/exp(d.scale.logC_l);
        physicalWindow=physicalX>=p.fixedPhysicalWindow(1)&physicalX<=p.fixedPhysicalWindow(2)&physicalY>=0&physicalY<=p.fixedPhysicalWindow(4);
        tail=physicalWindow & (abs(physicalX)>=p.tailThreshold | physicalY>=p.tailThreshold);assert(any(tail,'all'));
        defect=rhsDefects{j};normalizer=max(abs(nativeF),[],'all');
        defectMetric=struct('wholeBoxRelativeInf',max(abs(defect),[],'all')/normalizer, ...
            'wholeBoxRelativeL2',sqrt(sum(ops.integrationWeights.*defect.^2,'all')/sum(ops.integrationWeights.*nativeF.^2,'all')), ...
            'fixedPhysicalTailRelativeInf',max(abs(defect(tail)))/normalizer,'tailNodeCount',nnz(tail));
        audits{j}=struct('label',labels{j},'actualRates',rateRows(j,:),'parentRates', ...
            [rateRows(j,1:2)/d.cacheAudit.lineage.canonicalCovarianceFactor,rateRows(j,3)*d.cacheAudit.lineage.parentCx/d.cacheAudit.lineage.canonicalCovarianceFactor], ...
            'observers',{rows},'nativeFdModulatedResidual',nativeMetric,'coordinateRhsDifference',defectMetric, ...
            'rawDensityRhsNorms',raw_norms(allF{j},ops.integrationWeights), ...
            'rawOmegaRhsNorms',raw_norms(FX,ops.integrationWeights), ...
            'densityStateDifference',0,'physicalVelocityDifference',0,'scaleDifference',0, ...
            'stateAndVelocityDifferenceMeaning','Identical fixed arrays by construction; this is not a fresh changed-state field/tail validation.');
    end
    comparisons=cell(1,2);
    for j=2:3,comparisons{j-1}=compare_audits(audits{1},audits{j},p);end
    rootPassed=norm(r,inf)<=p.constraintTolerance;
    report=struct('status','completed_smooth_three_rate_fixed_flow_research','registration',p,'cacheAudit',d.cacheAudit, ...
        'sourceNativeRhsBitwiseReproduced',true,'sourceWallMode',ops.wallTransportMode,'sourceSymmetryMode',ops.symmetryMode, ...
        'nativeNormalizedConstraint',nativeR,'nativeSelector',nativeDetail.selector,'nativeOneSided',nativeSides,'nativeOneSidedHalfStep',nativeSidesHalf, ...
        'linearGeneratorMatrix',A,'linearGeneratorRcond',rcond(A),'seed',seed,'seedActualResidual',seedR, ...
        'records',{records},'stopReason',stopReason,'terminalQ',q,'terminalActualConstraint',r, ...
        'terminalActualRates',newRates,'terminalRawSmoothGeometry',detail.rates,'terminalSelector',detail.selector, ...
        'finalOneSided',finalSides,'finalOneSidedHalfStep',finalSidesHalf,'rootGatePassed',rootPassed, ...
        'audits',{audits},'comparisonsWithNative',{comparisons},'allNecessaryInstantaneousGatesPassed',rootPassed&&comparisons{2}.allPassed, ...
        'transportAssemblyCalls',calls,'newPoissonEvaluations',0,'pdeSteps',0,'productionQualification',false, ...
        'originalRejectedRateDecisionChanged',false,'oldRejectedReport',oldFailureReportFile, ...
        'spatialSensitivityReference','The same smooth observer retained 7--14 percent residual-field differences across the three registered late grids; this experiment does not remove that uncertainty.', ...
        'interpretation','No density, scale or physicalflow changed. Actual combined global-LF transport and complete smooth pullback residuals are evaluated at each rate. One-sided models are exploratory unless actual combined-direction checks support them; centered derivatives are diagnostic only. A root or instantaneous residual decrease does not establish acceleration, finite-step physical equivalence or production suitability.');
    save(fullfile(destination,'fields.mat'),'allF','seedF','residualFields','nativeFields','rhsDefects','rateRows','-v7.3');
    save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
        'records',{records},'lastQ',q,'lastResidual',r,'transportAssemblyCalls',calls);
    save(fullfile(destination,'failure.mat'),'failure','F','-v7.3');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('SMOOTH_THREE_RATES %s root=%d allNecessary=%d residual=%.3e calls=%d reason=%s\n', ...
    destination,rootPassed,report.allNecessaryInstantaneousGatesPassed,norm(r,inf),calls,stopReason);
end

function [v,F,d]=actual(q,rho,flow,ops,coordinate,a,w)
cl=flow.c_l+q(1);cw=flow.c_omega+q(2);cr=flow.c_r+q(3)*w-q(1)*a;
[F,~,u1,u2]=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,cl,cl,cw,cr,ops,ops.transportBoundaryMode);
[v,r]=coordinate(F);d=struct('rates',r,'selector',lf_selector(u1,u2,ops));
end
function [v,r]=functional(Omega,F,ops,window,w)
r=ipm_accellab_smooth_inner_rates(Omega,F*ops.Dx',ops.x,ops.y,ops.Dx,ops.Dy,window,.5);
assert(r.valid,'ipm:SmoothThreeGeometry','Smooth geometry invalid.');
v=[r.peakPrime/r.peak.value;r.translationRate/w;r.logScaleXRate];
end
function s=lf_selector(u1,u2,ops)
vx=u1./ops.metricX;vy=u2./ops.metricY;
[ax,ix]=max(abs(vx),[],2);[ay,iy]=max(abs(vy),[],1);
signX=sign(vx(sub2ind(size(vx),(1:size(vx,1))',ix)));
signY=sign(vy(sub2ind(size(vy),iy,1:size(vy,2))));
s=struct('xIndices',ix,'yIndices',iy,'xSigns',signX,'ySigns',signY,'xAlpha',ax,'yAlpha',ay, ...
    'wallAlpha',ax(1),'wallXIndex',ix(1),'wallXCoordinate',ops.x(ix(1)));
end
function [s,calls]=side_matrices(evaluate,q,base,h0)
plus=zeros(3);minus=zeros(3);selectors=cell(2,3);steps=zeros(3,1);
for k=1:3
    h=h0*max(1,abs(q(k)));e=zeros(3,1);e(k)=h;
    [vp,~,dp]=evaluate(q+e);[vm,~,dm]=evaluate(q-e);
    plus(:,k)=(vp-base)/h;minus(:,k)=(base-vm)/h;steps(k)=h;
    selectors{1,k}=dp.selector;selectors{2,k}=dm.selector;
end
s=struct('forward',plus,'backward',minus,'centeredDiagnosticOnly',(plus+minus)/2, ...
    'relativeOneSidedDifference',norm(plus-minus,'fro')/max([norm(plus,'fro'),norm(minus,'fro'),realmin]), ...
    'steps',steps,'probeSelectors',{selectors},'forwardRcond',rcond(plus),'backwardRcond',rcond(minus), ...
    'isCertifiedGeneralizedJacobian',false);calls=6;
end
function c=chart_coordinates(r,calibration)
c=struct('valid',true,'peak',r.peak,'peakPrime',r.peakPrime,'translationRate',r.translationRate, ...
    'wallCoreWidth',calibration.horizontalFactor*r.rawWallWidth,'verticalCoreWidth',calibration.verticalFactor*r.rawVerticalWidth, ...
    'wallCoreWidthPrime',calibration.horizontalFactor*r.rawWallWidthPrime, ...
    'verticalCoreWidthPrime',calibration.verticalFactor*r.rawVerticalWidthPrime, ...
    'logScaleXRate',r.logScaleXRate,'logScaleYRate',r.logScaleYRate,'signature',r.signature);
end
function [m,G]=native_metric(Omega,FX,ops,r)
ox=Omega*ops.Dx';oy=ops.Dy*Omega;P=r.peak.value;
G=(FX+(r.translationRate+r.logScaleXRate*(ops.X-r.peak.x)).*ox+r.logScaleYRate*ops.Y.*oy)/P-r.peakPrime/P*(Omega/P);
XI=(ops.X-r.peak.x)/r.wallCoreWidth;ETA=ops.Y/r.verticalCoreWidth;
core=abs(XI)<=1 & ETA<=1.5;hold=abs(XI)<=2 & ETA<=3 & ~core;
m=struct('wholeBox',raw_norms(G,ops.integrationWeights),'core',masked_norms(G,ops.integrationWeights,core), ...
    'holdout',masked_norms(G,ops.integrationWeights,hold));
end
function m=masked_norms(F,w,mask)
assert(any(mask,'all'));m=raw_norms(F(mask),w(mask));
end
function m=raw_norms(F,w)
m=struct('L2',sqrt(sum(w.*F.^2,'all')/sum(w,'all')),'infinity',max(abs(F),[],'all'));
end
function comparison=compare_audits(base,candidate,p)
rows={};checks=true(0,1);
for k=1:3
    for kind={'gauss','wallGauss'}
        field=kind{1};b=base.observers{k}.residual.(field).G;c=candidate.observers{k}.residual.(field).G;
        for name=fieldnames(b)'
            key=name{1};pass=c.(key)<=b.(key)*(1+p.nonincreaseRelativeAllowance);
            rows{end+1}=struct('theta',base.observers{k}.theta,'kind',field,'metric',key,'baseline',b.(key),'candidate',c.(key), ...
                'ratio',c.(key)/max(b.(key),realmin),'passed',pass);checks(end+1)=pass; %#ok<AGROW>
        end
    end
end
fdChecks={};
for region={'wholeBox','core','holdout'}
    reg=region{1};b=base.nativeFdModulatedResidual.(reg);c=candidate.nativeFdModulatedResidual.(reg);
    for name={'L2','infinity'}
        key=name{1};pass=c.(key)<=b.(key)*(1+p.nonincreaseRelativeAllowance);
        fdChecks{end+1}=struct('region',reg,'metric',key,'baseline',b.(key),'candidate',c.(key),'ratio',c.(key)/max(b.(key),realmin),'passed',pass);checks(end+1)=pass; %#ok<AGROW>
    end
end
rhsPassed=candidate.coordinateRhsDifference.wholeBoxRelativeInf<=p.coordinateRhsRelativeTolerance;
tailPassed=candidate.coordinateRhsDifference.fixedPhysicalTailRelativeInf<=p.coordinateRhsRelativeTolerance;
comparison=struct('label',candidate.label,'smoothResidualChecks',{rows},'nativeFdChecks',{fdChecks}, ...
    'allResidualNonincreasePassed',all(checks),'wholeBoxCoordinateRhsGatePassed',rhsPassed,'tailCoordinateRhsGatePassed',tailPassed, ...
    'allPassed',all(checks)&&rhsPassed&&tailPassed,'freshChangedStateTailVelocityTestPerformed',false);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
