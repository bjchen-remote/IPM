function report = ipm_accellab_test_jf_profile(fineDirectory,outputRoot,constraintMode)
%IPM_ACCELLAB_TEST_JF_PROFILE Projection-free full-field JF profile trial.
% Independent tiny profile candidates, never accepted native trajectory/CP.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
if nargin<3,constraintMode='unconstrained';end
assert(any(strcmp(constraintMode,{'unconstrained','tail_ellipsoid'})));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['jf_profile_',token]);mkdir(destination);
p=struct('grids',[49,25;65,33],'krylovDimension',3,'differenceSteps',[1e-5,5e-6], ...
    'trustRadius',.002,'svdRelativeFloor',1e-8,'minimumDecrease',.001, ...
    'maximumJacobianRelativeDifference',.02,'minimumPredictionAgreement',.1, ...
    'tailDensityTolerance',1e-3,'tailSourceTolerance',1e-3,'tailVelocityTolerance',1e-3, ...
    'maximumPhysicalPeakChange',1e-3,'maximumObservationWindowDensityChange',.003, ...
    'constraintMode',constraintMode,'linearTailBudgetFraction',.8, ...
    'physicalTarget',.10,'relaxationDt',[.001,.0005],'naturalBudgetDt',.001, ...
    'fineDirectory',fineDirectory,'outputDirectory',destination, ...
    'pde','Unchanged original quadratic/anchor PDE and original SSPRK54; no P projection or PPrime cancellation.', ...
    'construction','Whole-field F,JF,J2F Krylov directions with fixed-radius Gauss-Newton of G/candidateP; no local patch or field projection.', ...
    'acceptance','Every initial/final residual and tail gate, time refinement, both grids, and equal-RHS natural-progress gate must pass.', ...
    'matchingWindow','Fixed physical [-3,3] x [0,3], tail |x|>=2.5 or y>=2; not the entire computational box.', ...
    'relaxationStatus','Independent original-equation RK research states with finite/CFL checks, not a claimed trusted native continuation.');
write_json(fullfile(destination,'registration.json'),p);save(fullfile(destination,'registration.mat'),'p');
cases=cell(1,2);
for g=1:2
    caseDir=fullfile(destination,sprintf('n%d',p.grids(g,1)));mkdir(caseDir);
    try
        sourceFile=fullfile(fineDirectory,sprintf('native_quadratic_n%d_level3.mat',p.grids(g,1)));
        loaded=load(sourceFile,'run');run=loaded.run;
        assert(strcmp(run.mode,'native_quadratic') && strcmp(run.researchEquationId,'unmodified_native_quadratic_equation'));
        opts=flat(run.operatorSourceConfig);opts.saveResults=false;opts.storeSnapshots=false;
        opts.makePlots=false;opts.livePlot=false;opts.writeVideo=false;opts.verbose=false;
        state=ipm.evolve.initialize(opts);ops=state.ops;
        assert(isequal(run.x,ops.x) && isequal(run.y,ops.y) && isequal(run.Dx,ops.Dx) && isequal(run.Dy,ops.Dy));
        assert(ops.nx<=65 && ops.ny<=33 && ~state.config.remesh.adaptiveRemesh);
        rho=run.rho;z=run.z;clear run loaded state
        [base,baseData]=observe(rho,z,ops,[]);
        [candidate,trial,trialFields]=construct(rho,z,ops,base,baseData,p);
        save(fullfile(caseDir,'candidate_fields.mat'),'rho','z','candidate','trialFields','-v7.3');
        write_json(fullfile(caseDir,'instantaneous_trial.json'),trial);
        rhsCount=trial.constructionRhsEvaluations;
        naturalSteps=ceil(rhsCount/5);naturalRho=rho;naturalZ=z;naturalCache=[];
        for k=1:naturalSteps
            [naturalRho,naturalZ,naturalCache]=native_step(naturalRho,naturalZ,p.naturalBudgetDt,ops,naturalCache);
        end
        natural=observe(naturalRho,naturalZ,ops,[]);
        budget=struct('constructionRhsEvaluations',rhsCount,'naturalSteps',naturalSteps, ...
            'naturalStepRhsEvaluations',5*naturalSteps+double(naturalSteps>0), ...
            'naturalCanonicalIncrement',naturalZ(5)-z(5),'naturalPhysicalIncrement',naturalZ(3)-z(3), ...
            'candidateCoreOverNatural',trial.candidate.metrics.G.coreL2/natural.metrics.G.coreL2, ...
            'candidateIntrinsicCoreOverNatural',trial.candidate.metrics.intrinsicCore/natural.metrics.intrinsicCore, ...
            'passed',trial.candidate.metrics.G.coreL2<=natural.metrics.G.coreL2*(1-p.minimumDecrease) && ...
                trial.candidate.metrics.intrinsicCore<=natural.metrics.intrinsicCore*(1-p.minimumDecrease), ...
            'interpretation','Distinct-time profile-cost benchmark; no candidate physical clock is advanced or relabeled by this test.');
        relaxation=cell(1,2);endFields=cell(2,2);
        for level=1:2
            dt=p.relaxationDt(level);
            [baselineR,baselineZ,baselineTime]=physical_endpoint(rho,z,ops,p.physicalTarget,dt);
            [candidateR,candidateZ,candidateTime]=physical_endpoint(candidate,z,ops,p.physicalTarget,dt);
            [bo,bd]=observe(baselineR,baselineZ,ops,[]);
            [co,cd]=observe(candidateR,candidateZ,ops,[]);
            match=matching(cd,bd,p);gates=accept(co,bo,match,p);
            relaxation{level}=struct('dt',dt,'baseline',bo,'candidate',co,'matching',match,'gates',gates, ...
                'baselineTime',baselineTime,'candidateTime',candidateTime, ...
                'samePhysicalTarget',p.physicalTarget,'independentInitialValueBranch',true, ...
                'isOriginalInitialValueTrajectory',false);
            endFields{level,1}=struct('rho',baselineR,'z',baselineZ,'physical',bd.physical);
            endFields{level,2}=struct('rho',candidateR,'z',candidateZ,'physical',cd.physical);
        end
        refinement=struct('baselinePhysicalDensityInf',relative_inf(endFields{1,1}.physical.rho,endFields{2,1}.physical.rho), ...
            'candidatePhysicalDensityInf',relative_inf(endFields{1,2}.physical.rho,endFields{2,2}.physical.rho), ...
            'baselineCoreResidualRelative',abs(relaxation{1}.baseline.metrics.G.coreL2/relaxation{2}.baseline.metrics.G.coreL2-1), ...
            'candidateCoreResidualRelative',abs(relaxation{1}.candidate.metrics.G.coreL2/relaxation{2}.candidate.metrics.G.coreL2-1));
        refinement.passed=max(refinement.baselineCoreResidualRelative,refinement.candidateCoreResidualRelative)<p.minimumDecrease/10 && ...
            max(refinement.baselinePhysicalDensityInf,refinement.candidatePhysicalDensityInf)<1e-6;
        accepted=trial.gates.passed && trial.linearization.passed && budget.passed && ...
            all(cellfun(@(v)v.gates.passed,relaxation)) && refinement.passed;
        cases{g}=struct('status','completed','gridSize',p.grids(g,:),'sourceFile',sourceFile, ...
            'sourceTime',z,'base',base,'trial',trial,'naturalBudget',budget,'relaxation',{relaxation}, ...
            'timeRefinement',refinement,'acceptedIndependentProfileCandidate',accepted, ...
            'productionAccepted',false,'oldSecantDecisionChanged',false);
        save(fullfile(caseDir,'endpoint_fields.mat'),'endFields','-v7.3');
        value=cases{g};save(fullfile(caseDir,'report.mat'),'value');write_json(fullfile(caseDir,'report.json'),value);
        fprintf('JF_PROFILE n%d accepted=%d coreRatio=%.8g tailU=%.8g naturalRatio=%.8g relaxation=%.8g/%.8g\n', ...
            p.grids(g,1),accepted,trial.gates.ratios.coreL2,trial.matching.tailVelocityRelativeL2, ...
            budget.candidateCoreOverNatural,relaxation{1}.gates.ratios.coreL2,relaxation{2}.gates.ratios.coreL2);
        clear ops rho candidate baseData trialFields endFields
    catch exception
        cases{g}=struct('status','failed','identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
            'acceptedIndependentProfileCandidate',false,'productionAccepted',false);
        value=cases{g};save(fullfile(caseDir,'failure.mat'),'value');write_json(fullfile(caseDir,'failure.json'),value);
        fprintf('JF_PROFILE n%d FAILED %s: %s\n',p.grids(g,1),exception.identifier,exception.message);
    end
end
report=struct('registration',p,'cases',{cases}, ...
    'acceptedOnBothGrids',all(cellfun(@(c)c.acceptedIndependentProfileCandidate,cases)), ...
    'productionTrajectoryModified',false,'productionGaugeModified',false,'oldSecantDecisionChanged',false, ...
    'interpretation','Predeclared projection-free profile Newton trial with true forcing, matching physical-time relaxation, unchanged tail/holdout gates and a natural-progress cost control. Any single failure rejects the candidate.');
save(fullfile(destination,'report.mat'),'report');write_json(fullfile(destination,'report.json'),report);
fprintf('JF_PROFILE_REPORT %s\n',destination);
end

function [candidate,report,fields]=construct(rho,z,ops,base,bd,p)
w=ops.integrationWeights(:);rhoNorm=sqrt(sum(w.*rho(:).^2));m=p.krylovDimension;
V=zeros(numel(rho),m);v=bd.rhs(:);rhsCount=1;eps0=p.differenceSteps(1);
for k=1:m
    for pass=1:2
        for j=1:k-1,v=v-V(:,j)*sum(w.*V(:,j).*v)/rhoNorm^2;end
    end
    nv=sqrt(sum(w.*v.^2));assert(nv>1e-10*rhoNorm,'ipm:JFKrylovRank','PDE tangent Krylov direction is numerically dependent.');
    V(:,k)=rhoNorm*v/nv;
    if k<m
        [fp,~]=ipm.evolve.flow(rho+eps0*reshape(V(:,k),size(rho)),ops);
        [fm,~]=ipm.evolve.flow(rho-eps0*reshape(V(:,k),size(rho)),ops);rhsCount=rhsCount+2;
        v=(fp(:)-fm(:))/(2*eps0);
    end
end
obs=struct('xi',bd.detail.quadratureX,'eta',bd.detail.quadratureY,'weights',bd.detail.quadratureWeights, ...
    'collectPhysical',strcmp(p.constraintMode,'tail_ellipsoid'));
[XI,ETA]=meshgrid(obs.xi,obs.eta);obs.core=abs(XI)<1 & ETA<1.5;
obs.sqrtWeights=sqrt(obs.weights(obs.core)/sum(obs.weights(obs.core)));
r0=obs.sqrtWeights.*bd.detail.quadratureFields.G(obs.core)/base.P;
J=cell(1,2);signatures=true(2,m,2);
T=cell(1,m);peakDerivative=zeros(1,m);
for level=1:2
    h=p.differenceSteps(level);J{level}=zeros(numel(r0),m);
    for k=1:m
        [plus,pd]=observe(rho+h*reshape(V(:,k),size(rho)),z,ops,obs);
        [minus,md]=observe(rho-h*reshape(V(:,k),size(rho)),z,ops,obs);rhsCount=rhsCount+2;
        J{level}(:,k)=(pd.objective-md.objective)/(2*h);
        signatures(level,k,1)=isequaln(plus.signature,base.signature);
        signatures(level,k,2)=isequaln(minus.signature,base.signature);
        if level==2 && obs.collectPhysical
            b=bd.physical;mask=b.tail;tw=b.weights(mask);tw=tw/sum(tw);
            dr=(pd.physical.rho(mask)-md.physical.rho(mask))/(2*h);
            du=(pd.physical.u1(mask)-md.physical.u1(mask))/(2*h);
            dv=(pd.physical.u2(mask)-md.physical.u2(mask))/(2*h);
            densityScale=sqrt(sum(tw.*b.rho(mask).^2));
            velocityScale=sqrt(sum(tw.*(b.u1(mask).^2+b.u2(mask).^2)));
            T{k}=[sqrt(tw).*dr/(densityScale*p.tailDensityTolerance*p.linearTailBudgetFraction); ...
                sqrt(tw).*du/(velocityScale*p.tailVelocityTolerance*p.linearTailBudgetFraction); ...
                sqrt(tw).*dv/(velocityScale*p.tailVelocityTolerance*p.linearTailBudgetFraction)];
            peakDerivative(k)=(plus.physicalPeak-minus.physicalPeak)/(2*h*base.physicalPeak*p.maximumPhysicalPeakChange*p.linearTailBudgetFraction);
        end
    end
end
discrepancy=norm(J{1}-J{2},'fro')/max(norm(J{2},'fro'),realmin);
Q=eye(m)/p.trustRadius^2;
if obs.collectPhysical,tailMatrix=cat(2,T{:});Q=Q+tailMatrix'*tailMatrix+peakDerivative'*peakDerivative;end
L=chol(Q);transformedJ=J{2}/L;
[U,S,W]=svd(transformedJ,'econ');singular=diag(S);keep=singular>p.svdRelativeFloor*singular(1);
assert(any(keep));s=singular(keep);a=U(:,keep)'*r0;B=W(:,keep);
coefficient=@(mu)-B*((s.*a)./(s.^2+mu));mu=0;coeff=coefficient(mu);
if norm(coeff)>1
    hi=max(s.^2);while norm(coefficient(hi))>1,hi=2*hi;end
    lo=0;for k=1:80,mid=(lo+hi)/2;if norm(coefficient(mid))>1,lo=mid;else,hi=mid;end,end
    mu=hi;coeff=coefficient(mu);
end
constraintNorm=norm(coeff);coeff=L\coeff;
delta=reshape(V*coeff,size(rho));candidate=rho+delta;
[trial,td]=observe(candidate,z,ops,[]);rhsCount=rhsCount+1;
match=matching(td,bd,p);gates=accept(trial,base,match,p);
predicted=norm(r0)^2-norm(r0+J{2}*coeff)^2;
actual=norm(r0)^2-trial.metrics.intrinsicCore^2;
agreement=actual/max(predicted,realmin);
linearization=struct('jacobianRelativeDifference',discrepancy,'geometrySignaturesStable',signatures, ...
    'singularValues',singular,'retainedRank',nnz(keep),'regularization',mu, ...
    'trustRadius',p.trustRadius,'coefficientNorm',norm(coeff),'predictedSquaredDecrease',predicted, ...
    'constraintMode',p.constraintMode,'constraintQuadraticForm',Q,'constraintNorm',constraintNorm, ...
    'actualSquaredDecrease',actual,'predictionAgreement',agreement, ...
    'passed',discrepancy<=p.maximumJacobianRelativeDifference && all(signatures,'all') && ...
        predicted>0 && agreement>=p.minimumPredictionAgreement);
estimatedTimeShift=sum(w.*delta(:).*bd.rhs(:))/sum(w.*bd.rhs(:).^2);
% The independent scaling control measures numerical homogeneity; it is not
% a candidate, projection or adjustable normalization of the actual trial.
amplitudeControls=cell(1,2);
for k=1:2
    factor=[.999,1.001];control=observe(factor(k)*rho,z,ops,[]);rhsCount=rhsCount+1;
    amplitudeControls{k}=struct('factor',factor(k),'rawCoreRatio',control.metrics.G.coreL2/base.metrics.G.coreL2, ...
        'intrinsicCoreRatio',control.metrics.intrinsicCore/base.metrics.intrinsicCore);
end
report=struct('candidate',trial,'matching',match,'gates',gates,'linearization',linearization, ...
    'fullFieldRelativeL2Jump',sqrt(sum(w.*delta(:).^2))/rhoNorm, ...
    'estimatedNativeTangentTimeShift',estimatedTimeShift, ...
    'estimatedNativeTangentPhysicalTimeShift',estimatedTimeShift*exp(z(2)-z(1)), ...
    'candidateClockShiftApplied',false,'amplitudeControls',{amplitudeControls}, ...
    'constructionRhsEvaluations',rhsCount,'peakProjectionCount',0,'localPatchCount',0, ...
    'sourceFieldUnchanged',true,'independentProfileCandidate',true,'productionAccepted',false);
fields=struct('basis',V,'coefficients',coeff,'delta',delta,'baselineObjective',r0, ...
    'jacobians',{J},'fixedObservation',obs,'candidateRates',td.rates);
end

function [o,d]=observe(rho,z,ops,objective)
[rhs,flow]=ipm.evolve.flow(rho,ops);F=rhs*ops.Dx';Omega=flow.source;
r=ipm_accellab_continuous_inner_rates(Omega,F,ops.x,ops.y,ops.Dx,ops.Dy,[0,max(ops.x)]);
assert(r.valid,'ipm:JFCoordinates','True continuous coordinates are invalid.');
P=r.peak.value;clock=exp(z(1)-z(2));
d=struct('rhs',rhs,'flow',flow,'rho',rho,'z',z,'omega',Omega,'forcing',F,'rates',r);
[~,nativePeakNode]=max(Omega(1,:));
signature=struct('continuous',r.signature,'nativeQuadraticPeakNode',nativePeakNode);
o=struct('P',P,'physicalPeak',P*clock,'PPrime',r.peakPrime,'physicalTime',z(3), ...
    'canonicalTime',z(5),'signature',signature,'poissonResidual',flow.poissonResidual);
if ~isempty(objective)
    c=ipm_accellab_pullback_hermite(Omega,F,ops.x,ops.y,ops.Dx,ops.Dy,r,objective.xi,objective.eta);
    d.objective=objective.sqrtWeights.*c.G(objective.core)/P;
    if objective.collectPhysical,d.physical=physical_observation(rho,flow,z,ops);end
    return
end
[rs,detail]=ipm_accellab_continuous_residual(Omega,F,ops.x,ops.y,ops.Dx,ops.Dy,r);
metrics=struct('G',rs.gauss.G,'wall',rs.wallGauss.G,'intrinsicCore',rs.gauss.G.coreL2/P, ...
    'intrinsicHoldout',rs.gauss.G.holdoutL2/P,'intrinsicWall',rs.wallGauss.G.fullLocalL2/P);
for group={'G','wall'}
    for key=fieldnames(metrics.(group{1}))'
        metrics.(group{1}).(key{1})=clock*metrics.(group{1}).(key{1});
    end
end
OmegaX=Omega*ops.Dx';OmegaY=ops.Dy*Omega;
native=(F+(r.translationRate+r.logScaleXRate*(ops.X-r.peak.x)).*OmegaX+ ...
    r.logScaleYRate*ops.Y.*OmegaY-r.peakPrime/P*Omega)/P*clock;
nx=(ops.X-r.peak.x)/r.wallCoreWidth;ny=ops.Y/r.verticalCoreWidth;
core=abs(nx)<=1 & ny<=1.5;whole=abs(nx)<=2 & ny<=3;
metrics.nativeCore=weighted(native,ops.integrationWeights,core);
metrics.nativeHoldout=weighted(native,ops.integrationWeights,whole&~core);
o.metrics=metrics;o.rateSummary=struct('aPrime',r.translationRate,'beta',r.logScaleXRate,'gamma',r.logScaleYRate, ...
    'PPrime',r.peakPrime,'widthX',r.wallCoreWidth,'widthY',r.verticalCoreWidth);
d.detail=detail;d.physical=physical_observation(rho,flow,z,ops);
d.nativeResidual=native;d.nativeCoreMask=core;d.nativeHoldoutMask=whole&~core;d.nativeWeights=ops.integrationWeights;
end

function d=physical_observation(rho,flow,z,ops)
x=linspace(-3,3,49);y=linspace(0,3,25);[X,Y]=meshgrid(x,y);
Cx=exp(z(1));Cw=exp(z(2));qx=Cx*X+z(4);qy=Cx*Y;
assert(min(qx,[],'all')>=min(ops.x) && max(qx,[],'all')<=max(ops.x) && max(qy,[],'all')<=max(ops.y));
r=ipm_accellab_tensor_hermite(rho,ops.x,ops.y,ops.Dx,ops.Dy,qx,qy);
o=ipm_accellab_tensor_hermite(flow.source,ops.x,ops.y,ops.Dx,ops.Dy,qx,qy);
u=ipm_accellab_tensor_hermite(flow.u1,ops.x,ops.y,ops.Dx,ops.Dy,qx,qy);
v=ipm_accellab_tensor_hermite(flow.u2,ops.x,ops.y,ops.Dx,ops.Dy,qx,qy);
d=struct('x',x,'y',y,'rho',r.value/Cw,'omega',Cx/Cw*o.value,'u1',u.value/Cw,'u2',v.value/Cw, ...
    'tail',abs(X)>=2.5|Y>=2,'weights',trap_weights(y)'*trap_weights(x));
end

function m=matching(trial,base,p)
a=trial.physical;b=base.physical;mask=b.tail;w=b.weights(mask);w=w/sum(w);
assert(isequal(a.x,b.x)&&isequal(a.y,b.y));
m=struct('tailDensityRelativeL2',sqrt(sum(w.*(a.rho(mask)-b.rho(mask)).^2))/max(sqrt(sum(w.*b.rho(mask).^2)),realmin), ...
    'tailSourceAbsoluteInfOverPeak',max(abs(a.omega(mask)-b.omega(mask)))/(base.rates.peak.value*exp(base.z(1)-base.z(2))), ...
    'tailVelocityRelativeL2',sqrt(sum(w.*((a.u1(mask)-b.u1(mask)).^2+(a.u2(mask)-b.u2(mask)).^2)))/ ...
        max(sqrt(sum(w.*(b.u1(mask).^2+b.u2(mask).^2))),realmin), ...
    'observationWindowDensityRelativeInf',relative_inf(a.rho,b.rho), ...
    'physicalPeakRelativeChange',abs(trial.rates.peak.value*exp(trial.z(1)-trial.z(2))/(base.rates.peak.value*exp(base.z(1)-base.z(2)))-1), ...
    'physicalTimeDifference',abs(trial.z(3)-base.z(3)));
m.passed=m.tailDensityRelativeL2<=p.tailDensityTolerance && ...
    m.tailSourceAbsoluteInfOverPeak<=p.tailSourceTolerance && m.tailVelocityRelativeL2<=p.tailVelocityTolerance && ...
    m.observationWindowDensityRelativeInf<=p.maximumObservationWindowDensityChange && ...
    m.physicalPeakRelativeChange<=p.maximumPhysicalPeakChange && m.physicalTimeDifference<1e-12;
for region={'nativeCoreMask','nativeHoldoutMask'}
    name=region{1};mask=base.(name);aa=weighted(trial.nativeResidual,base.nativeWeights,mask);bb=weighted(base.nativeResidual,base.nativeWeights,mask);
    m.fixedNativeRatios.(name)=struct('L2',aa.L2/bb.L2,'infinity',aa.infinity/bb.infinity);
    m.passed=m.passed && aa.L2<=bb.L2*(1+1e-10) && aa.infinity<=bb.infinity*(1+1e-10);
end
end

function g=accept(a,b,matching,p)
failure={};ratio=struct();
for key={'coreL2','holdoutL2','coreSampledInfinity','holdoutSampledInfinity'}
    name=key{1};ratio.(name)=a.metrics.G.(name)/max(b.metrics.G.(name),realmin);
    threshold=1+1e-10;if strcmp(name,'coreL2'),threshold=1-p.minimumDecrease;end
    if ratio.(name)>threshold,failure{end+1}=name;end %#ok<AGROW>
end
ratio.intrinsicCore=a.metrics.intrinsicCore/b.metrics.intrinsicCore;
ratio.intrinsicHoldout=a.metrics.intrinsicHoldout/b.metrics.intrinsicHoldout;
if ratio.intrinsicCore>1-p.minimumDecrease,failure{end+1}='intrinsicCore';end
if ratio.intrinsicHoldout>1+1e-10,failure{end+1}='intrinsicHoldout';end
for key={'coreL2','holdoutL2','fullLocalSampledInfinity'}
    name=key{1};ratio.(['wall_',name])=a.metrics.wall.(name)/max(b.metrics.wall.(name),realmin);
    if ratio.(['wall_',name])>1+1e-10,failure{end+1}=['wall_',name];end %#ok<AGROW>
end
for group={'nativeCore','nativeHoldout'}
    for key={'L2','infinity'}
        name=[group{1},'_',key{1}];ratio.(name)=a.metrics.(group{1}).(key{1})/max(b.metrics.(group{1}).(key{1}),realmin);
        if ratio.(name)>1+1e-10,failure{end+1}=name;end %#ok<AGROW>
    end
end
if ~matching.passed,failure{end+1}='physical_tail_or_peak_matching';end
g=struct('passed',isempty(failure),'failures',{failure},'ratios',ratio);
end

function [rho,z,audit]=physical_endpoint(rho,z,ops,target,dt)
cache=[];steps=0;initialZ=z;
while z(3)<target
    leftRho=rho;leftZ=z;
    if isempty(cache),[f,flow]=ipm.evolve.flow(rho,ops);cache=ipm.evolve.makeRhsCache(rho,f,flow,ops,unpack(z));end
    leftCache=cache;
    [rho,z,cache]=native_step(rho,z,dt,ops,cache);steps=steps+1;assert(steps<1000);
end
h=z(3)-leftZ(3);s=(target-leftZ(3))/h;assert(s>=0&&s<=1);
leftRate=leftCache.rhoRate/exp(leftZ(2)-leftZ(1));rightRate=cache.rhoRate/exp(z(2)-z(1));
leftScaleRate=leftCache.scaleRate/exp(leftZ(2)-leftZ(1));
rightScaleRate=cache.scaleRate/exp(z(2)-z(1));
rho=hermite_time(leftRho,rho,leftRate,rightRate,h,s);
z=hermite_time(leftZ,z,leftScaleRate,rightScaleRate,h,s);
assert(abs(z(3)-target)<2e-13);
audit=struct('target',target,'actualInterpolatedPhysicalTime',z(3),'dt',dt,'steps',steps, ...
    'rhsEvaluations',1+5*steps,'leftPhysicalTime',leftZ(3),'rightPhysicalTime',leftZ(3)+h, ...
    'fraction',s,'initialZ',initialZ,'sameOriginalPde',true,'isNativeCheckpoint',false, ...
    'timeObservation','C1 temporal interpolation of original accepted RK bracketing states, including all five scale variables and physical-time derivatives; no field or P projection.');
end

function [rho,z,cache]=native_step(rho,z,dt,ops,cache)
if isempty(cache),cache=[];end
[rho,flow,scale,cache]=ipm.evolve.stepSsprk54(rho,dt,ops,unpack(z),cache);
z=[scale.logC_l;scale.logC_omega;scale.physicalTime;scale.X_shift;scale.canonicalTime];
rate=max(abs(flow.transportU1)./ops.hx+abs(flow.transportU2)./ops.hy,[],'all');
assert(dt*rate<.4 && all(isfinite(rho),'all') && all(isfinite(z)));
end
function z=hermite_time(a,b,da,db,h,s)
z=(2*s^3-3*s^2+1)*a+(s^3-2*s^2+s)*h*da+(-2*s^3+3*s^2)*b+(s^3-s^2)*h*db;
end
function s=unpack(z)
s=struct('logC_l',z(1),'logC_omega',z(2),'physicalTime',z(3),'X_shift',z(4),'canonicalTime',z(5));
end
function m=weighted(f,w,mask)
assert(any(mask,'all'));m=struct('L2',sqrt(sum(w(mask).*f(mask).^2)/sum(w(mask))), ...
    'infinity',max(abs(f(mask))),'nodeCount',nnz(mask));
end
function v=relative_inf(a,b)
v=max(abs(a-b),[],'all')/max(max(abs(b),[],'all'),realmin);
end
function w=trap_weights(x)
x=x(:)';h=diff(x);w=[h(1)/2,(h(1:end-1)+h(2:end))/2,h(end)/2];
end
function opts=flat(config)
schema=ipm.config.schema();opts=struct();
for k=1:numel(schema.domainNames)
    q=config.(schema.domainNames{k});for name=fieldnames(q)',opts.(name{1})=q.(name{1});end
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
