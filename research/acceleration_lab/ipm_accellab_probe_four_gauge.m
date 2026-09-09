function report = ipm_accellab_probe_four_gauge(frameFile,outputRoot)
%IPM_ACCELLAB_PROBE_FOUR_GAUGE Actual WENO rate constraints, no PDE step.
% One fresh half-plane physical flow per tiny grid; all rate probes reuse it.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['four_gauge_rhs_',token]);mkdir(destination);
p=struct('grids',[49,25;65,33],'constraintTolerance',5e-11,'minimumRcond',1e-8, ...
    'rateDifferenceStep',1e-5,'maximumNewtonIterations',12,'maximumBacktracks',10, ...
    'directionalSteps',[1e-3,5e-4,2.5e-4],'frameFile',frameFile,'outputDirectory',destination, ...
    'geometry','Independent fresh half_plane; anisotropic state and shift are allowed; no double-odd native lineage is inherited.', ...
    'constraints','Actual [PPrime/P,aPrime/wallWidth,beta,gamma] of the combined WENO RHS.', ...
    'equation','Actual maintained assembleRhs with cxX,cyY,cr transport velocities; four linear generators supply an initial guess only.', ...
    'variants',{{'four_constraints_anisotropic','three_constraints_isotropic_vertical_rate_measured'}}, ...
    'pdeSteps',0,'rateProbesBuildPoisson',false,'peakProjectionCount',0);
write_json(fullfile(destination,'registration.json'),p);save(fullfile(destination,'registration.mat'),'p');
source=load(frameFile,'config');opts=flat(source.config);clear source
opts.symmetryMode='half_plane';opts.dynamicScaleGeometry='anisotropic';opts.rescalingMode='dynamic';
% This isotropic-only configuration label is unused by anisotropic fixed
% rates, but must itself be legal for the independently initialized domain.
opts.lengthGauge='local_strain';
opts.anisotropicGaugeMode='fixed';opts.anisotropicFixedCX=0;opts.anisotropicFixedCY=0;
opts.anisotropicFixedCOmega=0;opts.anisotropicFixedCR=0;opts.anisotropicPoissonSolver='direct';
opts.saveResults=false;opts.storeSnapshots=false;opts.makePlots=false;opts.livePlot=false;
opts.writeVideo=false;opts.verbose=false;opts.adaptiveRemesh=false;opts.initialAnalyticRemesh=false;
cases=cell(1,2);
for g=1:2
    caseDir=fullfile(destination,sprintf('n%d',p.grids(g,1)));mkdir(caseDir);
    records={};q=[];residual=[];
    try
        opts.nx=p.grids(g,1);opts.ny=p.grids(g,2);state=ipm.evolve.initialize(opts);
        assert(strcmp(state.ops.symmetryMode,'half_plane') && strcmp(state.ops.dynamicScaleGeometry,'anisotropic') && ...
            strcmp(state.ops.transportScheme,'weno5_fd') && strcmp(state.ops.transportBoundaryMode,'open'));
        rho=state.rho;ops=state.ops;flow=state.flow;Omega=flow.source;
        assert(all([flow.c_x,flow.c_y,flow.c_omega,flow.c_r]==0) && ...
            state.scale.logC_l==0 && state.scale.logC_y==0 && state.scale.logC_omega==0 && state.scale.X_shift==0);
        F0=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,0,0,0,0,ops,ops.transportBoundaryMode);
        assert(isequal(F0,state.rhsCache.rhoRate));
        r0=ipm_accellab_continuous_inner_rates(Omega,F0*ops.Dx',ops.x,ops.y,ops.Dx,ops.Dy,[0,max(ops.x)]);
        assert(r0.valid);length=r0.wallCoreWidth;
        coordinate=@(F) functional(Omega,F,ops,length);
        b=coordinate(F0);rhoX=rho*ops.Dx';rhoY=ops.Dy*rho;
        generators=cat(3,-ops.X.*rhoX,-ops.Y.*rhoY,rho,-length*rhoX);
        A=zeros(4);
        for k=1:4,A(:,k)=coordinate(generators(:,:,k));end
        assert(rcond(A)>=p.minimumRcond,'ipm:FourGaugeLinearCondition','Generator gauge Jacobian is ill-conditioned.');
        initial=-A\b;generatorRhs=F0;
        for k=1:4,generatorRhs=generatorRhs+initial(k)*generators(:,:,k);end
        linearResidual=coordinate(generatorRhs);
        evaluate=@(q) actual(q,rho,flow,ops,coordinate,length);
        try
            [three,F3]=solve_three(rho,F0,flow,ops,coordinate,length,A,b,p);
        catch exception3
            three=struct('status','failed','passed',false,'identifier',exception3.identifier,'message',exception3.message);
            F3=[];
        end
        write_json(fullfile(caseDir,'three_isotropic_report.json'),three);
        save(fullfile(caseDir,'three_isotropic.mat'),'three','F3','-v7.3');
        [residual,initialActualRhs,~]=evaluate(initial);calls=1;
        initialActualResidual=residual;records=cell(1,0);q=initial;converged=false;
        for iteration=1:p.maximumNewtonIterations
            if norm(residual,inf)<=p.constraintTolerance,converged=true;break;end
            [J,jcalls]=jacobian(evaluate,q,p.rateDifferenceStep);calls=calls+jcalls;
            condition=rcond(J);assert(condition>=p.minimumRcond,'ipm:FourGaugeActualCondition','Actual WENO gauge Jacobian is ill-conditioned.');
            direction=-J\residual;attempts=cell(1,0);accepted=false;
            for backtrack=0:p.maximumBacktracks
                damping=2^-backtrack;candidate=q+damping*direction;
                [next,~,details]=evaluate(candidate);calls=calls+1;
                attempts{end+1}=struct('damping',damping,'rates',candidate,'residual',next,'norm',norm(next,inf), ...
                    'PPrime',details.rates.peakPrime); %#ok<AGROW>
                if norm(next,2)<norm(residual,2)*(1-1e-4*damping) || norm(next,inf)<=p.constraintTolerance
                    accepted=true;break
                end
            end
            records{end+1}=struct('iteration',iteration,'oldRates',q,'oldResidual',residual,'jacobian',J, ...
                'rcond',condition,'direction',direction,'attempts',{attempts},'accepted',accepted); %#ok<AGROW>
            assert(accepted,'ipm:FourGaugeLineSearch','All registered actual-WENO Newton backtracks failed.');
            q=candidate;residual=next;
        end
        [residual,F,details]=evaluate(q);calls=calls+1;
        converged=converged || norm(residual,inf)<=p.constraintTolerance;
        [J1,n1]=jacobian(evaluate,q,p.rateDifferenceStep);
        [J2,n2]=jacobian(evaluate,q,p.rateDifferenceStep/2);calls=calls+n1+n2;
        jacobianRelativeDifference=norm(J1-J2,'fro')/max(norm(J2,'fro'),realmin);
        derivativeErrors=zeros(size(p.directionalSteps));signatures=true(size(p.directionalSteps));
        for k=1:numel(p.directionalSteps)
            h=p.directionalSteps(k);
            plus=ipm_accellab_continuous_inner_rates((rho+h*F)*ops.Dx',zeros(size(rho)),ops.x,ops.y,ops.Dx,ops.Dy,[0,max(ops.x)]);
            minus=ipm_accellab_continuous_inner_rates((rho-h*F)*ops.Dx',zeros(size(rho)),ops.x,ops.y,ops.Dx,ops.Dy,[0,max(ops.x)]);
            assert(plus.valid && minus.valid);
            signatures(k)=isequaln(plus.signature,r0.signature) && isequaln(minus.signature,r0.signature);
            derivative=(values(plus,length)-values(minus,length))/(2*h);
            derivativeErrors(k)=norm(derivative-residual,inf);
        end
        transportRate=max(abs(details.transportU1)./ops.hx+abs(details.transportU2)./ops.hy,[],'all');
        cases{g}=struct('status','completed','gridSize',p.grids(g,:),'operatorSourceConfig',state.config, ...
            'sourceInitialScale',state.scale,'sourcePoissonResidual',flow.poissonResidual, ...
            'physicalBaseRhsBitwise',true,'geometry',r0,'baselineConstraints',b, ...
            'linearGeneratorMatrix',A,'linearGeneratorRcond',rcond(A),'linearSeed',initial, ...
            'linearGeneratorConstraintResidual',linearResidual,'seedActualWenoResidual',initialActualResidual, ...
            'generatorVsActualSeedRhsRelativeInf',max(abs(generatorRhs-initialActualRhs),[],'all')/max(abs(initialActualRhs),[],'all'), ...
            'newtonRecords',{records},'converged',converged,'normalizedTranslationUnit',length, ...
            'actualRates',struct('cx',q(1),'cy',q(2),'comega',q(3),'cr',q(4)*length), ...
            'actualConstraintResidual',residual,'actualContinuousRates',details.rates, ...
            'physicalRhsCoordinateConsistency',covariance_defect(F,F0,rho,ops,[q(1:3);q(4)*length]), ...
            'finalJacobian',J2,'finalJacobianRcond',rcond(J2),'jacobianRelativeDifference',jacobianRelativeDifference, ...
            'directionalSteps',p.directionalSteps,'directionalDerivativeErrors',derivativeErrors, ...
            'directionalGeometrySignaturesStable',signatures,'rateAssemblyCalls',calls, ...
            'initialPhysicalPoissonEvaluations',1,'additionalRateProbePoissonEvaluations',0, ...
            'threeConstraintIsotropic',three, ...
            'transportRate',transportRate,'nextStepCflBound',state.config.time.cfl/transportRate, ...
            'pdeSteps',0,'peakProjectionCount',0,'isNativeGauge',false, ...
            'passed',converged && rcond(J2)>=p.minimumRcond && jacobianRelativeDifference<1e-5 && ...
                all(signatures) && derivativeErrors(end)<1e-6, ...
            'interpretation','Only instantaneous constraints of a separately initialized half-plane case are tested. No time integration, finite-step conservation, physical-time agreement or production promotion is claimed.');
        value=cases{g};save(fullfile(caseDir,'report.mat'),'value','-v7.3');
        x=ops.x;y=ops.y;Dx=ops.Dx;Dy=ops.Dy;
        save(fullfile(caseDir,'fields.mat'),'rho','Omega','F0','F','generatorRhs','initialActualRhs','x','y','Dx','Dy','-v7.3');
        write_json(fullfile(caseDir,'report.json'),value);
        fprintf('FOUR_GAUGE n%d passed=%d seed=%.3e final=%.3e rateProbes=%d cx=%.8g cy=%.8g cw=%.8g cr=%.8g\n', ...
            p.grids(g,1),value.passed,norm(initialActualResidual,inf),norm(residual,inf),calls,q(1:3),q(4)*length);
        clear state ops flow
    catch exception
        failure=struct('status','failed','identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
            'passed',false,'pdeSteps',0,'newtonRecords',{records},'lastRates',q,'lastResidual',residual);
        save(fullfile(caseDir,'failure.mat'),'failure');write_json(fullfile(caseDir,'failure.json'),failure);cases{g}=failure;
        fprintf('FOUR_GAUGE n%d FAILED %s: %s\n',p.grids(g,1),exception.identifier,exception.message);
    end
end
fourPassed=all(cellfun(@(c)c.passed,cases));
threePassed=all(cellfun(@(c)isfield(c,'threeConstraintIsotropic') && c.threeConstraintIsotropic.passed,cases));
report=struct('registration',p,'cases',{cases},'passed',fourPassed&&threePassed, ...
    'fourConstraintPassed',fourPassed,'threeConstraintPassed',threePassed, ...
    'productionModified',false,'pdeSteps',0,'physicalTimeValidationPerformed',false);
save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
fprintf('FOUR_GAUGE_REPORT %s\n',destination);
end

function [r,F,d]=actual(q,rho,flow,ops,coordinate,length)
[F,~,u1,u2]=ipm.evolve.assembleRhs(rho,flow.u1,flow.u2,q(1),q(2),q(3),q(4)*length,ops,ops.transportBoundaryMode);
[r,rates]=coordinate(F);d=struct('rates',rates,'transportU1',u1,'transportU2',u2);
end
function [r,rates]=functional(Omega,F,ops,length)
rates=ipm_accellab_continuous_inner_rates(Omega,F*ops.Dx',ops.x,ops.y,ops.Dx,ops.Dy,[0,max(ops.x)]);
assert(rates.valid,'ipm:FourGaugeGeometry','The continuous geometric functional is not differentiable here.');
r=[rates.peakPrime/rates.peak.value;rates.translationRate/length;rates.logScaleXRate;rates.logScaleYRate];
end
function v=values(r,length)
v=[log(r.peak.value);r.peak.x/length;log(r.wallCoreWidth);log(r.verticalCoreWidth)];
end
function [J,calls]=jacobian(evaluate,q,relativeStep)
m=numel(q);J=zeros(m);calls=0;
for k=1:m
    h=relativeStep*max(1,abs(q(k)));unit=zeros(m,1);unit(k)=h;
    plus=evaluate(q+unit);minus=evaluate(q-unit);calls=calls+2;J(:,k)=(plus-minus)/(2*h);
end
end
function [report,F]=solve_three(rho,F0,flow,ops,coordinate,length,A,b,p)
B=[1,0,0;1,0,0;0,1,0;0,0,1];M=A(1:3,:)*B;
assert(rcond(M)>=p.minimumRcond);q=-M\b(1:3);seed=q;
evaluate=@(v) actual_three(v,rho,flow,ops,coordinate,length);
[residual,~,~]=evaluate(q);initialResidual=residual;calls=1;records={};
for iteration=1:p.maximumNewtonIterations
    if norm(residual,inf)<=p.constraintTolerance,break;end
    [J,n]=jacobian(evaluate,q,p.rateDifferenceStep);calls=calls+n;
    assert(rcond(J)>=p.minimumRcond,'ipm:ThreeGaugeCondition','Actual isotropic gauge Jacobian is ill-conditioned.');
    direction=-J\residual;attempts={};accepted=false;
    for backtrack=0:p.maximumBacktracks
        damping=2^-backtrack;candidate=q+damping*direction;[next,~,detail]=evaluate(candidate);calls=calls+1;
        attempts{end+1}=struct('damping',damping,'rates',candidate,'residual',next,'gamma',detail.rates.logScaleYRate); %#ok<AGROW>
        if norm(next,2)<norm(residual,2)*(1-1e-4*damping) || norm(next,inf)<=p.constraintTolerance
            accepted=true;break
        end
    end
    records{end+1}=struct('iteration',iteration,'oldRates',q,'oldResidual',residual,'jacobian',J, ...
        'rcond',rcond(J),'direction',direction,'attempts',{attempts},'accepted',accepted); %#ok<AGROW>
    if ~accepted
        report=struct('status','line_search_rejected','passed',false,'records',{records},'lastRates',q,'lastResidual',residual);F=[];return
    end
    q=candidate;residual=next;
end
[residual,F,details]=evaluate(q);calls=calls+1;
[J1,n1]=jacobian(evaluate,q,p.rateDifferenceStep);[J2,n2]=jacobian(evaluate,q,p.rateDifferenceStep/2);calls=calls+n1+n2;
discrepancy=norm(J1-J2,'fro')/max(norm(J2,'fro'),realmin);
report=struct('status','completed','linearGeneratorMatrix',M,'linearGeneratorRcond',rcond(M), ...
    'linearSeed',seed,'seedActualWenoResidual',initialResidual,'records',{records}, ...
    'actualRates',struct('cx',q(1),'cy',q(1),'comega',q(2),'cr',q(3)*length), ...
    'actualConstraintResidual',residual,'actualContinuousRates',details.rates, ...
    'physicalRhsCoordinateConsistency',covariance_defect(F,F0,rho,ops,[q(1);q(1);q(2);q(3)*length]), ...
    'verticalLogWidthRateMeasured',details.rates.logScaleYRate,'verticalConstraintApplied',false, ...
    'finalJacobian',J2,'finalJacobianRcond',rcond(J2),'jacobianRelativeDifference',discrepancy, ...
    'rateAssemblyCalls',calls,'additionalPoissonEvaluations',0,'pdeSteps',0, ...
    'futureIsotropicAspect',1,'fixedPoissonLuReuseSupported',true, ...
    'luEvidence','Maintained high_order poisson.m uses ops.poisson for kappa==1; unequal aspect uses a fresh decomposition per call.', ...
    'passed',norm(residual,inf)<=p.constraintTolerance && rcond(J2)>=p.minimumRcond && discrepancy<1e-5, ...
    'interpretation','Actual amplitude/position/wall-width constraints only; vertical motion remains measured. No integration or physical covariance certification.');
fprintf('THREE_GAUGE n%d passed=%d residual=%.3e gamma=%.8g rateProbes=%d\n',ops.nx,report.passed,norm(residual,inf),report.verticalLogWidthRateMeasured,calls);
end
function [r,F,d]=actual_three(q,rho,flow,ops,coordinate,length)
[r,F,d]=actual([q(1);q(1);q(2);q(3)],rho,flow,ops,coordinate,length);r=r(1:3);
end
function report=covariance_defect(F,F0,rho,ops,q)
% At the common initial scales/shift, remove the continuum coordinate
% generators using maintained derivatives. WENO assembly need not commute
% exactly with this discrete chain rule; this measures that defect.
rx=rho*ops.Dx';ry=ops.Dy*rho;
reconstructed=F+q(1)*ops.X.*rx+q(2)*ops.Y.*ry-q(3)*rho+q(4)*rx;
defect=reconstructed-F0;w=ops.integrationWeights;
report=struct('wholeBoxRelativeInf',max(abs(defect),[],'all')/max(abs(F0),[],'all'), ...
    'wholeBoxRelativeL2',sqrt(sum(w.*defect.^2,'all')/sum(w.*F0.^2,'all')), ...
    'wallRelativeInf',max(abs(defect(1,:)))/max(abs(F0(1,:))), ...
    'wholeBoxAbsoluteInf',max(abs(defect),[],'all'), ...
    'initialAspect',1,'wallNoPenetrationIdentity',0, ...
    'meaning','Discrete coordinate-chain-rule consistency defect, not an exact physical covariance or continuum error certificate. Nonzero values are retained although instantaneous gauge constraints pass.');
end
function opts=flat(config)
schema=ipm.config.schema();opts=struct();
for k=1:numel(schema.domainNames)
    q=config.(schema.domainNames{k});for name=fieldnames(q)',opts.(name{1})=q.(name{1});end
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
