function out=ipm_accellab_shadow_tap(action,varargin)
%IPM_ACCELLAB_SHADOW_TAP Observations only, in an isolated instrumented source.
% Neither hook has outputs used by the solver. Every oracle LU still runs.
persistent s
out=[];
switch action
    case 'reset'
        s=struct('directory',varargin{1},'window',0,'historyX',[], ...
            'historyB',[],'selector',NaN,'stopped',false,'attempts',0, ...
            'calibration',[],'flowCount',0,'poissonCount',0,'pending',[], ...
            'observations',{{}},'attemptReports',{{}},'windows',{{}}, ...
            'callbackSeconds',0,'poissonCaptureSeconds',0);
    case 'poisson'
        timer=tic;
        assert(~isempty(s));
        s.poissonCount=s.poissonCount+1;
        s.pending=struct('b',varargin{1},'psi',varargin{2},'A',varargin{3});
        s.poissonCaptureSeconds=s.poissonCaptureSeconds+toc(timer);
    case 'sample'
        timer=tic;
        rho=varargin{1};ops=varargin{2};scale=varargin{3};F=varargin{4};
        flow=varargin{5};nativeSeconds=varargin{6};
        assert(ops.nx<=65 && ops.ny<=33 && (ops.nx-2)*(ops.ny-2)<=1985);
        assert(~isempty(s.pending) && isequaln(s.pending.psi,flow.psi));
        s.flowCount=s.flowCount+1;
        tau=scale.canonicalTime; window=1+(tau>.032);
        selector=flow.omegaGaugeQuadraticStencilCenterX;
        change=~isequal(selector,s.selector);
        if window~=s.window
            if s.window>0
                s.windows{end+1}=struct('index',s.window,'attempts',s.attempts, ...
                    'stoppedOnFailure',s.stopped);
            end
            s.window=window;s.historyX=[];s.historyB=[];s.stopped=false;s.attempts=0;
        end
        if change, s.historyX=[];s.historyB=[]; end
        s.selector=selector;
        psi=flow.psi(2:end-1,2:end-1); b=s.pending.b;
        residual=norm(s.pending.A*psi(:)-b,inf)/max(norm(b,inf),eps);
        assert(isequaln(residual,flow.poissonResidual),'ipm:ShadowPoissonPair', ...
            'Captured original b and psi do not reproduce the original residual.');
        observation=struct('flowIndex',s.flowCount,'window',window, ...
            'canonicalTime',tau,'physicalTime',scale.physicalTime, ...
            'selector',selector,'selectorChanged',change, ...
            'nativeFlowSecondsIncludingCapture',nativeSeconds, ...
            'poissonResidualExact',true,'predicted',false,'windowStopped',s.stopped);
        if ~s.stopped && s.attempts<20 && size(s.historyX,2)==5
            s.attempts=s.attempts+1; observation.predicted=true;
            try
                attempt=compare(rho,ops,scale,F,flow,s.pending.A,b, ...
                    s.historyX,s.historyB,s.calibration);
                if isempty(s.calibration),s.calibration=attempt.calibration;end
                attempt.window=window;attempt.flowIndex=s.flowCount;
                attempt.canonicalTime=tau;attempt.physicalTime=scale.physicalTime;
                file=fullfile(s.directory,sprintf('attempt_%02d_flow_%04d.mat',window,s.flowCount));
                save(file,'attempt','rho','scale','F','flow','-v7.3');
                attempt.fields=[];s.attemptReports{end+1}=attempt;
                if ~attempt.passed,s.stopped=true;end
            catch exception
                attempt=struct('passed',false,'infrastructureOrObserverFailure',true, ...
                    'window',window,'flowIndex',s.flowCount,'canonicalTime',tau, ...
                    'identifier',exception.identifier,'message',exception.message);
                s.attemptReports{end+1}=attempt;s.stopped=true;
                save(fullfile(s.directory,sprintf('failed_%02d_flow_%04d.mat',window,s.flowCount)), ...
                    'attempt','rho','scale','F','flow','-v7.3');
            end
        end
        s.historyX=[s.historyX,psi(:)];s.historyB=[s.historyB,b];
        if size(s.historyX,2)>5,s.historyX=s.historyX(:,end-4:end);s.historyB=s.historyB(:,end-4:end);end
        s.observations{end+1}=observation;s.pending=[];
        s.callbackSeconds=s.callbackSeconds+toc(timer);
    case 'finish'
        s.windows{end+1}=struct('index',s.window,'attempts',s.attempts,'stoppedOnFailure',s.stopped);
        out=rmfield(s,{'historyX','historyB','pending'});s=[];
    otherwise
        error('ipm:ShadowAction','Unknown action.');
end
end

function a=compare(rho,ops,scale,F,flow,A,b,Xhistory,Bhistory,calibration)
timer=tic; nx=ops.nx;ny=ops.ny;
% Here the map is used only for the preliminary interior comparison; actual
% two-component maintained velocity derivatives below supply the hard gate.
shadow=ipm_accellab_recycle_shadow(A,[Bhistory,b], ...
    [Xhistory,reshape(flow.psi(2:end-1,2:end-1),[],1)],speye(size(A,1)));
predictionSeconds=toc(timer);e=shadow.entries(1);
psi=flow.psi;psi(2:end-1,2:end-1)=reshape(e.candidate,ny-2,nx-2);
u1=-(ops.Dy*psi);u2=psi*ops.Dx';u2(1,:)=0;
timer=tic;
[original,originalRates,originalDetails]=from_velocity(rho,flow.u1,flow.u2,flow,ops);
oracleGroupingExact=isequaln(original,F) && ...
    isequaln(originalRates,[flow.c_l,flow.c_omega,flow.c_r]);
assert(oracleGroupingExact,'ipm:ShadowRhsOracle','Native adapter parity failed.');
[trial,rates,trialDetails]=from_velocity(rho,u1,u2,flow,ops);
assemblyPairSeconds=toc(timer);
FX=F*ops.Dx';trialFX=trial*ops.Dx';Omega=flow.source;
window=[0,2*ops.rescaling.transportAnchorX];
coordinates=ipm_accellab_continuous_inner_rates(Omega,FX,ops.x,ops.y,ops.Dx,ops.Dy,window);
assert(coordinates.valid,'ipm:ShadowGeometry','Native geometry is unavailable.');
XI=(ops.X-coordinates.peak.x)/coordinates.wallCoreWidth;
ETA=ops.Y/coordinates.verticalCoreWidth;
Cx=exp(scale.logC_l);Cw=exp(scale.logC_omega);
px=(ops.X-scale.X_shift)/Cx;py=ops.Y/Cx;
mask=struct('wholeBox',true(size(rho)),'core',abs(XI)<=1 & ETA<=1.5, ...
    'holdout',abs(XI)<=2 & ETA<=3 & ~(abs(XI)<=1 & ETA<=1.5), ...
    'physicalTail',abs(px)<=3 & py<=3 & (abs(px)>=2.5 | py>=2));
names=fieldnames(mask);observations=struct();fieldPassed=true;
for k=1:numel(names)
    key=names{k};m=mask.(key);assert(any(m,'all'),'ipm:ShadowEmptyMask','Empty registered mask.');
    v1=pair_norms(u1/Cw,flow.u1/Cw,ops.integrationWeights,m);
    v2=pair_norms(u2/Cw,flow.u2/Cw,ops.integrationWeights,m);
    f=pair_norms(trial,F,ops.integrationWeights,m);
    fx=pair_norms(trialFX,FX,ops.integrationWeights,m);
    observations.(key)=struct('nodes',nnz(m),'u1',v1,'u2',v2,'F',f,'FX',fx);
    fieldPassed=fieldPassed && max([v1.relativeInf,v2.relativeInf])<=1e-10 && ...
        max([f.relativeInf,f.relativeL2,fx.relativeInf,fx.relativeL2])<=1e-9;
end
timer=tic;
other=ipm_accellab_continuous_inner_rates(Omega,trialFX,ops.x,ops.y,ops.Dx,ops.Dy,window);
assert(other.valid);
[r0,d0]=ipm_accellab_continuous_residual(Omega,FX,ops.x,ops.y,ops.Dx,ops.Dy,coordinates);
[r1,d1]=ipm_accellab_continuous_residual(Omega,trialFX,ops.x,ops.y,ops.Dx,ops.Dy,other);
continuous=shape_pair(r0,d0,r1,d1);
g0=native_modulated(Omega,FX,ops,coordinates);g1=native_modulated(Omega,trialFX,ops,other);
nativeModulated=struct();nativeModulatedPassed=true;
for k=1:numel(names)
    key=names{k};v=pair_norms(g1,g0,ops.integrationWeights,mask.(key));
    nativeModulated.(key)=v;
    nativeModulatedPassed=nativeModulatedPassed && max(v.relativeInf,v.relativeL2)<=1e-4;
end
smooth=cell(1,3);theta=[.5,.3,.7];
if isempty(calibration)
    calibration=repmat(struct('theta',0,'horizontalFactor',0,'verticalFactor',0),1,3);
    for k=1:3
        sr=ipm_accellab_smooth_inner_rates(Omega,FX,ops.x,ops.y,ops.Dx,ops.Dy,window,theta(k));
        assert(sr.valid);calibration(k)=struct('theta',theta(k), ...
            'horizontalFactor',coordinates.wallCoreWidth/sr.rawWallWidth, ...
            'verticalFactor',coordinates.verticalCoreWidth/sr.rawVerticalWidth);
    end
end
for k=1:3
    b0=ipm_accellab_smooth_inner_rates(Omega,FX,ops.x,ops.y,ops.Dx,ops.Dy,window,theta(k));
    b1=ipm_accellab_smooth_inner_rates(Omega,trialFX,ops.x,ops.y,ops.Dx,ops.Dy,window,theta(k));
    assert(b0.valid && b1.valid);
    c0=chart(b0,calibration(k));c1=chart(b1,calibration(k));
    [s0,v0]=ipm_accellab_continuous_residual(Omega,FX,ops.x,ops.y,ops.Dx,ops.Dy,c0);
    [s1,v1]=ipm_accellab_continuous_residual(Omega,trialFX,ops.x,ops.y,ops.Dx,ops.Dy,c1);
    smooth{k}=shape_pair(s0,v0,s1,v1);smooth{k}.theta=theta(k);
end
shapeSeconds=toc(timer);
rateError=abs(rates-originalRates)./max(1,abs(originalRates));
algebraPassed=e.rawCoefficientNorm<=32 && e.retainedRank>0 && ...
    e.relativeGlobalResidual<=5e-12 && e.componentwiseBackwardError<=5e-14;
shapePassed=continuous.passed && all(cellfun(@(v)v.passed,smooth)) && nativeModulatedPassed;
fields=struct('trialPsi',psi,'trialU1',u1,'trialU2',u2,'trialRhs',trial, ...
    'trialFX',trialFX,'originalFX',FX,'capturedB',b, ...
    'historyX',Xhistory,'historyB',Bhistory,'nativePsi',flow.psi, ...
    'originalContinuousRates',coordinates,'trialContinuousRates',other);
e.candidate=[];
a=struct('kind','tiny_IPM_oracle_recycle_shadow_v1','nodes',[nx,ny], ...
    'memory',4,'interiorProjection',e,'oracleGroupingExact',oracleGroupingExact, ...
    'nativeRates',originalRates,'trialRates',rates,'normalizedRateError',rateError, ...
    'nativeActualPPrime',originalDetails.omegaGaugeQuadraticRateResidual, ...
    'trialActualPPrime',trialDetails.omegaGaugeQuadraticRateResidual, ...
    'observations',observations,'continuous',continuous,'smooth',{smooth}, ...
    'nativeModulated',nativeModulated,'nativeModulatedPassed',nativeModulatedPassed, ...
    'calibration',calibration,'algebraPassed',algebraPassed,'fieldPassed',fieldPassed, ...
    'shapePassed',shapePassed,'ratesPassed',all(rateError<=1e-10), ...
    'passed',algebraPassed && fieldPassed && shapePassed && all(rateError<=1e-10), ...
    'predictionSeconds',predictionSeconds,'assemblyPairSeconds',assemblyPairSeconds, ...
    'shapeObservationSeconds',shapeSeconds,'fields',fields, ...
    'integrationApproved',false,'oracleLUSkipped',false);
end

function [F,rates,details]=from_velocity(rho,u1,u2,flow,ops)
p=struct('u1',u1,'u2',u2,'psi',flow.psi,'source',flow.source, ...
    'poissonResidual',flow.poissonResidual,'poissonSolveInfo',flow.poissonSolveInfo);
details=ipm.diagnostics.trackFeatures(rho,[],p,ops);
[cl,cw,cr,base,details]=ipm.evolve.isotropicGauge(rho,[],p,ops,details);
F=base+cw*rho;rates=[cl,cw,cr];
end

function o=pair_norms(trial,native,w,mask)
d=trial(mask)-native(mask);n=native(mask);w=w(mask);
scaleInf=max(abs(n));scaleL2=sqrt(sum(w.*n.^2)/sum(w));
assert(scaleInf>1e-14 && scaleL2>1e-14,'ipm:ShadowSmallDenominator','Relative observer denominator unavailable.');
o=struct('relativeInf',max(abs(d))/scaleInf, ...
    'relativeL2',sqrt(sum(w.*d.^2)/sum(w))/scaleL2, ...
    'absoluteInf',max(abs(d)),'nativeInf',scaleInf,'nativeL2',scaleL2);
end

function o=shape_pair(r0,d0,r1,d1)
assert(isequaln(d0.quadratureX,d1.quadratureX) && isequaln(d0.quadratureY,d1.quadratureY));
[X,Y]=meshgrid(d0.quadratureX,d0.quadratureY);core=abs(X)<1 & Y<1.5;
w=d0.quadratureWeights;
m0=pair_norms(d1.quadratureFields.G,d0.quadratureFields.G,w,core);
m1=pair_norms(d1.quadratureFields.G,d0.quadratureFields.G,w,~core);
o=struct('native',r0,'trial',r1,'coreDifference',m0,'holdoutDifference',m1, ...
    'passed',max([m0.relativeL2,m0.relativeInf,m1.relativeL2,m1.relativeInf])<=1e-4);
end

function c=chart(r,u)
c=struct('valid',true,'peak',r.peak,'peakPrime',r.peakPrime,'translationRate',r.translationRate, ...
    'wallCoreWidth',u.horizontalFactor*r.rawWallWidth,'verticalCoreWidth',u.verticalFactor*r.rawVerticalWidth, ...
    'wallCoreWidthPrime',u.horizontalFactor*r.rawWallWidthPrime,'verticalCoreWidthPrime',u.verticalFactor*r.rawVerticalWidthPrime);
end

function G=native_modulated(Omega,FX,ops,r)
P=r.peak.value;beta=r.wallCoreWidthPrime/r.wallCoreWidth;gamma=r.verticalCoreWidthPrime/r.verticalCoreWidth;
G=(FX+(r.translationRate+beta*(ops.X-r.peak.x)).*(Omega*ops.Dx')+ ...
    gamma*ops.Y.*(ops.Dy*Omega))/P-r.peakPrime/P*(Omega/P);
end
