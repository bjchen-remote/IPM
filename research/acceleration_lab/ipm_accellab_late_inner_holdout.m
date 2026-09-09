function report=ipm_accellab_late_inner_holdout(registrationFile,mode)
% Chronological holdout: training and future evaluation are separate calls.
reg=jsondecode(fileread(registrationFile));assert(maxNumCompThreads==10);
profile clear;profile on;timer=tic;
assert(any(strcmp(mode,{'train','evaluate'})));
out=fullfile(reg.outputDirectory,mode);assert(~isfolder(out));mkdir(out);
try
 if strcmp(mode,'train')
  files=reg.trainingFiles;frames=cell(1,numel(files));lineage=[];
  for k=1:numel(files)
   [frames{k},lineage]=read_frame(files{k},reg,lineage);fprintf('TRAIN frame=%d step=%d tau=%.12g\n',k,frames{k}.record.step,frames{k}.record.equivalentTau);
  end
  times=cellfun(@(f)f.record.equivalentTau,frames);dt=mean(diff(times));
  assert(numel(times)==8&&max(abs(diff(times)-dt))<1e-10&&abs(dt-.2)<1e-10);
  models=cell(2,2);methods={'hermite','linear'};
  for observer=1:2
   for direction=1:2
    z=linspace(reg.intervals(direction,1),reg.intervals(direction,2),reg.pointCounts(direction));
    values=zeros(8,numel(z));
    for k=1:8,values(k,:)=sample(frames{k}.traces{observer,direction},z,methods{observer});end
    weights=ones(size(z));weights([1,end])=.5;weights=weights/sum(weights);
    models{observer,direction}=fit_models(values,weights,reg.futureCount);
   end
  end
  geom=cell2mat(cellfun(@(f)geometry_values(f.record),frames,'UniformOutput',false)');
  scalarModels=cell(1,4);for j=1:4,scalarModels{j}=fit_models(geom(:,j),1,reg.futureCount);end
  frozen=struct('registrationFile',registrationFile,'trainingFrames',{frames},'trainingTimes',times,'dt',dt, ...
   'trainingLineage',lineage,'shapeModels',{models},'scalarModels',{scalarModels},'geometryTransforms', ...
   {{'log physical peak','physical peak position','log physical wall width','log physical vertical width'}}, ...
   'rank',1,'ridge',0,'futureRead',false,'futureFitOrParameterSelection',false);
  save(fullfile(out,'frozen_model.mat'),'frozen','-v7.3');
  report=struct('kind','late_inner_training_only_v1','registration',reg,'trainingTimes',times, ...
   'shapeDiagnostics',{{models{1,1}.diagnostics,models{1,2}.diagnostics,models{2,1}.diagnostics,models{2,2}.diagnostics}}, ...
   'geometryDiagnostics',{cellfun(@(m)m.diagnostics,scalarModels,'UniformOutput',false)}, ...
   'futureFieldsRead',false,'frozenModelFile',fullfile(out,'frozen_model.mat'));
 else
  d=load(fullfile(reg.outputDirectory,'train/frozen_model.mat'),'frozen');frozen=d.frozen;clear d;
  assert(~frozen.futureRead&&frozen.rank==1&&frozen.ridge==0&&numel(reg.futureFiles)==reg.futureCount);
  report=struct('kind','late_inner_chronological_holdout_v1','registration',reg,'futureCases',{{}}, ...
   'trainingTimes',frozen.trainingTimes,'shapeModelDiagnostics',{{}},'geometryModelDiagnostics',{{}}, ...
   'modelRefittedOnFuture',false,'nativeSignatureValidated',true,'allHistoryTrusted',true,'noLU',true,'noPDE',true, ...
   'limitEstimated',false,'acceptedTrajectoryWritten',false, ...
   'interpretation','Chronologically held-out existing native fields. Some later states had prior scientific plots; this is not a newly prospective experiment. Future arrays do not enter fit, rank, ridge, clocks or window selection.');
  methods={'hermite','linear'};modelNames={'persistence','linearTrend','rankOneIncrement'};
  futureFrames=cell(1,reg.futureCount);curveData=cell(reg.futureCount,2,2);
  for k=1:reg.futureCount
   [f,~]=read_frame(reg.futureFiles{k},reg,frozen.trainingLineage);futureFrames{k}=f;
   expected=frozen.trainingTimes(end)+k*frozen.dt;assert(abs(f.record.equivalentTau-expected)<1e-10);
   row=struct('record',f.record,'shape',{{}},'geometry',{{}},'observerDifferences',{{}});
   for direction=1:2
    interval=reg.intervals(direction,:);n=reg.pointCounts(direction);z=linspace(interval(1),interval(2),n);
    fine=linspace(interval(1),interval(2),2*n-1);[q,w]=gauss_nodes(frozen.trainingFrames,f,interval,direction);
    for observer=1:2
     name=methods{observer};m=frozen.shapeModels{observer,direction};
     actual=sample(f.traces{observer,direction},q,name);actualCoarse=sample(f.traces{observer,direction},z,name);
     actualFine=sample(f.traces{observer,direction},fine,name);
     trainQ=zeros(8,numel(q));trainZ=zeros(8,numel(z));trainFine=zeros(8,numel(fine));
     for j=1:8
      trace=frozen.trainingFrames{j}.traces{observer,direction};
      trainQ(j,:)=sample(trace,q,name);trainZ(j,:)=sample(trace,z,name);trainFine(j,:)=sample(trace,fine,name);
     end
     entry=struct('observer',name,'direction',direction,'modelErrors',{{}});
     cd=struct('coordinate',fine,'actual',actualFine,'predictions',{{}});
     for model=1:3
      coeff=m.coefficients{model}(k,:);pred=coeff*trainQ;predZ=coeff*trainZ;predFine=coeff*trainFine;
      metric=errors(pred,actual,w,predZ,actualCoarse,predFine,actualFine);metric.model=modelNames{model};
      entry.modelErrors{model}=metric;cd.predictions{model}=predFine;
     end
     row.shape{observer,direction}=entry;curveData{k,observer,direction}=cd;
    end
    h=sample(f.traces{1,direction},q,'hermite');l=sample(f.traces{2,direction},q,'linear');
    hf=sample(f.traces{1,direction},fine,'hermite');lf=sample(f.traces{2,direction},fine,'linear');
    row.observerDifferences{direction}=struct('direction',direction,'relativeL2',sqrt(sum(w.*(h-l).^2)/sum(w.*h.^2)), ...
     'sampledRelativeInf',max(abs(hf-lf))/max(abs(hf)));
   end
   actualGeometry=geometry_values(f.record);trainGeometry=cell2mat(cellfun(@(v)geometry_values(v.record),frozen.trainingFrames,'UniformOutput',false)');
   for channel=1:4
    entry=struct('channel',frozen.geometryTransforms{channel},'actualTransformed',actualGeometry(channel),'models',{{}});
    for model=1:3
     pred=frozen.scalarModels{channel}.coefficients{model}(k,:)*trainGeometry(:,channel);
     if channel==2
      a=struct('model',modelNames{model},'prediction',pred,'actual',actualGeometry(channel), ...
       'absoluteError',abs(pred-actualGeometry(channel)),'errorOverActualWallWidth',abs(pred-actualGeometry(channel))/f.record.physicalWidths(1));
     else
      a=struct('model',modelNames{model},'prediction',exp(pred),'actual',exp(actualGeometry(channel)), ...
       'absoluteError',abs(exp(pred)-exp(actualGeometry(channel))),'relativeError',abs(expm1(pred-actualGeometry(channel))));
     end
     entry.models{model}=a;
    end
    row.geometry{channel}=entry;
   end
   report.futureCases{k}=row;
   fprintf('HOLDOUT frame=%d step=%d tau=%.12g wallRank1=%.6g verticalRank1=%.6g\n', ...
    k,f.record.step,f.record.equivalentTau,row.shape{1,1}.modelErrors{3}.relativeL2,row.shape{1,2}.modelErrors{3}.relativeL2);
  end
  report.shapeModelDiagnostics={frozen.shapeModels{1,1}.diagnostics,frozen.shapeModels{1,2}.diagnostics};
  report.geometryModelDiagnostics=cellfun(@(m)m.diagnostics,frozen.scalarModels,'UniformOutput',false);
  report.decision=decision(report);save(fullfile(out,'observations.mat'),'report','futureFrames','curveData','-v7.3');
 end
 profile off;profileInfo=profile('info');files={profileInfo.FunctionTable.FileName};
 for suffix={'+mesh/build.m','+field/poisson.m','+field/velocity.m','+evolve/flow.m','+evolve/advance.m','+output/restoreCheckpoint.m','+remesh/transfer.m'}
  assert(~any(endsWith(files,suffix{1})),['Forbidden call: ',suffix{1}]);
 end
 report.callGraphNoLUVerified=true;report.wallSeconds=toc(timer);
 save(fullfile(out,'report.mat'),'report','profileInfo','-v7.3');write_json(fullfile(out,'report.json'),report);
catch e
 profile off;failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
 write_json(fullfile(out,'failure.json'),failure);rethrow(e)
end
end

function [f,lineage]=read_frame(file,reg,lineage)
cp=ipm.output.readCheckpoint(file);s=cp.payload.state;h=cp.payload.log.history;
assert(all(ipm.output.trustedMask(h,s.config))&&strcmp(s.runMetadata.caseId,reg.caseId)&& ...
 strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic')&&s.scale.X_shift==0);
l=s.runMetadata.caseMetadata.latePhysicalBoxBranch;
if isempty(lineage)
 parent=ipm.output.readCheckpoint(l.parentCheckpoint);p=parent.payload.state;
 assert(~isfield(p.runMetadata.caseMetadata,'latePhysicalBoxBranch')&&strcmp(p.runMetadata.caseId,l.parentCaseId)&& ...
  p.step==l.parentStep&&p.scale.canonicalTime==l.parentCanonicalTime&&p.scale.physicalTime==l.absolutePhysicalEpoch&& ...
  exp(p.scale.logC_l)==l.parentCx&&exp(p.scale.logC_omega)==l.parentComega&& ...
  l.canonicalCovarianceFactor==l.parentCx/l.parentComega&&~l.parentHistoryInherited&&~l.parentCheckpointModified);
 lineage=l;
else,assert(isequaln(l,lineage));end
Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);W=s.rho*Dx';
g=continuous_inner_geometry(W,s.x,s.y,Dx,Dy,s.config.scaling.transportAnchorX*reg.peakWindowAnchorFactors(:)');
x=(s.x-g.peakX)/g.wallWidth;y=s.y'/g.verticalWidth;P=g.peak;
wallJets=W(1,:)*Dx';v=ipm_accellab_tensor_hermite(W,s.x,s.y,Dx,Dy,g.peakX+zeros(size(s.y)),s.y);
linearVertical=interp2(s.x,s.y,W,g.peakX+zeros(size(s.y)),s.y,'linear');
traces=cell(2,2);
traces{1,1}=struct('x',x,'v',W(1,:)/P,'d',wallJets*g.wallWidth/P);
traces{1,2}=struct('x',y,'v',v.value(:)'/P,'d',v.derivativeY(:)'*g.verticalWidth/P);
traces{2,1}=struct('x',x,'v',W(1,:)/P,'d',[]);
traces{2,2}=struct('x',y,'v',linearVertical(:)'/P,'d',[]);
Cx=exp(s.scale.logC_l);Cw=exp(s.scale.logC_omega);
record=struct('checkpointFile',file,'step',s.step,'caseId',s.runMetadata.caseId, ...
 'canonicalTime',s.scale.canonicalTime,'equivalentTau',l.parentCanonicalTime+l.canonicalCovarianceFactor*s.scale.canonicalTime, ...
 'absolutePhysicalTime',l.absolutePhysicalEpoch+s.scale.physicalTime,'physicalPeak',P*Cx/Cw, ...
 'physicalPeakX',g.peakX/Cx,'physicalWidths',[g.wallWidth,g.verticalWidth]/Cx, ...
 'nodeCount',[numel(s.x),numel(s.y)],'remeshCount',s.remeshCount,'box',[s.x([1,end]),s.y([1,end])'], ...
 'coreCells',[h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)], ...
 'physicalRhoXInf',h.common.physicalRhoXInf(end),'geometry',g,'pairedNativeAxesAndField',true,'actualPPrimeEvaluated',false);
f=struct('record',record,'traces',{traces});
% Check the reduced 1D representation against the original tensor observer.
for direction=1:2
 z=linspace(reg.intervals(direction,1),reg.intervals(direction,2),reg.pointCounts(direction));
 if direction==1,xq=g.peakX+g.wallWidth*z;yq=zeros(size(z));else,xq=g.peakX+zeros(size(z));yq=g.verticalWidth*z;end
 full=ipm_accellab_tensor_hermite(W,s.x,s.y,Dx,Dy,xq,yq);
 err=max(abs(full.value/P-sample(traces{1,direction},z,'hermite')));
 assert(err<1e-10,'ipm:ResearchReducedTrace','Reduced and full tensor observations disagree.');
 f.record.reducedTraceMaximumDifference(direction)=err;
end
end
function m=fit_models(values,weights,count)
n=size(values,1);assert(n==8&&all(isfinite(values),'all'));
t=(0:n-1)';future=(n:n+count-1)';D=[ones(n,1),t-mean(t)];
linear=[ones(count,1),future-mean(t)]*(D\eye(n));persistence=zeros(count,n);persistence(:,end)=1;
W=values.*sqrt(weights(:)');delta=diff(W,1,1)';X=delta(:,1:end-1);Y=delta(:,2:end);
[U,S,V]=svd(X,'econ');sv=diag(S);
assert(~isempty(sv)&&sv(1)>1e-12*max(norm(W,'fro'),realmin),'ipm:ResearchLowModeSignal','Training increments below fixed signal floor.');
u=U(:,1);v=V(:,1);sigma=S(1,1);B=Y*v/sigma;lambda=u'*B;c=u'*delta(:,end);
coefB=zeros(1,n);for j=1:n-2,coefB(j+2)=coefB(j+2)+v(j)/sigma;coefB(j+1)=coefB(j+1)-v(j)/sigma;end
rankOne=persistence;gain=0;power=1;
for k=1:count,gain=gain+power;rankOne(k,:)=rankOne(k,:)+c*gain*coefB;power=power*lambda;end
fitResidual=norm(Y-B*(u'*X),'fro')/norm(Y,'fro');
assert(max(abs(sum(rankOne,2)-1))<1e-10&&all(isfinite(rankOne),'all'));
m=struct('coefficients',{{persistence,linear,rankOne}},'diagnostics',struct('rank',1,'ridge',0,'lambda',lambda, ...
 'singularValues',sv','trainingIncrementRelativeResidual',fitResidual,'fractionFirstSingularEnergy',sv(1)^2/sum(sv.^2), ...
 'stableDecayMode',abs(lambda)<1-1e-8,'negativeOrComplexDecayInterpretationExcluded',lambda<=0, ...
 'coefficientAbsoluteSumMaximum',max(sum(abs(rankOne),2))));
end
function val=sample(t,z,method)
assert(all(z>=t.x(1)&z<=t.x(end)));
if strcmp(method,'linear'),val=interp1(t.x,t.v,z,'linear');return;end
j=discretize(z,t.x);j(z==t.x(end))=numel(t.x)-1;assert(all(isfinite(j)));
h=t.x(j+1)-t.x(j);q=(z-t.x(j))./h;
val=(2*q.^3-3*q.^2+1).*t.v(j)+(q.^3-2*q.^2+q).*h.*t.d(j)+ ...
 (-2*q.^3+3*q.^2).*t.v(j+1)+(q.^3-q.^2).*h.*t.d(j+1);
end
function [q,w]=gauss_nodes(training,f,interval,direction)
knots=interval;
for frame=[training,{f}]
 x=frame{1}.traces{1,direction}.x;knots=[knots,x(x>interval(1)&x<interval(2))]; %#ok<AGROW>
end
knots=unique(knots);a=(knots(1:end-1)+knots(2:end))/2;b=diff(knots)/2;
z=[-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7),sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
v=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
q=reshape(a'+b'*z,1,[]);w=reshape(b'*v,1,[]);assert(all(w>0));
end
function e=errors(p,a,w,pc,ac,pf,af)
e=struct('relativeL2',sqrt(sum(w.*(p-a).^2)/sum(w.*a.^2)), ...
 'sampledRelativeInf',max(abs(pf-af))/max(abs(af)), ...
 'sampledInfRefinementChange',abs(max(abs(pf-af))/max(abs(af))-max(abs(pc-ac))/max(abs(ac))), ...
 'integration','Gauss4 on union of all training and heldout native pullback cells; exact for represented polynomial squared error, not PDE accuracy');
end
function a=geometry_values(r)
a=[log(r.physicalPeak),r.physicalPeakX,log(r.physicalWidths)];
end
function d=decision(r)
% Registered descriptive screening gates. Never a PDE acceptance criterion.
L=zeros(numel(r.futureCases),2,3);I=L;
for k=1:numel(r.futureCases)
 for j=1:2
  for m=1:3
   q=r.futureCases{k}.shape{1,j}.modelErrors{m};L(k,j,m)=q.relativeL2;I(k,j,m)=q.sampledRelativeInf;
  end
 end
end
better=all(L(:,:,3)<L(:,:,1)&L(:,:,3)<L(:,:,2),'all')&&all(I(:,:,3)<I(:,:,1)&I(:,:,3)<I(:,:,2),'all');
small=all(L(:,:,3)<=.005,'all')&&all(I(:,:,3)<=.01,'all');
stable=all(cellfun(@(v)v.stableDecayMode&&~v.negativeOrComplexDecayInterpretationExcluded,r.shapeModelDiagnostics));
d=struct('rankOneBeatsBothBaselinesAllTimesDirectionsNorms',better,'rankOneBelowRegisteredHalfPercentL2OnePercentInf',small, ...
 'positiveStableTrainingModes',stable,'finitePredictionScreenPassed',better&&small&&stable, ...
 'PDEAccelerationAccepted',false,'scientificSpaceErrorCertified',false,'limitEstablished',false);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
