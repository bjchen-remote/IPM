function report=ipm_accellab_test_full_tail_trace(out)
%IPM_ACCELLAB_TEST_FULL_TAIL_TRACE Registered complete-boundary pure MMS.
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
H=8;x=H*[-1,-.75,-.5,-.25,-.1,-1e-4,0,1e-4,.1,.25,.5,.75,1];
y=H*[0,1e-4,.01,.05,.125,.25,.375,.5]';
q=[x(1)+zeros(numel(y),1),y;x(end)+zeros(numel(y),1),y; ...
    x',y(end)+zeros(numel(x),1);x',zeros(numel(x),1)];
nb=size(q,1);holdout=[1,.008;2,.8;6,.8];q=[q;holdout];
reg=struct('kind','complete_Green_exterior_Dirichlet_trace_MMS_v1', ...
    'H',H,'x',x,'y',y,'query',q,'boundaryQueryCount',nb,'interiorHoldout',holdout, ...
    'orders',[32,64,128,256],'methods',{{'radial_analytic','target_polar'}}, ...
    'maximumTraceScaledError',1e-10,'maximumNearAxisQuotientScaledError',1e-9, ...
    'maximumSymmetryScaledError',1e-12,'maximumCovarianceScaledError',1e-10, ...
    'CxValues',[.7,1,1.7],'ComegaValues',[.6,1,1.8], ...
    'perturbation',.2,'perturbationHValues',[4,8,16,32], ...
    'perturbationOrders',[64,128,256],'reuseAmplitudeValues',[.125,.7,1,2], ...
    'reuseBatches',7,'reuseCallsPerBatch',1000,'rebuildRepeats',3, ...
    'source','Omega_leading=x^7/(x^8+y^8); perturbation=d_x[x^2/(x^2+y^2)^(3/2)] outside the box', ...
    'scope','Pure manufactured exterior integrals, not an unforced IPM solution or a production-admitted grid.', ...
    'LUCount',0,'PDECount',0,'injectedIntoSolver',false,'realEvolvedClosureValidated',false);
write_json(fullfile(out,'registration.json'),reg);timer=tic;
analytic=zeros(size(q,1),numel(reg.orders));numeric=analytic;timings=zeros(numel(reg.orders),2);
for k=1:numel(reg.orders)
    t=tic;analytic(:,k)=ipm_accellab_green_full_tail_trace(q,H,method='radial_analytic',order=reg.orders(k));timings(k,1)=toc(t);
    t=tic;numeric(:,k)=ipm_accellab_green_full_tail_trace(q,H,method='target_polar',order=reg.orders(k));timings(k,2)=toc(t);
end
report=reg;report.analyticRadial=analytic;report.targetCenteredPolar=numeric;report.constructionSeconds=timings;
scale=max(H,abs(analytic(:,end)));difference=analytic(:,end)-numeric(:,end);
report.maximumMethodDifference=max(abs(difference));report.methodScaledError=max(abs(difference)./scale);
report.lastAnalyticScaledChange=max(abs(analytic(:,end)-analytic(:,end-1))./scale);
report.lastTargetPolarScaledChange=max(abs(numeric(:,end)-numeric(:,end-1))./scale);
nonzero=q(:,1)~=0&q(:,2)~=0;
qa=H*analytic(nonzero,end)./(q(nonzero,1).*q(nonzero,2));
qn=H*numeric(nonzero,end)./(q(nonzero,1).*q(nonzero,2));
report.nearAxisQuotientAnalytic=qa;report.nearAxisQuotientTargetPolar=qn;
report.nearAxisQuotientScaledError=max(abs(qa-qn)./max(1,abs(qa)));
report.allAxisTracesExactlyZero=all(analytic(~nonzero,:)==0,'all')&&all(numeric(~nonzero,:)==0,'all');
ny=numel(y);nx=numel(x);
report.leftRightOddSymmetryScaledError=max(abs(analytic(1:ny,end)+analytic(ny+(1:ny),end)))/H;
report.cornerDuplicatesExact=isequaln(analytic(ny,end),analytic(2*ny+1,end))&& ...
    isequaln(analytic(2*ny,end),analytic(2*ny+nx,end));
report.quadraturePassed=max([report.methodScaledError,report.lastAnalyticScaledChange,report.lastTargetPolarScaledChange])<=reg.maximumTraceScaledError && ...
    report.nearAxisQuotientScaledError<=reg.maximumNearAxisQuotientScaledError;
template=analytic(1:nb,end);boundaryQuery=q(1:nb,:);
report.unitLeadingTemplate=pack_boundary(template,nx,ny);
report.notAnInnerPointApproximation=struct('innerPointXYAOverH',q(:,1).*q(:,2)*.42476990616722304/H, ...
    'maximumBoundaryScaledDifferenceFromFullTrace', ...
    max(abs(template-q(1:nb,1).*q(1:nb,2)*.42476990616722304/H))/H);
pert=zeros(nb,numel(reg.perturbationOrders));
for k=1:numel(reg.perturbationOrders)
    pert(:,k)=ipm_accellab_green_full_tail_trace(boundaryQuery,H,method='target_polar', ...
        order=reg.perturbationOrders(k),sourceOnly='perturbation');
end
report.manufacturedPerturbation=pert;
report.lastPerturbationScaledChange=max(abs(pert(:,end)-pert(:,end-1)))/max(1,max(abs(pert(:,end))));
baseSum=ipm_accellab_green_full_tail_trace(boundaryQuery,H,method='target_polar', ...
    order=256,sourceOnly='sum',perturbation=reg.perturbation);
report.manufacturedSourceLinearityScaledError=max(abs(baseSum-(numeric(1:nb,end)+reg.perturbation*pert(:,end))))/H;
covariance=struct([]);
for C=reg.CxValues
    for W=reg.ComegaValues
        v=ipm_accellab_green_full_tail_trace(C*boundaryQuery,C*H,method='target_polar', ...
            order=256,Cx=C,Comega=W,sourceOnly='sum',perturbation=reg.perturbation);
        recovered=v/(C*W);
        row=struct('Cx',C,'Comega',W,'canonicalTrace',v,'physicalTrace',recovered, ...
            'scaledError',max(abs(recovered-baseSum))/H);
        if isempty(covariance),covariance=row;else,covariance(end+1)=row;end
    end
end
report.covariance=covariance;report.maximumCovarianceScaledError=max([covariance.scaledError]);
report.perturbationBoxScaling=zeros(nb,numel(reg.perturbationHValues));
for k=1:numel(reg.perturbationHValues)
    hh=reg.perturbationHValues(k);
    report.perturbationBoxScaling(:,k)=ipm_accellab_green_full_tail_trace(hh/H*boundaryQuery,hh, ...
        method='target_polar',order=256,sourceOnly='perturbation');
end
report.perturbationBoundaryHZeroScalingError=max(abs(report.perturbationBoxScaling-pert(:,end)),[],'all');
reuse=struct([]);
for W=reg.reuseAmplitudeValues
    direct=ipm_accellab_green_full_tail_trace(boundaryQuery,H,method='target_polar',order=256,Cx=1.3,Comega=W);
    row=struct('Comega',W,'directFullIntegral',direct,'reusedTemplate',W*template, ...
        'scaledError',max(abs(direct-W*template))/(H*max(1,W)));
    if isempty(reuse),reuse=row;else,reuse(end+1)=row;end
end
report.scalarReuse=reuse;report.maximumScalarReuseScaledError=max([reuse.scaledError]);
% Synthetic existing-boundary data are only for timing scalar multiplication
% and addition; they are not native boundary values or scientific observations.
original=H*1e-3*sin((1:nb)');batchSeconds=zeros(1,reg.reuseBatches);checksum=0;
for b=1:reg.reuseBatches
    t=tic;
    for k=1:reg.reuseCallsPerBatch
        W=reg.reuseAmplitudeValues(1+mod(k-1,numel(reg.reuseAmplitudeValues)));
        combined=original+W*template;checksum=checksum+combined(1+mod(k-1,nb));
    end
    batchSeconds(b)=toc(t);
end
rebuild=zeros(1,reg.rebuildRepeats);
for k=1:reg.rebuildRepeats
    t=tic;built=ipm_accellab_green_full_tail_trace(boundaryQuery,H,method='radial_analytic',order=256);
    rebuild(k)=toc(t);assert(isequaln(built,template));
end
report.cost=struct('boundaryQueries',nb,'reuseBatchSeconds',batchSeconds, ...
    'reuseCallsPerBatch',reg.reuseCallsPerBatch,'medianSecondsPerScalarReuse',median(batchSeconds)/reg.reuseCallsPerBatch, ...
    'rebuildSeconds',rebuild,'medianSecondsPerTemplateRebuild',median(rebuild), ...
    'checksum',checksum,'existingBoundaryIsSynthetic',true,'nativeRHSRuntimeMeasured',false, ...
    'costMeaning','Small registered complete boundary only; includes MATLAB-call and rule overhead for rebuild. No production speedup claim.');
report.allPureChecksPassed=report.quadraturePassed&&report.allAxisTracesExactlyZero&&report.cornerDuplicatesExact&& ...
    report.leftRightOddSymmetryScaledError<=reg.maximumSymmetryScaledError&& ...
    report.lastPerturbationScaledChange<=reg.maximumTraceScaledError&& ...
    report.manufacturedSourceLinearityScaledError<=reg.maximumTraceScaledError&& ...
    report.maximumCovarianceScaledError<=reg.maximumCovarianceScaledError&& ...
    report.perturbationBoundaryHZeroScalingError<=reg.maximumTraceScaledError&& ...
    report.maximumScalarReuseScaledError<=reg.maximumTraceScaledError;
report.wallSeconds=toc(timer);
save(fullfile(out,'report.mat'),'report');write_json(fullfile(out,'report.json'),report);
disp(jsonencode(struct('allPureChecksPassed',report.allPureChecksPassed, ...
    'maximumMethodDifference',report.maximumMethodDifference,'methodScaledError',report.methodScaledError, ...
    'lastAnalyticScaledChange',report.lastAnalyticScaledChange,'lastTargetPolarScaledChange',report.lastTargetPolarScaledChange, ...
    'nearAxisQuotientScaledError',report.nearAxisQuotientScaledError, ...
    'maximumCovarianceScaledError',report.maximumCovarianceScaledError, ...
    'perturbationBoundaryHZeroScalingError',report.perturbationBoundaryHZeroScalingError, ...
    'maximumScalarReuseScaledError',report.maximumScalarReuseScaledError,'cost',report.cost,'wallSeconds',report.wallSeconds)));
end
function b=pack_boundary(values,nx,ny)
b=struct('left',values(1:ny),'right',values(ny+(1:ny)), ...
    'top',values(2*ny+(1:nx))','bottom',values(2*ny+nx+(1:nx))');
end
function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
