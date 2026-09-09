function report=ipm_accellab_green_initial_H_series(registrationFile)
%IPM_ACCELLAB_GREEN_INITIAL_H_SERIES Strict t0 history plus analytic tails.
% No evolution, Poisson, restore, or recomputation of the native rates.
reg=jsondecode(fileread(registrationFile));out=reg.outputDirectory;
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
profile clear;profile on;
original=jsondecode(fileread(reg.initialABReport));
originalMAT=load(reg.initialABMat,'data');originalConfig=originalMAT.data{1}.config;
integral=jsondecode(fileread(reg.shellIntegralReport));
assert(integral.quadratureAgreementPassed);
report=struct('kind','strict_t0_native_and_exact_source_Green_H_series_v1', ...
    'registrationFile',registrationFile,'cases',{{}},'orders',reg.orders, ...
    'infiniteCLComparator',integral.infiniteSourceCLFromInsidePlusTail(end,1), ...
    'infiniteCLComparatorMeaning','Analytic k8 t0 continuum Green integral; not native FD velocity truth.', ...
    'noPDE',true,'noLU',true,'noRestore',true,'endpointComparison',false, ...
    'pureBoundaryAttribution',false,'initialSelectedAxesDiffer',true);
for j=1:2
    old=original.cases(j);o=old.initialObservation;
    entry=struct('H',old.configuredBox(2),'sourceKind','previous_strict_native_t0_audit', ...
        'sourceFile',reg.initialABReport,'checkpointFile',old.checkpointFile, ...
        'initialCL',o.canonicalCL,'initialCOmega',o.canonicalCOmega, ...
        'physicalRateCL',o.physicalRateCL,'physicalRateCOmega',o.physicalRateCOmega, ...
        'physicalTime',o.physicalTime,'canonicalTime',o.canonicalTime, ...
        'initialScale',[o.Cx,o.Cy,o.Comega,o.Xshift], ...
        'physicalRhoXInf',o.physicalRhoXInf,'physicalGradInf',o.physicalGradInf, ...
        'physicalQuadraticPeak',o.physicalQuadraticPeak, ...
        'initialCore',o.historyCoreCells,'initializationSelectedAxes',old.initializationSelectedAxes, ...
        'strictNativeSignatureValidated',old.strictNativeSignatureValidated);
    report.cases{j}=entry;
end
data=cell(1,numel(reg.cases));
for j=1:numel(reg.cases)
    item=reg.cases(j);cp=ipm.output.readCheckpoint(item.checkpointFile);
    resultData=load(item.resultFile,'result');r=resultData.result;s=cp.payload.state;
    assert(isequaln(r.history,cp.payload.log.history) && isequaln(r.state.rho,s.rho));
    assert(isequaln(r.grid.x,s.x) && isequaln(r.grid.y,s.y));
    assert(isequaln(r.config,s.config));
    h=r.history.common;m=r.history.mesh;[trusted,~]=ipm.output.trustedMask(r.history,r.config);
    assert(all(trusted) && h.acceptedStep(1)==0 && h.canonicalTau(1)==0 && h.physicalTime(1)==0 && h.t(1)==0);
    assert(all(m.remeshCount==0) && s.step==4);
    assert(isequaln([h.C_l(1),h.C_y(1),h.C_omega(1),h.X_shift(1)],[1,1,1,0]));
    assert(isequaln(r.config.physics,originalConfig.physics));
    assert(isequaln(r.config.scaling,originalConfig.scaling));
    assert(isequaln(r.config.transport,originalConfig.transport));
    assert(isequaln(r.config.elliptic,originalConfig.elliptic));
    assert(isequaln([s.x(1),s.x(end),s.y(1),s.y(end)],[-item.H,item.H,0,item.H/2]));
    multiplier=h.C_l(1)/h.C_omega(1);
    entry=struct('H',item.H,'sourceKind','strict_native_checkpoint_and_result_first_history', ...
        'sourceFile',item.resultFile,'checkpointFile',item.checkpointFile, ...
        'initialCL',h.c_l(1),'initialCOmega',h.c_omega(1), ...
        'physicalRateCL',h.c_l(1)*multiplier,'physicalRateCOmega',h.c_omega(1)*multiplier, ...
        'physicalTime',h.physicalTime(1),'canonicalTime',h.canonicalTau(1), ...
        'initialScale',[h.C_l(1),h.C_y(1),h.C_omega(1),h.X_shift(1)], ...
        'physicalRhoXInf',h.physicalRhoXInf(1),'physicalGradInf',h.physicalGradInf(1), ...
        'physicalQuadraticPeak',h.physicalQuadraticPeak(1), ...
        'initialCore',[m.coreGridPoints(1),m.verticalCoreGridPoints(1)], ...
        'initialNodes',[numel(s.x),numel(s.y)],'strictNativeSignatureValidated',true, ...
        'completeNativeResultHistoryAndStatePairing',true,'allAvailableHistoryTrusted',true, ...
        'checkpointFinalRatesUsed',false);
    report.cases{j+2}=entry;
    data{j}=struct('x',s.x,'y',s.y,'config',s.config,'firstHistory',first_row(r.history));
end
report.tailCL=zeros(numel(report.cases),numel(reg.orders));
for j=1:numel(report.cases)
    for k=1:numel(reg.orders)
        report.tailCL(j,k)=tail(report.cases{j}.H,reg.orders(k));
    end
end
report.lastTailChange=abs(report.tailCL(:,end)-report.tailCL(:,end-1));
report.quadratureAgreementPassed=all(report.lastTailChange<=reg.absoluteTailAgreementTolerance);
report.analyticInsideCL=report.infiniteCLComparator-report.tailCL(:,end);
report.nativeInitialCL=cellfun(@(c)c.initialCL,report.cases)';
report.nativeInitialCOmega=cellfun(@(c)c.initialCOmega,report.cases)';
report.nativeMinusAnalyticInside=report.nativeInitialCL-report.analyticInsideCL;
report.analyticInsideAdjacentDifference=diff(report.analyticInsideCL);
report.nativeAdjacentDifference=diff(report.nativeInitialCL);
report.adjacentDifferenceMinusAnalyticShell=report.nativeAdjacentDifference-report.analyticInsideAdjacentDifference;
report.tailRelativeToInfinite=report.tailCL(:,end)/report.infiniteCLComparator;
report.HScaledTail=cellfun(@(c)c.H,report.cases)'.*report.tailCL(:,end);
report.claim='At t0 the exact-source exterior integral and native box/initial-grid joint sensitivity are compared. Remaining gaps are not pure Poisson or pure boundary errors.';
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for key={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.field.poisson', ...
        'ipm.output.restoreCheckpoint','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,key{1})),['Forbidden call: ',key{1}]);
end
report.noLUCallGraphChecked=true;
save(fullfile(out,'report.mat'),'report','data','profileInfo','-v7.3');
fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
disp(jsonencode(struct('H',[8,16,32,64,128],'nativeCL',report.nativeInitialCL, ...
    'nativeCOmega',report.nativeInitialCOmega,'analyticInside',report.analyticInsideCL, ...
    'tail',report.tailCL(:,end),'nativeMinusAnalyticInside',report.nativeMinusAnalyticInside, ...
    'quadratureAgreementPassed',report.quadratureAgreementPassed)));
end

function v=tail(H,n)
j=(1:n-1)';b=j./sqrt(4*j.^2-1);[Q,T]=eig(diag(b,1)+diag(b,-1),'vector');
[z,index]=sort(T);w=2*Q(1,index)'.^2;
edges=[0;atan(.5);pi/2];theta=(edges(1:end-1)+edges(2:end))/2+diff(edges)*z'/2;
wt=diff(edges)*w'/2;theta=theta(:)';wt=wt(:)';
cs=cos(theta);sn=sin(theta);r0=H./max(cs,2*sn);
invr=(.5+.5*z)./r0;winvr=.5*w./r0;r=1./invr;
x=r.*cs;y=r.*sn;
integrand=4*x.*y./(pi*((1-x).^2+y.^2).*((1+x).^2+y.^2)).*x.^7./(1+x.^8+y.^8);
v=sum(integrand.*r./invr.^2.*winvr.*wt,'all');
end
function first=first_row(h)
first=struct();
for g=fieldnames(h)'
    key=g{1};first.(key)=struct();
    for f=fieldnames(h.(key))'
        name=f{1};value=h.(key).(name);
        if isempty(value),first.(key).(name)=value;
        else,first.(key).(name)=value(1);end
    end
end
end
