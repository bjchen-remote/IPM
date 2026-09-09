function report=ipm_accellab_test_local_lf(outputRoot)
%IPM_ACCELLAB_TEST_LOCAL_LF No-LU research flux contracts and smooth MMS.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['local_lf_kernel_',token]);mkdir(destination);
p=struct('outputDirectory',destination,'grids',[49,65,97,129], ...
    'methods',{{'maintained_global','face_alpha_conditional_global_scale'}}, ...
    'testedLocalNormalization',false,'newPoissonEvaluations',0,'pdeSteps',0, ...
    'freeStreamRelativeTolerance',5e-11,'telescopingRelativeTolerance',5e-13, ...
    'interpretation','Only a different spatial research operator is tested; no original PDE trajectory or gauge is changed.');
save(fullfile(destination,'registration.mat'),'p');rows=cell(1,4);
try
    for j=1:4
        n=p.grids(j);xi=linspace(-1,1,n);h=xi(2)-xi(1);
        x=xi+.13*sin(2*pi*xi)/(2*pi);J=1+.13*cos(2*pi*xi);
        q=exp(-2*x.^2)+.2*sin(3*x);qx=-4*x.*exp(-2*x.^2)+.6*cos(3*x);
        a=.25+.7*x;exact=.7*q+a.*qx;
        native=ipm.field.weno5FluxDerivative(q,a,J,h);
        lineOptions=struct('alphaMode','line');[reference,ref]=ipm_accellab_local_lf_derivative(q,a,J,h,lineOptions);
        [local,detail]=ipm_accellab_local_lf_derivative(q,a,J,h);
        parity=isequal(native,reference);assert(parity,'ipm:LocalLfReference','The copied line-alpha reference is not bitwise native.');
        telescoping=abs(sum(J.*local)*h-(detail.faceFlux(end)-detail.faceFlux(1)))/max(1,max(abs(detail.faceFlux)));
        assert(telescoping<p.telescopingRelativeTolerance);
        constantErrors=[];
        for amplitude=[0,.7,-2,1e-12]
            dc=ipm_accellab_local_lf_derivative(amplitude*ones(size(q)),a,J,h);
            du=ipm_accellab_local_lf_derivative(ones(size(q)),a,J,h);
            constantErrors(end+1)=max(abs(dc-amplitude*du))/max(1,abs(amplitude)*max(abs(du))); %#ok<AGROW>
        end
        assert(max(constantErrors)<p.freeStreamRelativeTolerance);
        [constantVelocity,~]=ipm_accellab_local_lf_derivative(q,ones(size(q)),1,h);
        constantNative=ipm.field.weno5FluxDerivative(q,ones(size(q)),1,h);
        assert(isequal(constantVelocity,constantNative),'ipm:LocalLfConstantSpeed','Constant-speed conditional-global normalization lost native parity.');
        interior=4:n-3;errors=zeros(2,2);
        errors(1,:)=[max(abs(native-exact)),max(abs(native(interior)-exact(interior)))];
        errors(2,:)=[max(abs(local-exact)),max(abs(local(interior)-exact(interior)))];
        rows{j}=struct('nodes',n,'physicalX',x,'metric',J,'errorsFullAndInterior',errors, ...
            'nativeBitwiseReference',parity,'constantVelocityBitwiseNative',true,'freeStreamErrors',constantErrors, ...
            'weightedFluxBoundaryIdentityRelativeError',telescoping,'localFaceAlpha',detail.faceAlpha, ...
            'nativeLineAlpha',ref.faceAlpha(1),'ghostMetricMinimum',detail.ghostMetricMinimum);
    end
    n=129;x=linspace(-1,1,n);h=x(2)-x(1);q=exp(-3*x.^2)+.2*sin(4*x);
    velocity1=.3+.2*x;velocity2=velocity1+2e3*max(abs(x)-.55,0).^4;
    fields=cell(2,2);
    for j=1:2
        a=velocity1;if j==2,a=velocity2;end
        fields{1,j}=ipm.field.weno5FluxDerivative(q,a,1,h);
        fields{2,j}=ipm_accellab_local_lf_derivative(q,a,1,h);
    end
    core=abs(x)<.3;nativeResponse=max(abs(fields{1,2}(core)-fields{1,1}(core)));
    localResponse=max(abs(fields{2,2}(core)-fields{2,1}(core)));
    locality=struct('remoteVelocityChangeStartsAt',.55,'coreWindow',[-.3,.3], ...
        'nativeCoreResponse',nativeResponse,'localAlphaCoreResponse',localResponse, ...
        'ratio',localResponse/nativeResponse,'remainingGlobalWeightNormalization',true, ...
        'meaning','Same local samples and velocity, changed remote velocity on the same grid. This is operator locality only, not equivalence of two global PDEs.');
    reflect=struct('lowerBoundary','reflect','upperBoundary','reflect');
    velocity=sin(pi*x);velocity([1,end])=0;
    [dr,fr]=ipm_accellab_local_lf_derivative(exp(-x.^2),velocity,1,h,reflect);
    reflectionIdentity=abs(sum(dr)*h-(fr.faceFlux(end)-fr.faceFlux(1)));
    assert(reflectionIdentity<1e-11);
    metricRejected=false;
    try
        ipm_accellab_local_lf_derivative(ones(1,9),ones(1,9),linspace(.01,1,9),1/8);
    catch errorMetric
        metricRejected=strcmp(errorMetric.identifier,'ipm:ResearchLocalMetric');
    end
    assert(metricRejected);
    mmsErrors=zeros(4,2,2);
    for j=1:4,mmsErrors(j,:,:)=rows{j}.errorsFullAndInterior;end
    orders=zeros(3,2,2);
    for j=1:3,orders(j,:,:)=log(mmsErrors(j,:,:)./mmsErrors(j+1,:,:))/log((p.grids(j+1)-1)/(p.grids(j)-1));end
    report=struct('status','completed_no_LU_local_alpha_kernel_probe','registration',p,'mmsCases',{rows}, ...
        'observedMmsOrders',orders,'mmsErrors',mmsErrors,'remoteInfluence',locality, ...
        'reflectedBoundaryFluxIdentityAbsoluteError',reflectionIdentity,'nonpositiveGhostMetricRejected',metricRejected, ...
        'kernelContractChecksPassed',true,'pdeTimeValidationPerformed',false,'physicalGaugeCovarianceCertified',false, ...
        'interpretation','Native line-alpha arithmetic and constant-speed limit are bitwise verified. Free-stream and boundary flux contracts passed; MMS errors/orders are measured, not assumed. This alpha-only research variant still uses conditional whole-line weight scales and is not a production qualification or a proof of faster shape convergence.');
    save(fullfile(destination,'fields.mat'),'fields','x','velocity1','velocity2','q','-v7.3');
    save(fullfile(destination,'report.mat'),'report');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack,'completedCases',{rows});
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('LOCAL_LF_KERNEL %s contracts=%d remoteRatio=%.3e localFineFull=%.3e\n', ...
    destination,report.kernelContractChecksPassed,locality.ratio,mmsErrors(end,2,1));
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
