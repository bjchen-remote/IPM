function summary=ipm_perflab_rate_one_sided(directory)
%IPM_PERFLAB_RATE_ONE_SIDED Re-read saved fields; no transport, LU or PDE.
assert(maxNumCompThreads==10);
loaded=load(fullfile(directory,'report.mat'),'report');r=loaded.report;
d=load(r.registration.cacheFile,'Omega','Dx','Dy','x','y');
prior=load(r.registration.observationReportFile,'report');window=prior.report.localWindow;
width=r.nativeGeometry.wallCoreWidth;
forward=zeros(4,3,2,2);backward=forward;centered=forward;
for e=1:2
    for k=1:3
        record=r.records{e,k};fields=load(record.fieldFile,'pf','mf');h=record.actualStep;
        for mode=1:2
            plus=coordinates(fields.pf.F{mode},d,window,width);
            minus=coordinates(fields.mf.F{mode},d,window,width);
            forward(:,k,e,mode)=(plus-r.baseline.coordinates(:,mode))/h;
            backward(:,k,e,mode)=(r.baseline.coordinates(:,mode)-minus)/h;
            centered(:,k,e,mode)=(plus-minus)/(2*h);
            assert(isequaln(centered(:,k,e,mode),record.modes{mode}.geometryResponse));
        end
    end
end
rows=cell(2,2);
for e=1:2
    for mode=1:2
        jp=forward(:,:,e,mode);jm=backward(:,:,e,mode);
        rows{e,mode}=struct('relativeStep',r.registration.relativeSteps(e),'mode',r.registration.modes{mode}, ...
            'forwardJacobian',jp,'backwardJacobian',jm,'centeredJacobian',centered(:,:,e,mode), ...
            'oneSidedGap',jp-jm,'gapRelativeToCentered',norm(jp-jm,'fro')/norm(centered(:,:,e,mode),'fro'), ...
            'forwardThreeConstraintRcond',rcond(jp(1:3,:)),'backwardThreeConstraintRcond',rcond(jm(1:3,:)));
    end
end
summary=struct('kind','one_sided_fixed_native_geometry_response_from_saved_fields', ...
    'sourceReport',fullfile(directory,'report.mat'),'rows',{rows},'savedCentralDifferencesExactlyReproduced',true, ...
    'poissonEvaluations',0,'transportEvaluations',0,'pdeSteps',0, ...
    'interpretation','Finite one-sided differences at two registered steps. A persistent gap is evidence against a smooth Jacobian at this tied global-LF maximum; no asymptotic differentiability claim is made.');
save(fullfile(directory,'one_sided_geometry.mat'),'summary');
fid=fopen(fullfile(directory,'one_sided_geometry.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(summary,PrettyPrint=true),'char');
fprintf('SAVED_RATE_ONE_SIDED_COMPLETE noTransport=1 centerExact=1\n');
end
function v=coordinates(F,d,window,width)
r=ipm_accellab_continuous_inner_rates(d.Omega,F*d.Dx',d.x,d.y,d.Dx,d.Dy,window);
assert(r.valid);v=[r.peakPrime/r.peak.value;r.translationRate/width;r.logScaleXRate;r.logScaleYRate];
end
