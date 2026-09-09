function report = ipm_perflab_remainder_limits(sourceDirectory,outputDirectory)
%IPM_PERFLAB_REMAINDER_LIMITS Necessary RK2 bounds of the actual frozen F-B.
% Reads previously measured WENO Jacobians; no new PDE, LU or large-grid work.
% The result is a local linearized restriction, not nonlinear stability.
assert(~exist(outputDirectory,'dir'),'Pass a new directory.');
mkdir(outputDirectory);
items=cell(1,3);
sizes=[17,33,65];
for k=1:numel(sizes)
    data=load(fullfile(sourceDirectory,sprintf('linearized_row_%d.mat',sizes(k))));
    eigenvalues=eig(data.N);
    scale=max(abs(eigenvalues));
    decaying=real(eigenvalues)<-1e-8*scale;
    positive=real(eigenvalues)>1e-8*scale;
    limits=Inf(size(eigenvalues));
    for j=find(decaying).'
        z=eigenvalues(j); s=real(z); a=abs(z)^2;
        % |1+h*z+(h*z)^2/2|^2-1 divided by positive h.
        rootsH=roots([0.25*a^2,s*a,2*s^2,2*s]);
        rootsH=real(rootsH(abs(imag(rootsH))<1e-7 & real(rootsH)>0));
        limits(j)=min(rootsH);
    end
    [minimum,index]=min(limits);
    z=eigenvalues(index);
    h=[0.5,1,2]*minimum;
    amplification=abs(1+h*z+0.5*(h*z).^2);
    exact=abs(exp(h*z));
    items{k}=struct('nx',sizes(k),'zeroResidualVelocity',true, ...
        'decayingModeCount',nnz(decaying),'growingModeCount',nnz(positive), ...
        'nearNeutralModesExcluded',nnz(~decaying & ~positive), ...
        'minimumDecayModeRk2Dt',minimum, ...
        'equivalentOriginalFullVelocityCfl',minimum*data.item.fullRate, ...
        'criticalEigenvalueReal',real(z),'criticalEigenvalueImag',imag(z), ...
        'probeDt',h,'rk2CriticalModeAmplification',amplification, ...
        'exactCriticalModeAmplification',exact, ...
        'qualification','Necessary frozen decay-mode bound only; no nonlinear admission.');
    fprintf('N %d nodes: RK2 decaying-mode dt<=%.6g, original full CFL<=%.6g; 2x limit amplification %.4g vs exact %.4g.\n', ...
        sizes(k),minimum,minimum*data.item.fullRate,amplification(3),exact(3));
end
report=struct('schemaVersion',1,'kind','frozen_remainder_rk2_necessary_limits', ...
    'items',{items},'nonlinearStableDtEstablished',false);
save(fullfile(outputDirectory,'remainder_limits.mat'),'report');
fid=fopen(fullfile(outputDirectory,'remainder_limits.json'),'w');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
end
