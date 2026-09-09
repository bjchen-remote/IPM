function report=ipm_accellab_benchmark_local_lf_cache(cacheFile,outputRoot)
%IPM_ACCELLAB_BENCHMARK_LOCAL_LF_CACHE Actual small-line cost, no full RHS/LU.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['local_lf_cost_',token]);mkdir(destination);
d=load(cacheFile,'rho','flow','transportOps','cacheAudit');ops=d.transportOps;
u1=d.flow.u1+d.flow.c_l*ops.X+d.flow.c_r;u2=d.flow.u2+d.flow.c_l*ops.Y;
p=struct('cacheFile',cacheFile,'outputDirectory',destination,'lineCounts',[1,16], ...
    'fullOperatorPredictedSecondsCeiling',180,'newPoissonEvaluations',0,'pdeSteps',0, ...
    'fullRhsEvaluations',0,'normalization','conditional_global');
save(fullfile(destination,'registration.mat'),'p');
options=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',ops.wenoEpsilon,'alphaMode','face','normalization','conditional_global');
rows=cell(2,2);
for axis=1:2
    q=d.rho;a=u1;J=ops.metricX;h=ops.computationalSpacingX;
    if axis==2,q=d.rho';a=u2';J=ops.metricY';h=ops.computationalSpacingY;end
    % Maintained metrics are usually one-dimensional broadcasting arrays.
    J=ones(size(q),'like',q).*J;
    for j=1:2
        n=p.lineCounts(j);timer=tic;
        [~,detail]=ipm_accellab_local_lf_derivative(q(1:n,:),a(1:n,:),J(1:n,:),h,options);
        elapsed=toc(timer);
        rows{axis,j}=struct('axis',axis,'lineCount',n,'nodeCount',size(q,2),'seconds',elapsed, ...
            'ghostMetricMinimum',detail.ghostMetricMinimum,'maximumFaceAlpha',max(detail.faceAlpha,[],'all'));
    end
end
% Each axis requires one density and one unit call. This linear extrapolation
% is deliberately labelled an estimate; larger arrays can change bandwidth.
predicted=2*(rows{1,2}.seconds*ops.ny/16+rows{2,2}.seconds*ops.nx/16);
report=struct('registration',p,'sourceCheckpoint',d.cacheAudit.sourceCheckpoint,'measurements',{rows}, ...
    'predictedFullFourKernelSeconds',predicted,'withinRegisteredCeiling',predicted<=p.fullOperatorPredictedSecondsCeiling, ...
    'interpretation','Cost estimate from the first 1/16 actual lines in each axis. No full density RHS or smooth residual evaluated. Conditional-global split scaling is quadratic in node count per line; extrapolation may underpredict large-array costs.');
save(fullfile(destination,'report.mat'),'report');
fid=fopen(fullfile(destination,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('LOCAL_LF_COST %s predicted=%.3f within=%d\n',destination,predicted,report.withinRegisteredCeiling);
end
