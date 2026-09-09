function report=ipm_accellab_observe_audited_cache(cacheFile,reviewFile,outputRoot)
%IPM_ACCELLAB_OBSERVE_AUDITED_CACHE Complete C1 derivative of audited arrays.
% No checkpoint restore, Poisson operator, RHS evaluation, or time step.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['continuous_audited_cache_',token]);mkdir(destination);
registration=struct('cacheFile',cacheFile,'reviewFile',reviewFile,'outputDirectory',destination, ...
    'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0, ...
    'windowRule','[0,2*lineage.physicalAnchorX] in fresh computational coordinates');
save(fullfile(destination,'registration.mat'),'registration');
try
    data=load(cacheFile,'Omega','FX','x','y','Dx','Dy','scale','nativeRates','cacheAudit');
    evidence=jsondecode(fileread(reviewFile));review=evidence;
    if isfield(evidence,'audit') && isfield(evidence.audit,'strictNativeReview')
        assert(evidence.audit.passed && evidence.audit.sourcePrefixExact && evidence.audit.lineageAndRuntimeReferencesPreserved);
        review=evidence.audit.strictNativeReview;
    elseif isfield(evidence,'review')
        review=evidence.review;
    end
    a=data.cacheAudit;lineage=a.lineage;
    assert(a.rhsCachePairedExactly && a.extraTimeSteps==0 && strcmp(a.sourceCheckpoint,review.checkpointFile));
    assert(review.nativeAndResultEntireHistoryTrusted && review.nativeResultStateClocksAndHistoryExactlyPaired && ...
        review.freshLineagePreserved && review.initialConditionCapturedSamplesExactlyPaired && ...
        review.nativeResultScalesPhysicalSamplesSnapshotsAndConfigExactlyPaired);
    assert(strcmp(lineage.kind,'fresh_late_physical_initial_value_branch') && ...
        ~lineage.parentHistoryInherited && ~lineage.parentCheckpointModified && data.nativeRates(3)==0);
    lambda=lineage.canonicalCovarianceFactor;cx0=lineage.parentCx;
    assert(lambda>0 && cx0>0 && abs(lambda-cx0/lineage.parentComega)<1e-12*lambda);
    assert(isequal(data.Dx,ipm.mesh.fdMatrix(data.x,1,7)) && isequal(data.Dy,ipm.mesh.fdMatrix(data.y,1,7)));
    localWindow=[0,2*lineage.physicalAnchorX];
    rates=ipm_accellab_continuous_inner_rates(data.Omega,data.FX,data.x,data.y,data.Dx,data.Dy,localWindow);
    assert(rates.valid,'ipm:AuditedContinuousGeometry','Continuous geometry is not differentiable on this state.');
    peak=ipm_accellab_hermite_peak(data.Omega(1,:),data.FX(1,:),data.x,data.Dx,[0,max(data.x)]);
    assert(peak.valid && peak.x==rates.peak.x && peak.value==rates.peak.value);
    [residual,details]=ipm_accellab_continuous_residual(data.Omega,data.FX,data.x,data.y,data.Dx,data.Dy,rates);
    absoluteTime=lineage.absolutePhysicalEpoch+data.scale.physicalTime;
    parentTau=lineage.parentCanonicalTime+lambda*data.scale.canonicalTime;
    assert(abs(absoluteTime-review.absolutePhysicalTime)<1e-12 && abs(parentTau-review.parentEquivalentTau)<1e-12);
    physicalClock=exp(data.scale.logC_l-data.scale.logC_omega);
    parent=struct('peakPosition',cx0*rates.peak.x,'wallCoreWidth',cx0*rates.wallCoreWidth, ...
        'verticalCoreWidth',cx0*rates.verticalCoreWidth,'peakValue',rates.peak.value/lambda, ...
        'peakPrime',rates.peakPrime/lambda^2,'peakRelativeRate',rates.peakPrime/rates.peak.value/lambda, ...
        'translationRate',cx0*rates.translationRate/lambda,'logScaleXRate',rates.logScaleXRate/lambda, ...
        'logScaleYRate',rates.logScaleYRate/lambda,'nativeCL',data.nativeRates(1)/lambda, ...
        'nativeCOmega',data.nativeRates(2)/lambda);
    physical=struct('canonicalTimeRate',physicalClock,'peakValue',physicalClock*rates.peak.value, ...
        'wallCoreWidth',rates.wallCoreWidth/exp(data.scale.logC_l), ...
        'verticalCoreWidth',rates.verticalCoreWidth/exp(data.scale.logC_l), ...
        'logWidthXRate',(rates.logScaleXRate-data.nativeRates(1))*physicalClock, ...
        'logWidthYRate',(rates.logScaleYRate-data.nativeRates(1))*physicalClock, ...
        'logPeakRate',(rates.peakPrime/rates.peak.value+data.nativeRates(1)-data.nativeRates(2))*physicalClock);
    for key={'coreL2','holdoutL2','fullLocalL2','coreSampledInfinity','holdoutSampledInfinity','fullLocalSampledInfinity'}
        name=key{1};parent.UTimeDerivative.(name)=residual.gauss.G.(name)/lambda;
        physical.UTimeDerivative.(name)=residual.gauss.G.(name)*physicalClock;
    end
    report=struct('status','completed_audited_cache_no_LU','registration',registration,'cacheAudit',a, ...
        'strictNativeReview',review,'gridSize',size(data.Omega),'absolutePhysicalTime',absoluteTime, ...
        'parentEquivalentTau',parentTau,'freshScale',data.scale,'freshNativeRates',data.nativeRates, ...
        'localWindow',localWindow,'coordinates',rates,'residual',residual,'parentUnits',parent,'physicalUnits',physical, ...
        'interpretation','Complete instantaneous normalized inner-coordinate derivative of the true cached native RHS. FullLocal is the inner observation rectangle only. No accelerator, new gauge, finite-step invariance, convergence rate or singularity claim.');
    save(fullfile(destination,'residual_fields.mat'),'details','-v7.3');save(fullfile(destination,'report.mat'),'report');
    write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',registration,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('AUDITED_C1_CACHE %s coreParent=%.10g holdParent=%.10g corePhysical=%.10g holdPhysical=%.10g\n', ...
    destination,parent.UTimeDerivative.coreL2,parent.UTimeDerivative.holdoutL2,physical.UTimeDerivative.coreL2,physical.UTimeDerivative.holdoutL2);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
