function report = ipm_accellab_probe_continuous_fresh_cache(cacheFile,continuationReportFile,outputRoot)
%IPM_ACCELLAB_PROBE_CONTINUOUS_FRESH_CACHE Actual fresh case, zero new RHS/LU.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['continuous_inner_fresh_',token]);mkdir(destination);
registration=struct('cacheFile',cacheFile,'continuationReportFile',continuationReportFile, ...
    'outputDirectory',destination,'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0, ...
    'windowRule','[0,2*lineage.physicalAnchorX] in this fresh case computational coordinates');
save(fullfile(destination,'registration.mat'),'registration');
try
    data=load(cacheFile,'Omega','FX','x','y','Dx','Dy','scale','nativeRates','cacheAudit');
    source=jsondecode(fileread(continuationReportFile));a=data.cacheAudit;lineage=a.lineage;
    assert(a.rhsCachePairedExactly && a.extraTimeSteps==0 && source.continuationQualified && ...
        all(structfun(@(v)v,source.checks)) && strcmp(a.sourceCheckpoint,source.checkpointFile));
    assert(strcmp(lineage.kind,'fresh_late_physical_initial_value_branch') && ...
        ~lineage.parentHistoryInherited && ~lineage.parentCheckpointModified && ...
        a.absolutePhysicalEpoch==lineage.absolutePhysicalEpoch && data.nativeRates(3)==0);
    lambda=lineage.canonicalCovarianceFactor;cx0=lineage.parentCx;
    assert(lambda>0 && cx0>0 && abs(lambda-cx0/lineage.parentComega)<1e-12*lambda);
    assert(isequal(data.Dx,ipm.mesh.fdMatrix(data.x,1,7)) && isequal(data.Dy,ipm.mesh.fdMatrix(data.y,1,7)));
    localWindow=[0,2*lineage.physicalAnchorX];assert(localWindow(2)<=max(data.x));
    rates=ipm_accellab_continuous_inner_rates(data.Omega,data.FX,data.x,data.y,data.Dx,data.Dy,localWindow);
    assert(rates.valid,'ipm:FreshContinuousGeometry','Fresh cached continuous geometry is invalid.');
    globalPeak=ipm_accellab_hermite_peak(data.Omega(1,:),data.FX(1,:),data.x,data.Dx,[0,max(data.x)]);
    sameGlobalPeak=globalPeak.valid && globalPeak.x==rates.peak.x && globalPeak.value==rates.peak.value;
    assert(sameGlobalPeak,'ipm:FreshContinuousWindow','The registered local peak differs from the global positive-half maximum.');
    [residual,details]=ipm_accellab_continuous_residual(data.Omega,data.FX,data.x,data.y,data.Dx,data.Dy,rates);
    parent=struct('peakPosition',cx0*rates.peak.x,'wallCoreWidth',cx0*rates.wallCoreWidth, ...
        'verticalCoreWidth',cx0*rates.verticalCoreWidth,'peakValue',rates.peak.value/lambda, ...
        'peakPrime',rates.peakPrime/lambda^2,'peakRelativeRate',rates.peakPrime/rates.peak.value/lambda, ...
        'translationRate',cx0*rates.translationRate/lambda,'logScaleXRate',rates.logScaleXRate/lambda, ...
        'logScaleYRate',rates.logScaleYRate/lambda,'nativeCL',data.nativeRates(1)/lambda, ...
        'nativeCOmega',data.nativeRates(2)/lambda, ...
        'physicalInnerCompressionX',(data.nativeRates(1)-rates.logScaleXRate)/lambda, ...
        'physicalInnerCompressionY',(data.nativeRates(1)-rates.logScaleYRate)/lambda);
    for key={'coreL2','holdoutL2','fullLocalL2'}
        name=key{1};parent.UTimeDerivativeL2.(name)=residual.gauss.G.(name)/lambda;
        parent.UTimeDerivativeOverProfile.(name)=residual.gauss.GOverProfile.(name)/lambda;
        parent.UTimeDerivativeOverNormalizedRhs.(name)=residual.gauss.GOverAmplitudeNormalizedRhs.(name);
    end
    absoluteTime=lineage.absolutePhysicalEpoch+data.scale.physicalTime;
    parentTau=lineage.parentCanonicalTime+lambda*data.scale.canonicalTime;
    assert(abs(absoluteTime-source.absolutePhysicalTime)<1e-12 && abs(parentTau-source.parentEquivalentTau)<1e-12);
    report=struct('status','completed_actual_fresh_cache_no_LU','registration',registration,'cacheAudit',a, ...
        'sourceContinuationChecks',source.checks,'gridSize',size(data.Omega),'freshStep',source.step, ...
        'freshScale',data.scale,'freshNativeRates',data.nativeRates,'absolutePhysicalTime',absoluteTime, ...
        'parentEquivalentTau',parentTau,'localWindow',localWindow,'sameGlobalPositiveHalfPeak',sameGlobalPeak, ...
        'coordinates',rates,'residual',residual,'parentUnits',parent, ...
        'unitConversion','Xparent=Cx0*Xfresh, OmegaParent=OmegaFresh/lambda0, tauParentIncrement=lambda0*tauFresh, Uparent=Ufresh, U_tauParent=U_tauFresh/lambda0.', ...
        'interpretation','Actual original quadratic-gauge RHS of a separately qualified fresh finite-box case. Clock and coordinate conversion is explicit; these are instantaneous observations, not an identity between discretized cases, a new gauge, an accepted accelerator or a singularity proof.');
    save(fullfile(destination,'residual_fields.mat'),'details','-v7.3');save(fullfile(destination,'report.mat'),'report');
    write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',registration,'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('CONTINUOUS_INNER_FRESH %s\n',destination);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
