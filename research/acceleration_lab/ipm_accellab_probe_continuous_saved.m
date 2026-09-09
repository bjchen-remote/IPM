function report = ipm_accellab_probe_continuous_saved(observerFile,sourceReportFile,outputRoot)
%IPM_ACCELLAB_PROBE_CONTINUOUS_SAVED Reuse only the old native baseline cache.
% The rejected second field is never proposed, evaluated or relabeled here.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['continuous_inner_cached_',token]); mkdir(destination);
registration = struct('observerFile',observerFile,'sourceReportFile',sourceReportFile,'outputDirectory',destination, ...
    'selectedSavedField',1,'fieldIdentity','original accepted native baseline at step2270, tau4.783201944539746', ...
    'originalCandidateDecision','rejected_unchanged','newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0);
save(fullfile(destination,'registration.mat'),'registration');
try
    data = load(observerFile,'fields','x','y','Dx','Dy'); source = load(sourceReportFile,'report');
    assert(~source.report.accepted && strcmp(source.report.originalDecision,'rejected_unchanged') && ...
        all(source.report.sourceAudit.nativeSignaturesPassed) && source.report.sourceAudit.nativeSteps(end) == 2270, ...
        'ipm:ContinuousCachedSource','The previously audited source is not the registered native baseline.');
    assert(isequal(data.Dx,ipm.mesh.fdMatrix(data.x,1,7)) && isequal(data.Dy,ipm.mesh.fdMatrix(data.y,1,7)));
    baseline = data.fields{1}; x = data.x; y = data.y; Dx = data.Dx; Dy = data.Dy;
    r = ipm_accellab_continuous_inner_rates(baseline.omega,baseline.forcing,x,y,Dx,Dy,[0,x(end)]);
    assert(r.valid,'ipm:ContinuousCachedGeometry','The saved native field has invalid continuous derivative coordinates.');
    [residual,details] = ipm_accellab_continuous_residual(baseline.omega,baseline.forcing,x,y,Dx,Dy,r);
    old = baseline.exactCoordinates;
    clX = old.physicalCompressionXPerCanonicalTime+old.logScaleXRate;
    clY = old.physicalCompressionYPerCanonicalTime+old.logScaleYRate;
    assert(abs(clX-clY) < 1e-12);
    xi = linspace(-2,2,65); eta = linspace(0,3,49);
    chart = ipm_accellab_pullback_hermite(baseline.omega,baseline.forcing,x,y,Dx,Dy,r,xi,eta);
    consistency = max(abs(chart.G-details.sampleFields{1}.G),[],'all');
    assert(consistency < 1e-12);
    report = struct('status','completed_old_native_cache_no_LU','registration',registration, ...
        'step',2270,'canonicalTime',4.783201944539746,'gridSize',size(baseline.omega), ...
        'coordinates',r,'residual',residual,'nativeCLFromPairedCachedDiagnostics',clX, ...
        'physicalInnerCompressionX',clX-r.logScaleXRate,'physicalInnerCompressionY',clX-r.logScaleYRate, ...
        'previousQuadraticLinearCoordinates',old,'sourcePoissonResidual',source.report.baseline.poissonResidual, ...
        'sourceWarningIdentifier',source.report.baseline.warningIdentifier,'sourceWarningMessage',source.report.baseline.warningMessage, ...
        'fullChainRuleImplementationAgreement',consistency,'originalCandidateDecision','rejected_unchanged', ...
        'interpretation','An old native q512 baseline and its saved genuine full FX are reused. This is not the latest tau7.093 state, a fresh RHS, a reaccepted secant candidate or a trajectory prediction.');
    save(fullfile(destination,'residual_fields.mat'),'details','-v7.3');
    save(fullfile(destination,'report.mat'),'report'); write_json(fullfile(destination,'report.json'),report);
catch exception
    failure = struct('registration',registration,'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,'failure.mat'),'failure'); write_json(fullfile(destination,'failure.json'),failure); rethrow(exception);
end
fprintf('CONTINUOUS_INNER_CACHED %s\n',destination);
end

function write_json(file,value)
fid = fopen(file,'w'); assert(fid >= 0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
