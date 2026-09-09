function report = ipm_accellab_compare_saved_candidate(experimentFile,rejectedFieldFile,outputRoot)
%IPM_ACCELLAB_COMPARE_SAVED_CANDIDATE Parent-window-only fixed q512 comparison.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
assert(maxNumCompThreads == 10,'ipm:ObserverNativeThreads','Native signature validation requires ten threads.');
loaded = load(experimentFile,'experiment'); source = loaded.experiment;
saved = load(rejectedFieldFile,'rhoTrial','x','y','report');
assert(~source.trial.accepted && ~saved.report.accepted && ~saved.report.isNativeCheckpoint, ...
    'ipm:ObserverOriginalDecision','The source must be the recorded rejected independent candidate.');
[~,token] = fileparts(tempname);
destination = fullfile(outputRoot,['fixed_observers_',token]); mkdir(destination);
registration = struct('experimentFile',experimentFile,'rejectedFieldFile',rejectedFieldFile, ...
    'sourceCheckpoint',source.checkpointFiles{end},'outputDirectory',destination, ...
    'fixedDamping',source.trial.attempts.damping,'originalDecision','rejected_unchanged', ...
    'parameterScan',false,'plannedRestoreCalls',1,'plannedExtraRhsEvaluations',2);
save(fullfile(destination,'registration.mat'),'registration');
[checkpoint,sourceAudit,expectedDx,expectedDy] = ipm_accellab_verify_rejected_source(source,saved);
base = ipm.output.restoreCheckpoint(checkpoint,struct('saveResults',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false));
assert(isequal(base.ops.Dx,expectedDx) && isequal(base.ops.Dy,expectedDy), ...
    'ipm:ObserverNativeDerivatives','Restored maintained derivatives differ from the source-array reconstruction.');
jump = norm(saved.rhoTrial(:)-base.rho(:))/max(norm(base.rho(:)),realmin);
assert(abs(jump-source.trial.attempts.relativeJump) < 1e-13, ...
    'ipm:ObserverRejectedField','The fixed candidate field does not reproduce its recorded jump.');
[report,fields] = ipm_accellab_compare_observers(base,saved.rhoTrial,source.trial);
report.registration = registration;
report.sourceAudit = sourceAudit;
report.maintainedDerivativeOperatorsExact = true;
x = base.ops.x; y = base.ops.y; Dx = base.ops.Dx; Dy = base.ops.Dy;
save(fullfile(destination,'observer_fields.mat'),'fields','x','y','Dx','Dy','-v7.3');
save(fullfile(destination,'report.mat'),'report','-v7.3');
fid = fopen(fullfile(destination,'report.json'),'w');
assert(fid >= 0,'ipm:ObserverOutput','Cannot write observer comparison report.');
cleaner = onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
clear cleaner;
fprintf('FIXED_OBSERVER_OUTPUT %s original_decision=rejected_unchanged\n',destination);
end
