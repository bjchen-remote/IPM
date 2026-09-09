function report = ipm_accellab_analyze_saved_observers(observerFieldsFile,comparisonReportFile,outputRoot)
%IPM_ACCELLAB_ANALYZE_SAVED_OBSERVERS Only saved arrays; no native CP or LU.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
data = load(observerFieldsFile,'fields','x','y','Dx','Dy');
prior = load(comparisonReportFile,'report');
assert(strcmp(prior.report.originalDecision,'rejected_unchanged') && ~prior.report.accepted, ...
    'ipm:SamplingOriginalDecision','The original rejected comparison must remain unchanged.');
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['observer_sampling_',token]);
mkdir(destination);
registration = struct('observerFieldsFile',observerFieldsFile, ...
    'comparisonReportFile',comparisonReportFile,'outputDirectory',destination, ...
    'observerFactors',[1,2,4],'phases',[0,0.5],'nativeStrides',[1,2,4], ...
    'quadratureOrderPerAxis',4,'originalDecision','rejected_unchanged');
save(fullfile(destination,'registration.mat'),'registration');
try
    [report,details] = ipm_accellab_observer_sampling(data.fields,data.x,data.y,data.Dx,data.Dy);
    report.registration = registration;
    report.originalComparison = prior.report;
    save(fullfile(destination,'quadrature_fields.mat'),'details','-v7.3');
    save(fullfile(destination,'report.mat'),'report','-v7.3');
catch exception
    failure = struct('status','failed_no_decision_change','registration',registration, ...
        'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,'failure.mat'),'failure');
    rethrow(exception);
end
fid = fopen(fullfile(destination,'report.json'),'w');
assert(fid >= 0,'ipm:SamplingOutput','Cannot write sampling report.');
cleanup = onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
clear cleanup;
fprintf('NO_LU_OBSERVER_SAMPLING %s\n',destination);
end
