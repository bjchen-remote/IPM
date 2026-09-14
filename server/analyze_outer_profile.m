function report=analyze_outer_profile(runDirectory,outputFile,maximumFiles)
%ANALYZE_OUTER_PROFILE Read-only cutoff norms for a new-C Profile run.
%   REPORT=ANALYZE_OUTER_PROFILE(RUNDIRECTORY,OUTPUTFILE,MAXIMUMFILES)
%   validates accepted checkpoints and writes a fresh JSON report. It does
%   not resume or change the solver trajectory. OUTPUTFILE is optional;
%   the default is a new timestamped file under RUNDIRECTORY/analysis.
serverDirectory=fileparts(mfilename('fullpath'));
projectRoot=fileparts(serverDirectory);
addpath(projectRoot,fullfile(projectRoot,'research','acceleration_lab'));
solverFile=which('ipm.output.restoreCheckpoint');
assert(startsWith(solverFile,projectRoot), ...
    'ipm:OuterProfileSolverPath', ...
    'The analysis must use the solver from this package.');
if nargin<1 || isempty(runDirectory)
    settings=profile_settings(projectRoot);
    runDirectory=settings.outputDirectory;
end
if isstring(runDirectory) && isscalar(runDirectory)
    runDirectory=char(runDirectory);
end
assert(ischar(runDirectory) && isrow(runDirectory) && ...
    isfolder(runDirectory),'ipm:OuterProfileRunDirectory', ...
    'Provide an existing new-C run directory.');
if nargin<2 || isempty(outputFile)
    analysisDirectory=fullfile(runDirectory,'analysis');
    if ~isfolder(analysisDirectory),mkdir(analysisDirectory);end
    stamp=char(datetime('now','TimeZone','UTC', ...
        'Format','yyyyMMdd''T''HHmmssSSS'));
    outputFile=fullfile(analysisDirectory, ...
        ['outer_cutoff_' stamp '.json']);
end
if isstring(outputFile) && isscalar(outputFile)
    outputFile=char(outputFile);
end
if nargin<3 || isempty(maximumFiles),maximumFiles=8;end
report=analyze_outer_cutoff_checkpoints( ...
    runDirectory,outputFile,maximumFiles);
fprintf('Outer Profile cutoff report: %s\n',outputFile);
end
