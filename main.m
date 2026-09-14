% High-resolution 4n-by-n server entry. Run this file by its absolute path.
% The folder containing main.m is discovered at runtime, so the package may
% be copied to any absolute server location without editing MATLAB paths.
mainPath = [mfilename('fullpath'),'.m'];
releaseRoot = fileparts(mainPath);
assert(isfolder(fullfile(releaseRoot,'+ipm')) && ...
    isfolder(fullfile(releaseRoot,'server')), ...
    'ipm:HighResolutionPackage','Run main.m from the complete release package.');
restoredefaultpath;
addpath(releaseRoot,fullfile(releaseRoot,'server'));

% 1280 X cells and 320 Y cells: an exact 4:1 cell-count family.
n = 320;
maximumTotalNodes = 1700000;
settings = high_res_settings(releaseRoot,n,maximumTotalNodes);
settings.outputDirectory = fullfile(releaseRoot,'runs', ...
    sprintf('profile_H8_4to1_n%d',n));

% To resume this same run, set an absolute signed checkpoint path here and
% retain the same numerical settings. Use a new output directory for a new run.
settings.restartCheckpoint = '';
fprintf('High-resolution main: %s\n',mainPath);
if strcmp(getenv('IPM_PREFLIGHT_ONLY'),'1')
    preview = ipm.config.resolve(ipm.config.longTimeProfile( ...
        rmfield(settings,'restartCheckpoint')));
    fprintf('Preflight only: initial grid %d x %d, node cap %d, output %s\n', ...
        preview.grid.nx,preview.grid.ny, ...
        preview.remesh.autonomousMesh.nodeFamily.maximumTotalNodes, ...
        settings.outputDirectory);
    return
end
if isempty(settings.restartCheckpoint) && isfolder(settings.outputDirectory)
    prior = dir(settings.outputDirectory);
    assert(numel(prior)<=2,'ipm:HighResolutionOutputExists', ...
        'A fresh run needs an empty output directory; choose a new name or resume a signed checkpoint.');
end
launch_profile(settings,releaseRoot);
