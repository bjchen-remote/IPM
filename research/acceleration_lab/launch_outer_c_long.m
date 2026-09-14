% Start an isolated original-t=0 R2.2 outer-density-C experiment.
% Run with: matlab -batch "run('research/acceleration_lab/launch_outer_c_long.m')"

worktreeRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
workspaceRoot = fileparts(worktreeRoot);
releaseRoot = fullfile(workspaceRoot,'releases', ...
    'ipm_long_time_server_r2_2_20260914');
runDirectory = fullfile(worktreeRoot,'result','verification', ...
    'outer_c_long_20260914','run');
if ~isfolder(releaseRoot)
    error('ipm:OuterCLongReleaseMissing','R2.2 release directory is missing.');
end
if isfolder(runDirectory) && ~isempty(dir(fullfile(runDirectory,'*')))
    error('ipm:OuterCLongRunExists', ...
        'The long-run output directory already exists; never overwrite it.');
end
mkdir(runDirectory);

restoredefaultpath;
addpath(releaseRoot,fullfile(releaseRoot,'server'));
solverFile = which('ipm.solve');
assert(startsWith(solverFile,releaseRoot), ...
    'ipm:OuterCLongSolverPath','The solver must come from immutable R2.2.');

settings = profile_settings(releaseRoot);
settings.outputDirectory = runDirectory;
settings.checkpointEvery = 0.5;
settings.restartCheckpoint = '';
restartCheckpoint = settings.restartCheckpoint;
settings = rmfield(settings,'restartCheckpoint');
opts = ipm.config.longTimeProfile(settings);
config = ipm.config.resolve(opts);
assert(config.remesh.autonomousMesh.enabled && ...
    config.remesh.autonomousMesh.version == 5 && ...
    strcmp(config.scaling.cOmegaGauge,'outer_wall_density_window_l2') && ...
    config.remesh.autonomousMesh.nodeFamily.maximumTotalNodes == 310000, ...
    'ipm:OuterCLongConfig','The long-run numerical contract changed.');
save(fullfile(runDirectory,'launch_configuration.mat'), ...
    'settings','opts','config','restartCheckpoint','solverFile','-v7.3');

diary(fullfile(runDirectory,'console.log'));
fprintf('IPM R2.2 outer-C long experiment\n');
fprintf('  release: %s\n',releaseRoot);
fprintf('  solver: %s\n',solverFile);
fprintf('  output: %s\n',runDirectory);
fprintf('  MATLAB: %s\n',version);
fprintf('  initial grid: %d x %d\n',config.grid.nx,config.grid.ny);
fprintf('  autonomous mesh: version %d, node cap %d\n', ...
    config.remesh.autonomousMesh.version, ...
    config.remesh.autonomousMesh.nodeFamily.maximumTotalNodes);
fprintf('  C gauge: %s, radius %.6g\n', ...
    config.scaling.cOmegaGauge,config.scaling.omegaGaugeWindowRadius);
fprintf('  canonical horizon: %.6g, checkpoint interval %.6g\n', ...
    config.time.finalTime,settings.checkpointEvery);
result = ipm.solve(opts);
fprintf('  final reason: %s\n',result.state.stopReason);
diary off;
