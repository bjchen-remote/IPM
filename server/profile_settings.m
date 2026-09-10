function settings = profile_settings(projectRoot)
% Edit only the values in this file before starting a new server run.

settings = struct();
settings.boxHalfWidth = 8;
settings.boxHeight = 4;
settings.initialNodeCount = [321,161];
settings.canonicalFinalTime = 16;
settings.maximumSteps = 60000;
settings.maximumTimeStep = 0.005;
settings.minimumTimeStep = 1e-12;
settings.outputEvery = 0.01;
settings.checkpointEvery = 0.1;

% The run remains admissible while both adjacent-cell ratios are <= 2.
% This controls the grid-smoothness stop. Other numerical safety stops stay on.
settings.maximumAdjacentGridRatio = 2;
settings.maximumTotalNodes = 310000;

settings.outputDirectory = fullfile(projectRoot,'runs','profile_H8_tau16');
settings.storeSnapshots = false;
settings.makePlots = false;
settings.verbose = true;

% Leave empty for a new physical t=0 run. Set an absolute native checkpoint
% path to resume the same frozen numerical problem with longer horizons.
settings.restartCheckpoint = '';
end
