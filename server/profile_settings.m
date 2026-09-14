function settings = profile_settings(projectRoot)
% Edit only the values in this file before starting a new server run.

settings = struct();
settings.boxHalfWidth = 8;
settings.boxHeight = 4;
% Version 5 uses the known H8 starting grid; on larger boxes it screens the
% analytic t=0 datum before any LU for an affordable qualified X/Y pair.
% To pin an initial grid explicitly, replace 'auto' with [Nx,Ny].
settings.initialNodeCount = 'auto';
% These are deliberately distant operational ceilings: a fresh server run
% should not stop at the earlier tau=16 milestone before the Profile study.
% Numerical safety and mesh-capacity gates remain active independently.
settings.canonicalFinalTime = 1000;
settings.maximumSteps = 10000000;
settings.maximumTimeStep = 0.005;
settings.minimumTimeStep = 1e-12;
settings.outputEvery = 0.01;
settings.checkpointEvery = 0.1;

% Ratios <= 2 do not trigger the outer grid-smoothness stop. Automatic
% candidate axes still obey the stricter 1.08 quality gate and all other
% numerical safety stops remain active.
settings.maximumAdjacentGridRatio = 2;
% The outer-density window is centered at (2,0), away from the shrinking
% gradient core. The previous quadratic-peak gauge is available for a new
% comparison run, but cannot replace the gauge of an existing checkpoint.
settings.amplitudeGauge = 'outer_wall_density_window_l2';
settings.omegaGaugeWindowRadius = 0.5;
% Version 5 generates additional directional node levels automatically up to
% the resource cap chosen before the t=0 run. Version 4 remains available.
settings.autonomousMeshVersion = 5;
settings.maximumTotalNodes = 310000;

settings.outputDirectory = fullfile(projectRoot,'runs','profile_H8_long');
settings.storeSnapshots = false;
settings.makePlots = false;
settings.verbose = true;

% Leave empty for a new physical t=0 run. Set an absolute native checkpoint
% path to resume the same frozen numerical problem with longer horizons.
settings.restartCheckpoint = '';
end
