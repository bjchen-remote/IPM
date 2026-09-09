function opts = fourth_order_blowup_case(name,overrides)
%FOURTH_ORDER_BLOWUP_CASE Fourth-order studies of the maintained active datum.
%   OPTS = FOURTH_ORDER_BLOWUP_CASE(NAME) returns explicit options for
%   ipm.solve.  Case names ending in QUADRANT_256 contain 256 cells in each
%   first-quadrant direction; denser named cases state their first-quadrant
%   cell count explicitly.  They retain the active datum, dynamic gauge,
%   Green artificial boundary, and open outer transport while replacing the
%   numerical tuple by high_order/WENO5-FD/SSPRK(5,4).
%
%   These dynamic/Green/open calculations are research runs.  They do not
%   inherit the fixed-grid physical fourth-order convergence certificate.

if nargin < 1 || isempty(name)
    name = 'pilot_256';
end
if nargin < 2
    overrides = struct();
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('ipm:BadFourthOrderBlowupOverrides', ...
        'overrides must be a scalar structure.');
end

opts = ipm.config.activeCase();
opts.nx = 257;
opts.ny = 257;
opts.spatialDiscretization = 'high_order';
opts.transportScheme = 'weno5_fd';
opts.timeIntegrator = 'ssprk54';
opts.remeshTransferScheme = 'high_order';
opts.cfl = 0.25;
opts.initialAnalyticRemesh = false;
opts.adaptiveRemesh = false;
opts.storeSnapshots = false;
opts.saveResults = true;
opts.makePlots = false;
opts.livePlot = false;
opts.writeVideo = false;
opts.verbose = true;

metadata = struct( ...
    'campaign','fourth_order_active_blowup_20260902', ...
    'gridInterpretation','257 nodes = 256 stored cells per direction', ...
    'methodScope',[ ...
        'Fourth-order components; dynamic/Green/open blow-up evolution ' ...
        'remains outside the global fixed-grid certificate.']);
smoothGridStrength = NaN;
smoothGridWidth = NaN;
smoothGridParameters = struct();
smoothGridQualityPolicy = [];
% Archived named cases used the original abs(x) monitor.  Keep that choice
% explicit so their saved results remain reproducible; new analytic cases
% opt into a genuinely smooth double-Gaussian monitor below.
smoothGridProfile = 'legacy_abs_gaussian';

switch lower(char(name))
    case 'pilot_256'
        opts.physicalFinalTime = 5e-3;
        opts.maxDt = 1e-4;
        opts.outputEvery = 1e-4;
        caseName = 'pilot_256';
    case 'pilot_256_dt_half'
        opts.physicalFinalTime = 5e-3;
        opts.cfl = 0.125;
        opts.maxDt = 5e-5;
        opts.outputEvery = 1e-4;
        caseName = 'pilot_256_dt_half';
    case 'fixed_256'
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 5e-4;
        opts.outputEvery = 2e-3;
        caseName = 'fixed_256';
    case 'premeshed_256'
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 5e-4;
        opts.outputEvery = 2e-3;
        opts.initialAnalyticRemesh = true;
        opts.initialAnalyticRemeshPasses = 2;
        opts.adaptiveRemesh = true;
        opts.maxRemeshes = 2;
        opts.remeshAxisProfile = 'lattice';
        opts.remeshOuterAnchorFactor = 2;
        opts.remeshReferenceMode = 'current';
        caseName = 'premeshed_256';
    case 'premeshed_quadrant_256'
        % The solver stores both horizontal symmetry halves.  Therefore
        % 513-by-257 nodes give 256 cells in x>0 and 256 cells in y>0.
        opts.nx = 513;
        opts.ny = 257;
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 5e-4;
        opts.outputEvery = 2e-3;
        opts.initialAnalyticRemesh = true;
        opts.initialAnalyticRemeshPasses = 2;
        opts.adaptiveRemesh = true;
        opts.maxRemeshes = 2;
        opts.remeshAxisProfile = 'lattice';
        opts.remeshOuterAnchorFactor = 2;
        opts.remeshReferenceMode = 'current';
        caseName = 'premeshed_quadrant_256';
    case 'smooth_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 5e-4;
        opts.outputEvery = 2e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 2;
        smoothGridWidth = 0.85;
        caseName = 'smooth_quadrant_256';
    case 'smooth_box_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.02];
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 5e-4;
        opts.outputEvery = 2e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 0.5;
        smoothGridWidth = 0.7;
        caseName = 'smooth_box_quadrant_256';
    case 'focused_box_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.02];
        opts.physicalFinalTime = 1;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 0.7;
        caseName = 'focused_box_quadrant_256';
    case 'balanced_box_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.01];
        opts.physicalFinalTime = 1.2;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 2;
        smoothGridWidth = 0.7;
        caseName = 'balanced_box_quadrant_256';
    case 'concentrated_box_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.005];
        opts.physicalFinalTime = 1.5;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 3;
        smoothGridWidth = 0.7;
        caseName = 'concentrated_box_quadrant_256';
    case 'near_terminal_box_quadrant_256'
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.85;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 12;
        smoothGridWidth = 0.7;
        caseName = 'near_terminal_box_quadrant_256';
    case 'stable_limit_box_quadrant_256'
        % Strongest same-width monitor retained after short-time screens.
        % A=24 and the narrower A=28 limit grid excite early WENO peaks.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.7;
        caseName = 'stable_limit_box_quadrant_256';
    case 'narrow_limit_box_quadrant_256'
        % Same stable strength with a narrower legacy peak monitor.
        % Short screens show improved core resolution and far-field ratios.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.35;
        caseName = 'narrow_limit_box_quadrant_256';
    case 'moderate_narrow_limit_box_quadrant_256'
        % Intermediate-width monitor retained for a long-time stability test.
        % Widths 0.35 and 0.45 excite wall oscillations before the broad case.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.55;
        caseName = 'moderate_narrow_limit_box_quadrant_256';
    case 'analytic_stable_limit_box_quadrant_256'
        % Smooth-monitor replacement for the strongest stable legacy map.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.7;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_stable_limit_box_quadrant_256';
    case 'analytic_intermediate_limit_box_quadrant_256'
        % Brackets the stability boundary between strengths 16 and 20.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 18;
        smoothGridWidth = 0.7;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_intermediate_limit_box_quadrant_256';
    case 'analytic_moderate_narrow_limit_box_quadrant_256'
        % Width sensitivity check at the strongest stable broad strength.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.55;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_moderate_narrow_limit_box_quadrant_256';
    case 'analytic_moderate_narrow_quadrant_384'
        % Matched 1.5x spatial refinement of the retained A=16, w=0.55
        % analytic monitor.  The stored domain includes both x halves, so
        % 769-by-385 nodes represent 384-by-384 first-quadrant cells.
        opts.nx = 769;
        opts.ny = 385;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0193174285333333, ...
            0.000833333333333333];
        opts.physicalFinalTime = 2.05;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 16;
        smoothGridWidth = 0.55;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_moderate_narrow_quadrant_384';
    case 'analytic_lorentzian_a24_w020_c114_quadrant_256'
        % Broad algebraic shoulders cover the measured wall-peak path in
        % solver/grid x=0.97--1.28, avoiding a sharp Gaussian lattice exit.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 24;
        smoothGridWidth = 0.20;
        smoothGridProfile = 'analytic_double_lorentzian';
        smoothGridParameters = struct('peakCenter',1.14);
        caseName = 'analytic_lorentzian_A24_w020_c114_quadrant_256';
    case 'analytic_gaussian_a14_w042_c114_quadrant_256'
        % Low-amplitude analytic control centered on the observed path.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 14;
        smoothGridWidth = 0.42;
        smoothGridProfile = 'analytic_double_gaussian';
        smoothGridParameters = struct('peakCenter',1.14);
        caseName = 'analytic_gaussian_A14_w042_c114_quadrant_256';
    case 'analytic_lorentzian_a28_w020_c114_quadrant_256'
        % Middle member of the smooth-shoulder resolution/cost screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 28;
        smoothGridWidth = 0.20;
        smoothGridProfile = 'analytic_double_lorentzian';
        smoothGridParameters = struct('peakCenter',1.14);
        caseName = 'analytic_lorentzian_A28_w020_c114_quadrant_256';
    case 'analytic_lorentzian_a32_w020_c114_quadrant_256'
        % Finest member admitted to the dynamic smooth-shoulder screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 32;
        smoothGridWidth = 0.20;
        smoothGridProfile = 'analytic_double_lorentzian';
        smoothGridParameters = struct('peakCenter',1.14);
        caseName = 'analytic_lorentzian_A32_w020_c114_quadrant_256';
    case 'analytic_lorentzian_a28_w020_c114_quadrant_512'
        % Strict two-times refinement of the q256 A=28 smooth-shoulder map.
        % maxSteps is a normal-save guard for the overnight wall-clock budget.
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0144880714,0.000625];
        opts.physicalFinalTime = 2.10;
        opts.maxDt = 2e-3;
        opts.maxSteps = 4800;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 28;
        smoothGridWidth = 0.20;
        smoothGridProfile = 'analytic_double_lorentzian';
        smoothGridParameters = struct('peakCenter',1.14);
        caseName = 'analytic_lorentzian_A28_w020_c114_quadrant_512';
    case 'coordinated_log_f80_r16_c114_quadrant_256'
        % Ratio-first rounded-log grid.  A uniform 80-cell platform covers
        % x=0.95--1.30; both outer density transitions are solved so their
        % asymptotic adjacent-density ratios are nearly equal.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(16/256);
        caseName = 'coordinated_log_F80_R16_c114_quadrant_256';
    case 'coordinated_log_f80_r20_c114_quadrant_256'
        % Intermediate ratio/curvature member of the coordinated screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256);
        caseName = 'coordinated_log_F80_R20_c114_quadrant_256';
    case 'coordinated_log_f80_r24_c114_quadrant_256'
        % Curvature-favouring member of the coordinated q256 screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(24/256);
        caseName = 'coordinated_log_F80_R24_c114_quadrant_256';
    case 'coordinated_log_f80_r28_c114_quadrant_256'
        % Smoothest-curvature member admitted to the coordinated q256 screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(28/256);
        caseName = 'coordinated_log_F80_R28_c114_quadrant_256';
    case 'coordinated_log_f90_r20_c114_quadrant_256'
        % Density-raised candidate after the rounded-log q256 R screen.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256,90/256);
        caseName = 'coordinated_log_F90_R20_c114_quadrant_256';
    case 'coordinated_log_f92_r20_c114_quadrant_256'
        % Highest platform density admitted by the q256 ratio/curvature gate.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256,92/256);
        caseName = 'coordinated_log_F92_R20_c114_quadrant_256';
    case 'coordinated_log_f80_r20_c114_quadrant_512'
        % Strict two-times refinement of the registered F80/R20 q256 grid.
        % This is the only coordinated q512 case that is spatially matched
        % to the schema-4 gauge tournament's common grid.
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0144880714,0.000625];
        opts.physicalFinalTime = 2.10;
        opts.maxDt = 2e-3;
        opts.maxSteps = 5300;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256);
        caseName = 'coordinated_log_F80_R20_c114_quadrant_512';
    case 'coordinated_log_anchor_f80_r20_large_box_quadrant_512'
        % Fresh-root large-box grid for schema-4 box-refinement studies.
        % The fine interval is shifted by 0.0025 so X=1 is an exact node;
        % its 160-cell platform retains the q512 peak resolution while the
        % rounded logarithmic tails reach 1e6 without exceeding ratio 1.08.
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-1e6,1e6];
        opts.ymax = 1e6;
        opts.targetCenterSpacing = [0.0144880714,0.000625];
        opts.physicalFinalTime = 2.10;
        opts.maxDt = 2e-3;
        opts.maxSteps = 5300;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        opts.remeshMaximumCellRatio = 1.08;
        opts.rescalingMode = 'dynamic';
        opts.dynamicScaleGeometry = 'isotropic';
        opts.lengthGauge = 'transport_anchor';
        opts.transportAnchorX = 1;
        opts.cOmegaGauge = 'gradient_energy';
        opts.lengthScaleGain = 1;
        opts.widthGaugeGain = 0;
        opts.omegaGaugeGain = 0;
        opts.travelingWaveGain = 0;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = struct( ...
            'fineInterval',[0.9475,1.2975], ...
            'fineCellFraction',80/256, ...
            'roundingFraction',20/256,'peakCenter',1.14);
        smoothGridQualityPolicy = struct( ...
            'minimumStencilRcond',1e-9, ...
            'minimumQuadratureWeightRatio',1e-8, ...
            'minimumLocalQuadratureWeightRatio',0.35, ...
            'maximumLocalQuadratureWeightRatio',1.65, ...
            'quadratureWarningRatio',1e-4);
        caseName = ...
            'coordinated_log_anchor_F80_R20_large_box_quadrant_512';
    case 'coordinated_log_f90_r20_c114_quadrant_512'
        % Two-times fractional refinement of the F90/R20 coordinated map.
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0144880714,0.000625];
        opts.physicalFinalTime = 2.10;
        opts.maxDt = 2e-3;
        opts.maxSteps = 4800;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256,90/256);
        caseName = 'coordinated_log_F90_R20_c114_quadrant_512';
    case 'coordinated_log_f92_r20_c114_quadrant_512'
        % Two-times fractional refinement of the F92/R20 coordinated map.
        % The F92 short screen takes about 7.5% more steps than A28; 5300
        % avoids reproducing the earlier q512 campaign's premature 4800-step
        % budget stop near the target late-time window.
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0144880714,0.000625];
        opts.physicalFinalTime = 2.10;
        opts.maxDt = 2e-3;
        opts.maxSteps = 5300;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 1;
        smoothGridWidth = 1;
        smoothGridProfile = 'coordinated_log_spacing';
        smoothGridParameters = coordinated_grid_parameters(20/256,92/256);
        caseName = 'coordinated_log_F92_R20_c114_quadrant_512';
    case 'analytic_strong_limit_box_quadrant_256'
        % Strong analytic monitor admitted for stability screening.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 20;
        smoothGridWidth = 0.7;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_strong_limit_box_quadrant_256';
    case 'analytic_ratio_limit_box_quadrant_256'
        % Near-ratio-limit analytic monitor used only after short screens.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 24;
        smoothGridWidth = 0.7;
        smoothGridProfile = 'analytic_double_gaussian';
        caseName = 'analytic_ratio_limit_box_quadrant_256';
    case 'limit_box_quadrant_256'
        % Last fixed-size mapping admitted by the 1.5 adjacent-cell gate.
        % The narrower monitor balances the initial x/y 90%-core counts.
        opts.nx = 513;
        opts.ny = 257;
        opts.xlim = [-281,281];
        opts.ymax = 132;
        opts.targetCenterSpacing = [0.0289761428,0.00125];
        opts.physicalFinalTime = 1.95;
        opts.maxDt = 2e-3;
        opts.outputEvery = 5e-3;
        opts.initialAnalyticRemesh = false;
        opts.adaptiveRemesh = false;
        smoothGridStrength = 28;
        smoothGridWidth = 0.35;
        caseName = 'limit_box_quadrant_256';
    case 'adaptive_256'
        opts.physicalFinalTime = 0.5;
        opts.maxDt = 2.5e-4;
        opts.outputEvery = 2e-3;
        opts.initialAnalyticRemesh = true;
        opts.adaptiveRemesh = true;
        caseName = 'adaptive_256';
    otherwise
        error('ipm:UnknownFourthOrderBlowupCase', ...
            ['Use ''pilot_256'', ''pilot_256_dt_half'', ' ...
            '''fixed_256'', ''premeshed_256'', ' ...
            '''premeshed_quadrant_256'', ''smooth_quadrant_256'', ' ...
            '''smooth_box_quadrant_256'', ' ...
            '''focused_box_quadrant_256'', ' ...
            '''balanced_box_quadrant_256'', ' ...
            '''concentrated_box_quadrant_256'', ' ...
            '''near_terminal_box_quadrant_256'', ' ...
            '''stable_limit_box_quadrant_256'', ' ...
            '''narrow_limit_box_quadrant_256'', ' ...
            '''moderate_narrow_limit_box_quadrant_256'', ' ...
            '''analytic_stable_limit_box_quadrant_256'', ' ...
            '''analytic_intermediate_limit_box_quadrant_256'', ' ...
            '''analytic_moderate_narrow_limit_box_quadrant_256'', ' ...
            '''analytic_moderate_narrow_quadrant_384'', ' ...
            '''analytic_gaussian_A14_w042_c114_quadrant_256'', ' ...
            '''analytic_lorentzian_A24_w020_c114_quadrant_256'', ' ...
            '''analytic_lorentzian_A28_w020_c114_quadrant_256'', ' ...
            '''analytic_lorentzian_A32_w020_c114_quadrant_256'', ' ...
            '''analytic_lorentzian_A28_w020_c114_quadrant_512'', ' ...
            '''coordinated_log_F80_R16_c114_quadrant_256'', ' ...
            '''coordinated_log_F80_R20_c114_quadrant_256'', ' ...
            '''coordinated_log_F80_R24_c114_quadrant_256'', ' ...
            '''coordinated_log_F80_R28_c114_quadrant_256'', ' ...
            '''coordinated_log_F90_R20_c114_quadrant_256'', ' ...
            '''coordinated_log_F92_R20_c114_quadrant_256'', ' ...
            '''coordinated_log_F80_R20_c114_quadrant_512'', ' ...
            ['''coordinated_log_anchor_F80_R20_large_box_' ...
            'quadrant_512'', '] ...
            '''coordinated_log_F90_R20_c114_quadrant_512'', ' ...
            '''coordinated_log_F92_R20_c114_quadrant_512'', ' ...
            '''analytic_strong_limit_box_quadrant_256'', ' ...
            '''analytic_ratio_limit_box_quadrant_256'', ' ...
            '''limit_box_quadrant_256'', or ''adaptive_256''.']);
end

opts.resultFile = fullfile('result','verification', ...
    ['fourth_order_active_',caseName,'.mat']);
names = fieldnames(overrides);
for index = 1:numel(names)
    opts.(names{index}) = overrides.(names{index});
end

if isfinite(smoothGridStrength)
    if isfield(overrides,'customX') || isfield(overrides,'customY')
        error('ipm:FourthOrderSmoothGridOverride', ...
            ['Smooth-grid cases generate customX/customY from the final ' ...
            'geometry options; do not override either axis directly.']);
    end
    [opts.customX,opts.customY,gridInfo] = ...
        fourth_order_smooth_peak_grid( ...
            opts,smoothGridStrength,smoothGridWidth,smoothGridProfile, ...
        smoothGridParameters,smoothGridQualityPolicy);
    if strcmp(gridInfo.profile,'coordinated_log_spacing')
        targetSpacing = opts.targetCenterSpacing(:)';
        if isscalar(targetSpacing)
            targetSpacing = [targetSpacing,targetSpacing];
        end
        targetSpacing(1) = gridInfo.coordinatedDesign.originSpacing;
        opts.targetCenterSpacing = targetSpacing;
    end
    metadata.smoothGrid = gridInfo;
end

metadata.caseName = caseName;
if mod(opts.nx,2) == 1 && opts.ny >= 2 && ...
        opts.nx == 2*opts.ny-1
    metadata.gridInterpretation = ...
        sprintf(['%d x %d stored nodes = %d x %d cells in the ' ...
        'first quadrant'],opts.nx,opts.ny,(opts.nx-1)/2,opts.ny-1);
elseif opts.nx == 257 && opts.ny == 257
    metadata.gridInterpretation = ...
        '257 nodes = 256 stored cells per direction';
else
    metadata.gridInterpretation = sprintf( ...
        '%d x %d stored nodes; no first-quadrant cell-count claim', ...
        opts.nx,opts.ny);
end
if isfield(overrides,'caseMetadata')
    metadata = merge_metadata(metadata,overrides.caseMetadata);
end
opts.caseMetadata = metadata;
end

function parameters = coordinated_grid_parameters( ...
        roundingFraction,fineCellFraction)
if nargin < 2
    fineCellFraction = 80/256;
end
parameters = struct('fineCellFraction',fineCellFraction, ...
    'roundingFraction',roundingFraction, ...
    'peakCenter',1.14);
end

function merged = merge_metadata(base,extra)
if ~isstruct(extra) || ~isscalar(extra)
    error('ipm:BadFourthOrderBlowupMetadata', ...
        'caseMetadata override must be a scalar structure.');
end
reserved = {'campaign','caseName','gridInterpretation','methodScope', ...
    'smoothGrid'};
collisions = intersect(fieldnames(extra),reserved,'stable');
if ~isempty(collisions)
    error('ipm:FourthOrderReservedMetadata', ...
        'caseMetadata cannot override reserved field(s): %s.', ...
        strjoin(collisions,', '));
end
merged = base;
names = fieldnames(extra);
for index = 1:numel(names)
    merged.(names{index}) = extra.(names{index});
end
end
