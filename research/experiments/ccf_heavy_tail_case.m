function opts = ccf_heavy_tail_case(name)
%CCF_HEAVY_TAIL_CASE Reproducible tests of a regularized CCF-like IPM tail.
%   The only simulation entry point remains ipm.solve:
%       addpath('research/experiments');
%       result = ipm.solve(ccf_heavy_tail_case('screen_physical'));

if nargin < 1
    name = 'screen_physical';
end

common = struct('initialCondition','ccf_heavy_tail', ...
    'ccfTailOnset',1,'ccfTailRegularization',0.08, ...
    'ccfTailAmplitude',1,'ccfTailCutoffScale',16, ...
    'ccfTailVerticalScale',2, ...
    'ccfTailVerticalPower',4,'nx',513,'ny',257, ...
    'xlim',[-64,64],'ymax',32,'targetCenterSpacing',[0.01,0.01], ...
    'physicalFinalTime',0.25,'maxDt',1e-4, ...
    'initialAnalyticRemesh',true,'remeshSafetyTrigger',0.65, ...
    'remeshTargetSafety',0.35, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'saveResults',true,'verbose',true);

switch lower(char(name))
    case 'screen_physical'
        opts = common;
        opts.rescalingMode = 'physical';
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_screen_physical.mat');
    case 'screen_dynamic'
        opts = common;
        opts.rescalingMode = 'dynamic';
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_screen_dynamic.mat');
    case 'box_large'
        opts = common;
        opts.rescalingMode = 'physical';
        opts.nx = 1025;
        opts.ny = 513;
        opts.xlim = [-128,128];
        opts.ymax = 64;
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_box_large.mat');
    case 'epsilon_fine'
        opts = common;
        opts.rescalingMode = 'physical';
        opts.nx = 1025;
        opts.ny = 513;
        opts.ccfTailRegularization = 0.04;
        opts.targetCenterSpacing = [0.005,0.005];
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_epsilon_fine.mat');
    case 'long_physical'
        opts = long_mesh_options(common);
        opts.rescalingMode = 'physical';
        opts.physicalFinalTime = 0.5;
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_long_physical.mat');
    case 'long_dynamic'
        % Freeze the historical CCF comparison explicitly.  The maintained
        % no-argument ipm.solve datum may change independently of this named
        % experiment.
        opts = long_mesh_options(common);
        opts.xlim = [-128,128];
        opts.ymax = 64;
        opts.rescalingMode = 'dynamic';
        opts.physicalFinalTime = 0.5;
        opts.transportScheme = 'weno5_nonuniform';
        opts.lengthGauge = 'transport_anchor';
        opts.transportAnchorX = 1;
        opts.cOmegaGauge = 'gradient_energy';
        opts.makePlots = false;
        opts.livePlot = false;
        opts.writeVideo = false;
        opts.saveResults = true;
        opts.verbose = true;
        opts.resultFile = fullfile('result','verification', ...
            'ccf_tail_long_dynamic.mat');
    otherwise
        error('ipm:UnknownCcfHeavyTailCase', ...
            ['Use ''screen_physical'', ''screen_dynamic'', ' ...
            '''box_large'', ''epsilon_fine'', ''long_physical'', ' ...
            'or ''long_dynamic''.']);
end
end

function opts = long_mesh_options(common)
opts = common;
opts.nx = 1025;
opts.ny = 513;
opts.xlim = [-128,128];
opts.ymax = 64;
opts.targetPeakPoints = 128;
opts.targetLevelPoints = [192,128,96];
opts.targetPeakExpansion = 8;
opts.minimumPeakPoints = 16;
opts.minimumLevelPoints = [24,14,10];
opts.remeshSafetyTrigger = 0.15;
opts.remeshTargetSafety = 0.08;
opts.remeshOuterAnchorFactor = 3;
opts.remeshMaximumCellRatio = 1.5;
opts.remeshReferenceMode = 'current';
opts.initialAnalyticRemeshPasses = 12;
opts.remeshAxisProfile = 'lattice';
opts.remeshLatticeTransitionFraction = 0.125;
opts.remeshBridgeFloor = 0.005;
opts.remeshBridgePower = 4;
end
