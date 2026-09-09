function opts = bessel_two_scale_case(name,userOpts)
%BESSEL_TWO_SCALE_CASE Reproducible smooth Bessel-tail IPM experiments.
%   The simulation entry point remains IPM.SOLVE.  For example,
%
%     addpath('research/experiments','research/analysis');
%     opts = bessel_two_scale_case('screen_perturbed');
%     result = ipm.solve(opts);
%     diagnostics = ipm_analyze_bessel_two_scale(result);
%
%   These cases use the existing impermeable upper-half-plane solver.  The
%   initial density is the restriction of a full-plane Bessel steady density
%   with its pole below the wall, but the solver recomputes the streamfunction
%   with psi=0 on y=0.  Consequently this is a dynamical-selection test, not
%   an exact stationary regression of the Bessel density/velocity pair.
%   The k1_regularized initializer has an exactly zero wall trace. Physical
%   mode does not require a scaling-gauge peak, so no artificial wall monitor
%   is added by the maintained cases.

if nargin < 1 || isempty(name)
    name = 'screen_perturbed';
end
if nargin < 2 || isempty(userOpts)
    userOpts = struct();
end

common = struct( ...
    'symmetryMode','half_plane', ...
    'rescalingMode','physical', ...
    'farBoundaryMode','green', ...
    'transportBoundaryMode','open', ...
    'transportScheme','weno5_nonuniform', ...
    'nx',257,'ny',129,'xlim',[-24,24],'ymax',12, ...
    'targetCenterSpacing',[0.075,0.075], ...
    'physicalFinalTime',0.05,'maxDt',2.5e-4,'maxSteps',20000, ...
    'outputEvery',0.01,'cfl',0.3, ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'adaptiveLevels',[0.05,0.3,0.6], ...
    'targetLevelPoints',[32,16,8], ...
    'minimumLevelPoints',[8,4,2], ...
    'rangeStopTolerance',2e-2, ...
    'positiveWallNegativeTolerance',0.95, ...
    'oscillationTVTolerance',1.2, ...
    'storeSnapshots',true,'makePlots',false,'livePlot',false, ...
    'writeVideo',false,'saveResults',true,'verbose',true);
metadata = default_metadata(name);

switch lower(char(name))
    case 'smoke'
        opts = common;
        opts.nx = 129;
        opts.ny = 65;
        opts.targetCenterSpacing = [0.12,0.12];
        opts.physicalFinalTime = 0.01;
        opts.maxDt = 5e-4;
        opts.outputEvery = 0.005;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_two_scale_smoke.mat');
    case 'screen_unperturbed'
        opts = common;
        metadata.parameters.perturbationAmplitude = 0;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_two_scale_screen_unperturbed.mat');
    case 'screen_perturbed'
        opts = common;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_two_scale_screen_perturbed.mat');
    case 'refined_perturbed'
        opts = common;
        opts.nx = 513;
        opts.ny = 257;
        opts.targetCenterSpacing = [0.04,0.04];
        opts.maxDt = 1.25e-4;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_two_scale_refined_perturbed.mat');
    case 'long_perturbed'
        opts = common;
        opts.physicalFinalTime = 0.25;
        opts.outputEvery = 0.025;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_two_scale_long_perturbed.mat');
    case 'smoke_k1'
        opts = common;
        metadata.family = 'k1_regularized';
        opts.nx = 129;
        opts.ny = 65;
        opts.targetCenterSpacing = [0.12,0.12];
        opts.physicalFinalTime = 0.01;
        opts.maxDt = 5e-4;
        opts.outputEvery = 0.005;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_k1_two_scale_smoke.mat');
    case 'screen_k1_perturbed'
        opts = common;
        metadata.family = 'k1_regularized';
        opts.resultFile = fullfile('result','verification', ...
            'bessel_k1_two_scale_screen_perturbed.mat');
    case 'refined_k1_perturbed'
        opts = common;
        metadata.family = 'k1_regularized';
        opts.nx = 513;
        opts.ny = 257;
        opts.targetCenterSpacing = [0.04,0.04];
        opts.maxDt = 1.25e-4;
        opts.resultFile = fullfile('result','verification', ...
            'bessel_k1_two_scale_refined_perturbed.mat');
    otherwise
        error('ipm:UnknownBesselTwoScaleCase', ...
            ['Use ''smoke'', ''screen_unperturbed'', ' ...
            '''screen_perturbed'', ''refined_perturbed'', or ' ...
            '''long_perturbed'', ''smoke_k1'', ' ...
            '''screen_k1_perturbed'', or ''refined_k1_perturbed''.']);
end

[userOpts,metadata] = extract_metadata_overrides(userOpts,metadata);
if strcmp(metadata.family,'k1_regularized') && ...
        ~isfield(userOpts,'positiveWallNegativeTolerance')
    % The homogeneous-wall K1 family has no prescribed positive wall-omega
    % branch; that sign-specific stop criterion is therefore inapplicable.
    opts.positiveWallNegativeTolerance = Inf;
end
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
opts.caseMetadata = metadata;

parameters = metadata.parameters;
if strcmp(metadata.family,'k1_regularized')
    opts.initialCondition = @(x,y)bessel_k1_two_scale_initial_data( ...
        x,y,parameters);
elseif strcmp(metadata.family,'k0_external')
    opts.initialCondition = @(x,y)bessel_two_scale_initial_data(x,y,parameters);
else
    error('ipm:UnknownBesselTwoScaleFamily', ...
        'caseMetadata.family must be k0_external or k1_regularized.');
end
end

function metadata = default_metadata(name)
parameters = struct('mu',2,'poleDepth',0.6,'coreRadius',0.45, ...
    'centerX',0,'amplitude',1,'perturbationAmplitude',0.03, ...
    'perturbationCenterX',5,'perturbationScaleX',2, ...
    'perturbationScaleY',1.5,'monitorAmplitude',0, ...
    'monitorScaleX',0.8,'monitorScaleY',2);
metadata = struct('name',lower(char(name)),'family','k0_external', ...
    'parameters',parameters);
end

function [solverOpts,metadata] = extract_metadata_overrides(userOpts,metadata)
solverOpts = userOpts;
if isfield(solverOpts,'caseMetadata')
    metadata = merge_metadata(metadata,solverOpts.caseMetadata);
    solverOpts = rmfield(solverOpts,'caseMetadata');
end
metadata = validate_metadata(metadata);
end

function metadata = merge_metadata(metadata,override)
if ~isstruct(override) || ~isscalar(override)
    error('ipm:BesselCaseMetadata', ...
        'caseMetadata must be a scalar structure.');
end
unknown = setdiff(fieldnames(override),{'name';'family';'parameters'});
if ~isempty(unknown)
    error('ipm:BesselCaseMetadataField', ...
        'Unknown caseMetadata field(s): %s.',strjoin(unknown,', '));
end
if isfield(override,'name')
    metadata.name = override.name;
end
if isfield(override,'family')
    metadata.family = override.family;
end
if isfield(override,'parameters')
    parameters = override.parameters;
    if ~isstruct(parameters) || ~isscalar(parameters)
        error('ipm:BesselCaseParameters', ...
            'caseMetadata.parameters must be a scalar structure.');
    end
    unknown = setdiff(fieldnames(parameters),fieldnames(metadata.parameters));
    if ~isempty(unknown)
        error('ipm:BesselCaseParameter', ...
            'Unknown Bessel parameter(s): %s.',strjoin(unknown,', '));
    end
    names = fieldnames(parameters);
    for index = 1:numel(names)
        metadata.parameters.(names{index}) = parameters.(names{index});
    end
end
end

function metadata = validate_metadata(metadata)
if ~(ischar(metadata.name) || ...
        (isstring(metadata.name) && isscalar(metadata.name)))
    error('ipm:BesselCaseName','caseMetadata.name must be text.');
end
if ~(ischar(metadata.family) || ...
        (isstring(metadata.family) && isscalar(metadata.family)))
    error('ipm:BesselCaseFamily','caseMetadata.family must be text.');
end
metadata.name = char(metadata.name);
metadata.family = lower(char(metadata.family));
if ~ismember(metadata.family,{'k0_external','k1_regularized'})
    error('ipm:UnknownBesselTwoScaleFamily', ...
        'caseMetadata.family must be k0_external or k1_regularized.');
end
end
