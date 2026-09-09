function [axis,info] = ipm_gridlab_axis(xlim,nodeCount,model,options)
%IPM_GRIDLAB_AXIS Equidistribute and admit one symmetric analytic grid.
%   The positive half monitor is integrated on a fine reference lattice,
%   inverted monotonically, and mirrored exactly.  By default the lattice
%   is uniform in physical x.  OPTIONS.referenceAxis instead interprets the
%   analytic density as a multiplier of an already admitted grid map.  A
%   constant multiplier then reproduces that reference grid while a local
%   feature only reallocates its existing computational cells.  This avoids
%   asking a small analytic contrast to recreate a several-decade global
%   spacing ratio.  INFO.admissible is false whenever a mesh-smoothness or
%   high-order conditioning gate fails.

if nargin < 4 || isempty(options)
    options = struct();
end
options = resolve_options(options,nodeCount);
validateattributes(xlim,{'numeric'}, ...
    {'vector','numel',2,'real','finite','increasing'},mfilename,'xlim');
xlim = double(xlim(:)');
validateattributes(nodeCount,{'numeric'}, ...
    {'scalar','integer','>=',9},mfilename,'nodeCount');
if mod(nodeCount,2) ~= 1 || abs(sum(xlim)) > ...
        100*eps(max(abs(xlim)))
    error('ipm:gridlab:AxisSymmetry', ...
        'The analytic grid requires odd nodeCount and symmetric xlim.');
end

halfWidth = xlim(2);
[referenceCoordinate,physicalReference,baselinePositive,mappingKind] = ...
    integration_lattice(xlim,nodeCount,options);
density = ipm_gridlab_density(physicalReference,model);
mass = cumtrapz(referenceCoordinate,density);
mass = mass/mass(end);
targets = linspace(0,1,(nodeCount+1)/2);
positiveAxis = interp1(mass,physicalReference,targets,'pchip');
positiveAxis(1) = 0;
positiveAxis(end) = halfWidth;
[positiveAxis,projection] = project_to_admissible( ...
    positiveAxis,baselinePositive,options);
axis = [-fliplr(positiveAxis(2:end)),positiveAxis];
axis((nodeCount+1)/2) = 0;

if any(~isfinite(axis)) || any(diff(axis) <= 0)
    error('ipm:gridlab:AxisMonotonicity', ...
        'Monitor inversion did not produce a strictly increasing axis.');
end
[weights,quality,reasons] = axis_diagnostics(axis,options);

cellCenters = 0.5*(positiveAxis(1:end-1)+positiveAxis(2:end));
cellDensity = ipm_gridlab_density(cellCenters,model);
info = struct('admissible',isempty(reasons),'rejectionReasons',{reasons}, ...
    'quality',quality,'quadratureWeights',weights, ...
    'positiveAxis',positiveAxis, ...
    'mappingKind',mappingKind, ...
    'projection',projection, ...
    'monitorMinimum',min(density),'monitorMaximum',max(density), ...
    'monitorContrast',max(density)/min(density), ...
    'cellMonitorDensity',cellDensity,'fineCount',options.fineCount, ...
    'limits',rmfield(options,{'fineCount','referenceAxis'}));
end

function [coordinate,physical,baselinePositive,mappingKind] = ...
        integration_lattice(xlim,nodeCount,options)
if isempty(options.referenceAxis)
    coordinate = linspace(0,xlim(2),options.fineCount);
    physical = coordinate;
    baselinePositive = linspace(0,xlim(2),(nodeCount+1)/2);
    mappingKind = 'uniform_physical_coordinate';
    return
end
referenceAxis = double(options.referenceAxis(:)');
if mod(numel(referenceAxis),2) ~= 1 || any(diff(referenceAxis) <= 0) || ...
        max(abs(referenceAxis+fliplr(referenceAxis))) > ...
        1e3*eps(max(1,max(abs(referenceAxis))))
    error('ipm:gridlab:AxisReference', ...
        'options.referenceAxis must be odd, increasing, and symmetric.');
end
domainTolerance = 1e3*eps(max(1,max(abs(xlim))));
if max(abs(referenceAxis([1,end])-xlim)) > domainTolerance
    error('ipm:gridlab:AxisReference', ...
        'options.referenceAxis must have the requested domain endpoints.');
end
referencePositive = referenceAxis((numel(referenceAxis)+1)/2:end);
coarseCoordinate = linspace(0,1,numel(referencePositive));
targetCoordinate = linspace(0,1,(nodeCount+1)/2);
baselinePositive = interp1( ...
    coarseCoordinate,referencePositive,targetCoordinate,'pchip');
baselinePositive([1,end]) = [0,xlim(2)];
coordinate = linspace(0,1,options.fineCount);
physical = interp1(coarseCoordinate,referencePositive,coordinate,'pchip');
physical([1,end]) = [0,xlim(2)];
if any(diff(physical) <= 0)
    error('ipm:gridlab:AxisReference', ...
        'The fine reference map must remain strictly increasing.');
end
if numel(referenceAxis) == nodeCount
    mappingKind = 'same_resolution_reference_multiplier';
else
    mappingKind = 'resampled_reference_multiplier';
end
end

function [positiveAxis,projection] = project_to_admissible( ...
        requestedPositive,baselinePositive,options)
requestedAxis = mirrored_axis(requestedPositive);
[~,requestedQuality,requestedReasons] = ...
    axis_diagnostics(requestedAxis,options);
projection = struct('applied',false,'alpha',1, ...
    'requestedQuality',requestedQuality, ...
    'requestedRejectionReasons',{requestedReasons});
if isempty(requestedReasons) || ~options.projectToAdmissible
    positiveAxis = requestedPositive;
    return
end

baselineAxis = mirrored_axis(baselinePositive);
[~,~,baselineReasons] = axis_diagnostics(baselineAxis,options);
if ~isempty(baselineReasons)
    positiveAxis = requestedPositive;
    projection.baselineRejectionReasons = baselineReasons;
    return
end

baselineLogWidths = log(diff(baselinePositive));
requestedLogWidths = log(diff(requestedPositive));
sampleCount = 65;
sampleAlpha = linspace(1,0,sampleCount);
validAlpha = 0;
invalidUpper = 1;
for index = 1:sampleCount
    alpha = sampleAlpha(index);
    trialPositive = blended_positive_axis( ...
        baselineLogWidths,requestedLogWidths,alpha,requestedPositive(end));
    [~,~,trialReasons] = axis_diagnostics( ...
        mirrored_axis(trialPositive),options);
    if isempty(trialReasons)
        validAlpha = alpha;
        if index > 1
            invalidUpper = sampleAlpha(index-1);
        else
            invalidUpper = alpha;
        end
        break
    end
end
for iteration = 1:40
    alpha = 0.5*(validAlpha+invalidUpper);
    trialPositive = blended_positive_axis( ...
        baselineLogWidths,requestedLogWidths,alpha,requestedPositive(end));
    [~,~,trialReasons] = axis_diagnostics( ...
        mirrored_axis(trialPositive),options);
    if isempty(trialReasons)
        validAlpha = alpha;
    else
        invalidUpper = alpha;
    end
end
positiveAxis = blended_positive_axis( ...
    baselineLogWidths,requestedLogWidths,validAlpha,requestedPositive(end));
projection.applied = true;
projection.alpha = validAlpha;
projection.baselineRejectionReasons = baselineReasons;
end

function positiveAxis = blended_positive_axis( ...
        baselineLogWidths,requestedLogWidths,alpha,halfWidth)
logWidths = baselineLogWidths+alpha*( ...
    requestedLogWidths-baselineLogWidths);
widths = exp(logWidths-max(logWidths));
widths = halfWidth*widths/sum(widths);
positiveAxis = [0,cumsum(widths)];
positiveAxis(end) = halfWidth;
end

function axis = mirrored_axis(positiveAxis)
axis = [-fliplr(positiveAxis(2:end)),positiveAxis];
axis((numel(axis)+1)/2) = 0;
end

function [weights,quality,reasons] = axis_diagnostics(axis,options)
weights = ipm.mesh.quadrature(axis);
quality = ipm.mesh.quality(axis,7,weights);
reasons = {};
if quality.maximumAdjacentCellRatio > options.maxAdjacentCellRatio
    reasons{end+1} = 'adjacent_cell_ratio';
end
if quality.maximumLogSpacingCurvature > options.maxLogSpacingCurvature
    reasons{end+1} = 'log_spacing_curvature';
end
if quality.minimumStencilRcond < options.minStencilRcond
    reasons{end+1} = 'stencil_conditioning';
end
if quality.minimumQuadratureWeightRatio < ...
        options.minQuadratureWeightRatio || any(weights <= 0)
    reasons{end+1} = 'quadrature_margin';
end
end

function options = resolve_options(options,nodeCount)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:AxisOptions', ...
        'options must be a scalar structure.');
end
defaults = struct( ...
    'fineCount',max(131073,512*(nodeCount-1)+1), ...
    'referenceAxis',[], ...
    'projectToAdmissible',true, ...
    'maxAdjacentCellRatio',1.08, ...
    'maxLogSpacingCurvature',0.03, ...
    'minStencilRcond',1e-9, ...
    'minQuadratureWeightRatio',1e-4);
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:AxisOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
end
validateattributes(options.fineCount,{'numeric'}, ...
    {'scalar','integer','>=',4097},mfilename,'options.fineCount');
if ~isempty(options.referenceAxis)
    validateattributes(options.referenceAxis,{'numeric'}, ...
        {'vector','real','finite','increasing'},mfilename, ...
        'options.referenceAxis');
elseif ~isnumeric(options.referenceAxis)
    error('ipm:gridlab:AxisOptions', ...
        'options.referenceAxis must be a numeric vector or empty.');
end
if ~(islogical(options.projectToAdmissible) && ...
        isscalar(options.projectToAdmissible))
    error('ipm:gridlab:AxisOptions', ...
        'options.projectToAdmissible must be a logical scalar.');
end
for name = {'maxAdjacentCellRatio','maxLogSpacingCurvature', ...
        'minStencilRcond','minQuadratureWeightRatio'}
    validateattributes(options.(name{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename, ...
        ['options.' name{1}]);
end
end
