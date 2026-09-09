function [axis,info] = ipm_gridlab_axis_positive(ylim,nodeCount,model,options)
%IPM_GRIDLAB_AXIS_POSITIVE Build and admit a one-sided analytic y grid.
%   This is the positive-half companion of IPM_GRIDLAB_AXIS.  The monitor
%   is mirrored only for equidistribution, then the retained [0,YMAX] axis
%   is diagnosed with its actual one-sided high-order stencils and
%   quadrature.  An inadmissible request is projected in log cell widths
%   toward the supplied one-sided reference grid.

if nargin < 4 || isempty(options)
    options = struct();
end
validateattributes(ylim,{'numeric'}, ...
    {'vector','numel',2,'real','finite','increasing'},mfilename,'ylim');
ylim = double(ylim(:)');
validateattributes(nodeCount,{'numeric'}, ...
    {'scalar','integer','>=',5},mfilename,'nodeCount');
if abs(ylim(1)) > 100*eps(max(1,max(abs(ylim)))) || ylim(2) <= 0
    error('ipm:gridlab:PositiveAxisDomain', ...
        'ylim must be a positive-half interval [0,YMAX].');
end
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:PositiveAxisOptions', ...
        'options must be a scalar structure.');
end

[referencePositive,referenceFull] = reference_axes( ...
    options,nodeCount,ylim(2));
axisOptions = options;
axisOptions.referenceAxis = referenceFull;
requestedProjection = true;
if isfield(options,'projectToAdmissible')
    requestedProjection = options.projectToAdmissible;
end
axisOptions.projectToAdmissible = false;
[fullAxis,fullInfo] = ipm_gridlab_axis( ...
    [-ylim(2),ylim(2)],2*nodeCount-1,model,axisOptions);
requested = fullAxis(nodeCount:end);
limits = fullInfo.limits;
[requestedWeights,requestedQuality,requestedReasons] = ...
    positive_diagnostics(requested,limits);

projection = struct('applied',false,'alpha',1, ...
    'requestedQuality',requestedQuality, ...
    'requestedRejectionReasons',{requestedReasons});
axis = requested;
if ~isempty(requestedReasons) && requestedProjection
    [~,~,baselineReasons] = positive_diagnostics( ...
        referencePositive,limits);
    projection.baselineRejectionReasons = baselineReasons;
    if isempty(baselineReasons)
        [axis,alpha] = project_axis( ...
            requested,referencePositive,limits,ylim(2));
        projection.applied = true;
        projection.alpha = alpha;
    end
end
[weights,quality,reasons] = positive_diagnostics(axis,limits);
info = struct('admissible',isempty(reasons), ...
    'rejectionReasons',{reasons},'quality',quality, ...
    'quadratureWeights',weights,'projection',projection, ...
    'fullAxisInfo',fullInfo,'requestedQuadratureWeights',requestedWeights, ...
    'limits',limits);
end

function [positive,full] = reference_axes(options,nodeCount,halfWidth)
if ~isfield(options,'referenceAxis') || isempty(options.referenceAxis)
    positive = linspace(0,halfWidth,nodeCount);
    full = [];
    return
end
reference = double(options.referenceAxis(:)');
tolerance = 1e3*eps(max(1,halfWidth));
if numel(reference) < 5 || any(diff(reference) <= 0) || ...
        abs(reference(1)) > tolerance || ...
        abs(reference(end)-halfWidth) > tolerance
    error('ipm:gridlab:PositiveAxisReference', ...
        ['options.referenceAxis must increase from 0 to the requested ' ...
        'positive endpoint.']);
end
sourceCoordinate = linspace(0,1,numel(reference));
targetCoordinate = linspace(0,1,nodeCount);
positive = interp1(sourceCoordinate,reference,targetCoordinate,'pchip');
positive([1,end]) = [0,halfWidth];
if any(diff(positive) <= 0)
    error('ipm:gridlab:PositiveAxisReference', ...
        'The resampled positive reference axis is not increasing.');
end
full = [-fliplr(reference(2:end)),reference];
end

function [axis,alpha] = project_axis(requested,baseline,limits,halfWidth)
baselineLogWidths = log(diff(baseline));
requestedLogWidths = log(diff(requested));
validAlpha = 0;
invalidUpper = 1;
sampleAlpha = linspace(1,0,65);
for index = 1:numel(sampleAlpha)
    trialAlpha = sampleAlpha(index);
    trial = blend_axis( ...
        baselineLogWidths,requestedLogWidths,trialAlpha,halfWidth);
    [~,~,reasons] = positive_diagnostics(trial,limits);
    if isempty(reasons)
        validAlpha = trialAlpha;
        if index > 1
            invalidUpper = sampleAlpha(index-1);
        else
            invalidUpper = trialAlpha;
        end
        break
    end
end
for iteration = 1:40
    trialAlpha = 0.5*(validAlpha+invalidUpper);
    trial = blend_axis( ...
        baselineLogWidths,requestedLogWidths,trialAlpha,halfWidth);
    [~,~,reasons] = positive_diagnostics(trial,limits);
    if isempty(reasons)
        validAlpha = trialAlpha;
    else
        invalidUpper = trialAlpha;
    end
end
alpha = validAlpha;
axis = blend_axis( ...
    baselineLogWidths,requestedLogWidths,alpha,halfWidth);
end

function axis = blend_axis(baselineLogWidths,requestedLogWidths,alpha,limit)
logWidths = baselineLogWidths+alpha*( ...
    requestedLogWidths-baselineLogWidths);
widths = exp(logWidths-max(logWidths));
widths = limit*widths/sum(widths);
axis = [0,cumsum(widths)];
axis(end) = limit;
end

function [weights,quality,reasons] = positive_diagnostics(axis,limits)
weights = ipm.mesh.quadrature(axis);
quality = ipm.mesh.quality(axis,7,weights);
reasons = {};
if quality.maximumAdjacentCellRatio > limits.maxAdjacentCellRatio
    reasons{end+1} = 'adjacent_cell_ratio';
end
if quality.maximumLogSpacingCurvature > limits.maxLogSpacingCurvature
    reasons{end+1} = 'log_spacing_curvature';
end
if quality.minimumStencilRcond < limits.minStencilRcond
    reasons{end+1} = 'stencil_conditioning';
end
if quality.minimumQuadratureWeightRatio < ...
        limits.minQuadratureWeightRatio || any(weights <= 0)
    reasons{end+1} = 'quadrature_margin';
end
end
