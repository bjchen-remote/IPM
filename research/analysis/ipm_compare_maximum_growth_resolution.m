function [figureHandle,comparison] = ipm_compare_maximum_growth_resolution( ...
        coarseSource,fineSource,outputFile)
%IPM_COMPARE_MAXIMUM_GROWTH_RESOLUTION Compare trusted growth across grids.
%   IPM_COMPARE_MAXIMUM_GROWTH_RESOLUTION(COARSESOURCE,FINESOURCE) accepts
%   scalar solver-result/assessment structures or MAT-file paths containing
%   variables named result or assessment.  Result histories are restricted
%   to IPM.OUTPUT.CONTINUOUSTRUSTEDPREFIX before the growth-rate diagnostic
%   is evaluated; assessment inputs use their already-restricted
%   GRADIENTGROWTHRATEFIT and WALLPEAKGROWTHRATEFIT series.
%
%   IPM_COMPARE_MAXIMUM_GROWTH_RESOLUTION(...,OUTPUTFILE) exports the figure
%   at 300 dpi.  An omitted extension is replaced by .png.  The four panels
%   compare the physical maximum gradient, normalized logarithmic growth,
%   the representative logarithmic rate with its three-bandwidth envelope,
%   and the physical wall peak.  Circles and vertical dashed lines mark the
%   last sample of each continuous trusted prefix.
%
%   [FIGUREHANDLE,COMPARISON] = ... also returns concise source, endpoint,
%   dynamic-range, and rate summaries.

if nargin < 3
    outputFile = '';
end
if ~(ischar(outputFile) || (isstring(outputFile) && isscalar(outputFile)))
    error('ipm:MaximumGrowthResolutionPlotOutput', ...
        'The optional output file must be a character vector or scalar string.');
end
outputFile = char(outputFile);

coarse = read_growth_source(coarseSource,'coarse');
fine = read_growth_source(fineSource,'fine');
datasets = {coarse,fine};
colors = [0,0.4470,0.7410;0.8500,0.3250,0.0980];

if isempty(outputFile) && usejava('desktop')
    visibility = 'on';
else
    visibility = 'off';
end
figureHandle = figure('Color','w','Visible',visibility, ...
    'Position',[80,80,1500,980], ...
    'Name','Maximum-growth resolution comparison');
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');

maximumAxes = nexttile(layout,1);
hold(maximumAxes,'on');
maximumHandles = gobjects(1,2);
for index = 1:2
    data = datasets{index};
    maximumHandles(index) = endpoint_curve(maximumAxes, ...
        data.gradientFit.time,data.gradientFit.maximum,colors(index,:), ...
        data.displayLabel);
end
finish_axes(maximumAxes,'Physical time', ...
    'M(t), physical maximum gradient','Physical maximum-gradient growth');
add_trusted_lines(maximumAxes,datasets,colors);
maximumLegend = legend(maximumAxes,maximumHandles, ...
    {coarse.displayLabel,fine.displayLabel},'Location','northwest');
style_legend(maximumLegend);

logAxes = nexttile(layout,2);
hold(logAxes,'on');
logHandles = gobjects(1,2);
for index = 1:2
    data = datasets{index};
    normalizedLogMaximum = data.gradientFit.logMaximum - ...
        data.gradientFit.logMaximum(1);
    logHandles(index) = endpoint_curve(logAxes,data.gradientFit.time, ...
        normalizedLogMaximum,colors(index,:),data.displayLabel);
end
finish_axes(logAxes,'Physical time','log(M(t)/M(0))', ...
    'Normalized logarithmic growth');
add_trusted_lines(logAxes,datasets,colors);
logLegend = legend(logAxes,logHandles, ...
    {coarse.displayLabel,fine.displayLabel},'Location','northwest');
style_legend(logLegend);

rateAxes = nexttile(layout,3);
hold(rateAxes,'on');
rateHandles = gobjects(1,4);
rateLabels = cell(1,4);
handleIndex = 0;
for index = 1:2
    data = datasets{index};
    fit = data.gradientFit;
    lowerRate = min(fit.logarithmicRate,[],2,'omitnan');
    upperRate = max(fit.logarithmicRate,[],2,'omitnan');
    validBand = isfinite(fit.time) & isfinite(lowerRate) & ...
        isfinite(upperRate);
    handleIndex = handleIndex+1;
    rateHandles(handleIndex) = fill(rateAxes, ...
        [fit.time(validBand);flipud(fit.time(validBand))], ...
        [lowerRate(validBand);flipud(upperRate(validBand))], ...
        colors(index,:),'FaceAlpha',0.14,'EdgeColor','none');
    rateLabels{handleIndex} = sprintf('%s: bandwidth sensitivity', ...
        data.displayLabel);
    representativeIndex = representative_bandwidth(fit);
    representativeRate = fit.logarithmicRate(:,representativeIndex);
    handleIndex = handleIndex+1;
    rateHandles(handleIndex) = endpoint_curve(rateAxes,fit.time, ...
        representativeRate,colors(index,:),data.displayLabel);
    rateLabels{handleIndex} = sprintf('%s: gamma, h = %.4g', ...
        data.displayLabel,fit.bandwidths(representativeIndex));
end
finish_axes(rateAxes,'Physical time', ...
    'gamma(t) = d log(M)/dt','Logarithmic growth rate');
add_trusted_lines(rateAxes,datasets,colors);
rateLegend = legend(rateAxes,rateHandles,rateLabels, ...
    'Location','northwest');
style_legend(rateLegend);

wallAxes = nexttile(layout,4);
hold(wallAxes,'on');
wallHandles = gobjects(1,2);
for index = 1:2
    data = datasets{index};
    wallHandles(index) = endpoint_curve(wallAxes,data.wallFit.time, ...
        data.wallFit.maximum,colors(index,:),data.displayLabel);
end
finish_axes(wallAxes,'Physical time','Physical wall peak', ...
    'Wall-peak growth');
add_trusted_lines(wallAxes,datasets,colors);
wallLegend = legend(wallAxes,wallHandles, ...
    {coarse.displayLabel,fine.displayLabel},'Location','northwest');
style_legend(wallLegend);

titleLine = sprintf('Maximum-growth resolution comparison: %s versus %s', ...
    coarse.resolutionLabel,fine.resolutionLabel);
endpointLine = sprintf( ...
    'continuous trusted endpoints: coarse t_* = %.9g; fine t_* = %.9g', ...
    coarse.trustedThroughTime,fine.trustedThroughTime);
sgtitle(layout,{titleLine,endpointLine},'Interpreter','none', ...
    'FontWeight','bold','Color',[0.12,0.12,0.12]);

if ~isempty(outputFile)
    outputFile = export_figure(figureHandle,outputFile);
end

comparison = struct();
comparison.coarse = comparison_summary(coarse);
comparison.fine = comparison_summary(fine);
comparison.commonTrustedThroughTime = min( ...
    coarse.trustedThroughTime,fine.trustedThroughTime);
comparison.outputFile = outputFile;
end

function data = read_growth_source(source,role)
[payload,sourceKind,sourceLabel] = read_payload(source);
if is_assessment(payload)
    validate_assessment_prefix(payload);
    gradientFit = payload.gradientGrowthRateFit;
    wallFit = payload.wallPeakGrowthRateFit;
    metadata = optional_struct(payload,'sourceMetadata');
else
    result = ipm.output.validate(payload);
    [trusted,~] = ipm.output.trustedMask(result);
    trustedPrefix = ipm.output.continuousTrustedPrefix(trusted);
    if ~any(trustedPrefix)
        error('ipm:MaximumGrowthResolutionTrustedPrefix', ...
            'The %s result has no continuous trusted-prefix records.',role);
    end
    common = result.history.common;
    time = common.physicalTime(trustedPrefix);
    gradient = common.physicalGradInf(trustedPrefix);
    if isfield(common,'physicalTrackedWallOmegaPeak')
        wallPeak = common.physicalTrackedWallOmegaPeak(trustedPrefix);
    else
        wallPeak = common.physicalWallOmegaPeak(trustedPrefix);
    end
    gradientFit = ipm.diagnostics.maximumGrowthRateFit(time,gradient);
    wallFit = ipm.diagnostics.maximumGrowthRateFit(time,wallPeak);
    metadata = result.metadata;
    payload = result;
end
validate_rate_fit(gradientFit,[role,' gradient']);
validate_rate_fit(wallFit,[role,' wall peak']);
if ~same_time_grid(gradientFit.time,wallFit.time)
    error('ipm:MaximumGrowthResolutionTimeGrid', ...
        'The %s gradient and wall fits do not share one trusted time grid.',role);
end

resolutionLabel = resolution_label(payload,metadata,sourceLabel);
data = struct();
data.role = role;
data.sourceKind = sourceKind;
data.sourceLabel = sourceLabel;
data.resolutionLabel = resolutionLabel;
data.displayLabel = sprintf('%s: %s',role,resolutionLabel);
data.gradientFit = gradientFit;
data.wallFit = wallFit;
data.trustedThroughTime = gradientFit.time(end);
end

function [payload,sourceKind,sourceLabel] = read_payload(source)
sourceLabel = 'unnamed source';
if ischar(source) || (isstring(source) && isscalar(source))
    fileName = char(source);
    if ~isfile(fileName)
        error('ipm:MaximumGrowthResolutionFileMissing', ...
            'Input MAT-file does not exist: %s',fileName);
    end
    loaded = load(fileName);
    [~,sourceLabel] = fileparts(fileName);
    if isfield(loaded,'assessment')
        payload = loaded.assessment;
        sourceKind = 'assessment file';
    elseif isfield(loaded,'result')
        payload = loaded.result;
        sourceKind = 'result file';
    else
        error('ipm:MaximumGrowthResolutionFileContract', ...
            'MAT-file must contain a variable named result or assessment.');
    end
elseif isstruct(source) && isscalar(source)
    payload = source;
    if is_assessment(payload)
        sourceKind = 'assessment struct';
    else
        sourceKind = 'result struct';
    end
else
    error('ipm:MaximumGrowthResolutionInput', ...
        'Pass scalar result/assessment structures or their MAT-file paths.');
end
end

function valid = is_assessment(value)
valid = isstruct(value) && isscalar(value) && ...
    isfield(value,'gradientGrowthRateFit') && ...
    isfield(value,'wallPeakGrowthRateFit');
end

function validate_assessment_prefix(assessment)
gradientFit = assessment.gradientGrowthRateFit;
wallFit = assessment.wallPeakGrowthRateFit;
if isfield(assessment,'trustedPrefixMask')
    prefix = logical(assessment.trustedPrefixMask(:));
    if ~isequal(prefix,ipm.output.continuousTrustedPrefix(prefix))
        error('ipm:MaximumGrowthResolutionAssessmentPrefix', ...
            'Assessment trustedPrefixMask is not a continuous prefix.');
    end
    if isfield(assessment,'trustedRecords') && ...
            nnz(prefix) ~= assessment.trustedRecords
        error('ipm:MaximumGrowthResolutionAssessmentCount', ...
            'Assessment trusted-prefix record counts are inconsistent.');
    end
end
if isfield(assessment,'trustedRecords') && ...
        (gradientFit.points ~= assessment.trustedRecords || ...
        wallFit.points ~= assessment.trustedRecords)
    error('ipm:MaximumGrowthResolutionAssessmentFitCount', ...
        'Assessment growth fits do not contain exactly the trusted prefix.');
end
if isfield(assessment,'trustedThroughTime') && gradientFit.valid
    tolerance = 64*eps(max(1,abs(assessment.trustedThroughTime)));
    if abs(gradientFit.time(end)-assessment.trustedThroughTime) > tolerance
        error('ipm:MaximumGrowthResolutionAssessmentEnd', ...
            'Assessment growth fit extends beyond its trusted endpoint.');
    end
end
end

function validate_rate_fit(fit,label)
required = {'valid','points','time','maximum','logMaximum', ...
    'logarithmicRate','bandwidths','eFoldCount'};
if ~isstruct(fit) || ~isscalar(fit) || ~all(isfield(fit,required)) || ...
        ~fit.valid
    error('ipm:MaximumGrowthResolutionRateFit', ...
        'The %s maximum-growth-rate fit is unavailable or invalid.',label);
end
pointCount = numel(fit.time);
if pointCount ~= fit.points || numel(fit.maximum) ~= pointCount || ...
        numel(fit.logMaximum) ~= pointCount || ...
        size(fit.logarithmicRate,1) ~= pointCount || ...
        size(fit.logarithmicRate,2) ~= numel(fit.bandwidths) || ...
        any(~isfinite(fit.time)) || any(diff(fit.time) <= 0) || ...
        any(~isfinite(fit.maximum)) || any(fit.maximum <= 0)
    error('ipm:MaximumGrowthResolutionRateFitSize', ...
        'The %s growth-fit fields are inconsistent.',label);
end
end

function same = same_time_grid(first,second)
first = first(:);
second = second(:);
if numel(first) ~= numel(second)
    same = false;
    return;
end
tolerance = 64*eps(max(1,max(abs([first;second]))));
same = all(abs(first-second) <= tolerance);
end

function label = resolution_label(payload,metadata,fallback)
label = '';
caseMetadata = optional_struct(metadata,'caseMetadata');
if isfield(caseMetadata,'gridInterpretation')
    gridText = char(caseMetadata.gridInterpretation);
    tokens = regexp(gridText, ...
        '(\d+)\s*x\s*(\d+)\s*cells','tokens','once');
    if ~isempty(tokens)
        label = sprintf('%s x %s cells',tokens{1},tokens{2});
    end
end
if isempty(label) && isfield(caseMetadata,'caseName')
    caseName = char(caseMetadata.caseName);
    tokens = regexp(caseName,'quadrant_(\d+)','tokens','once');
    if ~isempty(tokens)
        label = sprintf('%s x %s cells',tokens{1},tokens{1});
    end
end
if isempty(label) && isfield(payload,'config') && ...
        isfield(payload.config,'grid') && ...
        all(isfield(payload.config.grid,{'nx','ny'}))
    nx = payload.config.grid.nx;
    ny = payload.config.grid.ny;
    label = sprintf('%d x %d stored nodes',nx,ny);
end
if isempty(label)
    label = fallback;
end
end

function value = optional_struct(source,fieldName)
value = struct();
if isstruct(source) && isscalar(source) && isfield(source,fieldName) && ...
        isstruct(source.(fieldName)) && isscalar(source.(fieldName))
    value = source.(fieldName);
end
end

function lineHandle = endpoint_curve(axesHandle,time,value,color,label)
time = time(:);
value = value(:);
lineHandle = plot(axesHandle,time,value,'Color',color,'LineWidth',1.8, ...
    'Marker','o','MarkerIndices',numel(time),'MarkerSize',7, ...
    'MarkerFaceColor',color,'DisplayName',label);
end

function add_trusted_lines(axesHandle,datasets,colors)
for index = 1:2
    xline(axesHandle,datasets{index}.trustedThroughTime,'--', ...
        'Color',colors(index,:),'LineWidth',1.0, ...
        'HandleVisibility','off');
end
end

function finish_axes(axesHandle,xLabel,yLabel,titleText)
hold(axesHandle,'off');
grid(axesHandle,'on');
box(axesHandle,'on');
set(axesHandle,'Color','w','XColor',[0.15,0.15,0.15], ...
    'YColor',[0.15,0.15,0.15],'GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.35,'LineWidth',0.8);
xlabel(axesHandle,xLabel);
ylabel(axesHandle,yLabel);
title(axesHandle,titleText,'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);
end

function style_legend(legendHandle)
set(legendHandle,'Color','w','TextColor',[0.12,0.12,0.12], ...
    'EdgeColor',[0.55,0.55,0.55],'Interpreter','none');
end

function index = representative_bandwidth(fit)
if isfield(fit,'bandwidthFactors') && ...
        numel(fit.bandwidthFactors) == numel(fit.bandwidths)
    [~,index] = min(abs(fit.bandwidthFactors-1));
else
    index = ceil(numel(fit.bandwidths)/2);
end
end

function outputFile = export_figure(figureHandle,outputFile)
[outputDirectory,~,extension] = fileparts(outputFile);
if isempty(extension)
    outputFile = [outputFile,'.png'];
elseif ~strcmpi(extension,'.png')
    close(figureHandle);
    error('ipm:MaximumGrowthResolutionPlotFormat', ...
        'The output file must use the PNG extension.');
end
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    [created,message] = mkdir(outputDirectory);
    if ~created
        close(figureHandle);
        error('ipm:MaximumGrowthResolutionPlotDirectory', ...
            'Could not create output directory: %s',message);
    end
end
drawnow;
exportgraphics(figureHandle,outputFile,'Resolution',300, ...
    'BackgroundColor','white');
end

function summary = comparison_summary(data)
fit = data.gradientFit;
representativeIndex = representative_bandwidth(fit);
summary = struct( ...
    'sourceKind',data.sourceKind, ...
    'sourceLabel',data.sourceLabel, ...
    'resolutionLabel',data.resolutionLabel, ...
    'trustedThroughTime',data.trustedThroughTime, ...
    'trustedRecords',fit.points, ...
    'gradientGrowth',fit.maximum(end)/fit.maximum(1), ...
    'gradientEFolds',fit.eFoldCount, ...
    'finalRepresentativeRate', ...
        fit.logarithmicRate(end,representativeIndex), ...
    'finalRateBandwidthRange',[min(fit.logarithmicRate(end,:), ...
        [],'omitnan'),max(fit.logarithmicRate(end,:),[],'omitnan')], ...
    'wallPeakGrowth',data.wallFit.maximum(end)/data.wallFit.maximum(1));
end
