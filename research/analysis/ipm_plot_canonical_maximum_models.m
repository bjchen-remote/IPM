function [figureHandle,analysis] = ipm_plot_canonical_maximum_models( ...
        firstSource,secondSource,outputFile)
%IPM_PLOT_CANONICAL_MAXIMUM_MODELS Compare canonical-tau growth models.
%   IPM_PLOT_CANONICAL_MAXIMUM_MODELS(FIRST,SECOND) accepts two version-2
%   solver results, fourth-order blow-up assessments, or MAT-file paths
%   containing variables named result or assessment.  Result histories are
%   restricted to their continuous trusted prefixes.  Assessment inputs
%   must carry a continuous trusted-prefix mask; their saved canonical-tau
%   series are refitted with IPM.DIAGNOSTICS.CANONICALMAXIMUMGROWTHFIT.
%
%   The four panels show the physical maximum gradient and tracked wall
%   maximum versus canonical tau, followed by their canonical logarithmic
%   rates.  Both the linear-maximum and exponential-maximum candidates are
%   shown over the largest accumulated-growth window shared by both signals
%   and both sources.  The quadratic-log model remains a classification
%   guard and is not plotted as a third growth law.
%
%   IPM_PLOT_CANONICAL_MAXIMUM_MODELS(...,OUTPUTFILE) exports a 300 dpi PNG.
%   [FIGUREHANDLE,ANALYSIS] also returns the freshly computed fits, selected
%   windows, model decisions, AICc separations, validation improvements, and
%   input-sampling stability flags.  Exporting a PNG also saves ANALYSIS to
%   a sibling *_analysis.mat file.  If an early stop leaves no e-fold window
%   common to all four fits, the trusted data and rates are still plotted,
%   while modelFitsAvailable is false and candidate curves are omitted.

if nargin < 3
    outputFile = '';
end
validate_text(outputFile,'output file');
outputFile = char(outputFile);

first = read_canonical_source(firstSource,'first');
second = read_canonical_source(secondSource,'second');
if strcmp(first.displayLabel,second.displayLabel)
    first.displayLabel = [first.displayLabel,' (first)'];
    second.displayLabel = [second.displayLabel,' (second)'];
end
datasets = {first,second};
sharedLogGrowth = shared_growth_window(datasets);
modelFitsAvailable = isfinite(sharedLogGrowth);

visibility = 'off';
if isempty(outputFile) && usejava('desktop')
    visibility = 'on';
end
figureHandle = figure('Color','w','Visible',visibility, ...
    'Position',[60,60,1680,1080], ...
    'Name','Canonical-tau maximum-growth models');
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');
colors = [0,0.4470,0.7410;0.8500,0.3250,0.0980];

draw_maximum_panel(nexttile(layout,1),datasets,colors, ...
    'gradientFit',sharedLogGrowth,'Physical maximum gradient');
draw_maximum_panel(nexttile(layout,2),datasets,colors, ...
    'wallFit',sharedLogGrowth,'Physical tracked wall maximum');
draw_rate_panel(nexttile(layout,3),datasets,colors, ...
    'gradientFit',sharedLogGrowth, ...
    'Canonical logarithmic rate of gradient maximum');
draw_rate_panel(nexttile(layout,4),datasets,colors, ...
    'wallFit',sharedLogGrowth, ...
    'Canonical logarithmic rate of wall maximum');

titleLine = sprintf('Canonical-tau maximum growth: %s versus %s', ...
    first.displayLabel,second.displayLabel);
if modelFitsAvailable
    windowLine = sprintf([ ...
        'linear M=A(1+k(\\tau-\\tau_0)); exponential ' ...
        'M=A exp(\\lambda(\\tau-\\tau_0)); shared %.3g e-fold window'], ...
        sharedLogGrowth);
else
    windowLine = [ ...
        'no canonical growth window is shared by all four fits; ' ...
        'trusted data/rates only'];
end
sgtitle(layout,{titleLine,windowLine},'Interpreter','tex', ...
    'FontWeight','bold','Color',[0.12,0.12,0.12]);

if ~isempty(outputFile)
    outputFile = export_figure(figureHandle,outputFile);
end

analysis = struct();
analysis.definition = struct( ...
    'linearMaximum','M=A*(1+k*(canonicalTau-referenceCanonicalTau))', ...
    'exponentialMaximum', ...
        'M=A*exp(lambda*(canonicalTau-referenceCanonicalTau))', ...
    'rate','d(log(M))/d(canonicalTau)', ...
    'selection',[ ...
        'Largest requested e-fold window present for gradient and wall ' ...
        'fits on both sources.']);
analysis.sharedRequestedLogGrowth = sharedLogGrowth;
analysis.modelFitsAvailable = modelFitsAvailable;
analysis.first = analysis_summary(first,sharedLogGrowth);
analysis.second = analysis_summary(second,sharedLogGrowth);
analysis.outputFile = outputFile;
analysis.analysisFile = '';
if ~isempty(outputFile)
    analysis.analysisFile = analysis_file_name(outputFile);
    save(analysis.analysisFile,'analysis');
end
end

function data = read_canonical_source(source,role)
[payload,sourceKind,sourceLabel] = read_payload(source,role);
if is_assessment(payload)
    [gradientFit,wallFit,metadata] = assessment_fits(payload,role);
    resolutionPayload = payload;
else
    result = ipm.output.validate(payload);
    [gradientFit,wallFit] = result_fits(result,role);
    metadata = result.metadata;
    resolutionPayload = result;
end

data = struct();
data.role = role;
data.sourceKind = sourceKind;
data.sourceLabel = sourceLabel;
data.displayLabel = resolution_label( ...
    resolutionPayload,metadata,sourceLabel);
data.gradientFit = gradientFit;
data.wallFit = wallFit;
end

function [payload,sourceKind,sourceLabel] = read_payload(source,role)
sourceLabel = [role,' source'];
if ischar(source) || (isstring(source) && isscalar(source))
    fileName = char(source);
    if ~isfile(fileName)
        error('ipm:CanonicalMaximumModelsFileMissing', ...
            'The %s input MAT-file does not exist: %s.',role,fileName);
    end
    loaded = load(fileName);
    [~,sourceLabel] = fileparts(fileName);
    if isfield(loaded,'result')
        payload = loaded.result;
        sourceKind = 'result file';
    elseif isfield(loaded,'assessment')
        payload = loaded.assessment;
        sourceKind = 'assessment file';
    else
        error('ipm:CanonicalMaximumModelsFileContract', ...
            ['The %s MAT-file must contain a variable named result or ' ...
            'assessment.'],role);
    end
elseif isstruct(source) && isscalar(source)
    payload = source;
    if is_assessment(payload)
        sourceKind = 'assessment struct';
    else
        sourceKind = 'result struct';
    end
else
    error('ipm:CanonicalMaximumModelsInput', ...
        ['The %s source must be a scalar result/assessment structure ' ...
        'or a MAT-file path.'],role);
end
end

function valid = is_assessment(value)
valid = isstruct(value) && isscalar(value) && ...
    isfield(value,'trustedPrefixMask') && ...
    isfield(value,'gradientCanonicalGrowthFit') && ...
    isfield(value,'wallPeakCanonicalGrowthFit');
end

function [gradientFit,wallFit,metadata] = assessment_fits(assessment,role)
required = {'trustedPrefixMask','trustedRecords', ...
    'trustedThroughCanonicalTau','gradientCanonicalGrowthFit', ...
    'wallPeakCanonicalGrowthFit'};
require_fields(assessment,required,[role,' assessment']);
prefix = assessment.trustedPrefixMask;
if ~islogical(prefix) || ~isvector(prefix) || isempty(prefix)
    error('ipm:CanonicalMaximumModelsAssessmentPrefix', ...
        'The %s assessment trustedPrefixMask must be a logical vector.',role);
end
prefix = prefix(:);
if ~isequal(prefix,ipm.output.continuousTrustedPrefix(prefix)) || ...
        ~any(prefix)
    error('ipm:CanonicalMaximumModelsAssessmentPrefix', ...
        ['The %s assessment trustedPrefixMask must be one nonempty ' ...
        'continuous prefix.'],role);
end
if assessment.trustedRecords ~= nnz(prefix)
    error('ipm:CanonicalMaximumModelsAssessmentCount', ...
        'The %s assessment trusted-record count is inconsistent.',role);
end

savedGradient = assessment.gradientCanonicalGrowthFit;
savedWall = assessment.wallPeakCanonicalGrowthFit;
validate_canonical_fit(savedGradient,[role,' saved gradient']);
validate_canonical_fit(savedWall,[role,' saved wall']);
if savedGradient.points ~= assessment.trustedRecords || ...
        savedWall.points ~= assessment.trustedRecords
    error('ipm:CanonicalMaximumModelsAssessmentFitCount', ...
        ['The %s canonical fits do not contain exactly the saved trusted ' ...
        'prefix.'],role);
end
if ~same_coordinate(savedGradient.canonicalTau,savedWall.canonicalTau)
    error('ipm:CanonicalMaximumModelsAssessmentTau', ...
        'The %s gradient and wall fits use different canonical-tau grids.',role);
end
tolerance = 64*eps(max(1,abs(assessment.trustedThroughCanonicalTau)));
if abs(savedGradient.canonicalTau(end)- ...
        assessment.trustedThroughCanonicalTau) > tolerance
    error('ipm:CanonicalMaximumModelsAssessmentEnd', ...
        ['The %s saved canonical fit does not end at the assessment ' ...
        'trusted endpoint.'],role);
end

% Refit rather than trusting a possibly stale saved model decision.
gradientFit = ipm.diagnostics.canonicalMaximumGrowthFit( ...
    savedGradient.canonicalTau,savedGradient.maximum);
wallFit = ipm.diagnostics.canonicalMaximumGrowthFit( ...
    savedWall.canonicalTau,savedWall.maximum);
validate_canonical_fit(gradientFit,[role,' gradient']);
validate_canonical_fit(wallFit,[role,' wall']);
metadata = optional_struct(assessment,'sourceMetadata');
end

function [gradientFit,wallFit] = result_fits(result,role)
[rawTrusted,~] = ipm.output.trustedMask(result);
trustedPrefix = ipm.output.continuousTrustedPrefix(rawTrusted);
trustedPrefix = logical(trustedPrefix(:));
if ~any(trustedPrefix)
    error('ipm:CanonicalMaximumModelsTrustedPrefix', ...
        'The %s result has no continuous trusted-prefix records.',role);
end
common = result.history.common;
required = {'canonicalTau','physicalGradInf','physicalWallOmegaPeak'};
require_fields(common,required,[role,' history.common']);
tau = trusted_series(common.canonicalTau,trustedPrefix, ...
    [role,' canonical tau']);
gradient = trusted_series(common.physicalGradInf,trustedPrefix, ...
    [role,' physical gradient']);
if isfield(common,'physicalTrackedWallOmegaPeak')
    wall = trusted_series(common.physicalTrackedWallOmegaPeak, ...
        trustedPrefix,[role,' tracked wall maximum']);
else
    wall = trusted_series(common.physicalWallOmegaPeak,trustedPrefix, ...
        [role,' wall maximum']);
end
gradientFit = canonical_fit_with_plot_series(tau,gradient);
wallFit = canonical_fit_with_plot_series(tau,wall);
validate_canonical_fit(gradientFit,[role,' gradient']);
validate_canonical_fit(wallFit,[role,' wall']);
end

function fit = canonical_fit_with_plot_series(tau,maximum)
fit = ipm.diagnostics.canonicalMaximumGrowthFit(tau,maximum);
if fit.valid
    fit.plotDataOnly = false;
    return;
end
% canonicalMaximumGrowthFit deliberately withholds a fit below 16 samples.
% Preserve those trusted samples only so an exceptionally early solver stop
% can still produce a truthful data-only artifact.  No rate or model is
% synthesized, and FIT.VALID remains false.
fit.points = numel(tau);
fit.canonicalTau = tau(:);
fit.maximum = maximum(:);
fit.logMaximum = log(maximum(:));
fit.canonicalLogarithmicRate = NaN(numel(tau),1);
fit.canonicalLogarithmicRateDerivative = NaN(numel(tau),1);
fit.bandwidthsCanonicalTau = NaN;
fit.bandwidthFactors = 1;
fit.totalLogGrowth = fit.logMaximum(end)-fit.logMaximum(1);
fit.eFoldCount = fit.totalLogGrowth;
fit.finalCanonicalRate = NaN;
fit.finalCanonicalRateDerivative = NaN;
fit.plotDataOnly = true;
end

function values = trusted_series(values,trustedPrefix,label)
if ~isnumeric(values) || ~isreal(values) || ~isvector(values) || ...
        numel(values) ~= numel(trustedPrefix)
    error('ipm:CanonicalMaximumModelsSeries', ...
        'The %s series is not aligned with the trusted mask.',label);
end
values = values(:);
values = values(trustedPrefix);
end

function validate_canonical_fit(fit,label)
required = {'valid','points','coordinateName','canonicalTau','maximum', ...
    'logMaximum','canonicalLogarithmicRate', ...
    'bandwidthsCanonicalTau','windows','eFoldCount', ...
    'enoughDynamicRange','modelConsensus','consensusModel'};
if ~isstruct(fit) || ~isscalar(fit) || ~all(isfield(fit,required)) || ...
        ~strcmp(fit.coordinateName,'canonical_tau')
    error('ipm:CanonicalMaximumModelsFit', ...
        'The %s canonicalMaximumGrowthFit output is unavailable.',label);
end
points = numel(fit.canonicalTau);
if fit.points ~= points || numel(fit.maximum) ~= points || ...
        numel(fit.logMaximum) ~= points || ...
        size(fit.canonicalLogarithmicRate,1) ~= points || ...
        size(fit.canonicalLogarithmicRate,2) ~= ...
            numel(fit.bandwidthsCanonicalTau) || ...
        any(~isfinite(fit.canonicalTau)) || ...
        any(diff(fit.canonicalTau) <= 0) || ...
        any(~isfinite(fit.maximum)) || any(fit.maximum <= 0)
    error('ipm:CanonicalMaximumModelsFitSize', ...
        'The %s canonical fit fields are inconsistent.',label);
end
end

function shared = shared_growth_window(datasets)
fits = {datasets{1}.gradientFit,datasets{1}.wallFit, ...
    datasets{2}.gradientFit,datasets{2}.wallFit};
shared = [fits{1}.windows.requestedLogGrowth];
for index = 2:numel(fits)
    candidate = [fits{index}.windows.requestedLogGrowth];
    shared = shared(arrayfun(@(value)any( ...
        abs(candidate-value) <= 64*eps(max(1,abs(value)))),shared));
end
if isempty(shared)
    shared = NaN;
    return;
end
shared = max(shared);
end

function draw_maximum_panel( ...
        axesHandle,datasets,colors,fitName,requestedGrowth,titleText)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:numel(datasets)
    data = datasets{index};
    fit = data.(fitName);
    normalizedMaximum = fit.maximum/fit.maximum(1);
    handles(end+1) = semilogy(axesHandle,fit.canonicalTau, ...
        normalizedMaximum,'Color',colors(index,:),'LineWidth',1.8, ...
        'Marker','o','MarkerIndices',fit.points,'MarkerSize',6, ...
        'MarkerFaceColor',colors(index,:)); %#ok<AGROW>
    if isfinite(requestedGrowth)
        window = selected_window(fit,requestedGrowth);
        labels{end+1,1} = ...
            decision_label(data.displayLabel,window); %#ok<AGROW>
    else
        labels{end+1,1} = ...
            unavailable_label(data.displayLabel); %#ok<AGROW>
        continue;
    end

    modelTau = linspace(window.startCanonicalTau, ...
        window.endCanonicalTau,300)';
    if window.linearMaximum.valid
        linearMaximum = exp(model_log_value( ...
            window.linearMaximum,modelTau))/fit.maximum(1);
        handles(end+1) = semilogy(axesHandle,modelTau,linearMaximum,':', ...
            'Color',colors(index,:),'LineWidth',1.7); %#ok<AGROW>
        labels{end+1,1} = [data.displayLabel,' linear M']; %#ok<AGROW>
    end
    if window.exponentialMaximum.valid
        exponentialMaximum = exp(model_log_value( ...
            window.exponentialMaximum,modelTau))/fit.maximum(1);
        handles(end+1) = semilogy(axesHandle,modelTau, ...
            exponentialMaximum,'--','Color',colors(index,:), ...
            'LineWidth',1.55); %#ok<AGROW>
        labels{end+1,1} = [data.displayLabel,' exponential M']; %#ok<AGROW>
    end
end
finish_axes(axesHandle,'Canonical tau','M(\tau)/M(\tau_{initial})', ...
    titleText);
fitLegend = legend(axesHandle,handles,labels,'Location','northwest');
style_legend(fitLegend);
end

function draw_rate_panel( ...
        axesHandle,datasets,colors,fitName,requestedGrowth,titleText)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:numel(datasets)
    data = datasets{index};
    fit = data.(fitName);
    rate = fit.canonicalLogarithmicRate;
    lowerRate = min(rate,[],2,'omitnan');
    upperRate = max(rate,[],2,'omitnan');
    validBand = isfinite(fit.canonicalTau) & isfinite(lowerRate) & ...
        isfinite(upperRate);
    if any(validBand)
        fill(axesHandle, ...
            [fit.canonicalTau(validBand); ...
                flipud(fit.canonicalTau(validBand))], ...
            [lowerRate(validBand);flipud(upperRate(validBand))], ...
            colors(index,:),'FaceAlpha',0.13,'EdgeColor','none', ...
            'HandleVisibility','off');
    end
    representativeIndex = representative_bandwidth(fit);
    handles(end+1) = plot(axesHandle,fit.canonicalTau, ...
        rate(:,representativeIndex),'Color',colors(index,:), ...
        'LineWidth',1.8,'Marker','o','MarkerIndices',fit.points, ...
        'MarkerSize',6,'MarkerFaceColor',colors(index,:)); %#ok<AGROW>
    if isfinite(requestedGrowth)
        window = selected_window(fit,requestedGrowth);
        labels{end+1,1} = ...
            decision_label(data.displayLabel,window); %#ok<AGROW>
    else
        labels{end+1,1} = ...
            unavailable_label(data.displayLabel); %#ok<AGROW>
        continue;
    end

    modelTau = linspace(window.startCanonicalTau, ...
        window.endCanonicalTau,300)';
    if window.linearMaximum.valid
        handles(end+1) = plot(axesHandle,modelTau, ...
            model_rate(window.linearMaximum,modelTau),':', ...
            'Color',colors(index,:),'LineWidth',1.7); %#ok<AGROW>
        labels{end+1,1} = [data.displayLabel,' linear rate']; %#ok<AGROW>
    end
    if window.exponentialMaximum.valid
        handles(end+1) = plot(axesHandle,modelTau, ...
            model_rate(window.exponentialMaximum,modelTau),'--', ...
            'Color',colors(index,:),'LineWidth',1.55); %#ok<AGROW>
        labels{end+1,1} = [data.displayLabel,' exponential rate']; %#ok<AGROW>
    end
end
finish_axes(axesHandle,'Canonical tau', ...
    'd log(M)/d\tau',titleText);
rateLegend = legend(axesHandle,handles,labels,'Location','northeast');
style_legend(rateLegend);
end

function window = selected_window(fit,requestedGrowth)
values = [fit.windows.requestedLogGrowth];
[distance,index] = min(abs(values-requestedGrowth));
if isempty(distance) || distance > 64*eps(max(1,abs(requestedGrowth)))
    error('ipm:CanonicalMaximumModelsWindowLookup', ...
        'A previously identified shared growth window is missing.');
end
window = fit.windows(index);
end

function value = model_log_value(model,tau)
x = tau-model.referenceCanonicalTau;
switch model.name
    case 'linear_maximum'
        value = model.parameters.logA+log1p(model.parameters.k*x);
    case 'exponential_maximum'
        value = model.parameters.logA+model.parameters.lambda*x;
    otherwise
        value = NaN(size(tau));
end
end

function value = model_rate(model,tau)
x = tau-model.referenceCanonicalTau;
switch model.name
    case 'linear_maximum'
        value = model.parameters.k./(1+model.parameters.k*x);
    case 'exponential_maximum'
        value = repmat(model.parameters.lambda,size(tau));
    otherwise
        value = NaN(size(tau));
end
end

function index = representative_bandwidth(fit)
if isfield(fit,'bandwidthFactors') && ...
        numel(fit.bandwidthFactors) == ...
            numel(fit.bandwidthsCanonicalTau)
    [~,index] = min(abs(fit.bandwidthFactors-1));
else
    index = ceil(numel(fit.bandwidthsCanonicalTau)/2);
end
end

function value = decision_label(sourceLabel,window)
validation = numeric_label(window.validationImprovementOverAlternative);
value = sprintf([ ...
    '%s data | decision=%s; AICc separation=%.4g; ' ...
    'validation improvement=%s; sampling stable=%d'], ...
    sourceLabel,window.modelDecision,window.aiccSeparation,validation, ...
    window.inputSamplingStable);
end

function value = unavailable_label(sourceLabel)
value = [sourceLabel,' data | no shared canonical model window'];
end

function value = numeric_label(number)
if isfinite(number)
    value = sprintf('%.4g',number);
else
    value = 'n/a';
end
end

function summary = analysis_summary(data,requestedGrowth)
summary = struct();
summary.sourceKind = data.sourceKind;
summary.sourceLabel = data.sourceLabel;
summary.displayLabel = data.displayLabel;
summary.trustedRecords = data.gradientFit.points;
summary.trustedThroughCanonicalTau = ...
    data.gradientFit.canonicalTau(end);
summary.gradient = fit_summary(data.gradientFit,requestedGrowth);
summary.wall = fit_summary(data.wallFit,requestedGrowth);
summary.gradientCanonicalGrowthFit = data.gradientFit;
summary.wallCanonicalGrowthFit = data.wallFit;
end

function summary = fit_summary(fit,requestedGrowth)
if isfinite(requestedGrowth)
    selectedWindow = window_summary( ...
        selected_window(fit,requestedGrowth));
else
    selectedWindow = unavailable_window_summary();
end
summary = struct( ...
    'eFoldCount',fit.eFoldCount, ...
    'enoughDynamicRange',fit.enoughDynamicRange, ...
    'modelConsensus',fit.modelConsensus, ...
    'consensusModel',fit.consensusModel, ...
    'availableRequestedLogGrowth',[fit.windows.requestedLogGrowth], ...
    'selectedWindow',selectedWindow);
end

function summary = window_summary(window)
summary = struct( ...
    'requestedLogGrowth',window.requestedLogGrowth, ...
    'actualLogGrowth',window.actualLogGrowth, ...
    'startCanonicalTau',window.startCanonicalTau, ...
    'endCanonicalTau',window.endCanonicalTau, ...
    'modelDecision',window.modelDecision, ...
    'preferredByAicc',window.preferredByAicc, ...
    'preferredByValidation',window.preferredByValidation, ...
    'aiccSeparation',window.aiccSeparation, ...
    'validationImprovementOverAlternative', ...
        window.validationImprovementOverAlternative, ...
    'preSamplingDecision',window.preSamplingDecision, ...
    'uniformInputDecision',window.uniformInputDecision, ...
    'inputSamplingStable',window.inputSamplingStable, ...
    'linearMaximum',window.linearMaximum, ...
    'exponentialMaximum',window.exponentialMaximum);
end

function summary = unavailable_window_summary()
summary = struct( ...
    'requestedLogGrowth',NaN, ...
    'actualLogGrowth',NaN, ...
    'startCanonicalTau',NaN, ...
    'endCanonicalTau',NaN, ...
    'modelDecision','unavailable', ...
    'preferredByAicc','none', ...
    'preferredByValidation','none', ...
    'aiccSeparation',NaN, ...
    'validationImprovementOverAlternative',NaN, ...
    'preSamplingDecision','unavailable', ...
    'uniformInputDecision','unavailable', ...
    'inputSamplingStable',false, ...
    'linearMaximum',struct(), ...
    'exponentialMaximum',struct());
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
    label = sprintf('%d x %d stored nodes', ...
        payload.config.grid.nx,payload.config.grid.ny);
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

function same = same_coordinate(first,second)
first = first(:);
second = second(:);
if numel(first) ~= numel(second)
    same = false;
    return;
end
tolerance = 64*eps(max(1,max(abs([first;second]))));
same = all(abs(first-second) <= tolerance);
end

function require_fields(value,names,context)
missing = names(~isfield(value,names));
if ~isempty(missing)
    error('ipm:CanonicalMaximumModelsFields', ...
        'The %s value is missing fields: %s.', ...
        context,strjoin(missing,', '));
end
end

function finish_axes(axesHandle,xLabel,yLabel,titleText)
hold(axesHandle,'off');
grid(axesHandle,'on');
box(axesHandle,'on');
set(axesHandle,'Color','w','XColor',[0.15,0.15,0.15], ...
    'YColor',[0.15,0.15,0.15],'GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.35,'LineWidth',0.8);
xlabel(axesHandle,xLabel,'Interpreter','tex');
ylabel(axesHandle,yLabel,'Interpreter','tex');
title(axesHandle,titleText,'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);
end

function style_legend(legendHandle)
set(legendHandle,'Color','w','TextColor',[0.12,0.12,0.12], ...
    'EdgeColor',[0.55,0.55,0.55],'Interpreter','none', ...
    'FontSize',8);
end

function outputFile = export_figure(figureHandle,outputFile)
[outputDirectory,~,extension] = fileparts(outputFile);
if isempty(extension)
    outputFile = [outputFile,'.png'];
elseif ~strcmpi(extension,'.png')
    close(figureHandle);
    error('ipm:CanonicalMaximumModelsPlotFormat', ...
        'The output file must use the PNG extension.');
end
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    [created,message] = mkdir(outputDirectory);
    if ~created
        close(figureHandle);
        error('ipm:CanonicalMaximumModelsPlotDirectory', ...
            'Could not create output directory: %s',message);
    end
end
drawnow;
exportgraphics(figureHandle,outputFile,'Resolution',300, ...
    'BackgroundColor','white');
end

function fileName = analysis_file_name(outputFile)
[directory,name] = fileparts(outputFile);
fileName = fullfile(directory,[name,'_analysis.mat']);
end

function validate_text(value,label)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('ipm:CanonicalMaximumModelsText', ...
        'The %s must be a character vector or scalar string.',label);
end
end
