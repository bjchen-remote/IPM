function [figureHandle,analysis] = ipm_analyze_growth_campaign( ...
        source,outputFile,targetMaximum)
%IPM_ANALYZE_GROWTH_CAMPAIGN Analyze trusted rho_x1 and gradient growth.
%   [FIGUREHANDLE,ANALYSIS] = IPM_ANALYZE_GROWTH_CAMPAIGN(SOURCE) accepts a
%   version-2 solver result or a MAT-file containing a variable named
%   result.  Only the continuous trusted prefix is used.  The analysis
%   compares max|rho_x1| and max|grad rho| on physical time and canonical
%   tau, reusing IPM.DIAGNOSTICS.MAXIMUMGROWTHRATEFIT and
%   IPM.DIAGNOSTICS.CANONICALMAXIMUMGROWTHFIT.
%
%   ... = IPM_ANALYZE_GROWTH_CAMPAIGN(SOURCE,OUTPUTFILE) exports a 300 dpi
%   PNG and saves ANALYSIS to a sibling *_analysis.mat file.  The figure is
%   invisible when an output file is requested, so this call is suitable
%   for MATLAB -batch post-processing and does not advance the PDE.
%
%   ... = IPM_ANALYZE_GROWTH_CAMPAIGN(SOURCE,OUTPUTFILE,TARGETMAXIMUM)
%   changes the forecast threshold from its default value 1e3.  Forecasts
%   are reported separately in physical time and canonical tau.  They are
%   model-dependent extrapolations, not numerical evidence beyond the
%   trusted endpoint.  The quadratic-log fit remains a curvature guard and
%   is never promoted to an asymptotic growth law.

if nargin < 2
    outputFile = '';
end
if nargin < 3 || isempty(targetMaximum)
    targetMaximum = 1e3;
end
validate_text(outputFile,'output file');
validateattributes(targetMaximum,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'targetMaximum');
outputFile = char(outputFile);

result = ipm.output.validate(source);
[rawTrusted,trustChecks] = ipm.output.trustedMask(result);
trustedPrefix = ipm.output.continuousTrustedPrefix(rawTrusted);
trustedPrefix = logical(trustedPrefix(:));
if ~any(trustedPrefix)
    error('ipm:GrowthCampaignTrustedPrefix', ...
        'The result has no records in its continuous trusted prefix.');
end

common = result.history.common;
required = {'physicalTime','canonicalTau','physicalRhoXInf', ...
    'physicalGradInf'};
require_fields(common,required,'history.common');
physicalTime = trusted_series( ...
    common.physicalTime,trustedPrefix,'physical time');
canonicalTau = trusted_series( ...
    common.canonicalTau,trustedPrefix,'canonical tau');
rhoXMaximum = trusted_series( ...
    common.physicalRhoXInf,trustedPrefix,'physical rho_x1 maximum');
gradientMaximum = trusted_series( ...
    common.physicalGradInf,trustedPrefix,'physical gradient maximum');
validate_coordinate(physicalTime,'physical time');
validate_coordinate(canonicalTau,'canonical tau');
validate_maximum(rhoXMaximum,'physical rho_x1 maximum');
validate_maximum(gradientMaximum,'physical gradient maximum');

rhoX = analyze_signal('rho_x1','max |rho_{x_1}|', ...
    physicalTime,canonicalTau,rhoXMaximum,targetMaximum);
gradient = analyze_signal('gradient','max |grad rho|', ...
    physicalTime,canonicalTau,gradientMaximum,targetMaximum);

sourceLabel = source_label(source,result);
caseLabel = optional_text(result.metadata,'caseId',sourceLabel);
figureHandle = make_figure(rhoX,gradient,caseLabel,targetMaximum, ...
    outputFile);

analysis = struct();
analysis.schemaVersion = 1;
analysis.kind = 'trusted_maximum_growth_campaign';
analysis.sourceLabel = sourceLabel;
analysis.caseLabel = caseLabel;
analysis.targetMaximum = targetMaximum;
analysis.trustedPrefixMask = trustedPrefix;
analysis.rawTrustedMask = logical(rawTrusted(:));
analysis.trustChecks = trustChecks;
analysis.trustedRecords = nnz(trustedPrefix);
analysis.totalRecords = numel(trustedPrefix);
analysis.rawTrustedRecords = nnz(rawTrusted);
analysis.trustedThroughPhysicalTime = physicalTime(end);
analysis.trustedThroughCanonicalTau = canonicalTau(end);
analysis.rawTerminalTrusted = logical(rawTrusted(end));
analysis.terminalInContinuousTrustedPrefix = trustedPrefix(end);
analysis.firstRawTrustFailure = find(~rawTrusted,1,'first');
analysis.qualityEndpoint = quality_endpoint( ...
    result,find(trustedPrefix,1,'last'));
analysis.rhoX1 = rhoX;
analysis.fullGradient = gradient;
analysis.outputFile = '';
analysis.analysisFile = '';
analysis.classificationCredible = false;
analysis.interpretation = { ...
    ['All samples and fits stop at the continuous numerical-quality ' ...
    'trusted-prefix endpoint.'], ...
    ['Canonical linear/exponential decisions require stability across ' ...
    'e-fold windows, validation, rate bandwidths, and input resampling.'], ...
    ['Quadratic log(M) is only a curvature guard; its target crossing is ' ...
    'reported as a sensitivity bound, not an asymptotic law.'], ...
    ['Target crossings outside the observed interval are extrapolations ' ...
    'and do not establish grid, time-step, or box convergence.'], ...
    ['One campaign result, even if locally trusted and above the target, ' ...
    'cannot by itself establish credible blow-up.']};

if ~isempty(outputFile)
    outputFile = export_figure(figureHandle,outputFile);
    analysis.outputFile = outputFile;
    analysis.analysisFile = analysis_file_name(outputFile);
    save(analysis.analysisFile,'analysis');
end
print_summary(analysis);
end

function signal = analyze_signal( ...
        name,label,physicalTime,canonicalTau,maximum,targetMaximum)
physicalFit = ipm.diagnostics.maximumGrowthRateFit( ...
    physicalTime,maximum);
canonicalFit = ipm.diagnostics.canonicalMaximumGrowthFit( ...
    canonicalTau,maximum);

signal = struct();
signal.name = name;
signal.label = label;
signal.physicalTime = physicalTime;
signal.canonicalTau = canonicalTau;
signal.maximum = maximum;
signal.initialMaximum = maximum(1);
signal.finalMaximum = maximum(end);
signal.observedMaximum = max(maximum);
signal.growthFactor = maximum(end)/maximum(1);
signal.eFoldCount = log(signal.growthFactor);
signal.physicalGrowthRateFit = physicalFit;
signal.canonicalGrowthFit = canonicalFit;
signal.physicalRateSummary = rate_summary(physicalFit,'physical');
signal.canonicalRateSummary = rate_summary(canonicalFit,'canonical');
signal.windowDependence = window_dependence(canonicalFit);
signal.targetForecast = target_forecast(physicalFit,canonicalFit, ...
    physicalTime,canonicalTau,maximum,targetMaximum);
end

function summary = rate_summary(fit,coordinateKind)
summary = struct('available',false,'coordinateKind',coordinateKind, ...
    'terminalByBandwidth',[],'terminalRange',[NaN,NaN], ...
    'maximumByBandwidth',[],'maximumRange',[NaN,NaN], ...
    'absoluteMaximumRateByBandwidth',[], ...
    'absoluteMaximumRateRange',[NaN,NaN]);
if ~fit.valid
    return;
end
if strcmp(coordinateKind,'physical')
    rates = fit.logarithmicRate;
    terminal = fit.finalRate;
else
    rates = fit.canonicalLogarithmicRate;
    terminal = fit.finalCanonicalRate;
end
maximumByBandwidth = column_maximum(rates);
absoluteRates = rates.*fit.maximum;
absoluteMaximumByBandwidth = column_maximum(absoluteRates);
summary.available = true;
summary.terminalByBandwidth = terminal;
summary.terminalRange = finite_range(terminal);
summary.maximumByBandwidth = maximumByBandwidth;
summary.maximumRange = finite_range(maximumByBandwidth);
summary.absoluteMaximumRateByBandwidth = absoluteMaximumByBandwidth;
summary.absoluteMaximumRateRange = ...
    finite_range(absoluteMaximumByBandwidth);
end

function summary = window_dependence(fit)
summary = struct('available',false,'requestedLogGrowth',[], ...
    'actualLogGrowth',[],'startCanonicalTau',[], ...
    'modelDecision',{{}},'preferredByAicc',{{}}, ...
    'preferredByValidation',{{}},'inputSamplingStable',[], ...
    'curvatureGuardTriggered',[],'consensusModel','none', ...
    'modelConsensus',false);
if ~fit.valid
    return;
end
windows = fit.windows;
summary.available = true;
summary.requestedLogGrowth = [windows.requestedLogGrowth];
summary.actualLogGrowth = [windows.actualLogGrowth];
summary.startCanonicalTau = [windows.startCanonicalTau];
summary.modelDecision = {windows.modelDecision};
summary.preferredByAicc = {windows.preferredByAicc};
summary.preferredByValidation = {windows.preferredByValidation};
summary.inputSamplingStable = [windows.inputSamplingStable];
summary.curvatureGuardTriggered = [windows.curvatureGuardTriggered];
summary.consensusModel = fit.consensusModel;
summary.modelConsensus = fit.modelConsensus;
end

function forecast = target_forecast(physicalFit,canonicalFit, ...
        physicalTime,canonicalTau,maximum,targetMaximum)
[observedPhysicalTime,observedIndex] = observed_crossing( ...
    physicalTime,maximum,targetMaximum);
[observedCanonicalTau,~] = observed_crossing( ...
    canonicalTau,maximum,targetMaximum);
reached = isfinite(observedPhysicalTime);

forecast = struct();
forecast.targetMaximum = targetMaximum;
forecast.reachedInTrustedData = reached;
forecast.observedCrossingIndex = observedIndex;
forecast.observedPhysicalTime = observedPhysicalTime;
forecast.observedCanonicalTau = observedCanonicalTau;
forecast.extrapolationRequired = ~reached;
forecast.physicalByWindow = physical_predictions( ...
    physicalFit,targetMaximum,physicalTime(end),reached);
forecast.canonicalByWindow = canonical_predictions( ...
    canonicalFit,targetMaximum,canonicalTau(end),reached);
forecast.physicalPredictionRanges = struct( ...
    'exponential',prediction_range( ...
        forecast.physicalByWindow,'exponential'), ...
    'finiteTimePower',prediction_range( ...
        forecast.physicalByWindow,'finiteTimePower'), ...
    'doubleExponential',prediction_range( ...
        forecast.physicalByWindow,'doubleExponential'));
forecast.canonicalPredictionRanges = struct( ...
    'linearMaximum',prediction_range( ...
        forecast.canonicalByWindow,'linearMaximum'), ...
    'exponentialMaximum',prediction_range( ...
        forecast.canonicalByWindow,'exponentialMaximum'), ...
    'quadraticLogGuard',prediction_range( ...
        forecast.canonicalByWindow,'quadraticLogGuard'));
forecast.warning = [ ...
    'Finite future coordinates are model crossings beyond the trusted ' ...
    'endpoint; the quadratic-log crossing is a guard-only sensitivity.'];
end

function predictions = physical_predictions( ...
        fit,targetMaximum,endpoint,reached)
template = struct('requestedLogGrowth',NaN,'startTime',NaN, ...
    'endTime',NaN,'modelDecision','indistinguishable', ...
    'exponential',empty_prediction(), ...
    'finiteTimePower',empty_prediction(), ...
    'doubleExponential',empty_prediction());
if ~fit.valid
    predictions = repmat(template,0,1);
    return;
end
windows = fit.windows;
predictions = repmat(template,numel(windows),1);
for index = 1:numel(windows)
    window = windows(index);
    predictions(index).requestedLogGrowth = window.requestedLogGrowth;
    predictions(index).startTime = window.startTime;
    predictions(index).endTime = window.endTime;
    predictions(index).modelDecision = window.modelDecision;
    if reached
        continue;
    end
    predictions(index).exponential = model_crossing( ...
        window.exponential,targetMaximum,endpoint,'physical');
    predictions(index).finiteTimePower = model_crossing( ...
        window.powerLaw,targetMaximum,endpoint,'physical');
    predictions(index).doubleExponential = model_crossing( ...
        window.doubleExponential,targetMaximum,endpoint,'physical');
end
end

function predictions = canonical_predictions( ...
        fit,targetMaximum,endpoint,reached)
template = struct('requestedLogGrowth',NaN, ...
    'startCanonicalTau',NaN,'endCanonicalTau',NaN, ...
    'modelDecision','indistinguishable', ...
    'linearMaximum',empty_prediction(), ...
    'exponentialMaximum',empty_prediction(), ...
    'quadraticLogGuard',empty_prediction());
if ~fit.valid
    predictions = repmat(template,0,1);
    return;
end
windows = fit.windows;
predictions = repmat(template,numel(windows),1);
for index = 1:numel(windows)
    window = windows(index);
    predictions(index).requestedLogGrowth = window.requestedLogGrowth;
    predictions(index).startCanonicalTau = window.startCanonicalTau;
    predictions(index).endCanonicalTau = window.endCanonicalTau;
    predictions(index).modelDecision = window.modelDecision;
    if reached
        continue;
    end
    predictions(index).linearMaximum = model_crossing( ...
        window.linearMaximum,targetMaximum,endpoint,'canonical');
    predictions(index).exponentialMaximum = model_crossing( ...
        window.exponentialMaximum,targetMaximum,endpoint,'canonical');
    predictions(index).quadraticLogGuard = model_crossing( ...
        window.quadraticLogGuard,targetMaximum,endpoint,'canonical');
end
end

function prediction = model_crossing( ...
        model,targetMaximum,endpoint,coordinateKind)
prediction = empty_prediction();
if ~isstruct(model) || ~isfield(model,'valid') || ~model.valid
    prediction.reason = 'model_invalid';
    return;
end
logTarget = log(targetMaximum);
coordinate = NaN;
switch model.name
    case 'exponential'
        coordinate = model.referenceTime + ...
            (logTarget-model.parameters.logM0)/model.parameters.lambda;
    case 'finite_time_power'
        distance = exp((model.parameters.logAmplitude-logTarget) / ...
            model.parameters.p);
        coordinate = model.parameters.T-distance;
    case 'double_exponential'
        kappa = model.parameters.kappa;
        argument = 1+kappa*(logTarget-model.parameters.logM0) / ...
            model.parameters.gamma0;
        if argument > 0
            coordinate = model.referenceTime+log(argument)/kappa;
        end
    case 'linear_maximum'
        logRatio = logTarget-model.parameters.logA;
        if logRatio <= log(realmax)
            coordinate = model.referenceCanonicalTau + ...
                expm1(logRatio)/model.parameters.k;
        end
    case 'exponential_maximum'
        coordinate = model.referenceCanonicalTau + ...
            (logTarget-model.parameters.logA)/model.parameters.lambda;
    case 'quadratic_log_guard'
        coordinate = quadratic_guard_crossing( ...
            model,logTarget,endpoint);
end
if ~isfinite(coordinate)
    prediction.reason = 'no_finite_crossing';
    return;
end
tolerance = 64*eps(max(1,abs(endpoint)));
if coordinate < endpoint-tolerance
    prediction.reason = 'crossing_not_beyond_trusted_endpoint';
    return;
end
prediction.valid = true;
prediction.coordinateKind = coordinateKind;
prediction.coordinate = coordinate;
prediction.deltaFromTrustedEndpoint = coordinate-endpoint;
prediction.reason = 'future_model_crossing';
end

function coordinate = quadratic_guard_crossing(model,logTarget,endpoint)
parameters = model.parameters;
a = parameters.quadraticCoefficient;
b = parameters.lambda;
c = parameters.logA-logTarget;
scale = max([1,abs(a),abs(b),abs(c)]);
if abs(a) <= 128*eps(scale)
    if b > 0
        coordinate = model.referenceCanonicalTau-c/b;
    else
        coordinate = NaN;
    end
    return;
end
discriminant = b^2-4*a*c;
if discriminant < 0
    coordinate = NaN;
    return;
end
rootsX = [(-b-sqrt(discriminant))/(2*a), ...
    (-b+sqrt(discriminant))/(2*a)];
rootsTau = model.referenceCanonicalTau+rootsX;
rates = b+2*a*rootsX;
tolerance = 64*eps(max(1,abs(endpoint)));
valid = isfinite(rootsTau) & rootsTau >= endpoint-tolerance & rates > 0;
if any(valid)
    coordinate = min(rootsTau(valid));
else
    coordinate = NaN;
end
end

function prediction = empty_prediction()
prediction = struct('valid',false,'coordinateKind','', ...
    'coordinate',NaN,'deltaFromTrustedEndpoint',NaN, ...
    'reason','unavailable');
end

function range = prediction_range(predictions,fieldName)
range = struct('available',false,'minimum',NaN,'maximum',NaN, ...
    'span',NaN,'count',0,'requestedLogGrowth',[],'coordinates',[]);
if isempty(predictions)
    return;
end
valid = arrayfun(@(entry)entry.(fieldName).valid,predictions);
if ~any(valid)
    return;
end
coordinates = arrayfun( ...
    @(entry)entry.(fieldName).coordinate,predictions(valid));
range.available = true;
range.minimum = min(coordinates);
range.maximum = max(coordinates);
range.span = range.maximum-range.minimum;
range.count = numel(coordinates);
range.requestedLogGrowth = ...
    [predictions(valid).requestedLogGrowth];
range.coordinates = coordinates;
end

function [coordinate,index] = observed_crossing( ...
        coordinates,maximum,targetMaximum)
index = find(maximum >= targetMaximum,1,'first');
if isempty(index)
    coordinate = NaN;
    index = NaN;
    return;
end
if index == 1
    coordinate = coordinates(1);
    return;
end
leftValue = log(maximum(index-1));
rightValue = log(maximum(index));
if rightValue == leftValue
    fraction = 1;
else
    fraction = (log(targetMaximum)-leftValue)/(rightValue-leftValue);
end
fraction = min(max(fraction,0),1);
coordinate = coordinates(index-1) + ...
    fraction*(coordinates(index)-coordinates(index-1));
end

function figureHandle = make_figure( ...
        rhoX,gradient,caseLabel,targetMaximum,outputFile)
if isempty(outputFile) && usejava('desktop')
    visibility = 'on';
else
    visibility = 'off';
end
figureHandle = figure('Color','w','Visible',visibility, ...
    'Position',[60,50,1580,1120], ...
    'Name','Trusted maximum-growth campaign');
layout = tiledlayout(figureHandle,3,2,'TileSpacing','compact', ...
    'Padding','compact');

rhoColor = [0,0.4470,0.7410];
gradientColor = [0.8500,0.3250,0.0980];
draw_maxima_panel(nexttile(layout,1),rhoX,gradient, ...
    'physical',rhoColor,gradientColor,targetMaximum);
draw_maxima_panel(nexttile(layout,2),rhoX,gradient, ...
    'canonical',rhoColor,gradientColor,targetMaximum);
draw_tau_model_panel(nexttile(layout,3),rhoX,rhoColor,targetMaximum);
draw_tau_model_panel(nexttile(layout,4),gradient, ...
    gradientColor,targetMaximum);
draw_rate_panel(nexttile(layout,5),rhoX,gradient, ...
    'physical',rhoColor,gradientColor);
draw_rate_panel(nexttile(layout,6),rhoX,gradient, ...
    'canonical',rhoColor,gradientColor);

endpointLine = sprintf([ ...
    'trusted endpoint: t = %.9g, canonical tau = %.9g; ' ...
    'target M = %.6g'],rhoX.physicalTime(end), ...
    rhoX.canonicalTau(end),targetMaximum);
sgtitle(layout,{caseLabel,endpointLine},'Interpreter','none', ...
    'FontWeight','bold','Color',[0.12,0.12,0.12]);
end

function draw_maxima_panel(axesHandle,rhoX,gradient,coordinateKind, ...
        rhoColor,gradientColor,targetMaximum)
if strcmp(coordinateKind,'physical')
    coordinate = rhoX.physicalTime;
    xLabel = 'Physical time t';
    titleText = 'Trusted maxima versus physical time';
else
    coordinate = rhoX.canonicalTau;
    xLabel = 'Canonical time tau';
    titleText = 'Trusted maxima versus canonical time';
end
semilogy(axesHandle,coordinate,rhoX.maximum,'Color',rhoColor, ...
    'LineWidth',1.8,'DisplayName','max |rho_{x_1}|');
hold(axesHandle,'on');
semilogy(axesHandle,coordinate,gradient.maximum, ...
    'Color',gradientColor,'LineWidth',1.8, ...
    'DisplayName','max |grad rho|');
yline(axesHandle,targetMaximum,':','target', ...
    'Color',[0.25,0.25,0.25],'HandleVisibility','off');
hold(axesHandle,'off');
finish_axes(axesHandle,xLabel,'Physical maximum',titleText);
legendHandle = legend(axesHandle,'Location','northwest');
style_legend(legendHandle);
end

function draw_tau_model_panel( ...
        axesHandle,signal,signalColor,targetMaximum)
fit = signal.canonicalGrowthFit;
dataHandle = semilogy(axesHandle,signal.canonicalTau,signal.maximum, ...
    'Color',signalColor,'LineWidth',1.8,'DisplayName','trusted data');
hold(axesHandle,'on');
window = largest_window(fit.windows);
handles = gobjects(1,4);
labels = cell(1,4);
handles(1) = dataHandle;
labels{1} = 'trusted data';
handleCount = 1;
if ~isempty(window)
    modelTau = linspace(window.startCanonicalTau, ...
        window.endCanonicalTau,240)';
    models = {window.linearMaximum,window.exponentialMaximum, ...
        window.quadraticLogGuard};
    styles = {':','--','-.'};
    colors = [0.4660,0.6740,0.1880;0.6350,0.0780,0.1840; ...
        0.4940,0.1840,0.5560];
    modelLabels = {'linear maximum','exponential maximum', ...
        'quadratic-log guard'};
    for index = 1:numel(models)
        model = models{index};
        if ~model.valid
            continue;
        end
        logPrediction = canonical_model_log_value(model,modelTau);
        modelHandle = semilogy(axesHandle,modelTau,exp(logPrediction), ...
            'Color',colors(index,:),'LineStyle',styles{index}, ...
            'LineWidth',1.55,'DisplayName',modelLabels{index});
        handleCount = handleCount+1;
        handles(handleCount) = modelHandle;
        labels{handleCount} = modelLabels{index};
    end
    xline(axesHandle,window.startCanonicalTau,':', ...
        sprintf('%.3g e-fold window',window.requestedLogGrowth), ...
        'HandleVisibility','off','LabelVerticalAlignment','bottom');
end
hold(axesHandle,'off');
predictionText = selected_prediction_text( ...
    signal.targetForecast,window,targetMaximum);
titleText = {sprintf('%s: canonical models',signal.label),predictionText};
finish_axes(axesHandle,'Canonical time tau','Physical maximum',titleText);
legendHandle = legend(axesHandle,handles(1:handleCount), ...
    labels(1:handleCount),'Location','northwest');
style_legend(legendHandle);
end

function draw_rate_panel(axesHandle,rhoX,gradient,coordinateKind, ...
        rhoColor,gradientColor)
hold(axesHandle,'on');
if strcmp(coordinateKind,'physical')
    rhoFit = rhoX.physicalGrowthRateFit;
    gradientFit = gradient.physicalGrowthRateFit;
    coordinate = rhoX.physicalTime;
    rhoRates = rhoFit.logarithmicRate;
    gradientRates = gradientFit.logarithmicRate;
    xLabel = 'Physical time t';
    yLabel = 'd log(M) / dt';
    titleText = rate_title(rhoX.physicalRateSummary, ...
        gradient.physicalRateSummary,'Physical-time');
else
    rhoFit = rhoX.canonicalGrowthFit;
    gradientFit = gradient.canonicalGrowthFit;
    coordinate = rhoX.canonicalTau;
    rhoRates = rhoFit.canonicalLogarithmicRate;
    gradientRates = gradientFit.canonicalLogarithmicRate;
    xLabel = 'Canonical time tau';
    yLabel = 'd log(M) / d tau';
    titleText = rate_title(rhoX.canonicalRateSummary, ...
        gradient.canonicalRateSummary,'Canonical-time');
end
handles = gobjects(0);
labels = cell(0,1);
if ~isempty(rhoRates)
    handles(end+1) = draw_rate_envelope( ...
        axesHandle,coordinate,rhoRates,rhoColor);
    labels{end+1,1} = 'max |rho_{x_1}|';
end
if ~isempty(gradientRates)
    handles(end+1) = draw_rate_envelope( ...
        axesHandle,coordinate,gradientRates,gradientColor);
    labels{end+1,1} = 'max |grad rho|';
end
hold(axesHandle,'off');
finish_axes(axesHandle,xLabel,yLabel,titleText);
if ~isempty(handles)
    legendHandle = legend(axesHandle,handles,labels,'Location','northwest');
    style_legend(legendHandle);
end
end

function lineHandle = draw_rate_envelope( ...
        axesHandle,coordinate,rates,color)
[lower,upper] = row_range(rates);
valid = isfinite(lower) & isfinite(upper);
if any(valid)
    fill(axesHandle,[coordinate(valid);flipud(coordinate(valid))], ...
        [lower(valid);flipud(upper(valid))],color, ...
        'FaceAlpha',0.13,'EdgeColor','none','HandleVisibility','off');
end
centralIndex = min(2,size(rates,2));
lineHandle = plot(axesHandle,coordinate,rates(:,centralIndex), ...
    'Color',color,'LineWidth',1.7);
end

function titleText = rate_title(rhoSummary,gradientSummary,prefix)
if ~rhoSummary.available || ~gradientSummary.available
    titleText = [prefix,' logarithmic growth rate unavailable'];
    return;
end
titleText = {sprintf('%s logarithmic growth rates',prefix), ...
    sprintf(['endpoint rho_{x_1} [%.4g, %.4g], ' ...
    'grad [%.4g, %.4g]'],rhoSummary.terminalRange(1), ...
    rhoSummary.terminalRange(2),gradientSummary.terminalRange(1), ...
    gradientSummary.terminalRange(2))};
end

function value = selected_prediction_text(forecast,window,targetMaximum)
if forecast.reachedInTrustedData
    value = sprintf('M = %.3g observed at tau = %.6g', ...
        targetMaximum,forecast.observedCanonicalTau);
    return;
end
if isempty(window)
    value = sprintf('M = %.3g forecast unavailable',targetMaximum);
    return;
end
predictions = forecast.canonicalByWindow;
growth = [predictions.requestedLogGrowth];
[~,index] = min(abs(growth-window.requestedLogGrowth));
selected = predictions(index);
coordinates = [selected.linearMaximum.coordinate, ...
    selected.exponentialMaximum.coordinate, ...
    selected.quadraticLogGuard.coordinate];
value = sprintf('tau_{%.3g}: linear / exp / guard = %s / %s / %s', ...
    targetMaximum,numeric_label(coordinates(1)), ...
    numeric_label(coordinates(2)),numeric_label(coordinates(3)));
end

function values = canonical_model_log_value(model,tau)
x = tau-model.referenceCanonicalTau;
switch model.name
    case 'linear_maximum'
        values = model.parameters.logA+log1p(model.parameters.k*x);
    case 'exponential_maximum'
        values = model.parameters.logA+model.parameters.lambda*x;
    case 'quadratic_log_guard'
        values = model.parameters.logA+model.parameters.lambda*x + ...
            model.parameters.quadraticCoefficient*x.^2;
    otherwise
        values = NaN(size(tau));
end
end

function window = largest_window(windows)
if isempty(windows)
    window = [];
    return;
end
[~,index] = max([windows.requestedLogGrowth]);
window = windows(index);
end

function endpoint = quality_endpoint(result,index)
common = result.history.common;
mesh = result.history.mesh;
endpoint = struct( ...
    'safetyFactor',optional_index(mesh,'safetyFactor',index), ...
    'coreGridPoints',optional_index(mesh,'coreGridPoints',index), ...
    'verticalCoreGridPoints', ...
        optional_index(mesh,'verticalCoreGridPoints',index), ...
    'maximumCellRatioX',optional_index(mesh,'maximumCellRatioX',index), ...
    'maximumCellRatioY',optional_index(mesh,'maximumCellRatioY',index), ...
    'physicalMassDrift', ...
        optional_index(common,'physicalMassDrift',index), ...
    'physicalRangeViolation', ...
        optional_index(common,'physicalRangeViolation',index), ...
    'farBoundaryVelocityRatio', ...
        optional_index(common,'farBoundaryVelocityRatio',index), ...
    'farBoundarySourceRatio', ...
        optional_index(common,'farBoundarySourceRatio',index));
end

function value = optional_index(source,name,index)
value = NaN;
if isfield(source,name)
    values = source.(name);
    if isnumeric(values) && isvector(values) && numel(values) >= index
        value = values(index);
    end
end
end

function print_summary(analysis)
signals = {analysis.rhoX1,analysis.fullGradient};
fprintf(['Growth campaign trusted through t=%.9g, tau=%.9g; ' ...
    '%d/%d continuous-prefix records.\n'], ...
    analysis.trustedThroughPhysicalTime, ...
    analysis.trustedThroughCanonicalTau,analysis.trustedRecords, ...
    analysis.totalRecords);
for index = 1:numel(signals)
    signal = signals{index};
    fprintf(['%s: final %.9g, observed max %.9g, growth %.6g ' ...
        '(%.6g e-fold).\n'],signal.name,signal.finalMaximum, ...
        signal.observedMaximum,signal.growthFactor,signal.eFoldCount);
    if signal.physicalRateSummary.available
        fprintf(['  max log-rate physical [%.6g, %.6g], ' ...
            'endpoint [%.6g, %.6g].\n'], ...
            signal.physicalRateSummary.maximumRange(1), ...
            signal.physicalRateSummary.maximumRange(2), ...
            signal.physicalRateSummary.terminalRange(1), ...
            signal.physicalRateSummary.terminalRange(2));
    end
    if signal.canonicalRateSummary.available
        fprintf(['  max log-rate canonical [%.6g, %.6g], ' ...
            'endpoint [%.6g, %.6g], windows: %s.\n'], ...
            signal.canonicalRateSummary.maximumRange(1), ...
            signal.canonicalRateSummary.maximumRange(2), ...
            signal.canonicalRateSummary.terminalRange(1), ...
            signal.canonicalRateSummary.terminalRange(2), ...
            strjoin(signal.windowDependence.modelDecision,', '));
    end
    print_target_forecast(signal.targetForecast);
end
fprintf(['Classification credible = 0 (single-run, model-dependent ' ...
    'campaign diagnostic).\n']);
end

function print_target_forecast(forecast)
if forecast.reachedInTrustedData
    fprintf('  target %.6g observed at t=%.9g, tau=%.9g.\n', ...
        forecast.targetMaximum,forecast.observedPhysicalTime, ...
        forecast.observedCanonicalTau);
    return;
end
physical = forecast.physicalPredictionRanges;
canonical = forecast.canonicalPredictionRanges;
fprintf('  target %.6g physical-model ranges exp/power/double = %s / %s / %s.\n', ...
    forecast.targetMaximum,range_label(physical.exponential), ...
    range_label(physical.finiteTimePower), ...
    range_label(physical.doubleExponential));
fprintf('  target %.6g tau-model ranges linear/exp/guard = %s / %s / %s.\n', ...
    forecast.targetMaximum,range_label(canonical.linearMaximum), ...
    range_label(canonical.exponentialMaximum), ...
    range_label(canonical.quadraticLogGuard));
end

function value = range_label(range)
if range.available
    value = sprintf('[%.6g, %.6g]',range.minimum,range.maximum);
else
    value = 'unavailable';
end
end

function value = numeric_label(number)
if isfinite(number)
    value = sprintf('%.6g',number);
else
    value = 'n/a';
end
end

function maxima = column_maximum(values)
maxima = NaN(1,size(values,2));
for index = 1:size(values,2)
    finite = values(isfinite(values(:,index)),index);
    if ~isempty(finite)
        maxima(index) = max(finite);
    end
end
end

function range = finite_range(values)
finite = values(isfinite(values));
if isempty(finite)
    range = [NaN,NaN];
else
    range = [min(finite),max(finite)];
end
end

function [lower,upper] = row_range(values)
lower = NaN(size(values,1),1);
upper = lower;
for index = 1:size(values,1)
    finite = values(index,isfinite(values(index,:)));
    if ~isempty(finite)
        lower(index) = min(finite);
        upper(index) = max(finite);
    end
end
end

function values = trusted_series(values,trustedPrefix,label)
if ~isnumeric(values) || ~isreal(values) || ~isvector(values) || ...
        numel(values) ~= numel(trustedPrefix)
    error('ipm:GrowthCampaignSeries', ...
        'The %s series is not aligned with the trusted mask.',label);
end
values = values(:);
values = values(trustedPrefix);
end

function validate_coordinate(values,label)
if any(~isfinite(values)) || any(diff(values) <= 0)
    error('ipm:GrowthCampaignCoordinate', ...
        'Trusted %s values must be finite and strictly increasing.',label);
end
end

function validate_maximum(values,label)
if any(~isfinite(values)) || any(values <= 0)
    error('ipm:GrowthCampaignMaximum', ...
        'Trusted %s values must be finite and positive.',label);
end
end

function label = source_label(source,result)
if ischar(source) || (isstring(source) && isscalar(source))
    [~,label] = fileparts(char(source));
else
    label = optional_text(result.metadata,'caseId','solver result');
end
end

function value = optional_text(source,name,fallback)
value = fallback;
if isstruct(source) && isscalar(source) && isfield(source,name)
    candidate = source.(name);
    if ischar(candidate) || (isstring(candidate) && isscalar(candidate))
        value = char(candidate);
    end
end
end

function require_fields(value,names,context)
for index = 1:numel(names)
    if ~isfield(value,names{index})
        error('ipm:GrowthCampaignField', ...
            '%s is missing required field "%s".',context,names{index});
    end
end
end

function finish_axes(axesHandle,xLabel,yLabel,titleText)
grid(axesHandle,'on');
box(axesHandle,'on');
xlabel(axesHandle,xLabel,'Interpreter','tex');
ylabel(axesHandle,yLabel,'Interpreter','tex');
title(axesHandle,titleText,'Interpreter','tex', ...
    'Color',[0.12,0.12,0.12]);
set(axesHandle,'Color','w','XColor',[0.15,0.15,0.15], ...
    'YColor',[0.15,0.15,0.15],'GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.35,'LineWidth',0.8);
end

function style_legend(legendHandle)
set(legendHandle,'Color','w','TextColor',[0.12,0.12,0.12], ...
    'EdgeColor',[0.55,0.55,0.55],'Interpreter','tex');
end

function outputFile = export_figure(figureHandle,outputFile)
[outputDirectory,~,extension] = fileparts(outputFile);
if isempty(extension)
    outputFile = [outputFile,'.png'];
elseif ~strcmpi(extension,'.png')
    close(figureHandle);
    error('ipm:GrowthCampaignPlotFormat', ...
        'The output file must use the PNG extension.');
end
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    [created,message] = mkdir(outputDirectory);
    if ~created
        close(figureHandle);
        error('ipm:GrowthCampaignPlotDirectory', ...
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
    error('ipm:GrowthCampaignText', ...
        'The %s must be a character vector or scalar string.',label);
end
end
