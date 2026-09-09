function [figureHandle,plotData] = ipm_plot_maximum_growth(source,outputFile)
%IPM_PLOT_MAXIMUM_GROWTH Plot trusted maximum-growth diagnostics.
%   FIGUREHANDLE = IPM_PLOT_MAXIMUM_GROWTH(SOURCE) accepts a version-2
%   solver result, a fourth-order blow-up assessment, or a MAT-file path
%   containing either variable.  All plotted samples come from the
%   continuous trusted prefix.  The logarithmic-rate curves and model fits
%   are read from IPM.DIAGNOSTICS.MAXIMUMGROWTHRATEFIT output fields.
%
%   IPM_PLOT_MAXIMUM_GROWTH(SOURCE,OUTPUTFILE) exports a 300 dpi PNG.  The
%   figure is created invisibly when an output file is requested, so this
%   form is suitable for MATLAB -batch runs.
%
%   [FIGUREHANDLE,PLOTDATA] = ... also returns the two rate-fit structures,
%   selected trailing model windows, case label, and exported file path.

if nargin < 2
    outputFile = '';
end
if ~(ischar(outputFile) || (isstring(outputFile) && isscalar(outputFile)))
    error('ipm:MaximumGrowthPlotOutput', ...
        'The optional output file must be a character vector or scalar string.');
end
outputFile = char(outputFile);

[gradientFit,wallFit,caseLabel,sourceKind] = growth_fits(source);
validate_rate_fit(gradientFit,'gradient');
validate_rate_fit(wallFit,'wall peak');

if isempty(outputFile) && usejava('desktop')
    visibility = 'on';
else
    visibility = 'off';
end
figureHandle = figure('Color','w','Visible',visibility, ...
    'Position',[80,80,1500,1050],'Name','Maximum growth diagnostics');
layout = tiledlayout(figureHandle,3,2,'TileSpacing','compact', ...
    'Padding','compact');

gradientColor = [0,0.4470,0.7410];
wallColor = [0.8500,0.3250,0.0980];
gradientWindow = draw_signal_column(layout,[1,3,5],gradientFit, ...
    gradientColor,'Physical maximum gradient');
wallWindow = draw_signal_column(layout,[2,4,6],wallFit,wallColor, ...
    'Physical wall peak');

titleLine = sprintf('%s | continuous trusted prefix through t = %.9g', ...
    caseLabel,min(gradientFit.time(end),wallFit.time(end)));
rangeLine = sprintf(['gradient %.3f e-folds (%d points); ' ...
    'wall peak %.3f e-folds (%d points)'], ...
    gradientFit.eFoldCount,gradientFit.points, ...
    wallFit.eFoldCount,wallFit.points);
sgtitle(layout,{titleLine,rangeLine},'Interpreter','none', ...
    'FontWeight','bold','Color',[0.12,0.12,0.12]);

if ~isempty(outputFile)
    [outputDirectory,~,extension] = fileparts(outputFile);
    if isempty(extension)
        outputFile = [outputFile,'.png'];
    elseif ~strcmpi(extension,'.png')
        close(figureHandle);
        error('ipm:MaximumGrowthPlotFormat', ...
            'The output file must use the PNG extension.');
    end
    if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
        [created,message] = mkdir(outputDirectory);
        if ~created
            close(figureHandle);
            error('ipm:MaximumGrowthPlotDirectory', ...
                'Could not create output directory: %s',message);
        end
    end
    drawnow;
    exportgraphics(figureHandle,outputFile,'Resolution',300, ...
        'BackgroundColor','white');
end

plotData = struct();
plotData.caseLabel = caseLabel;
plotData.sourceKind = sourceKind;
plotData.trustedThroughTime = min( ...
    gradientFit.time(end),wallFit.time(end));
plotData.gradientGrowthRateFit = gradientFit;
plotData.wallPeakGrowthRateFit = wallFit;
plotData.gradientModelWindow = gradientWindow;
plotData.wallModelWindow = wallWindow;
plotData.outputFile = outputFile;
end

function selectedWindow = draw_signal_column( ...
        layout,tileNumbers,fit,signalColor,signalName)
time = fit.time(:);
maximum = fit.maximum(:);

growthAxes = nexttile(layout,tileNumbers(1));
plot(growthAxes,time,maximum,'Color',signalColor,'LineWidth',1.7);
style_axes(growthAxes);
grid(growthAxes,'on');
box(growthAxes,'on');
xlabel(growthAxes,'Physical time');
ylabel(growthAxes,'Physical maximum');
title(growthAxes,signalName,'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);

selectedWindow = largest_window(fit.windows);
logAxes = nexttile(layout,tileNumbers(2));
dataLogGrowth = fit.logMaximum(:)-fit.logMaximum(1);
dataHandle = plot(logAxes,time,dataLogGrowth,'Color',signalColor, ...
    'LineWidth',1.8,'DisplayName','trusted data');
style_axes(logAxes);
hold(logAxes,'on');
legendHandles = dataHandle;
legendLabels = {'trusted data'};
if ~isempty(selectedWindow)
    modelTime = linspace(selectedWindow.startTime, ...
        selectedWindow.endTime,240)';
    models = {selectedWindow.exponential,selectedWindow.powerLaw, ...
        selectedWindow.doubleExponential};
    modelColors = [0.25,0.25,0.25;0.9290,0.6940,0.1250; ...
        0.4940,0.1840,0.5560];
    modelStyles = {':','--','-.'};
    modelLabels = {'exponential','finite-time power','double exponential'};
    for modelIndex = 1:numel(models)
        model = models{modelIndex};
        if ~model.valid
            continue;
        end
        prediction = model_log_value(model,modelTime)-fit.logMaximum(1);
        label = modelLabels{modelIndex};
        if isfield(model,'identifiable') && ~model.identifiable
            label = [label,' (not identifiable)']; %#ok<AGROW>
        end
        modelHandle = plot(logAxes,modelTime,prediction, ...
            'Color',modelColors(modelIndex,:), ...
            'LineStyle',modelStyles{modelIndex},'LineWidth',1.45, ...
            'DisplayName',label);
        legendHandles(end+1) = modelHandle; %#ok<AGROW>
        legendLabels{end+1} = label; %#ok<AGROW>
    end
    xline(logAxes,selectedWindow.startTime,':', ...
        sprintf('%.3g e-fold window',selectedWindow.requestedLogGrowth), ...
        'HandleVisibility','off','LabelVerticalAlignment','bottom');
end
hold(logAxes,'off');
grid(logAxes,'on');
box(logAxes,'on');
xlabel(logAxes,'Physical time');
ylabel(logAxes,'log(M/M_0)');
title(logAxes,model_window_title(selectedWindow),'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);
logLegend = legend(logAxes,legendHandles,legendLabels, ...
    'Location','northwest');
style_legend(logLegend);

rateAxes = nexttile(layout,tileNumbers(3));
style_axes(rateAxes);
hold(rateAxes,'on');
rateColors = [0.3010,0.7450,0.9330;signalColor;0.6350,0.0780,0.1840];
rateHandles = gobjects(0);
rateLabels = cell(0,1);
for bandwidthIndex = 1:size(fit.logarithmicRate,2)
    rateHandle = plot(rateAxes,time, ...
        fit.logarithmicRate(:,bandwidthIndex), ...
        'Color',rateColors(bandwidthIndex,:), ...
        'LineWidth',1.35);
    rateHandles(end+1) = rateHandle; %#ok<AGROW>
    rateLabels{end+1,1} = sprintf('h = %.4g', ...
        fit.bandwidths(bandwidthIndex)); %#ok<AGROW>
end
decisionModel = decision_model(selectedWindow);
if ~isempty(decisionModel)
    modelTime = linspace(selectedWindow.startTime, ...
        selectedWindow.endTime,240)';
    rateHandle = plot(rateAxes,modelTime,model_rate(decisionModel,modelTime), ...
        'k--','LineWidth',1.8);
    rateHandles(end+1) = rateHandle;
    rateLabels{end+1,1} = sprintf('%.3g e-fold decision: %s', ...
        selectedWindow.requestedLogGrowth, ...
        short_model_name(decisionModel.name));
end
hold(rateAxes,'off');
grid(rateAxes,'on');
box(rateAxes,'on');
xlabel(rateAxes,'Physical time');
ylabel(rateAxes,'gamma = d log(M) / dt');
title(rateAxes,decision_summary(fit.windows),'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);
rateLegend = legend(rateAxes,rateHandles,rateLabels, ...
    'Location','northwest');
style_legend(rateLegend);
end

function style_axes(axesHandle)
set(axesHandle,'Color','w','XColor',[0.15,0.15,0.15], ...
    'YColor',[0.15,0.15,0.15],'GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.35,'LineWidth',0.8);
end

function style_legend(legendHandle)
set(legendHandle,'Color','w','TextColor',[0.12,0.12,0.12], ...
    'EdgeColor',[0.55,0.55,0.55]);
end

function [gradientFit,wallFit,caseLabel,sourceKind] = growth_fits(source)
[payload,sourceKind] = read_source(source);
if is_assessment(payload)
    gradientFit = payload.gradientGrowthRateFit;
    wallFit = payload.wallPeakGrowthRateFit;
    caseLabel = optional_case_label(payload,'caseId','assessment');
    return;
end

result = ipm.output.validate(payload);
[trusted,~] = ipm.output.trustedMask(result);
trustedPrefix = ipm.output.continuousTrustedPrefix(trusted);
if ~any(trustedPrefix)
    error('ipm:MaximumGrowthPlotTrustedPrefix', ...
        'The result has no records in its continuous trusted prefix.');
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
caseLabel = optional_case_label(result.metadata,'caseId','solver result');
end

function [payload,sourceKind] = read_source(source)
sourceKind = 'struct';
if ischar(source) || (isstring(source) && isscalar(source))
    fileName = char(source);
    if ~isfile(fileName)
        error('ipm:MaximumGrowthPlotFileMissing', ...
            'Input MAT-file does not exist: %s',fileName);
    end
    loaded = load(fileName);
    if isfield(loaded,'assessment')
        payload = loaded.assessment;
        sourceKind = 'assessment file';
    elseif isfield(loaded,'result')
        payload = loaded.result;
        sourceKind = 'result file';
    else
        error('ipm:MaximumGrowthPlotFileContract', ...
            'MAT-file must contain a variable named result or assessment.');
    end
elseif isstruct(source) && isscalar(source)
    payload = source;
else
    error('ipm:MaximumGrowthPlotInput', ...
        'Pass a scalar result/assessment structure or its MAT-file path.');
end
end

function valid = is_assessment(value)
valid = isstruct(value) && isscalar(value) && ...
    isfield(value,'gradientGrowthRateFit') && ...
    isfield(value,'wallPeakGrowthRateFit');
end

function validate_rate_fit(fit,label)
required = {'valid','points','time','maximum','logMaximum', ...
    'logarithmicRate','bandwidths','eFoldCount','windows'};
if ~isstruct(fit) || ~isscalar(fit) || ~all(isfield(fit,required)) || ...
        ~fit.valid
    error('ipm:MaximumGrowthPlotRateFit', ...
        'The %s maximumGrowthRateFit structure is unavailable or invalid.', ...
        label);
end
pointCount = numel(fit.time);
if pointCount ~= fit.points || numel(fit.maximum) ~= pointCount || ...
        numel(fit.logMaximum) ~= pointCount || ...
        size(fit.logarithmicRate,1) ~= pointCount || ...
        size(fit.logarithmicRate,2) ~= numel(fit.bandwidths)
    error('ipm:MaximumGrowthPlotRateFitSize', ...
        'The %s maximumGrowthRateFit fields have inconsistent sizes.',label);
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

function value = model_log_value(model,time)
switch model.name
    case 'exponential'
        value = model.parameters.logM0 + model.parameters.lambda* ...
            (time-model.referenceTime);
    case 'finite_time_power'
        value = model.parameters.logAmplitude - model.parameters.p* ...
            log(model.parameters.T-time);
    case 'double_exponential'
        kappa = model.parameters.kappa;
        value = model.parameters.logM0 + model.parameters.gamma0* ...
            expm1(kappa*(time-model.referenceTime))/kappa;
    otherwise
        value = NaN(size(time));
end
end

function value = model_rate(model,time)
switch model.name
    case 'exponential'
        value = repmat(model.parameters.lambda,size(time));
    case 'finite_time_power'
        value = model.parameters.p./(model.parameters.T-time);
    case 'double_exponential'
        value = model.parameters.gamma0*exp( ...
            model.parameters.kappa*(time-model.referenceTime));
    otherwise
        value = NaN(size(time));
end
end

function model = decision_model(window)
model = [];
if isempty(window) || strcmp(window.modelDecision,'indistinguishable')
    return;
end
switch window.modelDecision
    case 'exponential'
        model = window.exponential;
    case 'finite_time_power'
        model = window.powerLaw;
    case 'double_exponential'
        model = window.doubleExponential;
end
if ~isempty(model) && ~model.valid
    model = [];
end
end

function value = model_window_title(window)
if isempty(window)
    value = 'Log growth; no fitted trailing window';
else
    value = sprintf('Log growth and %.3g e-fold fits; decision: %s', ...
        window.requestedLogGrowth,short_model_name(window.modelDecision));
end
end

function value = decision_summary(windows)
if isempty(windows)
    value = 'Logarithmic rate; no model windows';
    return;
end
parts = cell(1,numel(windows));
for index = 1:numel(windows)
    parts{index} = sprintf('%.3g:%s', ...
        windows(index).requestedLogGrowth, ...
        short_model_name(windows(index).modelDecision));
end
value = ['Logarithmic rate; windows ',strjoin(parts,', ')];
end

function value = short_model_name(name)
switch name
    case 'exponential'
        value = 'exp';
    case 'finite_time_power'
        value = 'power';
    case 'double_exponential'
        value = 'double-exp';
    case 'indistinguishable'
        value = 'inconclusive';
    otherwise
        value = char(name);
end
end

function value = optional_case_label(source,fieldName,fallback)
if isfield(source,fieldName) && ...
        (ischar(source.(fieldName)) || ...
        (isstring(source.(fieldName)) && isscalar(source.(fieldName))))
    value = char(source.(fieldName));
else
    value = fallback;
end
end
