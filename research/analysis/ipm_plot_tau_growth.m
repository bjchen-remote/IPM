function [figureHandle,analysis] = ipm_plot_tau_growth( ...
        coarseSource,fineSource,outputFile)
%IPM_PLOT_TAU_GROWTH Compare canonical-tau growth on two resolutions.
%   IPM_PLOT_TAU_GROWTH() loads the retained 256-by-256 and 384-by-384
%   first-quadrant results, restricts each history to its own continuous
%   trusted prefix, and exports a six-panel canonical-time diagnostic.
%
%   IPM_PLOT_TAU_GROWTH(COARSESOURCE,FINESOURCE,OUTPUTFILE) accepts solver
%   result structures, fourth-order assessment structures, or MAT-file paths
%   containing variables named result or assessment.  An assessment must
%   identify its source result through sourceMetadata.resultFile (or through
%   the sibling *_assessment.mat naming convention), because the complete
%   scale history is intentionally not copied into an assessment.  Its saved
%   trusted prefix must exactly match a fresh result-based reconstruction.
%
%   The panels deliberately distinguish the physical maxima M_phys, the
%   rescaled maxima M_R, and q=C_l/C_omega.  The maximum panels compare
%   M=a+b*tau and M=a*exp(lambda*tau) over the last quarter of canonical tau
%   using canonical-time quadrature weights.  The returned ANALYSIS also contains
%   full/half/quarter-window linear and
%   exponential slopes, relative RMSE values, scale-rate identities, and
%   terminal estimates T_tail=t+1/(kappa*q).  Exporting a PNG also saves the
%   returned structure to a sibling *_analysis.mat file.

[repositoryRoot,defaultCoarse,defaultFine,defaultOutput] = default_paths();
if nargin < 1 || isempty(coarseSource)
    coarseSource = defaultCoarse;
end
if nargin < 2 || isempty(fineSource)
    fineSource = defaultFine;
end
if nargin < 3 || isempty(outputFile)
    outputFile = defaultOutput;
end
validate_text(outputFile,'output file');
outputFile = char(outputFile);

coarse = read_tau_source(coarseSource,'coarse',repositoryRoot);
fine = read_tau_source(fineSource,'fine',repositoryRoot);
datasets = {coarse,fine};
colors = [0,0.4470,0.7410;0.8500,0.3250,0.0980];

visibility = 'off';
if isempty(outputFile) && usejava('desktop')
    visibility = 'on';
end
figureHandle = figure('Color','w','Visible',visibility, ...
    'Position',[60,60,1600,1120], ...
    'Name','Canonical-tau two-scale growth');
layout = tiledlayout(figureHandle,3,2,'TileSpacing','compact', ...
    'Padding','compact');

draw_physical_fit_panel(nexttile(layout,1),datasets,colors, ...
    'physicalGradient','Physical maximum gradient M_phys', ...
    'log(M_phys,grad / M_phys,grad(0))');
draw_physical_fit_panel(nexttile(layout,2),datasets,colors, ...
    'physicalWall','Physical tracked wall peak M_phys', ...
    'log(M_phys,wall / M_phys,wall(0))');
draw_scale_panel(nexttile(layout,3),datasets,colors);
draw_rate_tail_panel(nexttile(layout,4),datasets,colors);
draw_rescaled_maximum_panel(nexttile(layout,5),datasets,colors);
draw_rescaled_width_panel(nexttile(layout,6),datasets,colors);

titleLine = sprintf('Canonical-tau growth: %s versus %s', ...
    coarse.resolutionLabel,fine.resolutionLabel);
endpointLine = sprintf([ ...
    'separate continuous trusted prefixes: tau_* %.6g / %.6g; ' ...
    'physical t_* %.6g / %.6g'], ...
    coarse.tau(end),fine.tau(end),coarse.physicalTime(end), ...
    fine.physicalTime(end));
sgtitle(layout,{titleLine,endpointLine},'Interpreter','none', ...
    'FontWeight','bold','Color',[0.12,0.12,0.12]);

if ~isempty(outputFile)
    outputFile = export_figure(figureHandle,outputFile);
end

analysis = struct();
analysis.definition = struct( ...
    'q','C_l/C_omega', ...
    'physicalClock','dt/dcanonicalTau = 1/q', ...
    'physicalMaximum','M_phys = q*M_R for an isotropic derivative', ...
    'localScaleRate','kappa = canonicalCL-canonicalCOmega', ...
    'tailEstimator','T_tail = physicalTime+1/(kappa*q)');
analysis.fitFractions = [1,0.5,0.25];
analysis.coarse = analysis_summary(coarse);
analysis.fine = analysis_summary(fine);
analysis.outputFile = outputFile;
analysis.analysisFile = '';
if ~isempty(outputFile)
    analysis.analysisFile = analysis_file_name(outputFile);
    save(analysis.analysisFile,'analysis');
end
end

function [repositoryRoot,coarseFile,fineFile,outputFile] = default_paths()
analysisDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(analysisDirectory));
verificationDirectory = fullfile(repositoryRoot,'result','verification');
coarseFile = fullfile(verificationDirectory, ...
    'fourth_order_active_analytic_A16_w055_t195_dt002.mat');
fineFile = fullfile(verificationDirectory, ...
    'fourth_order_active_analytic_A16_w055_q384_t205_dt002.mat');
outputFile = fullfile(verificationDirectory, ...
    'fourth_order_active_A16_w055_q256_vs_q384_tau_growth.png');
end

function data = read_tau_source(source,role,repositoryRoot)
[result,assessment,sourceKind,sourceLabel] = read_payload( ...
    source,role,repositoryRoot);
result = ipm.output.validate(result);
[rawTrusted,~] = ipm.output.trustedMask(result);
computedPrefix = ipm.output.continuousTrustedPrefix(rawTrusted);
computedPrefix = logical(computedPrefix(:));
if ~isempty(assessment)
    validate_assessment_prefix(assessment,computedPrefix,result,role);
    trustedPrefix = logical(assessment.trustedPrefixMask(:));
else
    trustedPrefix = computedPrefix;
end
if ~any(trustedPrefix)
    error('ipm:TauGrowthTrustedPrefix', ...
        'The %s result has no records in its continuous trusted prefix.',role);
end

common = result.history.common;
mesh = result.history.mesh;
required = {'canonicalTau','physicalTime','C_l','C_omega', ...
    'canonicalCL','canonicalCOmega','physicalGradInf','gradInf'};
if ~all(isfield(common,required))
    missing = required(~isfield(common,required));
    error('ipm:TauGrowthHistoryFields', ...
        'The %s result is missing history fields: %s.', ...
        role,strjoin(missing,', '));
end

data = struct();
data.role = role;
data.sourceKind = sourceKind;
data.sourceLabel = sourceLabel;
data.resolutionLabel = resolution_label(result,sourceLabel);
data.trustedPrefixMask = trustedPrefix;
data.tau = trusted_series(common.canonicalTau,trustedPrefix);
data.physicalTime = trusted_series(common.physicalTime,trustedPrefix);
data.C_l = trusted_series(common.C_l,trustedPrefix);
data.C_omega = trusted_series(common.C_omega,trustedPrefix);
data.q = data.C_l./data.C_omega;
data.canonicalCL = trusted_series(common.canonicalCL,trustedPrefix);
data.canonicalCOmega = trusted_series( ...
    common.canonicalCOmega,trustedPrefix);
data.localKappa = data.canonicalCL-data.canonicalCOmega;
data.physicalGradient = trusted_series( ...
    common.physicalGradInf,trustedPrefix);
data.rescaledGradient = trusted_series(common.gradInf,trustedPrefix);

[data.physicalWall,data.rescaledWall,data.wallSeriesKind] = ...
    wall_series(common,mesh,trustedPrefix);
[data.rescaledWallFwhm,data.rescaledWallCoreWidth] = ...
    width_series(mesh,trustedPrefix);
if isempty(data.rescaledWallFwhm)
    data.physicalWallFwhm = [];
else
    data.physicalWallFwhm = data.rescaledWallFwhm./data.C_l;
end
if isempty(data.rescaledWallCoreWidth)
    data.physicalWallCoreWidth = [];
else
    data.physicalWallCoreWidth = ...
        data.rescaledWallCoreWidth./data.C_l;
end

validate_series(data,role);
data.tailEstimate = NaN(size(data.tau));
validTail = isfinite(data.localKappa) & data.localKappa > 0;
data.tailEstimate(validTail) = data.physicalTime(validTail) + ...
    1./(data.localKappa(validTail).*data.q(validTail));
data.numericalQRate = gradient(log(data.q),data.tau);
data.clockDerivative = gradient(data.physicalTime,data.tau);

data.fits = struct();
data.fits.physicalGradient = fit_positive_series( ...
    data.tau,data.physicalGradient);
data.fits.physicalWall = fit_positive_series( ...
    data.tau,data.physicalWall);
data.fits.q = fit_positive_series(data.tau,data.q);
data.fits.C_l = fit_positive_series(data.tau,data.C_l);
data.fits.C_omega = fit_positive_series(data.tau,data.C_omega);
data.fits.rescaledGradient = fit_positive_series( ...
    data.tau,data.rescaledGradient);
data.fits.rescaledWall = fit_positive_series( ...
    data.tau,data.rescaledWall);
data.fits.rescaledWallFwhm = optional_positive_fit( ...
    data.tau,data.rescaledWallFwhm);
data.fits.rescaledWallCoreWidth = optional_positive_fit( ...
    data.tau,data.rescaledWallCoreWidth);
data.fits.physicalWallFwhm = optional_positive_fit( ...
    data.tau,data.physicalWallFwhm);
data.fits.physicalWallCoreWidth = optional_positive_fit( ...
    data.tau,data.physicalWallCoreWidth);
end

function [result,assessment,sourceKind,sourceLabel] = read_payload( ...
        source,role,repositoryRoot)
assessment = [];
sourceLabel = [role,' source'];
sourceFile = '';
if ischar(source) || (isstring(source) && isscalar(source))
    sourceFile = char(source);
    if ~isfile(sourceFile)
        candidate = fullfile(repositoryRoot,sourceFile);
        if isfile(candidate)
            sourceFile = candidate;
        else
            error('ipm:TauGrowthFileMissing', ...
                'The %s input MAT-file does not exist: %s.',role,char(source));
        end
    end
    loaded = load(sourceFile);
    [~,sourceLabel] = fileparts(sourceFile);
    if isfield(loaded,'result')
        result = loaded.result;
        sourceKind = 'result file';
        return;
    elseif isfield(loaded,'assessment')
        assessment = loaded.assessment;
        sourceKind = 'assessment file';
    else
        error('ipm:TauGrowthFileContract', ...
            'The %s MAT-file must contain result or assessment.',role);
    end
elseif isstruct(source) && isscalar(source)
    if is_assessment(source)
        assessment = source;
        sourceKind = 'assessment struct';
    else
        result = source;
        sourceKind = 'result struct';
        return;
    end
else
    error('ipm:TauGrowthInput', ...
        ['The %s source must be a scalar result/assessment structure ' ...
        'or a MAT-file path.'],role);
end

resultFile = assessment_result_file( ...
    assessment,sourceFile,repositoryRoot);
if isempty(resultFile) || ~isfile(resultFile)
    error('ipm:TauGrowthAssessmentResultMissing', ...
        ['The %s assessment does not resolve to an existing source result. ' ...
        'Pass the result directly or preserve sourceMetadata.resultFile.'],role);
end
loadedResult = load(resultFile,'result');
if ~isfield(loadedResult,'result')
    error('ipm:TauGrowthAssessmentResultContract', ...
        'The associated %s MAT-file does not contain result.',role);
end
result = loadedResult.result;
end

function fileName = assessment_result_file( ...
        assessment,assessmentFile,repositoryRoot)
fileName = '';
if isfield(assessment,'sourceMetadata') && ...
        isstruct(assessment.sourceMetadata) && ...
        isscalar(assessment.sourceMetadata) && ...
        isfield(assessment.sourceMetadata,'resultFile')
    candidate = char(assessment.sourceMetadata.resultFile);
    if isfile(candidate)
        fileName = candidate;
        return;
    end
    candidate = fullfile(repositoryRoot,candidate);
    if isfile(candidate)
        fileName = candidate;
        return;
    end
end
if ~isempty(assessmentFile)
    suffix = '_assessment.mat';
    if endsWith(assessmentFile,suffix)
        candidate = [extractBefore(assessmentFile, ...
            strlength(assessmentFile)-strlength(suffix)+1),'.mat'];
        candidate = char(candidate);
        if isfile(candidate)
            fileName = candidate;
        end
    end
end
end

function validate_assessment_prefix(assessment,computedPrefix,result,role)
if ~is_assessment(assessment) || ...
        ~isfield(assessment,'trustedPrefixMask')
    error('ipm:TauGrowthAssessmentContract', ...
        'The %s assessment lacks a trustedPrefixMask.',role);
end
savedPrefix = logical(assessment.trustedPrefixMask(:));
if numel(savedPrefix) ~= numel(computedPrefix) || ...
        ~isequal(savedPrefix,computedPrefix)
    error('ipm:TauGrowthAssessmentPrefixMismatch', ...
        ['The %s assessment trusted prefix does not exactly match a fresh ' ...
        'continuous-prefix reconstruction from its result.'],role);
end
if ~isequal(savedPrefix,ipm.output.continuousTrustedPrefix(savedPrefix))
    error('ipm:TauGrowthAssessmentPrefixShape', ...
        'The %s assessment mask is not a continuous prefix.',role);
end
if isfield(assessment,'trustedRecords') && ...
        assessment.trustedRecords ~= nnz(savedPrefix)
    error('ipm:TauGrowthAssessmentPrefixCount', ...
        'The %s assessment trusted-record count is inconsistent.',role);
end
if isfield(assessment,'caseId') && isfield(result.metadata,'caseId') && ...
        ~strcmp(char(assessment.caseId),char(result.metadata.caseId))
    error('ipm:TauGrowthAssessmentCase', ...
        'The %s assessment and associated result case IDs differ.',role);
end
end

function valid = is_assessment(value)
valid = isstruct(value) && isscalar(value) && ...
    isfield(value,'gradientGrowthRateFit') && ...
    isfield(value,'wallPeakGrowthRateFit');
end

function values = trusted_series(values,trustedPrefix)
values = values(:);
if numel(values) ~= numel(trustedPrefix)
    error('ipm:TauGrowthHistorySize', ...
        'A history series does not match the trusted-mask length.');
end
values = values(trustedPrefix);
end

function [physicalWall,rescaledWall,kind] = wall_series( ...
        common,mesh,trustedPrefix)
if isfield(common,'physicalTrackedWallOmegaPeak') && ...
        isfield(mesh,'trackedWallPeak')
    physicalWall = trusted_series( ...
        common.physicalTrackedWallOmegaPeak,trustedPrefix);
    rescaledWall = trusted_series(mesh.trackedWallPeak,trustedPrefix);
    kind = 'tracked wall derivative peak';
elseif isfield(common,'physicalWallOmegaPeak') && ...
        isfield(common,'wallPeak')
    physicalWall = trusted_series(common.physicalWallOmegaPeak,trustedPrefix);
    rescaledWall = trusted_series(common.wallPeak,trustedPrefix);
    kind = 'global wall derivative peak';
else
    error('ipm:TauGrowthWallSeries', ...
        'The result lacks matching physical and rescaled wall-peak histories.');
end
end

function [fwhm,coreWidth] = width_series(mesh,trustedPrefix)
fwhm = [];
coreWidth = [];
if isfield(mesh,'trackedWallPeakWidth')
    fwhm = trusted_series(mesh.trackedWallPeakWidth,trustedPrefix);
end
if isfield(mesh,'trackedWallCoreWidth')
    coreWidth = trusted_series(mesh.trackedWallCoreWidth,trustedPrefix);
end
end

function validate_series(data,role)
positiveNames = {'C_l','C_omega','q','physicalGradient', ...
    'rescaledGradient','physicalWall','rescaledWall'};
for index = 1:numel(positiveNames)
    name = positiveNames{index};
    values = data.(name);
    if any(~isfinite(values)) || any(values <= 0)
        error('ipm:TauGrowthPositiveSeries', ...
            'The %s %s series must be finite and positive.',role,name);
    end
end
if any(~isfinite(data.tau)) || any(diff(data.tau) <= 0)
    error('ipm:TauGrowthCanonicalTime', ...
        'The %s canonicalTau series must be finite and increasing.',role);
end
if any(~isfinite(data.physicalTime)) || ...
        any(diff(data.physicalTime) <= 0)
    error('ipm:TauGrowthPhysicalTime', ...
        'The %s physicalTime series must be finite and increasing.',role);
end
optionalNames = {'rescaledWallFwhm','rescaledWallCoreWidth'};
for index = 1:numel(optionalNames)
    values = data.(optionalNames{index});
    if ~isempty(values) && (any(~isfinite(values)) || any(values <= 0))
        error('ipm:TauGrowthWidthSeries', ...
            'The %s %s series must be finite and positive when present.', ...
            role,optionalNames{index});
    end
end
end

function family = optional_positive_fit(tau,values)
if isempty(values)
    family = empty_fit_family();
else
    family = fit_positive_series(tau,values);
end
end

function family = fit_positive_series(tau,values)
tau = tau(:);
values = values(:);
if numel(tau) ~= numel(values) || numel(tau) < 3 || ...
        any(~isfinite(tau)) || any(diff(tau) <= 0) || ...
        any(~isfinite(values)) || any(values <= 0)
    error('ipm:TauGrowthFitInput', ...
        'Tau fits require at least three increasing times and positive values.');
end
fractions = [1,0.5,0.25];
windows = repmat(empty_fit_window(),numel(fractions),1);
span = tau(end)-tau(1);
for index = 1:numel(fractions)
    fraction = fractions(index);
    selected = tau >= tau(end)-fraction*span;
    windows(index) = fit_window(tau(selected),values(selected),fraction);
end
family = struct('valid',true,'fractions',fractions,'windows',windows);
end

function family = empty_fit_family()
family = struct('valid',false,'fractions',[], ...
    'windows',repmat(empty_fit_window(),0,1));
end

function fit = fit_window(tau,values,fraction)
referenceTau = tau(1);
centeredTau = tau-referenceTau;
weights = trapezoid_weights(tau);
weights = weights/sum(weights);
design = [ones(numel(tau),1),centeredTau];

linearCoefficients = weighted_coefficients(design,values,weights);
linearPrediction = design*linearCoefficients;
logValues = log(values);
exponentialCoefficients = weighted_coefficients( ...
    design,logValues,weights);
exponentialPrediction = exp(design*exponentialCoefficients);

fit = empty_fit_window();
fit.fraction = fraction;
fit.startTau = tau(1);
fit.endTau = tau(end);
fit.referenceTau = referenceTau;
fit.points = numel(tau);
fit.linearIntercept = linearCoefficients(1);
fit.linearSlope = linearCoefficients(2);
fit.linearRelativeRmse = relative_rmse( ...
    linearPrediction,values,weights);
fit.linearLogRmse = log_rmse(linearPrediction,values,weights);
fit.linearR2 = weighted_r2(linearPrediction,values,weights);
fit.exponentialLogIntercept = exponentialCoefficients(1);
fit.exponentialSlope = exponentialCoefficients(2);
fit.exponentialRelativeRmse = relative_rmse( ...
    exponentialPrediction,values,weights);
fit.exponentialLogRmse = sqrt(sum(weights.* ...
    (log(exponentialPrediction)-logValues).^2));
fit.exponentialR2 = weighted_r2( ...
    exponentialPrediction,values,weights);
end

function fit = empty_fit_window()
fit = struct('fraction',NaN,'startTau',NaN,'endTau',NaN, ...
    'referenceTau',NaN,'points',0,'linearIntercept',NaN, ...
    'linearSlope',NaN,'linearRelativeRmse',NaN,'linearLogRmse',NaN, ...
    'linearR2',NaN,'exponentialLogIntercept',NaN, ...
    'exponentialSlope',NaN,'exponentialRelativeRmse',NaN, ...
    'exponentialLogRmse',NaN,'exponentialR2',NaN);
end

function coefficients = weighted_coefficients(design,values,weights)
weightedDesign = design.*sqrt(weights);
weightedValues = values.*sqrt(weights);
coefficients = weightedDesign\weightedValues;
end

function value = relative_rmse(prediction,values,weights)
value = sqrt(sum(weights.*((prediction-values)./values).^2));
end

function value = log_rmse(prediction,values,weights)
if any(prediction <= 0)
    value = NaN;
else
    value = sqrt(sum(weights.*(log(prediction)-log(values)).^2));
end
end

function value = weighted_r2(prediction,values,weights)
meanValue = sum(weights.*values);
denominator = sum(weights.*(values-meanValue).^2);
value = 1-sum(weights.*(prediction-values).^2)/max(denominator,eps);
end

function weights = trapezoid_weights(tau)
if isscalar(tau)
    weights = 1;
    return;
end
spacing = diff(tau);
weights = [spacing(1)/2; ...
    (spacing(1:end-1)+spacing(2:end))/2;spacing(end)/2];
end

function draw_physical_fit_panel(axesHandle,datasets,colors, ...
        seriesName,titleText,yLabelText)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:2
    data = datasets{index};
    values = data.(seriesName);
    normalizedLog = log(values/values(1));
    handles(end+1) = endpoint_curve(axesHandle,data.tau, ...
        normalizedLog,colors(index,:),'-'); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s data', ...
        data.resolutionLabel); %#ok<AGROW>
    fit = late_window(data.fits.(seriesName));
    selected = data.tau >= fit.startTau;
    linearMaximum = fit.linearIntercept + ...
        fit.linearSlope*(data.tau(selected)-fit.referenceTau);
    if all(linearMaximum > 0)
        linearPrediction = log(linearMaximum/values(1));
        handles(end+1) = plot(axesHandle,data.tau(selected), ...
            linearPrediction,':','Color',colors(index,:), ...
            'LineWidth',1.7); %#ok<AGROW>
        labels{end+1,1} = sprintf( ...
            '%s linear M, b=%.4g, err=%.2g', ...
            data.resolutionLabel,fit.linearSlope, ...
            fit.linearRelativeRmse); %#ok<AGROW>
    end
    prediction = fit.exponentialLogIntercept + ...
        fit.exponentialSlope*(data.tau(selected)-fit.referenceTau);
    prediction = prediction-log(values(1));
    handles(end+1) = plot(axesHandle,data.tau(selected),prediction,'--', ...
        'Color',colors(index,:),'LineWidth',1.5); %#ok<AGROW>
    labels{end+1,1} = sprintf( ...
        '%s exponential M, lambda=%.4g, err=%.2g', ...
        data.resolutionLabel,fit.exponentialSlope, ...
        fit.exponentialRelativeRmse); %#ok<AGROW>
end
finish_axes(axesHandle,'Canonical tau',yLabelText,titleText);
fitLegend = legend(axesHandle,handles,labels,'Location','northwest');
style_legend(fitLegend);
end

function draw_scale_panel(axesHandle,datasets,colors)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:2
    data = datasets{index};
    handles(end+1) = endpoint_curve(axesHandle,data.tau, ...
        log(data.q/data.q(1)),colors(index,:),'-'); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s log(q/q_0)', ...
        data.resolutionLabel); %#ok<AGROW>
    handles(end+1) = plot(axesHandle,data.tau, ...
        log(data.C_l/data.C_l(1)),'--','Color',colors(index,:), ...
        'LineWidth',1.35); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s log(C_l/C_l0)', ...
        data.resolutionLabel); %#ok<AGROW>
    handles(end+1) = plot(axesHandle,data.tau, ...
        log(data.C_omega/data.C_omega(1)),':', ...
        'Color',colors(index,:),'LineWidth',1.7); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s log(C_omega/C_omega0)', ...
        data.resolutionLabel); %#ok<AGROW>
end
finish_axes(axesHandle,'Canonical tau','Normalized logarithmic scale', ...
    'q=C_l/C_omega and its two components');
scaleLegend = legend(axesHandle,handles,labels,'Location','northwest');
style_legend(scaleLegend);
end

function draw_rate_tail_panel(axesHandle,datasets,colors)
rateHandles = gobjects(0);
rateLabels = cell(0,1);
yyaxis(axesHandle,'left');
hold(axesHandle,'on');
for index = 1:2
    data = datasets{index};
    rateHandles(end+1) = endpoint_curve(axesHandle,data.tau, ...
        data.localKappa,colors(index,:),'-'); %#ok<AGROW>
    rateLabels{end+1,1} = sprintf('%s kappa_q', ...
        data.resolutionLabel); %#ok<AGROW>
end
ylabel(axesHandle,'kappa_q = canonicalCL-canonicalCOmega');
axesHandle.YAxis(1).Color = [0.18,0.18,0.18];

yyaxis(axesHandle,'right');
hold(axesHandle,'on');
for index = 1:2
    data = datasets{index};
    tailStart = data.tau(1)+0.5*(data.tau(end)-data.tau(1));
    selected = data.tau >= tailStart & isfinite(data.tailEstimate);
    rateHandles(end+1) = plot(axesHandle,data.tau(selected), ...
        data.tailEstimate(selected),'--','Color',colors(index,:), ...
        'LineWidth',1.55); %#ok<AGROW>
    rateLabels{end+1,1} = sprintf('%s T_tail (late half)', ...
        data.resolutionLabel); %#ok<AGROW>
end
ylabel(axesHandle,'T_tail = t+1/(kappa_q q)');
axesHandle.YAxis(2).Color = [0.18,0.18,0.18];
finish_axes(axesHandle,'Canonical tau','', ...
    'Canonical scale-rate difference and terminal estimator');
rateLegend = legend(axesHandle,rateHandles,rateLabels, ...
    'Location','best');
style_legend(rateLegend);
end

function draw_rescaled_maximum_panel(axesHandle,datasets,colors)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:2
    data = datasets{index};
    handles(end+1) = endpoint_curve(axesHandle,data.tau, ...
        log(data.rescaledGradient/data.rescaledGradient(1)), ...
        colors(index,:),'-'); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s rescaled gradInf', ...
        data.resolutionLabel); %#ok<AGROW>
    handles(end+1) = plot(axesHandle,data.tau, ...
        log(data.rescaledWall/data.rescaledWall(1)),'--', ...
        'Color',colors(index,:),'LineWidth',1.5); %#ok<AGROW>
    labels{end+1,1} = sprintf('%s rescaled wallPeak', ...
        data.resolutionLabel); %#ok<AGROW>
end
finish_axes(axesHandle,'Canonical tau','log(M_R/M_R(0))', ...
    'Rescaled maxima M_R (not physical maxima)');
maximumLegend = legend(axesHandle,handles,labels,'Location','best');
style_legend(maximumLegend);
end

function draw_rescaled_width_panel(axesHandle,datasets,colors)
hold(axesHandle,'on');
handles = gobjects(0);
labels = cell(0,1);
for index = 1:2
    data = datasets{index};
    if ~isempty(data.rescaledWallFwhm)
        handles(end+1) = endpoint_curve(axesHandle,data.tau, ...
            log(data.rescaledWallFwhm/data.rescaledWallFwhm(1)), ...
            colors(index,:),'-'); %#ok<AGROW>
        labels{end+1,1} = sprintf('%s rescaled FWHM', ...
            data.resolutionLabel); %#ok<AGROW>
    end
    if ~isempty(data.rescaledWallCoreWidth)
        handles(end+1) = plot(axesHandle,data.tau, ...
            log(data.rescaledWallCoreWidth/ ...
            data.rescaledWallCoreWidth(1)),'--', ...
            'Color',colors(index,:),'LineWidth',1.5); %#ok<AGROW>
        labels{end+1,1} = sprintf('%s rescaled 90%% core', ...
            data.resolutionLabel); %#ok<AGROW>
    end
end
finish_axes(axesHandle,'Canonical tau','log(width_R/width_R(0))', ...
    'Rescaled wall widths: contraction inside the outer scale');
if isempty(handles)
    text(axesHandle,0.5,0.5,'No trusted rescaled-width history', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'Interpreter','none');
else
    widthLegend = legend(axesHandle,handles,labels,'Location','southwest');
    style_legend(widthLegend);
end
end

function fit = late_window(family)
if ~family.valid || isempty(family.windows)
    error('ipm:TauGrowthLateFit', ...
        'A required late canonical-tau fit is unavailable.');
end
[~,index] = min(abs([family.windows.fraction]-0.25));
fit = family.windows(index);
end

function handle = endpoint_curve(axesHandle,tau,values,color,lineStyle)
handle = plot(axesHandle,tau,values,'Color',color, ...
    'LineStyle',lineStyle,'LineWidth',1.8,'Marker','o', ...
    'MarkerIndices',numel(tau),'MarkerSize',6, ...
    'MarkerFaceColor',color);
end

function finish_axes(axesHandle,xLabel,yLabel,titleText)
grid(axesHandle,'on');
box(axesHandle,'on');
set(axesHandle,'Color','w','XColor',[0.15,0.15,0.15], ...
    'GridColor',[0.65,0.65,0.65],'GridAlpha',0.35,'LineWidth',0.8);
xlabel(axesHandle,xLabel,'Interpreter','none');
if ~isempty(yLabel)
    ylabel(axesHandle,yLabel,'Interpreter','none');
end
title(axesHandle,titleText,'Interpreter','none', ...
    'Color',[0.12,0.12,0.12]);
end

function style_legend(legendHandle)
set(legendHandle,'Color','w','TextColor',[0.12,0.12,0.12], ...
    'EdgeColor',[0.55,0.55,0.55],'Interpreter','none');
end

function outputFile = export_figure(figureHandle,outputFile)
[outputDirectory,~,extension] = fileparts(outputFile);
if isempty(extension)
    outputFile = [outputFile,'.png'];
elseif ~strcmpi(extension,'.png')
    close(figureHandle);
    error('ipm:TauGrowthPlotFormat', ...
        'The output file must use the PNG extension.');
end
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    [created,message] = mkdir(outputDirectory);
    if ~created
        close(figureHandle);
        error('ipm:TauGrowthPlotDirectory', ...
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

function summary = analysis_summary(data)
lateQ = late_window(data.fits.q);
lateGradient = late_window(data.fits.physicalGradient);
lateWall = late_window(data.fits.physicalWall);
lateCl = late_window(data.fits.C_l);
lateComega = late_window(data.fits.C_omega);
lateRescaledGradient = late_window(data.fits.rescaledGradient);
lateRescaledWall = late_window(data.fits.rescaledWall);

qKappa = lateQ.exponentialSlope;
fitTailTime = NaN;
if qKappa > 0
    fitTailTime = data.physicalTime(end)+1/(qKappa*data.q(end));
end
lateSummary = struct( ...
    'fraction',lateQ.fraction, ...
    'qKappa',qKappa, ...
    'qLinearRelativeRmse',lateQ.linearRelativeRmse, ...
    'qExponentialRelativeRmse',lateQ.exponentialRelativeRmse, ...
    'physicalGradientLambda',lateGradient.exponentialSlope, ...
    'physicalGradientLinearRelativeRmse', ...
        lateGradient.linearRelativeRmse, ...
    'physicalGradientExponentialRelativeRmse', ...
        lateGradient.exponentialRelativeRmse, ...
    'physicalWallLambda',lateWall.exponentialSlope, ...
    'physicalWallLinearRelativeRmse',lateWall.linearRelativeRmse, ...
    'physicalWallExponentialRelativeRmse', ...
        lateWall.exponentialRelativeRmse, ...
    'C_lRate',lateCl.exponentialSlope, ...
    'C_omegaRate',lateComega.exponentialSlope, ...
    'rescaledGradientRate',lateRescaledGradient.exponentialSlope, ...
    'rescaledWallRate',lateRescaledWall.exponentialSlope, ...
    'rescaledWallFwhmRate',optional_late_slope( ...
        data.fits.rescaledWallFwhm), ...
    'rescaledWallCoreRate',optional_late_slope( ...
        data.fits.rescaledWallCoreWidth), ...
    'physicalGradientPowerRatio',safe_ratio( ...
        lateGradient.exponentialSlope,qKappa), ...
    'physicalWallPowerRatio',safe_ratio( ...
        lateWall.exponentialSlope,qKappa), ...
    'fitTailTime',fitTailTime);

tailValues = data.tailEstimate(isfinite(data.tailEstimate));
lateHalf = data.tau >= data.tau(1)+0.5*(data.tau(end)-data.tau(1)) & ...
    isfinite(data.tailEstimate);
clockDefect = data.clockDerivative-1./data.q;
wallIdentity = data.physicalWall-data.q.*data.rescaledWall;
summary = struct();
summary.sourceKind = data.sourceKind;
summary.sourceLabel = data.sourceLabel;
summary.resolutionLabel = data.resolutionLabel;
summary.wallSeriesKind = data.wallSeriesKind;
summary.trustedRecords = numel(data.tau);
summary.trustedCanonicalTau = data.tau(end);
summary.trustedPhysicalTime = data.physicalTime(end);
summary.growth = struct( ...
    'q',data.q(end)/data.q(1), ...
    'physicalGradient',data.physicalGradient(end)/ ...
        data.physicalGradient(1), ...
    'physicalWall',data.physicalWall(end)/data.physicalWall(1), ...
    'rescaledGradient',data.rescaledGradient(end)/ ...
        data.rescaledGradient(1), ...
    'rescaledWall',data.rescaledWall(end)/data.rescaledWall(1));
summary.lateQuarter = lateSummary;
summary.tailEstimate = struct( ...
    'finalLocal',last_finite(data.tailEstimate), ...
    'lateHalfMedian',median(data.tailEstimate(lateHalf),'omitnan'), ...
    'minimum',min(tailValues,[],'omitnan'), ...
    'maximum',max(tailValues,[],'omitnan'), ...
    'fitBased',fitTailTime);
summary.identityChecks = struct( ...
    'maximumWallReconstructionRelativeDefect', ...
        max(abs(wallIdentity)./max(abs(data.physicalWall),eps)), ...
    'clockDerivativeRelativeRms', ...
        sqrt(mean((clockDefect./max(abs(1./data.q),eps)).^2)), ...
    'maximumNumericalQRateDefect', ...
        max(abs(data.numericalQRate-data.localKappa)));
summary.fits = data.fits;
summary.series = struct( ...
    'canonicalTau',data.tau, ...
    'physicalTime',data.physicalTime, ...
    'physicalGradient',data.physicalGradient, ...
    'physicalWall',data.physicalWall, ...
    'rescaledGradient',data.rescaledGradient, ...
    'rescaledWall',data.rescaledWall, ...
    'q',data.q,'C_l',data.C_l,'C_omega',data.C_omega, ...
    'canonicalCL',data.canonicalCL, ...
    'canonicalCOmega',data.canonicalCOmega, ...
    'localKappa',data.localKappa, ...
    'numericalQRate',data.numericalQRate, ...
    'tailEstimate',data.tailEstimate, ...
    'rescaledWallFwhm',data.rescaledWallFwhm, ...
    'rescaledWallCoreWidth',data.rescaledWallCoreWidth, ...
    'physicalWallFwhm',data.physicalWallFwhm, ...
    'physicalWallCoreWidth',data.physicalWallCoreWidth);
end

function value = optional_late_slope(family)
if family.valid
    window = late_window(family);
    value = window.exponentialSlope;
else
    value = NaN;
end
end

function value = safe_ratio(numerator,denominator)
if isfinite(denominator) && abs(denominator) > eps
    value = numerator/denominator;
else
    value = NaN;
end
end

function value = last_finite(values)
index = find(isfinite(values),1,'last');
if isempty(index)
    value = NaN;
else
    value = values(index);
end
end

function label = resolution_label(result,fallback)
label = '';
if isfield(result.metadata,'caseMetadata') && ...
        isstruct(result.metadata.caseMetadata)
    metadata = result.metadata.caseMetadata;
    if isfield(metadata,'gridInterpretation')
        tokens = regexp(char(metadata.gridInterpretation), ...
            '(\d+)\s*x\s*(\d+)\s*cells','tokens','once');
        if ~isempty(tokens)
            label = sprintf('%s x %s cells',tokens{1},tokens{2});
        end
    end
    if isempty(label) && isfield(metadata,'caseName')
        tokens = regexp(char(metadata.caseName), ...
            'quadrant_(\d+)','tokens','once');
        if ~isempty(tokens)
            label = sprintf('%s x %s cells',tokens{1},tokens{1});
        end
    end
end
if isempty(label) && isfield(result,'config') && ...
        isfield(result.config,'grid') && ...
        all(isfield(result.config.grid,{'nx','ny'}))
    nx = result.config.grid.nx;
    ny = result.config.grid.ny;
    if nx == 2*ny-1
        label = sprintf('%d x %d cells',ny-1,ny-1);
    else
        label = sprintf('%d x %d stored nodes',nx,ny);
    end
end
if isempty(label)
    label = fallback;
end
end

function validate_text(value,label)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('ipm:TauGrowthTextInput', ...
        'The %s must be a character vector or scalar string.',label);
end
end
