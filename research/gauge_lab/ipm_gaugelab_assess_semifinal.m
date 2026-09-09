function audit = ipm_gaugelab_assess_semifinal(result,targetPhysicalTime)
%IPM_GAUGELAB_ASSESS_SEMIFINAL Apply preregistered q256 semifinal gates.
%   AUDIT = IPM_GAUGELAB_ASSESS_SEMIFINAL(RESULT,TARGETPHYSICALTIME)
%   first applies the group-stage audit.  It then measures coefficient and
%   shape stability on the physical-time windows [target-W,target], with
%   W=0.25 and W=0.125.  Every signal, including canonical tau, is linearly
%   interpolated to the exact common physical-time boundaries before tau is
%   used as the regression coordinate.

if nargin ~= 2
    error('ipm:GaugeLabSemifinalArguments', ...
        'result and targetPhysicalTime are required.');
end
baseAudit = ipm_gaugelab_assess(result,targetPhysicalTime);
validated = ipm.output.validate(result);

limits = semifinal_limits();
audit = empty_audit(baseAudit,targetPhysicalTime,limits);
failure = strings(0,1);
if ~baseAudit.hardPassed
    failure(end+1) = "base_audit_failed";
end

common = validated.history.common;
trusted = logical(baseAudit.trustedPrefixMask(:));
requiredNames = {'physicalTime','canonicalTau','canonicalCL', ...
    'canonicalCOmega','wallPeak','gradInf'};
series = struct();
historyComplete = ~isempty(trusted) && any(trusted);
for index = 1:numel(requiredNames)
    name = requiredNames{index};
    [series.(name),valid] = trusted_series(common,name,trusted);
    historyComplete = historyComplete && valid;
end
audit.historyComplete = historyComplete;
if ~historyComplete
    failure(end+1) = "semifinal_history_missing_or_nonfinite";
    audit.failureReasons = cellstr(failure);
    audit.hardFailureReason = char(strjoin(failure,','));
    audit.hardPassed = false;
    return;
end

physicalTime = series.physicalTime;
canonicalTau = series.canonicalTau;
cl = series.canonicalCL;
comega = series.canonicalCOmega;
kappa = cl-comega;
if any(series.gradInf <= 0)
    wallToGradient = NaN(size(series.gradInf));
else
    wallToGradient = series.wallPeak./series.gradInf;
end
signals = [cl,comega,kappa,wallToGradient];

clocksValid = all(diff(physicalTime) > 0) && ...
    all(diff(canonicalTau) > 0);
audit.clocksValid = clocksValid;
if ~clocksValid || any(~isfinite(signals),'all')
    failure(end+1) = "semifinal_clock_or_signal_invalid";
    audit.failureReasons = cellstr(failure);
    audit.hardFailureReason = char(strjoin(failure,','));
    audit.hardPassed = false;
    return;
end

for index = 1:numel(limits.physicalWindowWidths)
    width = limits.physicalWindowWidths(index);
    startTime = targetPhysicalTime-width;
    [windowData,valid] = exact_physical_window(physicalTime, ...
        canonicalTau,signals,startTime,targetPhysicalTime, ...
        limits.minimumWindowRecords,limits.boundaryEpsMultiplier);
    window = empty_window(width,startTime,targetPhysicalTime);
    window.covered = valid;
    if ~valid
        reason = "window_"+window_tag(width)+ ...
            "_not_covered_or_insufficient";
        failure(end+1) = reason; %#ok<AGROW>
        window.failureReasons = cellstr(reason);
        audit.windows(index) = window;
        continue;
    end

    window.records = numel(windowData.physicalTime);
    window.physicalStart = windowData.physicalTime(1);
    window.physicalEnd = windowData.physicalTime(end);
    window.canonicalStart = windowData.canonicalTau(1);
    window.canonicalEnd = windowData.canonicalTau(end);
    window.canonicalWidth = window.canonicalEnd-window.canonicalStart;
    window.leftBoundaryInterpolated = windowData.leftBoundaryInterpolated;
    window.rightBoundaryInterpolated = ...
        windowData.rightBoundaryInterpolated;
    window.canonicalCL = tail_stats( ...
        windowData.canonicalTau,windowData.signals(:,1));
    window.canonicalCOmega = tail_stats( ...
        windowData.canonicalTau,windowData.signals(:,2));
    window.canonicalKappa = tail_stats( ...
        windowData.canonicalTau,windowData.signals(:,3));
    window.wallToGradientRatio = tail_stats( ...
        windowData.canonicalTau,windowData.signals(:,4));

    window.canonicalCLSignDefinite = ...
        strictly_one_signed(windowData.signals(:,1));
    window.canonicalCOmegaSignDefinite = ...
        strictly_one_signed(windowData.signals(:,2));
    window.canonicalCLMeanMargin = mean_margin( ...
        window.canonicalCL,windowData.signals(:,1));
    window.canonicalCOmegaMeanMargin = mean_margin( ...
        window.canonicalCOmega,windowData.signals(:,2));
    window.canonicalCLMeanNondegenerate = ...
        window.canonicalCLMeanMargin >= limits.minimumMeanToMaximumRatio;
    window.canonicalCOmegaMeanNondegenerate = ...
        window.canonicalCOmegaMeanMargin >= ...
        limits.minimumMeanToMaximumRatio;

    window.maximumRateRelativeTrend = max(abs([ ...
        window.canonicalCL.relativeTrendAcrossWindow, ...
        window.canonicalCOmega.relativeTrendAcrossWindow]));
    window.maximumRateDetrendedRelativeRms = max([ ...
        window.canonicalCL.detrendedRelativeRms, ...
        window.canonicalCOmega.detrendedRelativeRms]);
    window.absoluteKappaRelativeTrend = abs( ...
        window.canonicalKappa.relativeTrendAcrossWindow);
    window.absoluteWallToGradientRelativeTrend = abs( ...
        window.wallToGradientRatio.relativeTrendAcrossWindow);

    prefix = "window_"+window_tag(width)+"_";
    if ~window.canonicalCLSignDefinite
        failure(end+1) = prefix+"canonical_cl_sign"; %#ok<AGROW>
    end
    if ~window.canonicalCOmegaSignDefinite
        failure(end+1) = prefix+"canonical_comega_sign"; %#ok<AGROW>
    end
    if ~window.canonicalCLMeanNondegenerate
        failure(end+1) = prefix+"canonical_cl_mean_degenerate"; %#ok<AGROW>
    end
    if ~window.canonicalCOmegaMeanNondegenerate
        failure(end+1) = prefix+"canonical_comega_mean_degenerate"; %#ok<AGROW>
    end
    if window.maximumRateRelativeTrend > ...
            limits.maximumRateRelativeTrend
        failure(end+1) = prefix+"rate_relative_trend"; %#ok<AGROW>
    end
    if window.maximumRateDetrendedRelativeRms > ...
            limits.maximumRateDetrendedRelativeRms
        failure(end+1) = prefix+"rate_detrended_relative_rms"; %#ok<AGROW>
    end
    if window.absoluteKappaRelativeTrend > ...
            limits.maximumKappaRelativeTrend
        failure(end+1) = prefix+"kappa_relative_trend"; %#ok<AGROW>
    end
    if window.absoluteWallToGradientRelativeTrend > ...
            limits.maximumWallToGradientRelativeTrend
        failure(end+1) = prefix+"wall_to_gradient_relative_trend"; %#ok<AGROW>
    end
    windowFailures = failure_with_prefix(failure,prefix);
    window.failureReasons = cellstr(windowFailures);
    window.hardPassed = isempty(windowFailures);
    audit.windows(index) = window;
end

audit.failureReasons = cellstr(failure);
audit.hardFailureReason = char(strjoin(failure,','));
audit.hardPassed = isempty(failure);
end

function limits = semifinal_limits()
limits = struct( ...
    'physicalWindowWidths',[0.25,0.125], ...
    'minimumWindowRecords',3, ...
    'boundaryEpsMultiplier',100, ...
    'maximumRateRelativeTrend',5e-2, ...
    'maximumRateDetrendedRelativeRms',2e-3, ...
    'maximumKappaRelativeTrend',5e-2, ...
    'maximumWallToGradientRelativeTrend',2e-2, ...
    'minimumMeanToMaximumRatio',0.1);
end

function audit = empty_audit(baseAudit,targetPhysicalTime,limits)
audit = struct();
audit.schemaVersion = 1;
audit.kind = 'q256_dynamic_gauge_semifinal_audit';
audit.caseId = baseAudit.caseId;
audit.cOmegaGauge = baseAudit.cOmegaGauge;
audit.targetPhysicalTime = targetPhysicalTime;
audit.baseAudit = baseAudit;
audit.baseHardPassed = baseAudit.hardPassed;
audit.hardLimits = limits;
audit.boundaryInterpolation = 'linear_in_physical_time';
audit.fitCoordinate = 'canonical_tau';
audit.fitWeighting = 'canonical_tau_trapezoid';
audit.historyComplete = false;
audit.clocksValid = false;
prototype = empty_window(NaN,NaN,NaN);
audit.windows = repmat(prototype,numel(limits.physicalWindowWidths),1);
audit.failureReasons = cell(0,1);
audit.hardFailureReason = '';
audit.hardPassed = false;
end

function window = empty_window(width,startTime,endTime)
window = struct( ...
    'physicalWidth',width, ...
    'requestedPhysicalStart',startTime, ...
    'requestedPhysicalEnd',endTime, ...
    'covered',false, ...
    'records',0, ...
    'physicalStart',NaN, ...
    'physicalEnd',NaN, ...
    'canonicalStart',NaN, ...
    'canonicalEnd',NaN, ...
    'canonicalWidth',NaN, ...
    'leftBoundaryInterpolated',false, ...
    'rightBoundaryInterpolated',false, ...
    'canonicalCL',empty_stats(), ...
    'canonicalCOmega',empty_stats(), ...
    'canonicalKappa',empty_stats(), ...
    'wallToGradientRatio',empty_stats(), ...
    'canonicalCLSignDefinite',false, ...
    'canonicalCOmegaSignDefinite',false, ...
    'canonicalCLMeanMargin',NaN, ...
    'canonicalCOmegaMeanMargin',NaN, ...
    'canonicalCLMeanNondegenerate',false, ...
    'canonicalCOmegaMeanNondegenerate',false, ...
    'maximumRateRelativeTrend',NaN, ...
    'maximumRateDetrendedRelativeRms',NaN, ...
    'absoluteKappaRelativeTrend',NaN, ...
    'absoluteWallToGradientRelativeTrend',NaN, ...
    'failureReasons',{cell(0,1)}, ...
    'hardPassed',false);
end

function [values,valid] = trusted_series(group,name,mask)
values = [];
valid = isstruct(group) && isscalar(group) && isfield(group,name);
if ~valid
    return;
end
source = group.(name);
valid = isnumeric(source) && isreal(source) && isvector(source) && ...
    numel(source) == numel(mask);
if ~valid
    return;
end
source = source(:);
values = source(mask);
valid = ~isempty(values) && all(isfinite(values));
end

function [window,valid] = exact_physical_window(physicalTime, ...
        canonicalTau,signals,startTime,endTime,minimumRecords, ...
        boundaryEpsMultiplier)
window = struct('physicalTime',[],'canonicalTau',[],'signals',[], ...
    'leftBoundaryInterpolated',false, ...
    'rightBoundaryInterpolated',false);
valid = isfinite(startTime) && isfinite(endTime) && startTime < endTime;
if ~valid
    return;
end
tolerance = boundaryEpsMultiplier*eps(max([1;abs(physicalTime); ...
    abs(startTime);abs(endTime)]));
valid = startTime >= physicalTime(1)-tolerance && ...
    endTime <= physicalTime(end)+tolerance;
if ~valid
    return;
end
if abs(startTime-physicalTime(1)) <= tolerance
    startTime = physicalTime(1);
end
if abs(endTime-physicalTime(end)) <= tolerance
    endTime = physicalTime(end);
end
inside = physicalTime > startTime & physicalTime < endTime;
leftTau = interp1(physicalTime,canonicalTau,startTime,'linear');
rightTau = interp1(physicalTime,canonicalTau,endTime,'linear');
leftSignals = interp1(physicalTime,signals,startTime,'linear');
rightSignals = interp1(physicalTime,signals,endTime,'linear');
window.physicalTime = [startTime;physicalTime(inside);endTime];
window.canonicalTau = [leftTau;canonicalTau(inside);rightTau];
window.signals = [leftSignals;signals(inside,:);rightSignals];
window.leftBoundaryInterpolated = ...
    ~any(abs(physicalTime-startTime) <= tolerance);
window.rightBoundaryInterpolated = ...
    ~any(abs(physicalTime-endTime) <= tolerance);
valid = numel(window.physicalTime) >= minimumRecords && ...
    all(isfinite(window.canonicalTau)) && ...
    all(isfinite(window.signals),'all') && ...
    all(diff(window.physicalTime) > 0) && ...
    all(diff(window.canonicalTau) > 0);
end

function stats = tail_stats(t,y)
t = t(:);
y = y(:);
stats = empty_stats();
stats.records = numel(t);
if numel(t) < 2 || any(~isfinite(t)) || any(~isfinite(y)) || ...
        any(diff(t) <= 0)
    return;
end
weights = trapezoid_weights(t);
weightSum = sum(weights);
meanValue = sum(weights.*y)/weightSum;
centeredTime = t-sum(weights.*t)/weightSum;
design = [ones(size(t)),centeredTime];
weightedDesign = design.*sqrt(weights);
weightedValues = y.*sqrt(weights);
coefficients = weightedDesign\weightedValues;
fitted = design*coefficients;
residual = y-fitted;
scale = max(abs(meanValue),eps);
stats.mean = meanValue;
stats.standardDeviation = ...
    sqrt(sum(weights.*(y-meanValue).^2)/weightSum);
stats.coefficientOfVariation = stats.standardDeviation/scale;
stats.slope = coefficients(2);
stats.relativeTrendAcrossWindow = ...
    coefficients(2)*(t(end)-t(1))/scale;
stats.detrendedRelativeRms = ...
    sqrt(sum(weights.*residual.^2)/weightSum)/scale;
stats.minimum = min(y);
stats.maximum = max(y);
stats.start = y(1);
stats.finish = y(end);
stats.canonicalWindow = t(end)-t(1);
stats.window = stats.canonicalWindow;
end

function stats = empty_stats()
stats = struct('records',0,'mean',NaN,'standardDeviation',NaN, ...
    'coefficientOfVariation',NaN,'slope',NaN, ...
    'relativeTrendAcrossWindow',NaN, ...
    'detrendedRelativeRms',NaN,'minimum',NaN,'maximum',NaN, ...
    'start',NaN,'finish',NaN,'canonicalWindow',NaN,'window',NaN);
end

function weights = trapezoid_weights(t)
weights = zeros(size(t));
weights(1) = (t(2)-t(1))/2;
weights(end) = (t(end)-t(end-1))/2;
if numel(t) > 2
    weights(2:end-1) = (t(3:end)-t(1:end-2))/2;
end
end

function definite = strictly_one_signed(values)
definite = all(values > 0) || all(values < 0);
end

function margin = mean_margin(stats,values)
margin = abs(stats.mean)/max(max(abs(values)),eps);
end

function matching = failure_with_prefix(failure,prefix)
matching = failure(startsWith(failure,prefix));
end

function tag = window_tag(width)
tag = replace(string(sprintf('%.3g',width)),'.','p');
end
