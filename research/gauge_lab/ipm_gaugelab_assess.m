function audit = ipm_gaugelab_assess(source,targetPhysicalTime)
%IPM_GAUGELAB_ASSESS Audit coefficient stability on one trusted q256 run.
%   A common physical-time tail selects the samples compared across gauges;
%   canonical tau remains the independent variable for scaling-rate fits.
%   Gauge-fixed observables are reported only as normalization diagnostics,
%   not as independent shape evidence.

result = ipm.output.validate(source);
[rawTrusted,trustChecks] = ipm.output.trustedMask(result);
trusted = ipm.output.continuousTrustedPrefix(rawTrusted);
common = result.history.common;
gauge = result.history.gauge;
mesh = result.history.mesh;
lastTrusted = find(trusted,1,'last');
limits = hard_limits();

if nargin < 2 || isempty(targetPhysicalTime)
    targetPhysicalTime = result.config.time.physicalFinalTime;
    targetSource = 'result.config.time.physicalFinalTime';
else
    targetSource = 'argument';
end
validateattributes(targetPhysicalTime,{'numeric'}, ...
    {'scalar','real','nonnan','nonnegative'},mfilename, ...
    'targetPhysicalTime');

gaugeName = result.config.scaling.cOmegaGauge;
isAnchorGauge = any(strcmp(gaugeName, ...
    {'anchor_wall_slope','anchor_wall_strain'}));
isWindowGauge = any(strcmp(gaugeName, ...
    {'anchor_wall_window_l2','anchor_wall_window_l4'}));
isTemplateGauge = strcmp(gaugeName,'anchor_wall_template_projection');
isBulkGauge = strcmp(gaugeName,'anchor_bulk_gradient_l2');
isQuadraticPeakGauge = strcmp(gaugeName,'wall_omega_quadratic_peak');
[specification,specificationStatus] = ...
    selected_specification(result,gaugeName);
audit = empty_audit(result,rawTrusted,trusted,trustChecks,limits, ...
    targetPhysicalTime,targetSource,isAnchorGauge,isWindowGauge, ...
    isTemplateGauge,isBulkGauge,isQuadraticPeakGauge, ...
    specification,lastTrusted);
audit.candidateSpecificationSource = specificationStatus.source;
audit.candidateSpecificationTournament = ...
    specificationStatus.tournament;
audit.candidateSpecificationValid = specificationStatus.valid;
audit.candidateSpecificationFailure = specificationStatus.failure;
audit.contract = experiment_contract(result,mesh,limits);
failure = contract_failures(audit.contract);
if ~specificationStatus.valid
    failure(end+1) = "candidate_specification_invalid";
end

if isempty(lastTrusted)
    failure(end+1) = "empty_trusted_prefix";
    audit.hardFailureReason = char(strjoin(failure,','));
    return;
end

[tau,tauFinite] = required_masked_series( ...
    common,'canonicalTau',trusted);
[physicalTime,physicalTimeFinite] = required_masked_series( ...
    common,'physicalTime',trusted);
tailClocksValid = tauFinite && physicalTimeFinite && ...
    all(diff(tau) > 0) && all(diff(physicalTime) > 0);
if tailClocksValid
    availablePhysicalSpan = physicalTime(end)-physicalTime(1);
    requestedTailPhysicalWidth = min( ...
        limits.tailPhysicalWidthCap, ...
        limits.tailPhysicalSpanFraction*availablePhysicalSpan);
    if requestedTailPhysicalWidth <= 0
        tail = true(size(physicalTime));
    else
        tail = physicalTime >= ...
            physicalTime(end)-requestedTailPhysicalWidth;
    end
    firstTail = find(tail,1,'first');
    audit.availablePhysicalSpan = availablePhysicalSpan;
    audit.requestedTailPhysicalWidth = requestedTailPhysicalWidth;
    audit.tailPhysicalWidth = physicalTime(end)-physicalTime(firstTail);
    audit.tailCanonicalWidth = tau(end)-tau(firstTail);
else
    tail = false(size(tau));
end
audit.tailRecords = nnz(tail);

[canonicalCL,canonicalCLFinite] = required_masked_series( ...
    common,'canonicalCL',trusted);
[canonicalCOmega,canonicalCOmegaFinite] = required_masked_series( ...
    common,'canonicalCOmega',trusted);
canonicalRatesFinite = canonicalCLFinite && canonicalCOmegaFinite;
if canonicalRatesFinite
    canonicalKappa = canonicalCL-canonicalCOmega;
else
    canonicalKappa = [];
end
audit.canonicalCL = selected_tail_stats( ...
    tau,canonicalCL,tail,tailClocksValid && canonicalCLFinite);
audit.canonicalCOmega = selected_tail_stats( ...
    tau,canonicalCOmega,tail,tailClocksValid && canonicalCOmegaFinite);
audit.canonicalKappa = selected_tail_stats( ...
    tau,canonicalKappa,tail,tailClocksValid && canonicalRatesFinite);

[rescaledGradient,rescaledGradientFinite] = required_masked_series( ...
    common,'gradInf',trusted);
[wallPeak,wallPeakFinite] = required_masked_series( ...
    common,'wallPeak',trusted);
shapeSeriesFinite = rescaledGradientFinite && wallPeakFinite;
if shapeSeriesFinite
    wallToGradient = wallPeak./max(rescaledGradient,eps);
else
    wallToGradient = [];
end
audit.rescaledGradient = selected_tail_stats( ...
    tau,rescaledGradient,tail,tailClocksValid && rescaledGradientFinite);
audit.wallToGradientRatio = selected_tail_stats( ...
    tau,wallToGradient,tail,tailClocksValid && shapeSeriesFinite);

[timeSpeed,timeSpeedFinite] = required_masked_series( ...
    common,'timeSpeed',trusted);
timeSpeedPassed = timeSpeedFinite && all(timeSpeed == 1);
if timeSpeedFinite
    audit.minimumTimeSpeed = min(timeSpeed);
    audit.maximumTimeSpeedDeviation = max(abs(timeSpeed-1));
end

lengthNames = {'c_lNominal','canonicalCLNominal','lengthScaleGain', ...
    'widthControlMode','lengthGaugeResidual','transportAnchorX', ...
    'canonicalTransportAnchorVelocity','transportAnchorVelocity'};
lengthSeries = struct();
lengthHistoryComplete = canonicalCLFinite;
for index = 1:numel(lengthNames)
    name = lengthNames{index};
    [lengthSeries.(name),valid] = required_masked_series( ...
        gauge,name,trusted);
    lengthHistoryComplete = lengthHistoryComplete && valid;
end
audit.lengthGaugeHistoryComplete = lengthHistoryComplete;
if lengthHistoryComplete
    cLNominal = lengthSeries.c_lNominal;
    canonicalCLNominal = lengthSeries.canonicalCLNominal;
    lengthScaleGain = lengthSeries.lengthScaleGain;
    widthControlMode = lengthSeries.widthControlMode;
    lengthGaugeResidual = lengthSeries.lengthGaugeResidual;
    anchorX = lengthSeries.transportAnchorX;
    canonicalAnchorVelocity = ...
        lengthSeries.canonicalTransportAnchorVelocity;
    anchorVelocity = lengthSeries.transportAnchorVelocity;
    audit.minimumLengthScaleGain = min(lengthScaleGain);
    audit.maximumLengthScaleGain = max(lengthScaleGain);
    audit.maximumAbsoluteWidthControlMode = max(abs(widthControlMode));
    audit.lengthScaleGainExactlyOne = all(lengthScaleGain == 1);
    audit.widthControlModeExactlyZero = all(widthControlMode == 0);
    audit.maximumCanonicalCLNominalIdentityError = max( ...
        relative_identity_error(canonicalCL,canonicalCLNominal));
    audit.maximumCLNominalIdentityError = max( ...
        relative_identity_error(canonicalCL,cLNominal));
    audit.canonicalCLMatchesCanonicalCLNominal = ...
        all(canonicalCL == canonicalCLNominal);
    audit.canonicalCLMatchesCLNominal = ...
        all(canonicalCL == cLNominal);
    audit.maximumTransportAnchorXError = max(abs( ...
        anchorX-limits.requiredTransportAnchorX));
    audit.transportAnchorXFixedAtOne = all( ...
        anchorX == limits.requiredTransportAnchorX);
    [normalizedLengthResidual,lengthResidualScale] = ...
        normalized_anchor_residual( ...
        lengthGaugeResidual,cLNominal,anchorX);
    [normalizedCanonicalAnchorVelocity,canonicalAnchorScale] = ...
        normalized_anchor_residual( ...
        canonicalAnchorVelocity,canonicalCLNominal,anchorX);
    [normalizedAnchorVelocity,anchorVelocityScale] = ...
        normalized_anchor_residual( ...
        anchorVelocity,cLNominal,anchorX);
    audit.transportAnchorNormalizationFinite = all(isfinite([ ...
        normalizedLengthResidual;normalizedCanonicalAnchorVelocity; ...
        normalizedAnchorVelocity;lengthResidualScale; ...
        canonicalAnchorScale;anchorVelocityScale]));
    audit.maximumNormalizedLengthGaugeResidual = ...
        max(normalizedLengthResidual);
    audit.maximumRawLengthGaugeResidual = max(abs(lengthGaugeResidual));
    audit.maximumNormalizedCanonicalTransportAnchorVelocity = ...
        max(normalizedCanonicalAnchorVelocity);
    audit.maximumRawCanonicalTransportAnchorVelocity = ...
        max(abs(canonicalAnchorVelocity));
    audit.maximumNormalizedTransportAnchorVelocity = ...
        max(normalizedAnchorVelocity);
    audit.maximumRawTransportAnchorVelocity = max(abs(anchorVelocity));
    audit.maximumNormalizedTransportAnchorDefect = max([ ...
        audit.maximumNormalizedLengthGaugeResidual, ...
        audit.maximumNormalizedCanonicalTransportAnchorVelocity, ...
        audit.maximumNormalizedTransportAnchorVelocity]);
    audit.minimumLengthGaugeResidualScale = min(lengthResidualScale);
    audit.minimumCanonicalTransportAnchorScale = ...
        min(canonicalAnchorScale);
    audit.minimumTransportAnchorScale = min(anchorVelocityScale);
end

[canonicalWidthCorrection,canonicalWidthCorrectionFinite] = ...
    required_masked_series(gauge,'canonicalWidthCorrection',trusted);
[widthCorrection,widthCorrectionFinite] = required_masked_series( ...
    gauge,'widthCorrection',trusted);
widthCorrectionHistoryComplete = canonicalWidthCorrectionFinite && ...
    widthCorrectionFinite;
audit.widthCorrectionHistoryComplete = widthCorrectionHistoryComplete;
if widthCorrectionHistoryComplete
    audit.maximumCanonicalWidthCorrection = ...
        max(abs(canonicalWidthCorrection));
    audit.maximumWidthCorrection = max(abs(widthCorrection));
    audit.maximumWidthCorrectionMagnitude = max( ...
        audit.maximumCanonicalWidthCorrection, ...
        audit.maximumWidthCorrection);
end

[scaleLength,scaleLengthFinite] = required_masked_series( ...
    common,'C_l',trusted);
[scaleAmplitude,scaleAmplitudeFinite] = required_masked_series( ...
    common,'C_omega',trusted);
if timeSpeedPassed && scaleLengthFinite && scaleAmplitudeFinite && ...
        scaleLength(end) > 0 && scaleAmplitude(end) > 0
    audit.trustedEndpointPhysicalClockSpeed = ...
        scaleAmplitude(end)/scaleLength(end);
end
if isfinite(targetPhysicalTime) && ...
        isfinite(audit.trustedEndpointPhysicalClockSpeed) && ...
        audit.trustedEndpointPhysicalClockSpeed > 0
    roundoffTolerance = limits.targetRoundoffEpsMultiplier * ...
        eps(max(1,abs(targetPhysicalTime)));
    minimumStepTolerance = limits.targetMinimumStepMultiplier * ...
        result.config.time.minDt * ...
        audit.trustedEndpointPhysicalClockSpeed;
    audit.targetPhysicalTimeTolerance = ...
        roundoffTolerance+minimumStepTolerance;
    audit.reachedTarget = isfinite(audit.physicalTime) && ...
        audit.physicalTime >= ...
        targetPhysicalTime-audit.targetPhysicalTimeTolerance;
end

[gaugeResidual,gaugeResidualFinite] = required_masked_series( ...
    gauge,'omegaGaugeResidual',trusted);
if gaugeResidualFinite
    audit.maximumGaugeResidual = max(abs(gaugeResidual));
    if isAnchorGauge
        audit.maximumAnchorNormalizedRateResidual = ...
            audit.maximumGaugeResidual;
    end
end

if isAnchorGauge
    anchorNames = {'omegaGaugeAnchorX','omegaGaugeAnchorValue', ...
        'omegaGaugeAnchorCondition','omegaGaugeAnchorRateResidual', ...
        'omegaGaugeReferenceAnchorValue'};
    anchorSeries = struct();
    anchorHistoryComplete = true;
    for index = 1:numel(anchorNames)
        name = anchorNames{index};
        [anchorSeries.(name),valid] = required_masked_series( ...
            gauge,name,trusted);
        anchorHistoryComplete = anchorHistoryComplete && valid;
    end
    audit.anchorHistoryComplete = anchorHistoryComplete;
    if anchorHistoryComplete
        anchorX = anchorSeries.omegaGaugeAnchorX;
        anchorValue = anchorSeries.omegaGaugeAnchorValue;
        anchorCondition = anchorSeries.omegaGaugeAnchorCondition;
        boundedRateResidual = ...
            anchorSeries.omegaGaugeAnchorRateResidual;
        referenceValue = anchorSeries.omegaGaugeReferenceAnchorValue;
        audit.minimumAnchorWallCondition = min(anchorCondition);
        audit.maximumBoundedAnchorRateResidual = ...
            max(abs(boundedRateResidual));
        audit.anchorReferenceValue = referenceValue(1);
        audit.anchorCoordinateTolerance = ...
            limits.contractNumericEpsMultiplier * eps(max( ...
            1,abs(audit.transportAnchorX)));
        audit.maximumAnchorCoordinateError = ...
            max(abs(anchorX-audit.transportAnchorX));
        audit.anchorXMatchesConfig = ...
            audit.maximumAnchorCoordinateError <= ...
            audit.anchorCoordinateTolerance;
        audit.anchorReferenceTolerance = ...
            limits.contractNumericEpsMultiplier * ...
            eps(max(1,abs(referenceValue(1))));
        audit.maximumAnchorReferenceDrift = ...
            max(abs(referenceValue-referenceValue(1)));
        audit.anchorReferenceConstant = ...
            audit.maximumAnchorReferenceDrift <= ...
            audit.anchorReferenceTolerance;
        audit.anchorSignPreserved = ...
            all(anchorValue.*referenceValue > 0);
        if audit.anchorSignPreserved
            relativeDrift = anchorValue./referenceValue-1;
            audit.maximumAnchorRelativeDrift = ...
                max(abs(relativeDrift));
            audit.finalAnchorRelativeDrift = relativeDrift(end);
        end
    end
end

if isWindowGauge
    windowNames = {'omegaGaugeWindowCenter','omegaGaugeWindowRadius', ...
        'omegaGaugeWindowValue','omegaGaugeWindowCondition', ...
        'omegaGaugeWindowEnergy','omegaGaugeWindowRateResidual', ...
        'omegaGaugeReferenceWindowValue','omegaGaugeWindowMomentOrder', ...
        'omegaGaugeWindowMoment','omegaGaugeWindowSupportPoints', ...
        'omegaGaugeWindowRateScale'};
    windowSeries = struct();
    windowHistoryComplete = true;
    for index = 1:numel(windowNames)
        name = windowNames{index};
        [windowSeries.(name),valid] = required_masked_series( ...
            gauge,name,trusted);
        windowHistoryComplete = windowHistoryComplete && valid;
    end
    audit.windowHistoryComplete = windowHistoryComplete;
    if windowHistoryComplete
        windowCenter = windowSeries.omegaGaugeWindowCenter;
        windowRadius = windowSeries.omegaGaugeWindowRadius;
        windowValue = windowSeries.omegaGaugeWindowValue;
        windowCondition = windowSeries.omegaGaugeWindowCondition;
        windowEnergy = windowSeries.omegaGaugeWindowEnergy;
        windowRateResidual = ...
            windowSeries.omegaGaugeWindowRateResidual;
        windowRateScale = windowSeries.omegaGaugeWindowRateScale;
        referenceValue = ...
            windowSeries.omegaGaugeReferenceWindowValue;
        momentOrder = windowSeries.omegaGaugeWindowMomentOrder;
        moment = windowSeries.omegaGaugeWindowMoment;
        supportPoints = windowSeries.omegaGaugeWindowSupportPoints;
        audit.minimumWindowCondition = min(windowCondition);
        audit.maximumWindowRawRateResidual = max(abs(windowRateResidual));
        audit.maximumWindowNormalizedRateResidual = max( ...
            abs(windowRateResidual)./max(windowRateScale,realmin));
        audit.maximumWindowRateResidual = ...
            audit.maximumWindowNormalizedRateResidual;
        audit.minimumWindowRateScale = min(windowRateScale);
        audit.windowReferenceValue = referenceValue(1);
        audit.windowCenterTolerance = ...
            limits.contractNumericEpsMultiplier * ...
            eps(max(1,abs(specification.fixedSupport.center(1))));
        audit.maximumWindowCenterError = ...
            max(abs(windowCenter-specification.fixedSupport.center(1)));
        audit.windowCenterMatchesConfig = ...
            audit.maximumWindowCenterError <= ...
            audit.windowCenterTolerance;
        audit.windowCenterMatchesSpecification = ...
            audit.windowCenterMatchesConfig;
        audit.windowRadiusTolerance = ...
            limits.contractNumericEpsMultiplier * ...
            eps(max(1,abs(specification.fixedSupport.radius)));
        audit.maximumWindowRadiusError = ...
            max(abs(windowRadius-specification.fixedSupport.radius));
        audit.windowRadiusMatchesConfig = ...
            audit.maximumWindowRadiusError <= ...
            audit.windowRadiusTolerance;
        audit.windowRadiusMatchesSpecification = ...
            audit.windowRadiusMatchesConfig;
        audit.windowReferenceTolerance = ...
            limits.contractNumericEpsMultiplier * ...
            eps(max(1,abs(referenceValue(1))));
        audit.maximumWindowReferenceDrift = ...
            max(abs(referenceValue-referenceValue(1)));
        audit.windowReferenceConstant = ...
            audit.maximumWindowReferenceDrift <= ...
            audit.windowReferenceTolerance;
        audit.windowReferencePositive = all(referenceValue > 0);
        audit.windowValuePositive = all(windowValue > 0);
        audit.windowSignPreserved = all(windowValue.*referenceValue > 0);
        audit.windowEnergyPositive = all(windowEnergy > 0);
        audit.windowMomentPositive = all(moment > 0);
        audit.windowRateScalePositive = all(windowRateScale > 0);
        audit.windowSupportPositive = all(supportPoints >= 3);
        audit.windowSupportFixed = integer_constant_series(supportPoints);
        audit.windowSupportPoints = supportPoints(1);
        [audit.windowSupportGeometryMatchesSpecification, ...
            audit.maximumWindowSupportGeometryError] = ...
            fixed_support_geometry(windowCenter,zeros(size(windowCenter)), ...
            windowRadius,specification.fixedSupport,limits,'wall');
        audit.windowMomentOrderMatchesSpecification = ...
            all(momentOrder == specification.momentOrder);
        audit.maximumWindowMomentIdentityError = max( ...
            relative_identity_error(windowValue.^momentOrder,moment));
        if specification.momentOrder == 2
            audit.maximumWindowEnergyIdentityError = max( ...
                relative_identity_error(windowEnergy,moment));
        end
        if audit.windowReferencePositive && audit.windowValuePositive
            relativeDrift = windowValue./referenceValue-1;
            audit.maximumWindowRelativeDrift = max(abs(relativeDrift));
            audit.finalWindowRelativeDrift = relativeDrift(end);
        end
    end
end

if isTemplateGauge
    templateNames = {'omegaGaugeTemplateCenter', ...
        'omegaGaugeTemplateRadius','omegaGaugeTemplateSupportPoints', ...
        'omegaGaugeTemplateValue','omegaGaugeTemplateCondition', ...
        'omegaGaugeTemplateRateResidual','omegaGaugeTemplateRateScale', ...
        'omegaGaugeReferenceTemplateValue'};
    templateSeries = struct();
    templateHistoryComplete = true;
    for index = 1:numel(templateNames)
        name = templateNames{index};
        [templateSeries.(name),valid] = required_masked_series( ...
            gauge,name,trusted);
        templateHistoryComplete = templateHistoryComplete && valid;
    end
    audit.templateHistoryComplete = templateHistoryComplete;
    if templateHistoryComplete
        center = templateSeries.omegaGaugeTemplateCenter;
        radius = templateSeries.omegaGaugeTemplateRadius;
        supportPoints = templateSeries.omegaGaugeTemplateSupportPoints;
        value = templateSeries.omegaGaugeTemplateValue;
        condition = templateSeries.omegaGaugeTemplateCondition;
        rateResidual = templateSeries.omegaGaugeTemplateRateResidual;
        rateScale = templateSeries.omegaGaugeTemplateRateScale;
        reference = templateSeries.omegaGaugeReferenceTemplateValue;
        audit.minimumTemplateCondition = min(condition);
        audit.maximumTemplateRawRateResidual = max(abs(rateResidual));
        audit.maximumTemplateNormalizedRateResidual = max( ...
            abs(rateResidual)./max(rateScale,realmin));
        audit.minimumTemplateRateScale = min(rateScale);
        [audit.templateCenterMatchesSpecification, ...
            audit.maximumTemplateCenterError, ...
            audit.templateCenterTolerance] = fixed_numeric_series( ...
            center,specification.fixedSupport.center(1),limits);
        [audit.templateRadiusMatchesSpecification, ...
            audit.maximumTemplateRadiusError, ...
            audit.templateRadiusTolerance] = fixed_numeric_series( ...
            radius,specification.fixedSupport.radius,limits);
        audit.templateSupportPositive = all(supportPoints >= 3);
        audit.templateSupportFixed = integer_constant_series(supportPoints);
        audit.templateSupportPoints = supportPoints(1);
        [audit.templateSupportGeometryMatchesSpecification, ...
            audit.maximumTemplateSupportGeometryError] = ...
            fixed_support_geometry(center,zeros(size(center)),radius, ...
            specification.fixedSupport,limits,'wall');
        audit.templateReferencePositive = all(reference > 0);
        audit.templateReferenceConstant = constant_series(reference,limits);
        audit.templateValuePositive = all(value > 0);
        audit.templateSignPreserved = all(value.*reference > 0);
        if audit.templateReferencePositive && audit.templateSignPreserved
            relativeDrift = value./reference-1;
            audit.maximumTemplateRelativeDrift = max(abs(relativeDrift));
            audit.finalTemplateRelativeDrift = relativeDrift(end);
        end
    end
end

if isBulkGauge
    bulkNames = {'omegaGaugeBulkCenterX','omegaGaugeBulkCenterY', ...
        'omegaGaugeBulkRadius','omegaGaugeBulkSupportPointsX', ...
        'omegaGaugeBulkSupportPointsY','omegaGaugeBulkValue', ...
        'omegaGaugeBulkCondition','omegaGaugeBulkEnergy', ...
        'omegaGaugeBulkRateResidual','omegaGaugeBulkRateScale', ...
        'omegaGaugeReferenceBulkValue'};
    bulkSeries = struct();
    bulkHistoryComplete = true;
    for index = 1:numel(bulkNames)
        name = bulkNames{index};
        [bulkSeries.(name),valid] = required_masked_series( ...
            gauge,name,trusted);
        bulkHistoryComplete = bulkHistoryComplete && valid;
    end
    audit.bulkHistoryComplete = bulkHistoryComplete;
    if bulkHistoryComplete
        centerX = bulkSeries.omegaGaugeBulkCenterX;
        centerY = bulkSeries.omegaGaugeBulkCenterY;
        radius = bulkSeries.omegaGaugeBulkRadius;
        supportX = bulkSeries.omegaGaugeBulkSupportPointsX;
        supportY = bulkSeries.omegaGaugeBulkSupportPointsY;
        value = bulkSeries.omegaGaugeBulkValue;
        condition = bulkSeries.omegaGaugeBulkCondition;
        energy = bulkSeries.omegaGaugeBulkEnergy;
        rateResidual = bulkSeries.omegaGaugeBulkRateResidual;
        rateScale = bulkSeries.omegaGaugeBulkRateScale;
        reference = bulkSeries.omegaGaugeReferenceBulkValue;
        audit.minimumBulkCondition = min(condition);
        audit.maximumBulkRawRateResidual = max(abs(rateResidual));
        audit.maximumBulkNormalizedRateResidual = max( ...
            abs(rateResidual)./max(rateScale,realmin));
        audit.minimumBulkRateScale = min(rateScale);
        [audit.bulkCenterXMatchesSpecification, ...
            audit.maximumBulkCenterXError,audit.bulkCenterXTolerance] = ...
            fixed_numeric_series(centerX, ...
            specification.fixedSupport.center(1),limits);
        [audit.bulkCenterYMatchesSpecification, ...
            audit.maximumBulkCenterYError,audit.bulkCenterYTolerance] = ...
            fixed_numeric_series(centerY, ...
            specification.fixedSupport.center(2),limits);
        [audit.bulkRadiusMatchesSpecification, ...
            audit.maximumBulkRadiusError,audit.bulkRadiusTolerance] = ...
            fixed_numeric_series(radius,specification.fixedSupport.radius, ...
            limits);
        audit.bulkSupportPositive = ...
            all(supportX >= 3) && all(supportY >= 3);
        audit.bulkSupportFixed = integer_constant_series(supportX) && ...
            integer_constant_series(supportY);
        audit.bulkSupportPoints = [supportX(1),supportY(1)];
        [audit.bulkSupportGeometryMatchesSpecification, ...
            audit.maximumBulkSupportGeometryError] = ...
            fixed_support_geometry(centerX,centerY,radius, ...
            specification.fixedSupport,limits,'bulk');
        audit.bulkReferencePositive = all(reference > 0);
        audit.bulkReferenceConstant = constant_series(reference,limits);
        audit.bulkValuePositive = all(value > 0);
        audit.bulkEnergyPositive = all(energy > 0);
        audit.bulkSignPreserved = all(value.*reference > 0);
        audit.maximumBulkEnergyIdentityError = max( ...
            relative_identity_error(value.^2,energy));
        if audit.bulkReferencePositive && audit.bulkSignPreserved
            relativeDrift = value./reference-1;
            audit.maximumBulkRelativeDrift = max(abs(relativeDrift));
            audit.finalBulkRelativeDrift = relativeDrift(end);
        end
    end
end

if isQuadraticPeakGauge
    peakNames = {'omegaGaugePeakX','omegaGaugePeakCondition', ...
        'omegaGaugeQuadraticPeakX','omegaGaugeQuadraticPeakValue', ...
        'omegaGaugeQuadraticForcing', ...
        'omegaGaugeQuadraticCurvature', ...
        'omegaGaugeQuadraticNormalizedStationarityResidual', ...
        'omegaGaugeQuadraticCurvatureCondition', ...
        'omegaGaugeQuadraticWeightCondition', ...
        'omegaGaugeQuadraticActiveNodeMargin', ...
        'omegaGaugeQuadraticActiveSelectionMargin', ...
        'omegaGaugeQuadraticStencilCenterX', ...
        'omegaGaugeQuadraticVertexInteriorFraction', ...
        'omegaGaugeQuadraticRelativeSubgridLift', ...
        'omegaGaugeQuadraticRateResidual', ...
        'omegaGaugeQuadraticRateScale'};
    peakSeries = struct();
    peakHistoryComplete = true;
    for index = 1:numel(peakNames)
        name = peakNames{index};
        [peakSeries.(name),valid] = required_masked_series( ...
            gauge,name,trusted);
        peakHistoryComplete = peakHistoryComplete && valid;
    end
    audit.quadraticPeakHistoryComplete = peakHistoryComplete;
    if peakHistoryComplete
        audit.minimumQuadraticPeakValue = min( ...
            peakSeries.omegaGaugeQuadraticPeakValue);
        audit.maximumQuadraticPeakCurvature = max( ...
            peakSeries.omegaGaugeQuadraticCurvature);
        audit.minimumQuadraticPeakCurvatureCondition = min( ...
            peakSeries.omegaGaugeQuadraticCurvatureCondition);
        audit.minimumQuadraticPeakWeightCondition = min( ...
            peakSeries.omegaGaugeQuadraticWeightCondition);
        audit.minimumQuadraticPeakVertexInteriorFraction = min( ...
            peakSeries.omegaGaugeQuadraticVertexInteriorFraction);
        audit.minimumQuadraticPeakActiveNodeMargin = min( ...
            peakSeries.omegaGaugeQuadraticActiveNodeMargin);
        audit.minimumQuadraticPeakActiveSelectionMargin = min( ...
            peakSeries.omegaGaugeQuadraticActiveSelectionMargin);
        audit.quadraticPeakStencilSwitchCount = nnz(diff( ...
            peakSeries.omegaGaugeQuadraticStencilCenterX) ~= 0);
        audit.minimumQuadraticPeakGlobalCondition = min( ...
            peakSeries.omegaGaugePeakCondition);
        audit.minimumQuadraticPeakRateScale = min( ...
            peakSeries.omegaGaugeQuadraticRateScale);
        audit.maximumQuadraticPeakNormalizedRateResidual = max( ...
            abs(peakSeries.omegaGaugeQuadraticRateResidual) ./ ...
            max(peakSeries.omegaGaugeQuadraticRateScale,realmin));
        reconstructedResidual = ...
            peakSeries.omegaGaugeQuadraticForcing + canonicalCOmega.* ...
            peakSeries.omegaGaugeQuadraticPeakValue;
        reconstructedScale = max([ ...
            abs(peakSeries.omegaGaugeQuadraticForcing), ...
            abs(canonicalCOmega.* ...
                peakSeries.omegaGaugeQuadraticPeakValue), ...
            realmin*ones(size(canonicalCOmega))],[],2);
        audit.maximumQuadraticPeakReconstructedNormalizedResidual = ...
            max(abs(reconstructedResidual)./reconstructedScale);
        audit.maximumQuadraticPeakReportedResidualMismatch = max( ...
            abs(reconstructedResidual - ...
                peakSeries.omegaGaugeQuadraticRateResidual) ./ ...
            reconstructedScale);
        audit.maximumQuadraticPeakCoordinateIdentityError = max(abs( ...
            peakSeries.omegaGaugePeakX - ...
            peakSeries.omegaGaugeQuadraticPeakX));
        audit.minimumQuadraticPeakRelativeSubgridLift = min( ...
            peakSeries.omegaGaugeQuadraticRelativeSubgridLift);
        audit.maximumQuadraticPeakStationarityResidual = max(abs( ...
            peakSeries.omegaGaugeQuadraticNormalizedStationarityResidual));
    end
end

[horizontalCore,horizontalCoreFinite] = required_masked_series( ...
    mesh,'coreGridPoints',trusted);
[verticalCore,verticalCoreFinite] = required_masked_series( ...
    mesh,'verticalCoreGridPoints',trusted);
if horizontalCoreFinite && verticalCoreFinite
    audit.finalCoreCells = [horizontalCore(end),verticalCore(end)];
end
[safetyFactor,safetyFactorFinite] = required_masked_series( ...
    mesh,'safetyFactor',trusted);
if safetyFactorFinite
    audit.finalSafetyFactor = safetyFactor(end);
    audit.maximumSafetyFactor = max(safetyFactor);
end
audit.finalPhysicalRhoX = optional_indexed( ...
    common,'physicalRhoXInf',lastTrusted);
audit.finalPhysicalGradient = optional_indexed( ...
    common,'physicalGradInf',lastTrusted);
audit.finalFarVelocityRatio = optional_indexed( ...
    common,'farBoundaryVelocityRatio',lastTrusted);
audit.finalFarSourceRatio = optional_indexed( ...
    common,'farBoundarySourceRatio',lastTrusted);

if ~audit.terminalTrusted
    failure(end+1) = "terminal_not_in_continuous_trusted_prefix";
end
if ~audit.reachedTarget
    failure(end+1) = "trusted_physical_target_not_reached";
end
if ~canonicalRatesFinite
    failure(end+1) = "missing_or_nonfinite_canonical_rate";
end
if ~timeSpeedPassed
    failure(end+1) = "noncanonical_time_reparameterization";
end
if ~audit.lengthGaugeHistoryComplete
    failure(end+1) = ...
        "length_gauge_history_missing_or_nonfinite";
else
    if ~audit.lengthScaleGainExactlyOne
        failure(end+1) = "length_scale_gain_not_exactly_one";
    end
    if ~audit.widthControlModeExactlyZero
        failure(end+1) = "width_control_mode_not_exactly_zero";
    end
    if ~audit.canonicalCLMatchesCanonicalCLNominal || ...
            ~audit.canonicalCLMatchesCLNominal
        failure(end+1) = "canonical_cl_nominal_identity";
    end
    if ~audit.transportAnchorXFixedAtOne
        failure(end+1) = "transport_anchor_history_not_fixed_at_one";
    end
    if ~audit.transportAnchorNormalizationFinite || ...
            audit.minimumLengthGaugeResidualScale <= 0 || ...
            audit.minimumCanonicalTransportAnchorScale <= 0 || ...
            audit.minimumTransportAnchorScale <= 0 || ...
            audit.maximumNormalizedTransportAnchorDefect > ...
            limits.maximumExactGaugeResidual
        failure(end+1) = "transport_anchor_residual_or_velocity";
    end
end
if ~widthCorrectionHistoryComplete || ...
        audit.maximumWidthCorrectionMagnitude > ...
        limits.maximumWidthCorrectionMagnitude
    failure(end+1) = "nonzero_or_nonfinite_width_feedback";
end
if ~tailClocksValid
    failure(end+1) = "nonfinite_or_nonmonotone_tail_clock";
end
tailStatistics = [audit.canonicalCL,audit.canonicalCOmega, ...
    audit.canonicalKappa,audit.rescaledGradient, ...
    audit.wallToGradientRatio];
if audit.tailRecords < limits.minimumTailRecords || ...
        audit.tailPhysicalWidth <= 0 || ...
        audit.tailCanonicalWidth <= 0 || ...
        ~all(arrayfun(@stats_finite,tailStatistics))
    failure(end+1) = "insufficient_or_nonfinite_tail_statistics";
end
if any(~isfinite(audit.finalCoreCells)) || ...
        any(audit.finalCoreCells < limits.minimumTerminalCoreCells)
    failure(end+1) = "terminal_core_below_14";
end
if ~isfinite(audit.finalSafetyFactor) || ...
        audit.finalSafetyFactor >= ...
        limits.terminalSafetyFactorExclusiveUpperBound
    failure(end+1) = "terminal_safety_not_below_0p70";
end
if isAnchorGauge
    if ~audit.anchorHistoryComplete
        failure(end+1) = "anchor_history_missing_or_nonfinite";
    else
        if ~audit.anchorXMatchesConfig
            failure(end+1) = "anchor_history_coordinate_mismatch";
        end
        if ~audit.anchorReferenceConstant
            failure(end+1) = "anchor_history_reference_not_constant";
        end
        if ~audit.anchorSignPreserved
            failure(end+1) = "anchor_history_value_reference_sign";
        end
        if ~isfinite(specification.minimumCondition) || ...
                audit.minimumAnchorWallCondition < ...
                specification.minimumCondition
            failure(end+1) = "anchor_wall_condition";
        end
        if ~isfinite(specification.maximumNormalizedResidual) || ...
                audit.maximumAnchorNormalizedRateResidual > ...
                specification.maximumNormalizedResidual
            failure(end+1) = "anchor_wall_rate_residual";
        end
    end
end
if isWindowGauge
    if ~audit.windowHistoryComplete
        failure(end+1) = "window_history_missing_or_nonfinite";
    else
        if ~audit.windowCenterMatchesConfig
            failure(end+1) = "window_history_center_mismatch";
        end
        if ~audit.windowRadiusMatchesConfig
            failure(end+1) = "window_history_radius_mismatch";
        end
        if ~audit.windowSupportFixed || ~audit.windowSupportPositive || ...
                ~audit.windowSupportGeometryMatchesSpecification
            failure(end+1) = "window_history_support_not_fixed_positive";
        end
        if ~audit.windowMomentOrderMatchesSpecification
            failure(end+1) = "window_moment_order_mismatch";
        end
        if ~audit.windowReferenceConstant || ...
                ~audit.windowReferencePositive || ...
                ~audit.windowSignPreserved
            failure(end+1) = ...
                "window_history_reference_not_constant_positive";
        end
        if ~audit.windowValuePositive || ~audit.windowEnergyPositive || ...
                ~audit.windowMomentPositive || ...
                ~audit.windowRateScalePositive
            failure(end+1) = "window_history_value_not_positive";
        end
        if ~isfinite(specification.minimumCondition) || ...
                audit.minimumWindowCondition < ...
                specification.minimumCondition
            failure(end+1) = "window_condition";
        end
        if ~isfinite(specification.maximumNormalizedResidual) || ...
                audit.maximumWindowNormalizedRateResidual > ...
                specification.maximumNormalizedResidual
            failure(end+1) = "window_rate_residual";
        end
        if ~isfinite(specification.identityTolerance) || ...
                audit.maximumWindowMomentIdentityError > ...
                specification.identityTolerance
            failure(end+1) = "window_moment_identity";
        end
        if specification.momentOrder == 2 && ...
                audit.maximumWindowEnergyIdentityError > ...
                specification.identityTolerance
            failure(end+1) = "window_l2_energy_identity";
        end
    end
end
if isTemplateGauge
    if ~audit.templateHistoryComplete
        failure(end+1) = "template_history_missing_or_nonfinite";
    else
        if ~audit.templateCenterMatchesSpecification || ...
                ~audit.templateRadiusMatchesSpecification
            failure(end+1) = "template_fixed_geometry_mismatch";
        end
        if ~audit.templateSupportFixed || ...
                ~audit.templateSupportPositive || ...
                ~audit.templateSupportGeometryMatchesSpecification
            failure(end+1) = "template_support_not_fixed_positive";
        end
        if ~audit.templateReferenceConstant || ...
                ~audit.templateReferencePositive || ...
                ~audit.templateValuePositive || ...
                ~audit.templateSignPreserved
            failure(end+1) = "template_reference_sign_or_positivity";
        end
        if ~isfinite(specification.minimumCondition) || ...
                audit.minimumTemplateCondition < ...
                specification.minimumCondition
            failure(end+1) = "template_condition";
        end
        if ~isfinite(specification.maximumNormalizedResidual) || ...
                audit.minimumTemplateRateScale <= 0 || ...
                audit.maximumTemplateNormalizedRateResidual > ...
                specification.maximumNormalizedResidual
            failure(end+1) = "template_rate_residual";
        end
    end
end
if isBulkGauge
    if ~audit.bulkHistoryComplete
        failure(end+1) = "bulk_history_missing_or_nonfinite";
    else
        if ~audit.bulkCenterXMatchesSpecification || ...
                ~audit.bulkCenterYMatchesSpecification || ...
                ~audit.bulkRadiusMatchesSpecification
            failure(end+1) = "bulk_fixed_geometry_mismatch";
        end
        if ~audit.bulkSupportFixed || ~audit.bulkSupportPositive || ...
                ~audit.bulkSupportGeometryMatchesSpecification
            failure(end+1) = "bulk_support_not_fixed_positive";
        end
        if ~audit.bulkReferenceConstant || ...
                ~audit.bulkReferencePositive || ...
                ~audit.bulkValuePositive || ~audit.bulkEnergyPositive || ...
                ~audit.bulkSignPreserved
            failure(end+1) = "bulk_reference_sign_or_positivity";
        end
        if ~isfinite(specification.minimumCondition) || ...
                audit.minimumBulkCondition < ...
                specification.minimumCondition
            failure(end+1) = "bulk_condition";
        end
        if ~isfinite(specification.maximumNormalizedResidual) || ...
                audit.minimumBulkRateScale <= 0 || ...
                audit.maximumBulkNormalizedRateResidual > ...
                specification.maximumNormalizedResidual
            failure(end+1) = "bulk_rate_residual";
        end
        if ~isfinite(specification.identityTolerance) || ...
                audit.maximumBulkEnergyIdentityError > ...
                specification.identityTolerance
            failure(end+1) = "bulk_energy_identity";
        end
    end
end
if isQuadraticPeakGauge
    if ~audit.quadraticPeakHistoryComplete
        failure(end+1) = ...
            "quadratic_peak_history_missing_or_nonfinite";
    else
        minimumCondition = specification.minimumCondition;
        conditionValues = [ ...
            audit.minimumQuadraticPeakWeightCondition, ...
            audit.minimumQuadraticPeakVertexInteriorFraction, ...
            audit.minimumQuadraticPeakGlobalCondition];
        if ~isfinite(minimumCondition) || ...
                any(conditionValues < minimumCondition) || ...
                audit.minimumQuadraticPeakCurvatureCondition <= 0
            failure(end+1) = "quadratic_peak_condition";
        end
        if audit.minimumQuadraticPeakValue <= 0 || ...
                audit.maximumQuadraticPeakCurvature >= 0 || ...
                audit.minimumQuadraticPeakActiveNodeMargin < 0 || ...
                audit.minimumQuadraticPeakActiveSelectionMargin < ...
                -limits.contractNumericEpsMultiplier*eps(1) || ...
                audit.minimumQuadraticPeakRelativeSubgridLift < ...
                -limits.contractNumericEpsMultiplier*eps(1) || ...
                audit.maximumQuadraticPeakStationarityResidual > ...
                specification.maximumNormalizedResidual || ...
                audit.minimumQuadraticPeakRateScale <= 0 || ...
                audit.maximumQuadraticPeakCoordinateIdentityError > ...
                limits.contractNumericEpsMultiplier*eps(1)
            failure(end+1) = "quadratic_peak_geometry_or_sign";
        end
        if ~isfinite(specification.maximumNormalizedResidual) || ...
                audit.maximumQuadraticPeakNormalizedRateResidual > ...
                specification.maximumNormalizedResidual || ...
                audit.maximumQuadraticPeakReconstructedNormalizedResidual > ...
                specification.maximumNormalizedResidual || ...
                audit.maximumQuadraticPeakReportedResidualMismatch > ...
                specification.maximumNormalizedResidual
            failure(end+1) = "quadratic_peak_rate_residual";
        end
    end
end
if ~gaugeResidualFinite || ...
        audit.maximumGaugeResidual > audit.appliedGaugeResidualLimit
    failure(end+1) = "amplitude_gauge_residual";
end
audit.hardPassed = isempty(failure);
audit.hardFailureReason = char(strjoin(failure,','));
audit.productionEligible = audit.hardPassed && ...
    specification.productionCandidate;
if audit.productionEligible
    audit.artifactDisposition = 'production_candidate';
elseif audit.contract.legacySchema3ReadOnly
    audit.artifactDisposition = 'read_only_history';
elseif ~specification.productionCandidate
    audit.artifactDisposition = 'diagnostic_only';
else
    audit.artifactDisposition = 'failed_schema4_candidate';
end
end

function limits = hard_limits()
limits = struct( ...
    'targetRoundoffEpsMultiplier',100, ...
    'targetMinimumStepMultiplier',1, ...
    'tailPhysicalWidthCap',0.25, ...
    'tailPhysicalSpanFraction',0.40, ...
    'minimumTailRecords',8, ...
    'maximumExactGaugeResidual',1e-10, ...
    'maximumWidthCorrectionMagnitude',0, ...
    'minimumTerminalCoreCells',[14,14], ...
    'terminalSafetyFactorExclusiveUpperBound',0.70, ...
    'contractNumericEpsMultiplier',100, ...
    'requiredNodeCounts',[513,257], ...
    'requiredTransportAnchorX',1, ...
    'requiredOmegaGaugeWindowRadius',0.75, ...
    'requiredLocalizedConditionFloor',0.25, ...
    'requiredLengthScaleGain',1);
end

function audit = empty_audit(result,rawTrusted,trusted,trustChecks,limits, ...
        targetPhysicalTime,targetSource,isAnchorGauge,isWindowGauge, ...
        isTemplateGauge,isBulkGauge,isQuadraticPeakGauge, ...
        specification,lastTrusted)
audit = struct();
audit.schemaVersion = 4;
audit.kind = 'q256_exact_no_feedback_gauge_audit';
audit.caseId = result.metadata.caseId;
audit.cOmegaGauge = result.config.scaling.cOmegaGauge;
audit.lengthGauge = result.config.scaling.lengthGauge;
audit.transportAnchorX = result.config.scaling.transportAnchorX;
if isfield(result.config.scaling,'omegaGaugeWindowRadius')
    audit.omegaGaugeWindowRadius = ...
        result.config.scaling.omegaGaugeWindowRadius;
else
audit.omegaGaugeWindowRadius = NaN;
end
audit.candidateSpecification = specification;
audit.appliedMinimumCondition = specification.minimumCondition;
audit.appliedMaximumNormalizedResidual = ...
    specification.maximumNormalizedResidual;
audit.appliedIdentityTolerance = specification.identityTolerance;
audit.hardLimits = limits;
audit.targetPhysicalTime = targetPhysicalTime;
audit.targetPhysicalTimeSource = targetSource;
audit.solverMinimumTimeStep = result.config.time.minDt;
audit.targetToleranceRule = [ ...
    'targetRoundoffEpsMultiplier*eps(target) + ' ...
    'targetMinimumStepMultiplier*minDt*physicalClockSpeed'];
audit.targetPhysicalTimeTolerance = NaN;
audit.trustedEndpointPhysicalClockSpeed = NaN;
audit.rawTrustedMask = rawTrusted;
audit.trustedPrefixMask = trusted;
audit.trustChecks = trustChecks;
audit.trustedRecords = nnz(trusted);
audit.totalRecords = numel(trusted);
audit.terminalTrusted = ~isempty(trusted) && trusted(end);
audit.reachedTarget = false;
audit.stopReason = result.state.stopReason;
audit.steps = result.state.steps;
audit.rawPhysicalTime = result.state.physicalTime;
audit.rawCanonicalTau = result.state.canonicalTime;
audit.physicalTime = optional_indexed( ...
    result.history.common,'physicalTime',lastTrusted);
audit.canonicalTau = optional_indexed( ...
    result.history.common,'canonicalTau',lastTrusted);
audit.availablePhysicalSpan = NaN;
audit.requestedTailPhysicalWidth = NaN;
audit.tailPhysicalWidth = NaN;
audit.tailCanonicalWidth = NaN;
audit.tailRecords = 0;
audit.canonicalCL = empty_stats();
audit.canonicalCOmega = empty_stats();
audit.canonicalKappa = empty_stats();
audit.rescaledGradient = empty_stats();
audit.wallToGradientRatio = empty_stats();
audit.minimumTimeSpeed = NaN;
audit.maximumTimeSpeedDeviation = Inf;
audit.widthCorrectionHistoryComplete = false;
audit.maximumCanonicalWidthCorrection = Inf;
audit.maximumWidthCorrection = Inf;
audit.maximumWidthCorrectionMagnitude = Inf;
audit.lengthGaugeHistoryComplete = false;
audit.minimumLengthScaleGain = NaN;
audit.maximumLengthScaleGain = NaN;
audit.maximumAbsoluteWidthControlMode = Inf;
audit.lengthScaleGainExactlyOne = false;
audit.widthControlModeExactlyZero = false;
audit.maximumCanonicalCLNominalIdentityError = Inf;
audit.maximumCLNominalIdentityError = Inf;
audit.canonicalCLMatchesCanonicalCLNominal = false;
audit.canonicalCLMatchesCLNominal = false;
audit.maximumTransportAnchorXError = Inf;
audit.transportAnchorXFixedAtOne = false;
audit.maximumRawLengthGaugeResidual = Inf;
audit.maximumNormalizedLengthGaugeResidual = Inf;
audit.maximumRawCanonicalTransportAnchorVelocity = Inf;
audit.maximumNormalizedCanonicalTransportAnchorVelocity = Inf;
audit.maximumRawTransportAnchorVelocity = Inf;
audit.maximumNormalizedTransportAnchorVelocity = Inf;
audit.maximumNormalizedTransportAnchorDefect = Inf;
audit.transportAnchorNormalizationFinite = false;
audit.minimumLengthGaugeResidualScale = -Inf;
audit.minimumCanonicalTransportAnchorScale = -Inf;
audit.minimumTransportAnchorScale = -Inf;
audit.transportAnchorNormalizationRule = ...
    'abs(residual)/max(abs(rate*X),abs(residual-rate*X),realmin)';
audit.maximumGaugeResidual = Inf;
audit.appliedGaugeResidualLimit = limits.maximumExactGaugeResidual;
audit.anchorHistoryRequired = isAnchorGauge;
audit.anchorHistoryComplete = ~isAnchorGauge;
audit.minimumAnchorWallCondition = NaN;
audit.maximumBoundedAnchorRateResidual = NaN;
audit.maximumAnchorNormalizedRateResidual = Inf;
audit.anchorReferenceValue = NaN;
audit.anchorCoordinateTolerance = NaN;
audit.maximumAnchorCoordinateError = NaN;
audit.anchorXMatchesConfig = ~isAnchorGauge;
audit.anchorReferenceTolerance = NaN;
audit.maximumAnchorReferenceDrift = NaN;
audit.anchorReferenceConstant = ~isAnchorGauge;
audit.anchorSignPreserved = ~isAnchorGauge;
audit.maximumAnchorRelativeDrift = NaN;
audit.finalAnchorRelativeDrift = NaN;
audit.windowHistoryRequired = isWindowGauge;
audit.windowHistoryComplete = ~isWindowGauge;
audit.minimumWindowCondition = NaN;
audit.maximumWindowRawRateResidual = Inf;
audit.maximumWindowNormalizedRateResidual = Inf;
audit.maximumWindowRateResidual = Inf;
audit.minimumWindowRateScale = -Inf;
audit.windowReferenceValue = NaN;
audit.windowCenterTolerance = NaN;
audit.maximumWindowCenterError = NaN;
audit.windowCenterMatchesConfig = ~isWindowGauge;
audit.windowCenterMatchesSpecification = ~isWindowGauge;
audit.windowRadiusTolerance = NaN;
audit.maximumWindowRadiusError = NaN;
audit.windowRadiusMatchesConfig = ~isWindowGauge;
audit.windowRadiusMatchesSpecification = ~isWindowGauge;
audit.windowReferenceTolerance = NaN;
audit.maximumWindowReferenceDrift = NaN;
audit.windowReferenceConstant = ~isWindowGauge;
audit.windowReferencePositive = ~isWindowGauge;
audit.windowValuePositive = ~isWindowGauge;
audit.windowSignPreserved = ~isWindowGauge;
audit.windowEnergyPositive = ~isWindowGauge;
audit.windowMomentPositive = ~isWindowGauge;
audit.windowRateScalePositive = ~isWindowGauge;
audit.windowSupportPositive = ~isWindowGauge;
audit.windowSupportFixed = ~isWindowGauge;
audit.windowSupportPoints = NaN;
audit.windowSupportGeometryMatchesSpecification = ~isWindowGauge;
audit.maximumWindowSupportGeometryError = NaN;
audit.windowMomentOrderMatchesSpecification = ~isWindowGauge;
audit.maximumWindowMomentIdentityError = Inf;
audit.maximumWindowEnergyIdentityError = Inf;
audit.maximumWindowRelativeDrift = NaN;
audit.finalWindowRelativeDrift = NaN;
audit.templateHistoryRequired = isTemplateGauge;
audit.templateHistoryComplete = ~isTemplateGauge;
audit.minimumTemplateCondition = NaN;
audit.maximumTemplateRawRateResidual = Inf;
audit.maximumTemplateNormalizedRateResidual = Inf;
audit.minimumTemplateRateScale = -Inf;
audit.templateCenterMatchesSpecification = ~isTemplateGauge;
audit.maximumTemplateCenterError = NaN;
audit.templateCenterTolerance = NaN;
audit.templateRadiusMatchesSpecification = ~isTemplateGauge;
audit.maximumTemplateRadiusError = NaN;
audit.templateRadiusTolerance = NaN;
audit.templateSupportPositive = ~isTemplateGauge;
audit.templateSupportFixed = ~isTemplateGauge;
audit.templateSupportPoints = NaN;
audit.templateSupportGeometryMatchesSpecification = ~isTemplateGauge;
audit.maximumTemplateSupportGeometryError = NaN;
audit.templateReferencePositive = ~isTemplateGauge;
audit.templateReferenceConstant = ~isTemplateGauge;
audit.templateValuePositive = ~isTemplateGauge;
audit.templateSignPreserved = ~isTemplateGauge;
audit.maximumTemplateRelativeDrift = NaN;
audit.finalTemplateRelativeDrift = NaN;
audit.bulkHistoryRequired = isBulkGauge;
audit.bulkHistoryComplete = ~isBulkGauge;
audit.minimumBulkCondition = NaN;
audit.maximumBulkRawRateResidual = Inf;
audit.maximumBulkNormalizedRateResidual = Inf;
audit.minimumBulkRateScale = -Inf;
audit.bulkCenterXMatchesSpecification = ~isBulkGauge;
audit.maximumBulkCenterXError = NaN;
audit.bulkCenterXTolerance = NaN;
audit.bulkCenterYMatchesSpecification = ~isBulkGauge;
audit.maximumBulkCenterYError = NaN;
audit.bulkCenterYTolerance = NaN;
audit.bulkRadiusMatchesSpecification = ~isBulkGauge;
audit.maximumBulkRadiusError = NaN;
audit.bulkRadiusTolerance = NaN;
audit.bulkSupportPositive = ~isBulkGauge;
audit.bulkSupportFixed = ~isBulkGauge;
audit.bulkSupportPoints = [NaN,NaN];
audit.bulkSupportGeometryMatchesSpecification = ~isBulkGauge;
audit.maximumBulkSupportGeometryError = NaN;
audit.bulkReferencePositive = ~isBulkGauge;
audit.bulkReferenceConstant = ~isBulkGauge;
audit.bulkValuePositive = ~isBulkGauge;
audit.bulkEnergyPositive = ~isBulkGauge;
audit.bulkSignPreserved = ~isBulkGauge;
audit.maximumBulkEnergyIdentityError = Inf;
audit.maximumBulkRelativeDrift = NaN;
audit.finalBulkRelativeDrift = NaN;
audit.quadraticPeakHistoryRequired = isQuadraticPeakGauge;
audit.quadraticPeakHistoryComplete = ~isQuadraticPeakGauge;
audit.minimumQuadraticPeakValue = NaN;
audit.maximumQuadraticPeakCurvature = Inf;
audit.minimumQuadraticPeakCurvatureCondition = NaN;
audit.minimumQuadraticPeakWeightCondition = NaN;
audit.minimumQuadraticPeakVertexInteriorFraction = NaN;
audit.minimumQuadraticPeakActiveNodeMargin = NaN;
audit.minimumQuadraticPeakActiveSelectionMargin = NaN;
audit.quadraticPeakStencilSwitchCount = NaN;
audit.minimumQuadraticPeakGlobalCondition = NaN;
audit.minimumQuadraticPeakRateScale = -Inf;
audit.maximumQuadraticPeakNormalizedRateResidual = Inf;
audit.maximumQuadraticPeakReconstructedNormalizedResidual = Inf;
audit.maximumQuadraticPeakReportedResidualMismatch = Inf;
audit.maximumQuadraticPeakCoordinateIdentityError = Inf;
audit.minimumQuadraticPeakRelativeSubgridLift = NaN;
audit.maximumQuadraticPeakStationarityResidual = Inf;
audit.finalCoreCells = [NaN,NaN];
audit.finalSafetyFactor = NaN;
audit.maximumSafetyFactor = Inf;
audit.finalPhysicalRhoX = NaN;
audit.finalPhysicalGradient = NaN;
audit.finalFarVelocityRatio = NaN;
audit.finalFarSourceRatio = NaN;
audit.contract = struct();
audit.hardPassed = false;
audit.productionEligible = false;
if result.config.schemaVersion == 3
    audit.artifactDisposition = 'read_only_history';
else
    audit.artifactDisposition = 'failed_schema4_candidate';
end
audit.hardFailureReason = '';
end

function contract = experiment_contract(result,mesh,limits)
config = result.config;
allRecords = true(numel(result.history.common.t),1);
[remeshCount,remeshCountFinite] = required_masked_series( ...
    mesh,'remeshCount',allRecords);
if remeshCountFinite
    historyMaximumRemeshCount = max(remeshCount);
else
    historyMaximumRemeshCount = NaN;
end
numericTolerance = limits.contractNumericEpsMultiplier*eps(1);
contract = struct();
contract.expected = struct( ...
    'configSchemaVersion',4, ...
    'nodeCounts',limits.requiredNodeCounts, ...
    'spatialDiscretization','high_order', ...
    'transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54', ...
    'remeshTransferScheme','high_order', ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'remeshCount',0,'rescalingMode','dynamic', ...
    'dynamicScaleGeometry','isotropic', ...
    'lengthGauge','transport_anchor', ...
    'transportAnchorX',limits.requiredTransportAnchorX, ...
    'omegaGaugeWindowRadius', ...
        limits.requiredOmegaGaugeWindowRadius, ...
    'pointStrainConditionFloor', ...
        limits.requiredLocalizedConditionFloor, ...
    'lengthScaleGain',limits.requiredLengthScaleGain, ...
    'widthGaugeGain',0,'omegaGaugeGain',0,'travelingWaveGain',0, ...
    'scalingContract','exact_gauge_no_feedback_v1', ...
    'retiredAdaptiveFields',{retired_adaptive_fields()});
contract.actual = struct( ...
    'configSchemaVersion',config.schemaVersion, ...
    'nodeCounts',[numel(result.grid.x),numel(result.grid.y)], ...
    'configuredNodeCounts',[config.grid.nx,config.grid.ny], ...
    'spatialDiscretization',config.transport.spatialDiscretization, ...
    'transportScheme',config.transport.transportScheme, ...
    'timeIntegrator',config.time.timeIntegrator, ...
    'remeshTransferScheme',config.remesh.remeshTransferScheme, ...
    'adaptiveRemesh',config.remesh.adaptiveRemesh, ...
    'initialAnalyticRemesh',config.remesh.initialAnalyticRemesh, ...
    'terminalRemeshCount',result.grid.remeshCount, ...
    'historyMaximumRemeshCount',historyMaximumRemeshCount, ...
    'rescalingMode',config.scaling.rescalingMode, ...
    'dynamicScaleGeometry',config.scaling.dynamicScaleGeometry, ...
    'lengthGauge',config.scaling.lengthGauge, ...
    'transportAnchorX',config.scaling.transportAnchorX, ...
    'omegaGaugeWindowRadius',optional_scalar( ...
        config.scaling,'omegaGaugeWindowRadius'), ...
    'pointStrainConditionFloor',optional_scalar( ...
        config.scaling,'pointStrainConditionFloor'), ...
    'lengthScaleGain',config.scaling.lengthScaleGain, ...
    'widthGaugeGain',optional_scalar(config.scaling,'widthGaugeGain'), ...
    'omegaGaugeGain',optional_scalar(config.scaling,'omegaGaugeGain'), ...
    'travelingWaveGain',optional_scalar( ...
        config.scaling,'travelingWaveGain'), ...
    'scalingContract',optional_text(config.scaling,'scalingContract'), ...
    'retiredAdaptiveFieldsPresent',{retired_fields_present( ...
        config.scaling)});
contract.q256Grid = ...
    isequal(contract.actual.nodeCounts,limits.requiredNodeCounts) && ...
    isequal(contract.actual.configuredNodeCounts,limits.requiredNodeCounts);
contract.fixedGrid = ~config.remesh.adaptiveRemesh && ...
    ~config.remesh.initialAnalyticRemesh && ...
    result.grid.remeshCount == 0 && remeshCountFinite && ...
    all(remeshCount == 0);
contract.fourthOrderTuple = ...
    strcmp(config.transport.spatialDiscretization,'high_order') && ...
    strcmp(config.transport.transportScheme,'weno5_fd') && ...
    strcmp(config.time.timeIntegrator,'ssprk54') && ...
    strcmp(config.remesh.remeshTransferScheme,'high_order');
contract.dynamicIsotropic = ...
    strcmp(config.scaling.rescalingMode,'dynamic') && ...
    strcmp(config.scaling.dynamicScaleGeometry,'isotropic');
contract.transportAnchor = ...
    strcmp(config.scaling.lengthGauge,'transport_anchor') && ...
    abs(config.scaling.transportAnchorX- ...
    limits.requiredTransportAnchorX) <= numericTolerance && ...
    abs(config.scaling.lengthScaleGain- ...
    limits.requiredLengthScaleGain) <= numericTolerance;
contract.fixedLocalizationGeometry = abs(optional_scalar( ...
    config.scaling,'omegaGaugeWindowRadius') - ...
    limits.requiredOmegaGaugeWindowRadius) <= numericTolerance && ...
    abs(optional_scalar(config.scaling,'pointStrainConditionFloor') - ...
    limits.requiredLocalizedConditionFloor) <= numericTolerance;
contract.exactNoFeedbackGains = ...
    config.scaling.lengthScaleGain == 1 && ...
    optional_scalar(config.scaling,'widthGaugeGain') == 0 && ...
    optional_scalar(config.scaling,'omegaGaugeGain') == 0 && ...
    optional_scalar(config.scaling,'travelingWaveGain') == 0;
contract.noRetiredAdaptiveOptions = ...
    isempty(contract.actual.retiredAdaptiveFieldsPresent);
contract.exactScalingContract = config.schemaVersion == 4 && ...
    strcmp(contract.actual.scalingContract, ...
    'exact_gauge_no_feedback_v1');
contract.directCanonicalClock = config.schemaVersion == 4 && ...
    ~isfield(config.scaling,'maxDynamicRate');
contract.legacySchema3ReadOnly = config.schemaVersion == 3;
contract.productionEligibleSchema = config.schemaVersion == 4 && ...
    contract.exactScalingContract && contract.exactNoFeedbackGains && ...
    contract.noRetiredAdaptiveOptions && contract.directCanonicalClock;
contract.passed = contract.q256Grid && contract.fixedGrid && ...
    contract.fourthOrderTuple && contract.dynamicIsotropic && ...
    contract.transportAnchor && contract.fixedLocalizationGeometry && ...
    contract.productionEligibleSchema;
end

function failure = contract_failures(contract)
failure = strings(0,1);
if ~contract.q256Grid
    failure(end+1) = "q256_grid_contract";
end
if ~contract.fixedGrid
    failure(end+1) = "fixed_grid_contract";
end
if ~contract.fourthOrderTuple
    failure(end+1) = "fourth_order_tuple_contract";
end
if ~contract.dynamicIsotropic
    failure(end+1) = "dynamic_isotropic_contract";
end
if ~contract.transportAnchor
    failure(end+1) = "transport_anchor_x1_contract";
end
if ~contract.fixedLocalizationGeometry
    failure(end+1) = "fixed_localization_radius_contract";
end
if contract.legacySchema3ReadOnly
    failure(end+1) = "legacy_schema3_read_only_history";
elseif ~contract.productionEligibleSchema
    failure(end+1) = "schema4_exact_no_feedback_contract";
end
if ~contract.exactNoFeedbackGains
    failure(end+1) = "exact_no_feedback_gain_contract";
end
if ~contract.noRetiredAdaptiveOptions
    failure(end+1) = "retired_adaptive_option_present";
end
if ~contract.directCanonicalClock
    failure(end+1) = "direct_canonical_clock_contract";
end
end

function [specification,status] = selected_specification(result,gaugeName)
status = struct('source','','tournament',false,'valid',false, ...
    'failure','');
caseMetadata = struct();
if isstruct(result.metadata) && isscalar(result.metadata) && ...
        isfield(result.metadata,'caseMetadata')
    caseMetadata = result.metadata.caseMetadata;
end
if isstruct(caseMetadata) && isscalar(caseMetadata) && ...
        isfield(caseMetadata,'gaugeTournament')
    % A tournament artifact is judged by its frozen preregistration.  If
    % that record is malformed, fail it instead of substituting live policy.
    status.source = ...
        'metadata.caseMetadata.gaugeTournament.candidateSpecification';
    status.tournament = true;
    [specification,status.valid,status.failure] = ...
        tournament_specification(caseMetadata.gaugeTournament,gaugeName);
    return;
end

status.source = 'live_registry_non_tournament_fallback';
registry = ipm_gaugelab_candidate_specifications(true);
index = find(strcmp({registry.cOmegaGauge},gaugeName));
if numel(index) ~= 1
    error('ipm:GaugeLabUnregisteredGauge', ...
        'Gauge "%s" is not registered for gauge-lab assessment.',gaugeName);
end
specification = registry(index);
[status.valid,status.failure] = ...
    candidate_specification_format(specification,gaugeName);
end

function [specification,valid,failureText] = ...
        tournament_specification(tournament,gaugeName)
specification = blank_candidate_specification(gaugeName);
failure = strings(0,1);
required = {'schemaVersion','stage','candidate','baseCaseName', ...
    'onlyAmplitudeGaugeVaried','candidateSpecification', ...
    'configurationContract','physicalCovarianceMethod'};
if ~isstruct(tournament) || ~isscalar(tournament)
    failure(end+1) = "tournament_metadata_not_scalar_struct";
    valid = false;
    failureText = char(strjoin(failure,','));
    return;
end
missing = required(~isfield(tournament,required));
if ~isempty(missing)
    failure(end+1) = "tournament_metadata_missing_fields:" + ...
        strjoin(string(missing),'+');
end
if ~isfield(tournament,'schemaVersion') || ...
        ~isnumeric(tournament.schemaVersion) || ...
        ~isreal(tournament.schemaVersion) || ...
        ~isscalar(tournament.schemaVersion) || ...
        ~isfinite(tournament.schemaVersion) || ...
        tournament.schemaVersion ~= 4
    failure(end+1) = "tournament_metadata_schema";
end
[candidateName,candidateValid] = required_text(tournament,'candidate');
if ~candidateValid || ~strcmp(candidateName,gaugeName)
    failure(end+1) = "tournament_candidate_gauge_mismatch";
end
textFields = {'stage','baseCaseName','physicalCovarianceMethod'};
textFailures = strings(numel(textFields),1);
textFailureCount = 0;
for index = 1:numel(textFields)
    [~,textValid] = required_text(tournament,textFields{index});
    if ~textValid
        textFailureCount = textFailureCount+1;
        textFailures(textFailureCount) = ...
            "tournament_metadata_text_format:" + string(textFields{index});
    end
end
failure = [failure;textFailures(1:textFailureCount)];
[configurationContract,contractValid] = ...
    required_text(tournament,'configurationContract');
if ~contractValid || ...
        ~strcmp(configurationContract,'exact_gauge_no_feedback_v1')
    failure(end+1) = "tournament_scaling_contract";
end
if ~isfield(tournament,'onlyAmplitudeGaugeVaried') || ...
        ~islogical(tournament.onlyAmplitudeGaugeVaried) || ...
        ~isscalar(tournament.onlyAmplitudeGaugeVaried) || ...
        ~tournament.onlyAmplitudeGaugeVaried
    failure(end+1) = "tournament_variation_contract";
end
if isfield(tournament,'candidateSpecification')
    candidate = tournament.candidateSpecification;
    [formatValid,formatFailure] = ...
        candidate_specification_format(candidate,gaugeName);
    if formatValid
        specification = candidate;
    else
        failure(end+1) = "candidate_specification_format:" + ...
            string(formatFailure);
    end
end
valid = isempty(failure);
failureText = char(strjoin(failure,','));
end

function [valid,failureText] = ...
        candidate_specification_format(specification,gaugeName)
failure = strings(0,1);
required = {'cOmegaGauge','role','family','momentOrder', ...
    'fixedSupport','minimumCondition','maximumNormalizedResidual', ...
    'identityTolerance','productionCandidate'};
if ~isstruct(specification) || ~isscalar(specification)
    valid = false;
    failureText = 'not_scalar_struct';
    return;
end
if ~isequal(sort(fieldnames(specification)),sort(required(:)))
    valid = false;
    failureText = 'field_set';
    return;
end
[specifiedGauge,gaugeValid] = ...
    required_text(specification,'cOmegaGauge');
if ~gaugeValid || ~strcmp(specifiedGauge,gaugeName)
    failure(end+1) = "gauge_mismatch";
end
for name = {'role','family'}
    [~,textValid] = required_text(specification,name{1});
    if ~textValid
        failure(end+1) = "text_format:" + string(name{1}); %#ok<AGROW>
    end
end
if ~islogical(specification.productionCandidate) || ...
        ~isscalar(specification.productionCandidate)
    failure(end+1) = "production_candidate_format";
end
numericFormat = ...
    valid_optional_real_scalar(specification.momentOrder) && ...
    valid_optional_real_scalar(specification.minimumCondition) && ...
    valid_optional_real_scalar(specification.identityTolerance) && ...
    isnumeric(specification.maximumNormalizedResidual) && ...
    isreal(specification.maximumNormalizedResidual) && ...
    isscalar(specification.maximumNormalizedResidual) && ...
    isfinite(specification.maximumNormalizedResidual) && ...
    specification.maximumNormalizedResidual > 0;
if ~numericFormat
    failure(end+1) = "numeric_threshold_format";
end
supportFormat = candidate_support_format(specification.fixedSupport);
if ~supportFormat
    failure(end+1) = "fixed_support_format";
end
if ~(numericFormat && supportFormat)
    valid = false;
    failureText = char(strjoin(failure,','));
    return;
end

localized = any(strcmp(gaugeName,{'anchor_wall_slope', ...
    'anchor_wall_strain','anchor_wall_window_l2', ...
    'anchor_wall_window_l4','anchor_wall_template_projection', ...
    'anchor_bulk_gradient_l2'}));
if localized
    condition = specification.minimumCondition;
    support = specification.fixedSupport;
    if ~isfinite(condition) || condition <= 0 || condition > 1 || ...
            ~all(isfinite([support.center,support.radius, ...
            support.supportX,support.supportY])) || support.radius <= 0
        failure(end+1) = "localized_contract_format";
    end
end
expectedOrder = expected_moment_order(gaugeName);
if isfinite(expectedOrder)
    if specification.momentOrder ~= expectedOrder
        failure(end+1) = "moment_order";
    end
elseif ~isnan(specification.momentOrder)
    failure(end+1) = "unexpected_moment_order";
end
identityRequired = any(strcmp(gaugeName, ...
    {'anchor_wall_window_l2','anchor_wall_window_l4', ...
    'anchor_bulk_gradient_l2'}));
if identityRequired
    if ~isfinite(specification.identityTolerance) || ...
            specification.identityTolerance <= 0
        failure(end+1) = "identity_tolerance_format";
    end
elseif ~isnan(specification.identityTolerance)
    failure(end+1) = "unexpected_identity_tolerance";
end
valid = isempty(failure);
failureText = char(strjoin(failure,','));
end

function specification = blank_candidate_specification(gaugeName)
support = struct('center',[NaN,NaN],'radius',NaN, ...
    'supportX',[NaN,NaN],'supportY',[NaN,NaN]);
specification = struct('cOmegaGauge',gaugeName,'role','invalid', ...
    'family','invalid','momentOrder',NaN,'fixedSupport',support, ...
    'minimumCondition',NaN,'maximumNormalizedResidual',NaN, ...
    'identityTolerance',NaN,'productionCandidate',false);
end

function [value,valid] = required_text(group,name)
value = optional_text(group,name);
valid = ~isempty(value);
end

function valid = valid_optional_real_scalar(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    ~isinf(value);
end

function valid = candidate_support_format(support)
required = {'center','radius','supportX','supportY'};
valid = isstruct(support) && isscalar(support) && ...
    isequal(sort(fieldnames(support)),sort(required(:)));
if ~valid
    return;
end
vectors = {'center','supportX','supportY'};
for index = 1:numel(vectors)
    value = support.(vectors{index});
    valid = valid && isnumeric(value) && isreal(value) && ...
        isvector(value) && numel(value) == 2 && all(~isinf(value));
end
valid = valid && isnumeric(support.radius) && ...
    isreal(support.radius) && isscalar(support.radius) && ...
    ~isinf(support.radius);
end

function order = expected_moment_order(gaugeName)
switch gaugeName
    case {'anchor_wall_window_l2','anchor_bulk_gradient_l2'}
        order = 2;
    case 'anchor_wall_window_l4'
        order = 4;
    otherwise
        order = NaN;
end
end

function names = retired_adaptive_fields()
names = {'maxDynamicRate','adaptiveLengthScaling', ...
    'widthExpansionStrength', ...
    'widthContractionOnset','widthContractionStrength', ...
    'maxWidthRateCorrection'};
end

function present = retired_fields_present(scaling)
names = retired_adaptive_fields();
present = names(isfield(scaling,names));
end

function value = optional_scalar(group,name)
value = NaN;
if isstruct(group) && isscalar(group) && isfield(group,name)
    candidate = group.(name);
    if isnumeric(candidate) && isreal(candidate) && isscalar(candidate) && ...
            isfinite(candidate)
        value = double(candidate);
    end
end
end

function value = optional_text(group,name)
value = '';
if ~isstruct(group) || ~isscalar(group) || ~isfield(group,name)
    return;
end
candidate = group.(name);
if isstring(candidate) && isscalar(candidate) && ~ismissing(candidate)
    value = char(candidate);
elseif ischar(candidate) && isrow(candidate)
    value = candidate;
end
end

function [matches,maximumError,tolerance] = fixed_numeric_series( ...
        values,expected,limits)
tolerance = limits.contractNumericEpsMultiplier * ...
    eps(max(1,abs(expected)));
maximumError = max(abs(values-expected));
matches = isfinite(maximumError) && maximumError <= tolerance;
end

function fixed = integer_constant_series(values)
fixed = ~isempty(values) && all(values >= 0) && ...
    all(values == round(values)) && all(values == values(1));
end

function fixed = constant_series(values,limits)
if isempty(values)
    fixed = false;
    return;
end
tolerance = limits.contractNumericEpsMultiplier * ...
    eps(max(1,abs(values(1))));
fixed = max(abs(values-values(1))) <= tolerance;
end

function [matches,maximumError] = fixed_support_geometry( ...
        centerX,centerY,radius,support,limits,kind)
switch kind
    case 'wall'
        actual = [centerX-radius,centerX+radius, ...
            zeros(size(centerX)),zeros(size(centerX))];
    case 'bulk'
        actual = [centerX-radius,centerX+radius, ...
            max(0,centerY),centerY+radius];
    otherwise
        error('ipm:GaugeLabSupportKind', ...
            'Unknown fixed support geometry "%s".',kind);
end
expected = [support.supportX,support.supportY];
maximumError = max(abs(actual-expected),[],'all');
tolerance = limits.contractNumericEpsMultiplier * ...
    eps(max([1,abs(expected)]));
matches = isfinite(maximumError) && maximumError <= tolerance;
end

function errorValue = relative_identity_error(left,right)
scale = max(abs(left),abs(right));
scale = max(scale,realmin);
errorValue = abs(left-right)./scale;
end

function [normalized,scale] = normalized_anchor_residual( ...
        residual,rate,anchorX)
% RESIDUAL = anchorU1 + rate*X, so both denominator terms are recorded.
anchorU1 = residual-rate.*anchorX;
scale = max(abs(rate.*anchorX),abs(anchorU1));
scale = max(scale,realmin);
normalized = abs(residual)./scale;
end

function stats = selected_tail_stats(t,y,tail,valid)
if ~valid || numel(t) ~= numel(tail) || numel(y) ~= numel(tail)
    stats = empty_stats();
    return;
end
stats = tail_stats(t(tail),y(tail));
end

function stats = tail_stats(t,y)
t = column(t);
y = column(y);
stats = empty_stats();
stats.records = numel(t);
if numel(t) < 2 || any(~isfinite(t)) || any(~isfinite(y)) || ...
        any(diff(t) <= 0)
    return;
end
w = trapezoid_weights(t);
weightSum = sum(w);
meanValue = sum(w.*y)/weightSum;
centeredTime = t-sum(w.*t)/weightSum;
design = [ones(size(t)),centeredTime];
weightedDesign = design.*sqrt(w);
weightedValues = y.*sqrt(w);
coefficients = weightedDesign\weightedValues;
fitted = design*coefficients;
residual = y-fitted;
scale = max(abs(meanValue),eps);
stats.mean = meanValue;
stats.standardDeviation = sqrt(sum(w.*(y-meanValue).^2)/weightSum);
stats.coefficientOfVariation = stats.standardDeviation/scale;
stats.slope = coefficients(2);
stats.relativeTrendAcrossWindow = ...
    coefficients(2)*(t(end)-t(1))/scale;
stats.detrendedRelativeRms = ...
    sqrt(sum(w.*residual.^2)/weightSum)/scale;
stats.minimum = min(y);
stats.maximum = max(y);
stats.start = y(1);
stats.finish = y(end);
stats.window = t(end)-t(1);
end

function finite = stats_finite(stats)
names = setdiff(fieldnames(stats),{'records'},'stable');
finite = stats.records > 0;
for index = 1:numel(names)
    finite = finite && isfinite(stats.(names{index}));
end
end

function stats = empty_stats()
stats = struct('records',0,'mean',NaN,'standardDeviation',NaN, ...
    'coefficientOfVariation',NaN,'slope',NaN, ...
    'relativeTrendAcrossWindow',NaN,'detrendedRelativeRms',NaN, ...
    'minimum',NaN,'maximum',NaN,'start',NaN,'finish',NaN, ...
    'window',NaN);
end

function weights = trapezoid_weights(t)
weights = zeros(size(t));
weights(1) = (t(2)-t(1))/2;
weights(end) = (t(end)-t(end-1))/2;
if numel(t) > 2
    weights(2:end-1) = (t(3:end)-t(1:end-2))/2;
end
end

function [values,valid] = required_masked_series(group,name,mask)
values = [];
valid = isstruct(group) && isscalar(group) && isfield(group,name);
if ~valid
    return;
end
series = group.(name);
valid = isnumeric(series) && isreal(series) && isvector(series) && ...
    numel(series) == numel(mask);
if ~valid
    return;
end
values = column(series(mask));
valid = ~isempty(values) && all(isfinite(values));
end

function value = optional_indexed(group,name,index)
value = NaN;
if isempty(index) || ~isfield(group,name)
    return;
end
series = group.(name);
if ~isnumeric(series) || ~isreal(series) || ...
        ~isvector(series) || numel(series) < index
    return;
end
value = series(index);
end

function value = column(value)
value = double(value(:));
end
