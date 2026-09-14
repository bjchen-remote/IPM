function [c_l,c_omega,c_r,baseRhs,details,transportU1,transportU2] = ...
    isotropicGauge(rho,physicalRhs,flow,ops,details)
%IPM.EVOLVE.ISOTROPICGAUGE Compute exact isotropic gauge rates and residuals.

if ~ops.rescaling.enabled
    c_l = 0;
    c_omega = 0;
    c_r = 0;
    baseRhs = physicalRhs;
    transportU1 = flow.u1;
    transportU2 = flow.u2;
    return;
end

r = ops.rescaling;
trackedPeakX = details.trackedPeakX;
trackedPeakU1 = interp1(ops.x,flow.u1(1,:),trackedPeakX,'linear');
needsStrain = strcmp(r.lengthGauge,'local_strain') || ...
    any(strcmp(r.cOmegaGauge,{'strain_point','anchor_wall_slope', ...
        'anchor_wall_window_l2','anchor_wall_window_l4', ...
        'anchor_wall_template_projection','anchor_bulk_gradient_l2'}));
if needsStrain
    u1X = flow.u1*ipm.mesh.oddDx(ops)';
else
    u1X = [];
end
peakLocationRhoXX = NaN;
peakLocationRhoXXScale = NaN;
peakLocationPhaseDefect = NaN;
peakLocationResponseCondition = NaN;
peakLocationFullPhaseResidual = NaN;
peakLocationNormalizedFullPhaseResidual = NaN;
switch r.lengthGauge
    case 'omega_peak_location'
        rhoXX = (rho*ops.Dx')*ipm.mesh.oddDx(ops)';
        rhoXXX = rhoXX*ops.Dx';
        physicalRhsXX = (physicalRhs*ops.Dx')*ipm.mesh.oddDx(ops)';
        peakRhoXXX = interp1( ...
            ops.x,rhoXXX(1,:),trackedPeakX,'linear');
        peakRhoXX = interp1( ...
            ops.x,rhoXX(1,:),trackedPeakX,'linear');
        peakPhysicalRhsXX = interp1( ...
            ops.x,physicalRhsXX(1,:),trackedPeakX,'linear');
        peakLocationResponse = trackedPeakX*peakRhoXXX;
        responseScale = max(abs(ops.x.*rhoXXX(1,:)));
        if abs(peakLocationResponse) <= 1e-8*max(responseScale,eps)
            error('ipm:DegenerateLengthGauge', ...
                ['The wall-omega peak-location gauge has negligible ' ...
                'curvature response.']);
        end
        c_lNominal = peakPhysicalRhsXX/peakLocationResponse;
        peakLocationRhoXX = peakRhoXX;
        peakLocationRhoXXScale = max(abs(rhoXX(1,:)));
        peakLocationPhaseDefect = abs(peakRhoXX)/ ...
            max(peakLocationRhoXXScale,realmin);
        peakLocationResponseCondition = abs(peakLocationResponse)/ ...
            max(responseScale,realmin);
        lengthGaugeStrain = peakLocationResponse;
        lengthGaugeResidual = peakPhysicalRhsXX- ...
            c_lNominal*peakLocationResponse;
    case 'wall_omega_width'
        currentWidth = details.wallWidthCurrent;
        lengthGaugeStrain = currentWidth;
        if lengthGaugeStrain <= 100*eps(max(abs(ops.x)))
            error('ipm:DegenerateLengthGauge', ...
                ['The wall-omega half-height width is too small for ' ...
                'normalization.']);
        end
        c_lNominal = -details.wallWidthRate/lengthGaugeStrain;
        lengthGaugeResidual = details.wallWidthRate + ...
            c_lNominal*lengthGaugeStrain;
    case 'symmetry_peak'
        originX = ops.x(r.originIndex);
        originU1 = flow.u1(1,r.originIndex);
        peakSeparation = trackedPeakX-originX;
        if peakSeparation <= 100*eps(max(abs(ops.x)))
            error('ipm:DegenerateLengthGauge', ...
                ['The tracked positive omega peak is too close to the ' ...
                'symmetry origin.']);
        end
        lengthGaugeStrain = (trackedPeakU1-originU1)/peakSeparation;
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = trackedPeakU1-originU1 + ...
            c_lNominal*peakSeparation;
    case 'transport_anchor'
        anchorX = r.transportAnchorX;
        anchorU1 = interp1(ops.x,flow.u1(1,:),anchorX,'linear');
        c_lNominal = -anchorU1/anchorX;
        lengthGaugeStrain = anchorX;
        lengthGaugeResidual = anchorU1+c_lNominal*anchorX;
    case 'wall_density_width'
        lowerX = details.wallDensityLowerX;
        upperX = details.wallDensityUpperX;
        lengthGaugeWidth = details.wallDensityWidth;
        lowerU1 = interp1(ops.x,flow.u1(1,:),lowerX,'linear');
        upperU1 = interp1(ops.x,flow.u1(1,:),upperX,'linear');
        lengthGaugeStrain = (upperU1-lowerU1)/lengthGaugeWidth;
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = upperU1-lowerU1 + ...
            c_lNominal*lengthGaugeWidth;
    case 'local_strain'
        lengthGaugeStrain = interp1( ...
            ops.x,u1X(1,:),trackedPeakX,'linear');
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = lengthGaugeStrain+c_lNominal;
    otherwise
        error('ipm:BadLengthGaugeRuntime', ...
            'Unknown runtime length gauge: %s.',r.lengthGauge);
end

c_l = c_lNominal;
widthCorrection = 0;
widthControlMode = 0;
resolutionFactor = details.resolutionFactor;
safetyFactor = details.safetyFactor;
peakExpansionRatio = details.peakExpansionRatio;
if strcmp(ops.symmetryMode,'double_odd_omega')
    c_r = 0;
else
    c_r = -trackedPeakU1-c_l*trackedPeakX;
end

rhoX = flow.source;
wallOmegaPeak = details.trackedWallPeak;
wallOmegaPeakIndex = details.wallOmegaPeakIndex;
[~,baseRhs,transportU1,transportU2] = ipm.evolve.assembleRhs( ...
    rho,flow.u1,flow.u2,c_l,c_l,0,c_r,ops,'open');
if strcmp(r.lengthGauge,'transport_anchor')
    transportAnchorX = r.transportAnchorX;
    transportAnchorVelocity = interp1( ...
        ops.x,transportU1(1,:),transportAnchorX,'linear');
else
    transportAnchorX = NaN;
    transportAnchorVelocity = NaN;
end

if needsStrain
    strain = u1X(1,r.originIndex);
    strainScale = max(abs(u1X),[],'all');
    strainCondition = abs(strain)/max(strainScale,realmin);
else
    strain = NaN;
    strainScale = NaN;
    strainCondition = NaN;
end
omegaGaugeAnchorX = NaN;
omegaGaugeAnchorValue = NaN;
omegaGaugeAnchorCondition = NaN;
omegaGaugeAnchorRateResidual = NaN;
omegaGaugeReferenceAnchorValue = NaN;
omegaGaugeWindowCenter = NaN;
omegaGaugeWindowRadius = NaN;
omegaGaugeWindowValue = NaN;
omegaGaugeWindowCondition = NaN;
omegaGaugeWindowEnergy = NaN;
omegaGaugeWindowRateResidual = NaN;
omegaGaugeReferenceWindowValue = NaN;
omegaGaugeWindowMomentOrder = NaN;
omegaGaugeWindowMoment = NaN;
omegaGaugeWindowSupportPoints = NaN;
omegaGaugeWindowRateScale = NaN;
omegaGaugeTemplateCenter = NaN;
omegaGaugeTemplateRadius = NaN;
omegaGaugeTemplateSupportPoints = NaN;
omegaGaugeTemplateValue = NaN;
omegaGaugeTemplateCondition = NaN;
omegaGaugeTemplateRateResidual = NaN;
omegaGaugeTemplateRateScale = NaN;
omegaGaugeReferenceTemplateValue = NaN;
omegaGaugeBulkCenterX = NaN;
omegaGaugeBulkCenterY = NaN;
omegaGaugeBulkRadius = NaN;
omegaGaugeBulkSupportPointsX = NaN;
omegaGaugeBulkSupportPointsY = NaN;
omegaGaugeBulkValue = NaN;
omegaGaugeBulkCondition = NaN;
omegaGaugeBulkEnergy = NaN;
omegaGaugeBulkRateResidual = NaN;
omegaGaugeBulkRateScale = NaN;
omegaGaugeReferenceBulkValue = NaN;
omegaGaugeEnergyX = NaN;
omegaGaugeEnergyY = NaN;
omegaGaugeForcingX = NaN;
omegaGaugeForcingY = NaN;
omegaGaugeForcing = NaN;
omegaGaugeAbsoluteForcing = NaN;
omegaGaugeForcingCancellationRatio = NaN;
omegaGaugeQuadraticPeakX = NaN;
omegaGaugeQuadraticPeakValue = NaN;
omegaGaugeQuadraticForcing = NaN;
omegaGaugeQuadraticCurvature = NaN;
omegaGaugeQuadraticStationarityResidual = NaN;
omegaGaugeQuadraticNormalizedStationarityResidual = NaN;
omegaGaugeQuadraticCurvatureCondition = NaN;
omegaGaugeQuadraticWeightCondition = NaN;
omegaGaugeQuadraticActiveNodeMargin = NaN;
omegaGaugeQuadraticActiveSelectionMargin = NaN;
omegaGaugeQuadraticStencilCenterX = NaN;
omegaGaugeQuadraticVertexInteriorFraction = NaN;
omegaGaugeQuadraticRelativeSubgridLift = NaN;
omegaGaugeQuadraticRateResidual = NaN;
omegaGaugeQuadraticRateScale = NaN;
switch r.cOmegaGauge
    case 'wall_omega_peak'
        baseRhsX = baseRhs*ops.Dx';
        if wallOmegaPeak < r.cOmegaGaugeFloor
            error('ipm:DegenerateNormalization', ...
                ['The wall-omega peak c_omega gauge denominator became ' ...
                'too small.']);
        end
        gaugeForcing = baseRhsX(1,wallOmegaPeakIndex);
        gaugeEnergy = wallOmegaPeak;
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaPeakRateResidual = gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeResidual = omegaPeakRateResidual / ...
            max([abs(gaugeForcing),abs(c_omega*gaugeEnergy),realmin]);
        strainRateResidual = NaN;
    case 'wall_omega_quadratic_peak'
        peak = ipm.evolve.quadraticPeakFunctional( ...
            rhoX(1,:),ops.x,wallOmegaPeakIndex);
        if any(abs(ops.x(peak.indices)-r.pinX) > ...
                r.peakTrackingHalfWidth)
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega peak stencil crosses the ' ...
                'registered tracking-window boundary.']);
        end
        sourceScale = max(abs(rhoX),[],'all');
        peakCondition = peak.value/max(sourceScale,realmin);
        if peak.value < r.cOmegaGaugeFloor || ...
                peakCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega peak became too small ' ...
                'relative to the full R_X field.']);
        end
        stencilForcing = baseRhs(1,:)*ops.Dx(peak.indices,:)';
        gaugeForcing = sum(peak.weights.*stencilForcing);
        gaugeEnergy = peak.value;
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The quadratic wall-omega peak forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeQuadraticRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeQuadraticRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),realmin]);
        omegaGaugeResidual = omegaGaugeQuadraticRateResidual / ...
            omegaGaugeQuadraticRateScale;
        omegaPeakRateResidual = omegaGaugeQuadraticRateResidual;
        omegaGaugeQuadraticPeakX = peak.x;
        omegaGaugeQuadraticPeakValue = peak.value;
        omegaGaugeQuadraticForcing = gaugeForcing;
        omegaGaugeQuadraticCurvature = peak.curvature;
        omegaGaugeQuadraticStationarityResidual = ...
            peak.stationarityResidual;
        omegaGaugeQuadraticNormalizedStationarityResidual = ...
            peak.normalizedStationarityResidual;
        omegaGaugeQuadraticCurvatureCondition = ...
            peak.curvatureCondition;
        omegaGaugeQuadraticWeightCondition = peak.weightCondition;
        omegaGaugeQuadraticActiveNodeMargin = peak.activeNodeMargin;
        trackingMask = abs(ops.x-r.pinX) <= r.peakTrackingHalfWidth;
        competitorMask = trackingMask;
        competitorMask(wallOmegaPeakIndex) = false;
        selectionCompetitors = rhoX(1,competitorMask);
        if isempty(selectionCompetitors)
            error('ipm:DegenerateNormalization', ...
                'The quadratic wall-omega tracking window has no peer node.');
        end
        omegaGaugeQuadraticActiveSelectionMargin = ( ...
            rhoX(1,wallOmegaPeakIndex)-max(selectionCompetitors)) / ...
            peak.value;
        selectionRoundoffTolerance = 512*eps(1)*max( ...
            max(abs(rhoX(1,trackingMask))),realmin)/peak.value;
        if omegaGaugeQuadraticActiveSelectionMargin < ...
                -selectionRoundoffTolerance
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega stencil center is not the ' ...
                'active maximum in the registered tracking window.']);
        end
        omegaGaugeQuadraticStencilCenterX = ...
            ops.x(wallOmegaPeakIndex);
        omegaGaugeQuadraticVertexInteriorFraction = ...
            peak.vertexInteriorFraction;
        omegaGaugeQuadraticRelativeSubgridLift = ...
            peak.relativeSubgridLift;
        details.omegaGaugePeakX = peak.x;
        details.omegaGaugePeakCondition = peakCondition;
        strainRateResidual = NaN;
    case 'gradient_energy'
        integrationWeights = ops.integrationWeights;
        rhoY = ops.Dy*rho;
        baseRhsX = baseRhs*ops.Dx';
        baseRhsY = ops.Dy*baseRhs;
        omegaGaugeEnergyX = sum( ...
            rhoX.^2.*integrationWeights,'all');
        omegaGaugeEnergyY = sum( ...
            rhoY.^2.*integrationWeights,'all');
        omegaGaugeForcingX = sum( ...
            rhoX.*baseRhsX.*integrationWeights,'all');
        omegaGaugeForcingY = sum( ...
            rhoY.*baseRhsY.*integrationWeights,'all');
        gaugeEnergy = sum((rhoX.^2+rhoY.^2).*integrationWeights,'all');
        gaugeForcing = sum((rhoX.*baseRhsX+rhoY.*baseRhsY) .* ...
            integrationWeights,'all');
        omegaGaugeAbsoluteForcing = sum(abs( ...
            rhoX.*baseRhsX+rhoY.*baseRhsY).*integrationWeights,'all');
        finiteGaugeValues = [omegaGaugeEnergyX,omegaGaugeEnergyY, ...
            gaugeEnergy,omegaGaugeForcingX,omegaGaugeForcingY, ...
            gaugeForcing,omegaGaugeAbsoluteForcing];
        if any(~isfinite(finiteGaugeValues))
            error('ipm:DegenerateNormalization', ...
                ['The full-gradient c_omega gauge energy or forcing ' ...
                'diagnostic is not finite.']);
        end
        if gaugeEnergy < r.cOmegaGaugeFloor
            error('ipm:DegenerateNormalization', ...
                ['The full-gradient c_omega gauge energy became too ' ...
                'small.']);
        end
        omegaGaugeForcing = gaugeForcing;
        omegaGaugeForcingCancellationRatio = abs(gaugeForcing) / ...
            max(omegaGaugeAbsoluteForcing,realmin);
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeResidual = (gaugeForcing+c_omega*gaugeEnergy) / ...
            max([abs(gaugeForcing),abs(c_omega*gaugeEnergy),realmin]);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'strain_point'
        if strainCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The point-strain c_omega gauge is singular: ' ...
                'd_x u_1(0,0) is too small relative to the strain field.']);
        end
        [baseU1,~] = ipm.field.velocity(baseRhs,ops);
        baseStrain = baseU1*ipm.mesh.oddDx(ops)';
        gaugeForcing = baseStrain(1,r.originIndex);
        gaugeEnergy = strain;
        c_omega = -gaugeForcing/strain;
        strainRateResidual = gaugeForcing+c_omega*strain;
        omegaGaugeResidual = strainRateResidual / ...
            max([abs(gaugeForcing),abs(c_omega*strain),realmin]);
        omegaPeakRateResidual = NaN;
    case {'anchor_wall_window_l2','anchor_wall_window_l4'}
        baseRhsX = baseRhs*ops.Dx';
        if strcmp(r.cOmegaGauge,'anchor_wall_window_l4')
            momentOrder = 4;
        else
            momentOrder = 2;
        end
        window = ipm.evolve.wallWindowMoment( ...
            rhoX(1,:),ops,r,momentOrder);
        if ~isfield(r,'referenceAnchorWallWindowValue') || ...
                ~isnumeric(r.referenceAnchorWallWindowValue) || ...
                ~isreal(r.referenceAnchorWallWindowValue) || ...
                ~isscalar(r.referenceAnchorWallWindowValue) || ...
                ~isfinite(r.referenceAnchorWallWindowValue) || ...
                r.referenceAnchorWallWindowValue <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall window gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        omegaGaugeWindowCenter = window.center;
        omegaGaugeWindowRadius = window.radius;
        omegaGaugeWindowValue = window.value;
        omegaGaugeWindowCondition = window.condition;
        omegaGaugeWindowEnergy = window.energy;
        omegaGaugeWindowMomentOrder = window.order;
        omegaGaugeWindowMoment = window.moment;
        omegaGaugeWindowSupportPoints = window.supportPoints;
        omegaGaugeReferenceWindowValue = ...
            r.referenceAnchorWallWindowValue;
        gaugeEnergy = window.moment;
        gaugeForcing = sum(window.normalizedWeights .* ...
            rhoX(1,:).^(momentOrder-1).*baseRhsX(1,:));
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall moment forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeWindowRateResidual = gaugeForcing + ...
            c_omega*gaugeEnergy;
        omegaGaugeWindowRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy), ...
            strainScale*gaugeEnergy]);
        omegaGaugeResidual = omegaGaugeWindowRateResidual / ...
            max(omegaGaugeWindowRateScale,realmin);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'outer_wall_density_window_l2'
        window = ipm.evolve.outerWallDensityWindow( ...
            rho(1,:),ops,r);
        if ~isfield(r,'referenceOuterWallDensityWindowValue') || ...
                ~isnumeric(r.referenceOuterWallDensityWindowValue) || ...
                ~isreal(r.referenceOuterWallDensityWindowValue) || ...
                ~isscalar(r.referenceOuterWallDensityWindowValue) || ...
                ~isfinite(r.referenceOuterWallDensityWindowValue) || ...
                r.referenceOuterWallDensityWindowValue <= 0
            error('ipm:DegenerateNormalization', ...
                ['The outer wall-density gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        gaugeEnergy = window.moment;
        gaugeForcing = sum(window.normalizedWeights .* ...
            rho(1,:).*baseRhs(1,:));
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The outer wall-density forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeWindowRateResidual = gaugeForcing + ...
            c_omega*gaugeEnergy;
        omegaGaugeWindowRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),realmin]);
        omegaGaugeResidual = omegaGaugeWindowRateResidual / ...
            omegaGaugeWindowRateScale;
        omegaGaugeWindowCenter = window.center;
        omegaGaugeWindowRadius = window.radius;
        omegaGaugeWindowValue = window.value;
        omegaGaugeWindowCondition = window.condition;
        omegaGaugeWindowEnergy = window.energy;
        omegaGaugeWindowMomentOrder = window.order;
        omegaGaugeWindowMoment = window.moment;
        omegaGaugeWindowSupportPoints = window.supportPoints;
        omegaGaugeReferenceWindowValue = ...
            r.referenceOuterWallDensityWindowValue;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_wall_template_projection'
        baseRhsX = baseRhs*ops.Dx';
        geometry = ipm.evolve.wallWindowGeometry(ops,r);
        if ~isfield(r,'referenceAnchorWallTemplate') || ...
                ~isnumeric(r.referenceAnchorWallTemplate) || ...
                ~isreal(r.referenceAnchorWallTemplate) || ...
                ~isequal(size(r.referenceAnchorWallTemplate),size(ops.x)) || ...
                any(~isfinite(r.referenceAnchorWallTemplate),'all') || ...
                ~isfield(r,'referenceAnchorWallTemplateProjection') || ...
                ~isnumeric(r.referenceAnchorWallTemplateProjection) || ...
                ~isreal(r.referenceAnchorWallTemplateProjection) || ...
                ~isscalar(r.referenceAnchorWallTemplateProjection) || ...
                ~isfinite(r.referenceAnchorWallTemplateProjection) || ...
                r.referenceAnchorWallTemplateProjection <= 0
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall template gauge lacks a valid frozen template.');
        end
        template = r.referenceAnchorWallTemplate;
        weights = geometry.normalizedWeights;
        templateNorm = sqrt(sum(weights.*template.^2));
        currentNorm = sqrt(sum(weights.*rhoX(1,:).^2));
        projection = sum(weights.*rhoX(1,:).*template);
        condition = abs(projection)/max(templateNorm*currentNorm,realmin);
        if ~all(isfinite([templateNorm,currentNorm,projection,condition])) || ...
                abs(projection) < r.cOmegaGaugeFloor || ...
                projection*r.referenceAnchorWallTemplateProjection <= 0 || ...
                condition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall template projection became singular or ' ...
                'lost the sign of its frozen reference.']);
        end
        gaugeEnergy = projection;
        gaugeForcing = sum(weights.*baseRhsX(1,:).*template);
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall template forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeTemplateRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeTemplateRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),strainScale*abs(gaugeEnergy)]);
        omegaGaugeResidual = omegaGaugeTemplateRateResidual / ...
            max(omegaGaugeTemplateRateScale,realmin);
        omegaGaugeTemplateCenter = geometry.center;
        omegaGaugeTemplateRadius = geometry.radius;
        omegaGaugeTemplateSupportPoints = geometry.supportPoints;
        omegaGaugeTemplateValue = projection;
        omegaGaugeTemplateCondition = condition;
        omegaGaugeReferenceTemplateValue = ...
            r.referenceAnchorWallTemplateProjection;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_bulk_gradient_l2'
        rhoY = ops.Dy*rho;
        baseRhsX = baseRhs*ops.Dx';
        baseRhsY = ops.Dy*baseRhs;
        bulk = ipm.evolve.bulkGradientWindow(rhoX,rhoY,ops,r);
        if ~isfield(r,'referenceAnchorBulkGradientL2') || ...
                ~isnumeric(r.referenceAnchorBulkGradientL2) || ...
                ~isreal(r.referenceAnchorBulkGradientL2) || ...
                ~isscalar(r.referenceAnchorBulkGradientL2) || ...
                ~isfinite(r.referenceAnchorBulkGradientL2) || ...
                r.referenceAnchorBulkGradientL2 <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-bulk gradient L2 gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        gaugeEnergy = bulk.energy;
        gaugeForcing = sum(bulk.normalizedWeights .* ...
            (rhoX.*baseRhsX+rhoY.*baseRhsY),'all');
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-bulk gradient forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeBulkRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeBulkRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),strainScale*gaugeEnergy]);
        omegaGaugeResidual = omegaGaugeBulkRateResidual / ...
            max(omegaGaugeBulkRateScale,realmin);
        omegaGaugeBulkCenterX = bulk.centerX;
        omegaGaugeBulkCenterY = bulk.centerY;
        omegaGaugeBulkRadius = bulk.radius;
        omegaGaugeBulkSupportPointsX = bulk.supportPointsX;
        omegaGaugeBulkSupportPointsY = bulk.supportPointsY;
        omegaGaugeBulkValue = bulk.value;
        omegaGaugeBulkCondition = bulk.condition;
        omegaGaugeBulkEnergy = bulk.energy;
        omegaGaugeReferenceBulkValue = ...
            r.referenceAnchorBulkGradientL2;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_wall_slope'
        baseRhsX = baseRhs*ops.Dx';
        omegaGaugeAnchorX = r.transportAnchorX;
        % Evaluate every member of the anchor constraint with one common
        % fixed-coordinate linear interpolation functional.
        anchorSamples = interp1(ops.x,[rhoX(1,:);baseRhsX(1,:)]', ...
            omegaGaugeAnchorX,'linear');
        omegaGaugeAnchorValue = anchorSamples(1);
        gaugeForcing = anchorSamples(2);
        wallSlopeScale = max(abs(rhoX(1,:)));
        omegaGaugeAnchorCondition = abs(omegaGaugeAnchorValue) / ...
            max(wallSlopeScale,realmin);
        if abs(omegaGaugeAnchorValue) < r.cOmegaGaugeFloor || ...
                omegaGaugeAnchorCondition < ...
                r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall c_omega gauge slope became too small ' ...
                'relative to the wall R_X profile.']);
        end
        if ~isfield(r,'referenceAnchorWallSlope') || ...
                omegaGaugeAnchorValue*r.referenceAnchorWallSlope <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall c_omega gauge lost the sign of its ' ...
                'frozen reference slope.']);
        end
        omegaGaugeReferenceAnchorValue = r.referenceAnchorWallSlope;
        gaugeEnergy = omegaGaugeAnchorValue;
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeAnchorRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        anchorSlopeRateScale = max([ ...
            abs(gaugeForcing),abs(c_omega*gaugeEnergy), ...
            strainScale*abs(gaugeEnergy)]);
        omegaGaugeResidual = omegaGaugeAnchorRateResidual / ...
            max(anchorSlopeRateScale,realmin);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'none'
        c_omega = 0;
        gaugeEnergy = sum(rho.^2.*ops.integrationWeights,'all');
        omegaGaugeResidual = 0;
        % C_omega=1 is prescribed by this identity/material-amplitude
        % gauge, so its selected zero rate is not empirical convergence
        % evidence.  Independently evaluate the local quadratic wall-peak
        % functional.  Here F/P is the actual peak log-rate and -F/P is
        % the counterfactual c_omega that would freeze the peak.  Neither
        % quantity feeds back into the selected zero rate.
        omegaPeakRateResidual = NaN;
        try
        peak = ipm.evolve.quadraticPeakFunctional( ...
            rhoX(1,:),ops.x,wallOmegaPeakIndex);
        if any(abs(ops.x(peak.indices)-r.pinX) > ...
                r.peakTrackingHalfWidth)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic stencil crosses ' ...
                'the registered tracking-window boundary.']);
        end
        sourceScale = max(abs(rhoX),[],'all');
        peakCondition = peak.value/max(sourceScale,realmin);
        if peak.value < r.cOmegaGaugeFloor || ...
                peakCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic peak became too ' ...
                'small relative to the full R_X field.']);
        end
        stencilForcing = baseRhs(1,:)*ops.Dx(peak.indices,:)';
        gaugeForcing = sum(peak.weights.*stencilForcing);
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic forcing is not ' ...
                'finite.']);
        end
        omegaGaugeQuadraticRateResidual = gaugeForcing;
        omegaGaugeQuadraticRateScale = max(abs(gaugeForcing),realmin);
        omegaGaugeQuadraticPeakX = peak.x;
        omegaGaugeQuadraticPeakValue = peak.value;
        omegaGaugeQuadraticForcing = gaugeForcing;
        omegaGaugeQuadraticCurvature = peak.curvature;
        omegaGaugeQuadraticStationarityResidual = ...
            peak.stationarityResidual;
        omegaGaugeQuadraticNormalizedStationarityResidual = ...
            peak.normalizedStationarityResidual;
        omegaGaugeQuadraticCurvatureCondition = peak.curvatureCondition;
        omegaGaugeQuadraticWeightCondition = peak.weightCondition;
        omegaGaugeQuadraticActiveNodeMargin = peak.activeNodeMargin;
        trackingMask = abs(ops.x-r.pinX) <= r.peakTrackingHalfWidth;
        competitorMask = trackingMask;
        competitorMask(wallOmegaPeakIndex) = false;
        selectionCompetitors = rhoX(1,competitorMask);
        if isempty(selectionCompetitors)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic window has no ' ...
                'peer node.']);
        end
        omegaGaugeQuadraticActiveSelectionMargin = ( ...
            rhoX(1,wallOmegaPeakIndex)-max(selectionCompetitors)) / ...
            peak.value;
        selectionRoundoffTolerance = 512*eps(1)*max( ...
            max(abs(rhoX(1,trackingMask))),realmin)/peak.value;
        if omegaGaugeQuadraticActiveSelectionMargin < ...
                -selectionRoundoffTolerance
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic center is not ' ...
                'the active maximum in the registered window.']);
        end
        omegaGaugeQuadraticStencilCenterX = ...
            ops.x(wallOmegaPeakIndex);
        omegaGaugeQuadraticVertexInteriorFraction = ...
            peak.vertexInteriorFraction;
        omegaGaugeQuadraticRelativeSubgridLift = peak.relativeSubgridLift;
        details.omegaGaugePeakX = peak.x;
        details.omegaGaugePeakCondition = peakCondition;
        omegaPeakRateResidual = gaugeForcing;
        catch exception
            % A diagnostic peak may disappear in a generic identity-gauge
            % calculation.  That is recorded by absent/NaN telemetry and
            % must not turn an otherwise valid prescribed gauge into a
            % hidden normalization failure.
            if ~strcmp(exception.identifier,'ipm:DegenerateNormalization')
                rethrow(exception)
            end
        end
        strainRateResidual = NaN;
    case 'anchor_wall_strain'
        error('ipm:RetiredCOmegaGauge', ...
            ['anchor_wall_strain is not an exact discrete normalization ' ...
            'and is retired from the current runtime contract.']);
    otherwise
        error('ipm:BadCOmegaGaugeRuntime', ...
            'Unknown runtime c_omega gauge: %s.',r.cOmegaGauge);
end

if strcmp(r.lengthGauge,'omega_peak_location')
    fullLengthContribution = c_l*( ...
        2*peakLocationRhoXX+peakLocationResponse);
    fullTranslationContribution = c_r*peakRhoXXX;
    fullAmplitudeContribution = c_omega*peakLocationRhoXX;
    peakLocationFullPhaseResidual = peakPhysicalRhsXX - ...
        fullLengthContribution-fullTranslationContribution + ...
        fullAmplitudeContribution;
    fullPhaseScale = max([abs(peakPhysicalRhsXX), ...
        abs(fullLengthContribution),abs(fullTranslationContribution), ...
        abs(fullAmplitudeContribution),realmin]);
    peakLocationNormalizedFullPhaseResidual = ...
        peakLocationFullPhaseResidual/fullPhaseScale;
end

details.strain = strain;
details.strainCondition = strainCondition;
details.c_lNominal = c_lNominal;
details.lengthScaleGain = 1;
details.widthCorrection = widthCorrection;
details.widthRateScale = abs(c_lNominal);
details.widthControlMode = widthControlMode;
details.peakExpansionRatio = peakExpansionRatio;
details.resolutionFactor = resolutionFactor;
details.safetyFactor = safetyFactor;
details.strainError = strain-r.strainTarget;
details.strainRateResidual = strainRateResidual;
details.omegaGaugeEnergy = gaugeEnergy;
details.omegaGaugeResidual = omegaGaugeResidual;
details.omegaGaugeEnergyX = omegaGaugeEnergyX;
details.omegaGaugeEnergyY = omegaGaugeEnergyY;
details.omegaGaugeForcingX = omegaGaugeForcingX;
details.omegaGaugeForcingY = omegaGaugeForcingY;
details.omegaGaugeForcing = omegaGaugeForcing;
details.omegaGaugeAbsoluteForcing = omegaGaugeAbsoluteForcing;
details.omegaGaugeForcingCancellationRatio = ...
    omegaGaugeForcingCancellationRatio;
details.omegaGaugeQuadraticPeakX = omegaGaugeQuadraticPeakX;
details.omegaGaugeQuadraticPeakValue = omegaGaugeQuadraticPeakValue;
details.omegaGaugeQuadraticForcing = omegaGaugeQuadraticForcing;
details.omegaGaugeQuadraticCurvature = omegaGaugeQuadraticCurvature;
details.omegaGaugeQuadraticStationarityResidual = ...
    omegaGaugeQuadraticStationarityResidual;
details.omegaGaugeQuadraticNormalizedStationarityResidual = ...
    omegaGaugeQuadraticNormalizedStationarityResidual;
details.omegaGaugeQuadraticCurvatureCondition = ...
    omegaGaugeQuadraticCurvatureCondition;
details.omegaGaugeQuadraticWeightCondition = ...
    omegaGaugeQuadraticWeightCondition;
details.omegaGaugeQuadraticActiveNodeMargin = ...
    omegaGaugeQuadraticActiveNodeMargin;
details.omegaGaugeQuadraticActiveSelectionMargin = ...
    omegaGaugeQuadraticActiveSelectionMargin;
details.omegaGaugeQuadraticStencilCenterX = ...
    omegaGaugeQuadraticStencilCenterX;
details.omegaGaugeQuadraticVertexInteriorFraction = ...
    omegaGaugeQuadraticVertexInteriorFraction;
details.omegaGaugeQuadraticRelativeSubgridLift = ...
    omegaGaugeQuadraticRelativeSubgridLift;
details.omegaGaugeQuadraticRateResidual = ...
    omegaGaugeQuadraticRateResidual;
details.omegaGaugeQuadraticRateScale = omegaGaugeQuadraticRateScale;
details.omegaPeakRateResidual = omegaPeakRateResidual;
details.omegaGaugeAnchorX = omegaGaugeAnchorX;
details.omegaGaugeAnchorValue = omegaGaugeAnchorValue;
details.omegaGaugeAnchorCondition = omegaGaugeAnchorCondition;
details.omegaGaugeAnchorRateResidual = omegaGaugeAnchorRateResidual;
details.omegaGaugeReferenceAnchorValue = ...
    omegaGaugeReferenceAnchorValue;
details.omegaGaugeWindowCenter = omegaGaugeWindowCenter;
details.omegaGaugeWindowRadius = omegaGaugeWindowRadius;
details.omegaGaugeWindowValue = omegaGaugeWindowValue;
details.omegaGaugeWindowCondition = omegaGaugeWindowCondition;
details.omegaGaugeWindowEnergy = omegaGaugeWindowEnergy;
details.omegaGaugeWindowRateResidual = omegaGaugeWindowRateResidual;
details.omegaGaugeReferenceWindowValue = ...
    omegaGaugeReferenceWindowValue;
details.omegaGaugeWindowMomentOrder = omegaGaugeWindowMomentOrder;
details.omegaGaugeWindowMoment = omegaGaugeWindowMoment;
details.omegaGaugeWindowSupportPoints = omegaGaugeWindowSupportPoints;
details.omegaGaugeWindowRateScale = omegaGaugeWindowRateScale;
details.omegaGaugeTemplateCenter = omegaGaugeTemplateCenter;
details.omegaGaugeTemplateRadius = omegaGaugeTemplateRadius;
details.omegaGaugeTemplateSupportPoints = ...
    omegaGaugeTemplateSupportPoints;
details.omegaGaugeTemplateValue = omegaGaugeTemplateValue;
details.omegaGaugeTemplateCondition = omegaGaugeTemplateCondition;
details.omegaGaugeTemplateRateResidual = ...
    omegaGaugeTemplateRateResidual;
details.omegaGaugeTemplateRateScale = omegaGaugeTemplateRateScale;
details.omegaGaugeReferenceTemplateValue = ...
    omegaGaugeReferenceTemplateValue;
details.omegaGaugeBulkCenterX = omegaGaugeBulkCenterX;
details.omegaGaugeBulkCenterY = omegaGaugeBulkCenterY;
details.omegaGaugeBulkRadius = omegaGaugeBulkRadius;
details.omegaGaugeBulkSupportPointsX = omegaGaugeBulkSupportPointsX;
details.omegaGaugeBulkSupportPointsY = omegaGaugeBulkSupportPointsY;
details.omegaGaugeBulkValue = omegaGaugeBulkValue;
details.omegaGaugeBulkCondition = omegaGaugeBulkCondition;
details.omegaGaugeBulkEnergy = omegaGaugeBulkEnergy;
details.omegaGaugeBulkRateResidual = omegaGaugeBulkRateResidual;
details.omegaGaugeBulkRateScale = omegaGaugeBulkRateScale;
details.omegaGaugeReferenceBulkValue = omegaGaugeReferenceBulkValue;
details.lengthGaugeStrain = lengthGaugeStrain;
details.lengthGaugeResidual = lengthGaugeResidual;
details.peakLocationRhoXX = peakLocationRhoXX;
details.peakLocationRhoXXScale = peakLocationRhoXXScale;
details.peakLocationPhaseDefect = peakLocationPhaseDefect;
details.peakLocationResponseCondition = peakLocationResponseCondition;
details.peakLocationFullPhaseResidual = ...
    peakLocationFullPhaseResidual;
details.peakLocationNormalizedFullPhaseResidual = ...
    peakLocationNormalizedFullPhaseResidual;
details.pinGaugeResidual = details.lengthGaugeResidual;
details.transportAnchorX = transportAnchorX;
details.transportAnchorVelocity = transportAnchorVelocity;
if strcmp(ops.symmetryMode,'double_odd_omega')
    details.pinTransportVelocity = transportU1(1,r.originIndex);
else
    details.pinTransportVelocity = ...
        interp1(ops.x,transportU1(1,:),r.pinX,'linear');
end
details.peakTransportVelocity = interp1( ...
    ops.x,transportU1(1,:),trackedPeakX,'linear');
if strcmp(ops.symmetryMode,'double_odd_omega')
    details.travelingWaveResidual = transportU1(1,r.originIndex);
else
    details.travelingWaveResidual = details.peakTransportVelocity;
end
end
