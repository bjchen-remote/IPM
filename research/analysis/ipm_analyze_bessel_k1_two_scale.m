function diagnostics = ipm_analyze_bessel_k1_two_scale(inputResult,userOpts)
%IPM_ANALYZE_BESSEL_K1_TWO_SCALE Analyze a regularized K_1*sin(theta) run.
%   The parabolic far-tail model is
%
%     rho = C*y*x^(-3/2)*exp(-mu*y^2/(4*x)),
%
%   whose transverse maximum and rho_x zero curve are respectively
%
%     y_peak=sqrt(2*x/mu),   y_zero=sqrt(6*x/mu).

if nargin < 2 || isempty(userOpts)
    userOpts = struct();
end
canonicalResult = ipm.output.validate(inputResult);
result = canonicalResult.grid;
caseMetadata = ipm_bessel_case_metadata(canonicalResult, ...
    {'mu','centerX','perturbationAmplitude'});
parameters = caseMetadata.parameters;
opts = analysis_options(result,userOpts,parameters);
if ~strcmp(canonicalResult.config.scaling.rescalingMode,'physical')
    error('ipm:BesselK1PhysicalOnly', ...
        'This diagnostic currently expects physical-mode snapshots.');
end
if ~strcmpi(string(caseMetadata.family),'k1_regularized')
    error('ipm:BesselK1Family', ...
        'Use this diagnostic with a k1_regularized experiment result.');
end

snapshots = ipm.output.snapshotAt(canonicalResult,'all');
fields = {snapshots.physicalRho}';
times = [snapshots.physicalTime]';
template = struct('time',NaN,'mu',NaN,'profileError',NaN, ...
    'fullKernelError',NaN,'peakAmplitudeVariation',NaN, ...
    'peakWidthExponent',NaN,'rhoXZeroRelativeError',NaN, ...
    'rhoXZeroCoverage',0,'wallLeakage',NaN,'tailX',[], ...
    'observedPeakY',[],'predictedPeakY',[],'observedZeroY',[], ...
    'predictedZeroY',[],'rhoX',[]);
fits = repmat(template,numel(fields),1);
for index = 1:numel(fields)
    fits(index) = analyze_snapshot(fields{index},times(index), ...
        result.x,result.y,opts);
end

diagnostics = struct('times',times(:),'mu',[fits.mu]', ...
    'profileError',[fits.profileError]', ...
    'fullKernelError',[fits.fullKernelError]', ...
    'peakAmplitudeVariation',[fits.peakAmplitudeVariation]', ...
    'peakWidthExponent',[fits.peakWidthExponent]', ...
    'rhoXZeroRelativeError',[fits.rhoXZeroRelativeError]', ...
    'rhoXZeroCoverage',[fits.rhoXZeroCoverage]', ...
    'wallLeakage',[fits.wallLeakage]','fits',fits,'options',opts, ...
    'caseMetadata',caseMetadata, ...
    'initialParameters',struct('mu',parameters.mu, ...
        'perturbationAmplitude',parameters.perturbationAmplitude));

diagnostics.fullRelativeChangeFromInitial = zeros(numel(fields),1);
diagnostics.tailRelativeChangeFromInitial = zeros(numel(fields),1);
centeredX = result.x-opts.centerX;
tailMaskX = centeredX >= opts.tailXRange(1) & ...
    centeredX <= opts.tailXRange(2);
tailMaskY = result.y <= opts.tailYMax;
for index = 1:numel(fields)
    diagnostics.fullRelativeChangeFromInitial(index) = ...
        weighted_relative_change(fields{index},fields{1},result.x,result.y);
    diagnostics.tailRelativeChangeFromInitial(index) = ...
        weighted_relative_change(fields{index}(tailMaskY,tailMaskX), ...
        fields{1}(tailMaskY,tailMaskX),result.x(tailMaskX), ...
        result.y(tailMaskY));
end
gradient = canonicalResult.history.common.physicalGradInf;
diagnostics.gradientInfGrowth = gradient(end)/max(gradient(1),eps);

if opts.makePlot
    make_diagnostic_figure(diagnostics,fields,result,opts);
end
end

function fit = analyze_snapshot(rho,time,x,y,opts)
centeredX = x-opts.centerX;
tailColumns = find(centeredX >= opts.tailXRange(1) & ...
    centeredX <= opts.tailXRange(2));
tailX = centeredX(tailColumns);
tailRho = rho(:,tailColumns);
observedPeakY = nan(size(tailX));
peakValues = nan(size(tailX));
for column = 1:numel(tailX)
    [peakValues(column),peakIndex] = max(tailRho(:,column));
    observedPeakY(column) = ipm.diagnostics.peakLocation( ...
        tailRho(:,column)',y',peakIndex);
end
muSamples = 2*tailX./observedPeakY.^2;
muMask = isfinite(muSamples) & muSamples > 0;
mu = median(muSamples(muMask));
predictedPeakY = sqrt(2*tailX/mu);

[xGrid,yGrid] = meshgrid(tailX,y);
eta = yGrid.*sqrt(mu./xGrid);
predictedRatio = (eta/sqrt(2)).*exp((2-eta.^2)/4);
ratio = tailRho./peakValues;
comparisonMask = yGrid <= opts.tailYMax & ratio >= 0 & ...
    ratio <= opts.maximumComparisonRatio & isfinite(ratio);
profileError = relative_error(ratio(comparisonMask), ...
    predictedRatio(comparisonMask));

kernelShape = yGrid.*xGrid.^(-3/2).*exp(-mu*yGrid.^2./(4*xGrid));
fullMask = comparisonMask & isfinite(kernelShape);
coefficient = dot(kernelShape(fullMask),tailRho(fullMask)) / ...
    max(dot(kernelShape(fullMask),kernelShape(fullMask)),eps);
fullKernelError = relative_error(tailRho(fullMask), ...
    coefficient*kernelShape(fullMask));
peakAmplitude = tailX.*peakValues;
peakAmplitudeVariation = std(peakAmplitude)/max(abs(mean(peakAmplitude)),eps);

widthMask = isfinite(observedPeakY) & observedPeakY > 0;
if nnz(widthMask) >= 3
    exponentFit = polyfit(log(tailX(widthMask)), ...
        log(observedPeakY(widthMask)),1);
    peakWidthExponent = exponentFit(1);
else
    peakWidthExponent = NaN;
end

rhoX = differentiate_x(rho,x);
observedZeroY = nan(size(tailX));
predictedZeroY = sqrt(6*tailX/mu);
for column = 1:numel(tailX)
    observedZeroY(column) = first_negative_to_positive_crossing( ...
        y,rhoX(:,tailColumns(column)));
end
zeroMask = isfinite(observedZeroY) & predictedZeroY <= opts.tailYMax;
rhoXZeroRelativeError = relative_error(observedZeroY(zeroMask), ...
    predictedZeroY(zeroMask));
rhoXZeroCoverage = nnz(zeroMask)/numel(tailX);
wallLeakage = max(abs(rho(1,:)))/max(abs(rho),[],'all');

fit = struct('time',time,'mu',mu,'profileError',profileError, ...
    'fullKernelError',fullKernelError, ...
    'peakAmplitudeVariation',peakAmplitudeVariation, ...
    'peakWidthExponent',peakWidthExponent, ...
    'rhoXZeroRelativeError',rhoXZeroRelativeError, ...
    'rhoXZeroCoverage',rhoXZeroCoverage,'wallLeakage',wallLeakage, ...
    'tailX',tailX,'observedPeakY',observedPeakY, ...
    'predictedPeakY',predictedPeakY,'observedZeroY',observedZeroY, ...
    'predictedZeroY',predictedZeroY,'rhoX',rhoX);
end

function crossing = first_negative_to_positive_crossing(y,trace)
crossing = NaN;
index = find(trace(1:end-1) <= 0 & trace(2:end) > 0,1,'first');
if isempty(index)
    return;
end
weight = -trace(index)/(trace(index+1)-trace(index));
crossing = y(index)+weight*(y(index+1)-y(index));
end

function derivative = differentiate_x(field,x)
derivative = zeros(size(field));
for row = 1:size(field,1)
    derivative(row,:) = gradient(field(row,:),x);
end
end

function errorValue = relative_error(observed,predicted)
if isempty(observed) || any(~isfinite(predicted))
    errorValue = NaN;
else
    errorValue = norm(observed-predicted)/max(norm(observed),eps);
end
end

function change = weighted_relative_change(field,reference,x,y)
defectSquared = abs(field-reference).^2;
referenceSquared = abs(reference).^2;
change = sqrt(trapz(y,trapz(x,defectSquared,2))) / ...
    max(sqrt(trapz(y,trapz(x,referenceSquared,2))),eps);
end

function opts = analysis_options(result,userOpts,parameters)
positiveExtent = max(result.x-parameters.centerX);
opts = struct('centerX',parameters.centerX, ...
    'tailXRange',[max(2,0.15*positiveExtent),0.75*positiveExtent], ...
    'tailYMax',min(8,max(result.y)), ...
    'maximumComparisonRatio',1.25,'makePlot',true, ...
    'outputFile',fullfile('result','verification', ...
    'bessel_k1_two_scale_diagnostics.png'));
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
validateattributes(opts.tailXRange,{'numeric'}, ...
    {'vector','numel',2,'positive','finite','increasing'});
validateattributes(opts.tailYMax,{'numeric'},{'scalar','positive','finite'});
validateattributes(opts.makePlot,{'logical','numeric'},{'scalar'});
opts.makePlot = logical(opts.makePlot);
end

function make_diagnostic_figure(diagnostics,fields,result,opts)
initialRho = fields{1};
finalRho = fields{end};
finalFit = diagnostics.fits(end);
positive = result.x >= opts.centerX;

figureHandle = figure('Color','w','Position',[50,50,1500,840]);
layout = tiledlayout(figureHandle,2,3,'TileSpacing','compact', ...
    'Padding','compact');

axisHandle = nexttile(layout);
imagesc(axisHandle,result.x(positive)-opts.centerX,result.y, ...
    initialRho(:,positive));
set(axisHandle,'YDir','normal'); axis(axisHandle,'tight'); colorbar(axisHandle);
xlabel(axisHandle,'x-x_0'); ylabel(axisHandle,'y'); title(axisHandle,'initial \rho');

axisHandle = nexttile(layout);
imagesc(axisHandle,result.x(positive)-opts.centerX,result.y, ...
    finalRho(:,positive));
set(axisHandle,'YDir','normal'); axis(axisHandle,'tight'); colorbar(axisHandle);
xlabel(axisHandle,'x-x_0'); ylabel(axisHandle,'y');
title(axisHandle,sprintf('final rho, t=%.3g',diagnostics.times(end)));

axisHandle = nexttile(layout);
finalRhoX = finalFit.rhoX(:,positive);
limit = max(abs(finalRhoX),[],'all');
imagesc(axisHandle,result.x(positive)-opts.centerX,result.y,finalRhoX);
set(axisHandle,'YDir','normal'); axis(axisHandle,'tight'); colorbar(axisHandle);
clim(axisHandle,[-limit,limit]); hold(axisHandle,'on');
plot(axisHandle,finalFit.tailX,finalFit.observedZeroY,'w.','MarkerSize',9);
plot(axisHandle,finalFit.tailX,finalFit.predictedZeroY,'k--','LineWidth',1.4);
xlabel(axisHandle,'x-x_0'); ylabel(axisHandle,'y');
title(axisHandle,'final \rho_x and zero curve');

axisHandle = nexttile(layout);
plot(axisHandle,diagnostics.times,diagnostics.mu,'o-','LineWidth',1.3);
hold(axisHandle,'on');
yline(axisHandle,diagnostics.initialParameters.mu,'--');
yyaxis(axisHandle,'right');
plot(axisHandle,diagnostics.times,diagnostics.wallLeakage,'s-','LineWidth',1.3);
ylabel(axisHandle,'wall leakage'); yyaxis(axisHandle,'left');
ylabel(axisHandle,'\mu_{fit}'); xlabel(axisHandle,'physical time');
grid(axisHandle,'on'); title(axisHandle,'tail parameter and wall monitor');

axisHandle = nexttile(layout);
semilogy(axisHandle,diagnostics.times,diagnostics.profileError, ...
    'o-','LineWidth',1.3); hold(axisHandle,'on');
semilogy(axisHandle,diagnostics.times,diagnostics.fullKernelError, ...
    's-','LineWidth',1.3);
semilogy(axisHandle,diagnostics.times,diagnostics.rhoXZeroRelativeError, ...
    'd-','LineWidth',1.3);
grid(axisHandle,'on'); xlabel(axisHandle,'physical time'); ylabel(axisHandle,'relative error');
legend(axisHandle,'normalized profile','full K_1 heat kernel', ...
    '\rho_x zero line','Location','best'); title(axisHandle,'two-scale defects');

axisHandle = nexttile(layout);
selected = unique(round(linspace(1,numel(finalFit.tailX),5)));
for index = 1:numel(selected)
    [~,column] = min(abs((result.x-opts.centerX)- ...
        finalFit.tailX(selected(index))));
    xValue = finalFit.tailX(selected(index));
    profile = finalRho(:,column)/max(finalRho(:,column));
    eta = result.y*sqrt(finalFit.mu/xValue);
    plot(axisHandle,eta,profile,'LineWidth',1.1); hold(axisHandle,'on');
end
etaReference = linspace(0,5,240);
shape = (etaReference/sqrt(2)).*exp((2-etaReference.^2)/4);
plot(axisHandle,etaReference,shape,'k--','LineWidth',1.8);
xlim(axisHandle,[0,5]); ylim(axisHandle,[0,1.1]); grid(axisHandle,'on');
xlabel(axisHandle,'eta = y sqrt(mu/x)'); ylabel(axisHandle,'normalized rho');
title(axisHandle,'final K_1 parabolic collapse');

title(layout,sprintf(['regularized K_1 evolution; Poisson solve enforces ' ...
    'psi=0 on y=0; peak exponent %.3f'], ...
    diagnostics.peakWidthExponent(end)));
outputDirectory = fileparts(opts.outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle,opts.outputFile,'Resolution',180);
end
