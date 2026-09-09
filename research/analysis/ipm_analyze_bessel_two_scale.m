function diagnostics = ipm_analyze_bessel_two_scale(inputResult,userOpts)
%IPM_ANALYZE_BESSEL_TWO_SCALE Test a parabolic two-scale Bessel tail.
%   DIAGNOSTICS = IPM_ANALYZE_BESSEL_TWO_SCALE(RESULT) fits each stored
%   physical snapshot to
%
%     rho(x,y) = C*x^(-1/2)*exp(-mu*(y+h)^2/(4*x)),   x>0,
%
%   and independently checks the predicted rho_x zero curve
%
%     y_zero(x) = sqrt(2*x/mu)-h.
%
%   INPUTRESULT may be a result structure or the path to a MAT-file holding
%   one.  The experiment cases deliberately disable remeshing, so all stored
%   fields share RESULT.x and RESULT.y.

if nargin < 2 || isempty(userOpts)
    userOpts = struct();
end
canonicalResult = ipm.output.validate(inputResult);
result = canonicalResult.grid;
caseMetadata = ipm_bessel_case_metadata(canonicalResult, ...
    {'mu','poleDepth','centerX','perturbationAmplitude'});
parameters = caseMetadata.parameters;
opts = analysis_options(result,userOpts,parameters);
if ~strcmp(canonicalResult.config.scaling.rescalingMode,'physical')
    error('ipm:BesselTwoScalePhysicalOnly', ...
        'This diagnostic currently expects physical-mode snapshots.');
end

snapshots = ipm.output.snapshotAt(canonicalResult,'all');
fields = {snapshots.physicalRho}';
times = [snapshots.physicalTime]';
template = struct('time',NaN,'mu',NaN,'poleDepth',NaN, ...
    'normalizedProfileError',NaN,'fullKernelError',NaN, ...
    'wallAmplitudeVariation',NaN,'widthExponent',NaN, ...
    'rhoXZeroRelativeError',NaN,'rhoXZeroCoverage',0, ...
    'tailX',[],'eFoldWidth',[],'observedZeroY',[], ...
    'predictedZeroY',[],'rhoX',[]);
fits = repmat(template,numel(fields),1);
for index = 1:numel(fields)
    fits(index) = analyze_snapshot(fields{index},times(index), ...
        result.x,result.y,opts);
end

diagnostics = struct();
diagnostics.times = times(:);
diagnostics.mu = [fits.mu]';
diagnostics.poleDepth = [fits.poleDepth]';
diagnostics.normalizedProfileError = [fits.normalizedProfileError]';
diagnostics.fullKernelError = [fits.fullKernelError]';
diagnostics.wallAmplitudeVariation = [fits.wallAmplitudeVariation]';
diagnostics.widthExponent = [fits.widthExponent]';
diagnostics.rhoXZeroRelativeError = [fits.rhoXZeroRelativeError]';
diagnostics.rhoXZeroCoverage = [fits.rhoXZeroCoverage]';
diagnostics.fits = fits;
diagnostics.options = opts;
diagnostics.caseMetadata = caseMetadata;
diagnostics.initialParameters = struct('mu',parameters.mu, ...
    'poleDepth',parameters.poleDepth, ...
    'perturbationAmplitude',parameters.perturbationAmplitude);
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
wallValues = tailRho(1,:);
[xGrid,yGrid] = meshgrid(tailX,y);
ratio = tailRho./wallValues;

fitMask = yGrid > 0 & yGrid <= opts.tailYMax & ...
    ratio >= opts.fitMinimumFraction & ratio <= opts.fitMaximumFraction & ...
    isfinite(ratio) & wallValues > 0;
design = [yGrid(fitMask).^2,yGrid(fitMask)];
response = -xGrid(fitMask).*log(ratio(fitMask));
if size(design,1) < opts.minimumFitPoints
    coefficients = [NaN;NaN];
else
    coefficients = design\response;
end
mu = 4*coefficients(1);
poleDepth = coefficients(2)/(2*coefficients(1));
if ~isfinite(mu) || mu <= 0 || ~isfinite(poleDepth)
    mu = NaN;
    poleDepth = NaN;
end

comparisonMask = yGrid <= opts.tailYMax & ratio > 0 & ...
    ratio <= opts.maximumComparisonRatio & isfinite(ratio);
predictedRatio = exp(-mu*(yGrid.^2+2*poleDepth*yGrid)./(4*xGrid));
normalizedProfileError = relative_error(ratio(comparisonMask), ...
    predictedRatio(comparisonMask));

kernelShape = xGrid.^(-1/2).* ...
    exp(-mu*(yGrid+poleDepth).^2./(4*xGrid));
fullMask = comparisonMask & isfinite(kernelShape);
coefficient = dot(kernelShape(fullMask),tailRho(fullMask)) / ...
    max(dot(kernelShape(fullMask),kernelShape(fullMask)),eps);
fullKernelError = relative_error(tailRho(fullMask), ...
    coefficient*kernelShape(fullMask));

correctedWallAmplitude = sqrt(tailX).*wallValues.* ...
    exp(mu*poleDepth^2./(4*tailX));
wallAmplitudeVariation = std(correctedWallAmplitude) / ...
    max(abs(mean(correctedWallAmplitude)),eps);

eFoldWidth = nan(size(tailX));
for column = 1:numel(tailX)
    eFoldWidth(column) = first_level_crossing(y,ratio(:,column),exp(-1));
end
correctedWidth = sqrt(max((eFoldWidth+poleDepth).^2-poleDepth^2,0));
widthMask = isfinite(correctedWidth) & correctedWidth > 0;
if nnz(widthMask) >= 3
    exponentFit = polyfit(log(tailX(widthMask)), ...
        log(correctedWidth(widthMask)),1);
    widthExponent = exponentFit(1);
else
    widthExponent = NaN;
end

rhoX = differentiate_x(rho,x);
observedZeroY = nan(size(tailX));
predictedZeroY = sqrt(2*tailX/mu)-poleDepth;
for column = 1:numel(tailX)
    trace = rhoX(:,tailColumns(column));
    observedZeroY(column) = first_negative_to_positive_crossing(y,trace);
end
zeroMask = isfinite(observedZeroY) & isfinite(predictedZeroY) & ...
    predictedZeroY > 0 & predictedZeroY <= opts.tailYMax;
rhoXZeroRelativeError = relative_error(observedZeroY(zeroMask), ...
    predictedZeroY(zeroMask));
rhoXZeroCoverage = nnz(zeroMask)/numel(tailX);

fit = struct('time',time,'mu',mu,'poleDepth',poleDepth, ...
    'normalizedProfileError',normalizedProfileError, ...
    'fullKernelError',fullKernelError, ...
    'wallAmplitudeVariation',wallAmplitudeVariation, ...
    'widthExponent',widthExponent, ...
    'rhoXZeroRelativeError',rhoXZeroRelativeError, ...
    'rhoXZeroCoverage',rhoXZeroCoverage,'tailX',tailX, ...
    'eFoldWidth',eFoldWidth,'observedZeroY',observedZeroY, ...
    'predictedZeroY',predictedZeroY,'rhoX',rhoX);
end

function crossing = first_level_crossing(y,trace,level)
crossing = NaN;
valid = isfinite(trace);
if ~all(valid) || trace(1) < level
    return;
end
index = find(trace(1:end-1) >= level & trace(2:end) < level,1,'first');
if isempty(index)
    return;
end
weight = (level-trace(index)) / ...
    (trace(index+1)-trace(index));
crossing = y(index)+weight*(y(index+1)-y(index));
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
defectNorm = sqrt(trapz(y,trapz(x,defectSquared,2)));
referenceNorm = sqrt(trapz(y,trapz(x,referenceSquared,2)));
change = defectNorm/max(referenceNorm,eps);
end

function opts = analysis_options(result,userOpts,parameters)
opts = struct('centerX',parameters.centerX, ...
    'tailXRange',[],'tailYMax',min(8,max(result.y)), ...
    'fitMinimumFraction',0.03,'fitMaximumFraction',0.97, ...
    'maximumComparisonRatio',1.25,'minimumFitPoints',40, ...
    'makePlot',true,'outputFile',fullfile('result','verification', ...
    'bessel_two_scale_diagnostics.png'));
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
centeredX = result.x-opts.centerX;
positiveExtent = max(centeredX);
if isempty(opts.tailXRange)
    opts.tailXRange = [max(2,0.15*positiveExtent),0.75*positiveExtent];
end
validateattributes(opts.centerX,{'numeric'},{'scalar','real','finite'});
validateattributes(opts.tailXRange,{'numeric'}, ...
    {'vector','numel',2,'real','finite','positive','increasing'});
validateattributes(opts.tailYMax,{'numeric'},{'scalar','real','positive'});
validateattributes(opts.fitMinimumFraction,{'numeric'}, ...
    {'scalar','real','positive','<',1});
validateattributes(opts.fitMaximumFraction,{'numeric'}, ...
    {'scalar','real','>',opts.fitMinimumFraction});
validateattributes(opts.maximumComparisonRatio,{'numeric'}, ...
    {'scalar','real','>',1});
validateattributes(opts.minimumFitPoints,{'numeric'}, ...
    {'scalar','integer','>=',10});
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
yyaxis(axisHandle,'left');
plot(axisHandle,diagnostics.times,diagnostics.mu,'o-','LineWidth',1.3);
hold(axisHandle,'on');
yline(axisHandle,diagnostics.initialParameters.mu,'--');
ylabel(axisHandle,'\mu_{fit}');
yyaxis(axisHandle,'right');
plot(axisHandle,diagnostics.times,diagnostics.poleDepth,'s-','LineWidth',1.3);
yline(axisHandle,diagnostics.initialParameters.poleDepth,'--');
ylabel(axisHandle,'h_{fit}');
xlabel(axisHandle,'physical time'); grid(axisHandle,'on');
title(axisHandle,'tail parameters');

axisHandle = nexttile(layout);
semilogy(axisHandle,diagnostics.times, ...
    diagnostics.normalizedProfileError,'o-','LineWidth',1.3);
hold(axisHandle,'on');
semilogy(axisHandle,diagnostics.times,diagnostics.fullKernelError, ...
    's-','LineWidth',1.3);
semilogy(axisHandle,diagnostics.times,diagnostics.rhoXZeroRelativeError, ...
    'd-','LineWidth',1.3);
grid(axisHandle,'on'); xlabel(axisHandle,'physical time'); ylabel(axisHandle,'relative error');
legend(axisHandle,'normalized profile','full heat kernel','\rho_x zero line', ...
    'Location','best'); title(axisHandle,'two-scale defects');

axisHandle = nexttile(layout);
selectedColumns = unique(round(linspace(1,numel(finalFit.tailX),5)));
tailColumns = zeros(size(selectedColumns));
for index = 1:numel(selectedColumns)
    [~,tailColumns(index)] = min(abs((result.x-opts.centerX)- ...
        finalFit.tailX(selectedColumns(index))));
    xValue = finalFit.tailX(selectedColumns(index));
    ratio = finalRho(:,tailColumns(index))/finalRho(1,tailColumns(index));
    chi = finalFit.mu*(result.y.^2+2*finalFit.poleDepth*result.y) / ...
        (4*xValue);
    plot(axisHandle,chi,ratio,'LineWidth',1.1); hold(axisHandle,'on');
end
chiReference = linspace(0,4,200);
plot(axisHandle,chiReference,exp(-chiReference),'k--','LineWidth',1.8);
xlim(axisHandle,[0,4]); ylim(axisHandle,[0,1.1]); grid(axisHandle,'on');
xlabel(axisHandle,'\chi=\mu(y^2+2hy)/(4x)'); ylabel(axisHandle,'\rho(x,y)/\rho(x,0)');
title(axisHandle,'final parabolic collapse');

title(layout,sprintf(['smooth Bessel-tail evolution: perturbation %.3g, ' ...
    'width exponent %.3f'], ...
    diagnostics.initialParameters.perturbationAmplitude, ...
    diagnostics.widthExponent(end)));
outputDirectory = fileparts(opts.outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle,opts.outputFile,'Resolution',180);
end
