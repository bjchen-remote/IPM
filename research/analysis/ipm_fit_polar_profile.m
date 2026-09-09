function report = ipm_fit_polar_profile(source,userOpts)
%IPM_FIT_POLAR_PROFILE Fit a one-sided IPM profile about a wall point.
%   REPORT = IPM_FIT_POLAR_PROFILE(SOURCE,OPTS) accepts a ipm.solve result
%   structure or a MAT filename. The stored solver coordinates are treated
%   as coordinates relative to opts.sourceCenter. opts.displayCenter only
%   translates labels, so the same symmetric computation can be displayed
%   with its blow-up point at (1,0).

if nargin < 2 || isempty(userOpts)
    userOpts = struct();
end
profile = load_profile(source);
opts = profile_options(userOpts);
[x,y,omega,fieldName] = select_field(profile,opts.coordinate);

validateattributes(x,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(y,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(omega,{'numeric'},{'2d','real','finite'});
x = x(:)';
y = y(:);
if ~isequal(size(omega),[numel(y),numel(x)])
    error('ipm:PolarProfileSize', ...
        'omega must have size numel(y)-by-numel(x).');
end

[~,wallRow] = min(abs(y-opts.sourceCenter(2)));
rWall = x-opts.sourceCenter(1);
right = rWall > 0;
if nnz(right) < opts.minimumFitPoints
    error('ipm:PolarProfileWall', ...
        'Too few wall nodes to the right of sourceCenter.');
end
rWall = rWall(right);
wallOmega = omega(wallRow,right);
[~,peakIndex] = max(abs(wallOmega));
orientation = sign(wallOmega(peakIndex));
if orientation == 0
    orientation = 1;
end
wallOmega = orientation*wallOmega;
wallFit = fit_tail(rWall,wallOmega,opts);

polar = sample_polar(x,y,orientation*omega,wallFit,opts);
report = struct();
report.coordinate = opts.coordinate;
report.field = fieldName;
report.sourceCenter = opts.sourceCenter;
report.displayCenter = opts.displayCenter;
report.coordinateTranslation = opts.displayCenter-opts.sourceCenter;
report.orientation = orientation;
report.wall = wallFit;
report.polar = polar;
report.hypothesis = hypothesis_summary(wallFit.alpha);
report.acceptedSingleScale = wallFit.acceptedPlateau && ...
    abs(wallFit.alpha-polar.angularL2Fit.alpha) <= 0.08 && ...
    polar.separabilityError <= 0.15;
report.options = opts;

if opts.makePlot
    make_profile_plot(report,opts);
end
end

function profile = load_profile(source)
if ischar(source) || (isstring(source) && isscalar(source))
    loaded = load(char(source));
    names = fieldnames(loaded);
    isStructure = cellfun(@(name) isstruct(loaded.(name)),names);
    structures = names(isStructure);
    if numel(structures) ~= 1
        error('ipm:PolarProfileFile', ...
            'The MAT file must contain exactly one result structure.');
    end
    profile = loaded.(structures{1});
elseif isstruct(source)
    profile = source;
else
    error('ipm:PolarProfileInput', ...
        'source must be a ipm.solve result structure or MAT filename.');
end
isSolverResult = isfield(profile,'schemaVersion') || ...
    all(isfield(profile,{'rho','opts','history','physicalTime'}));
if isSolverResult
    result = ipm.output.validate(profile);
    profile = struct('x',result.grid.x,'y',result.grid.y, ...
        'omega',result.state.omega, ...
        'physicalX',result.grid.physicalX, ...
        'physicalY',result.grid.physicalY, ...
        'physicalOmega',result.physical.omega);
end
end

function opts = profile_options(userOpts)
opts = struct();
opts.coordinate = 'rescaled';
opts.sourceCenter = [0,0];
opts.displayCenter = [1,0];
opts.radialWindow = [];
opts.tailStartFactor = 1.5;
opts.boundaryFraction = 0.65;
opts.minimumAmplitudeFraction = 1e-7;
opts.minimumFitPoints = 24;
opts.minimumTailDecades = 0.45;
opts.radialSamples = 160;
opts.angularSamples = 121;
opts.maximumAngle = pi/2;
opts.fanThresholds = [0.1,0.25,0.5];
% Neutral diagnostic grid only; no entry is a theoretically selected exponent.
opts.candidateExponents = [1/3,1/2,2/3,3/4];
opts.makePlot = false;
opts.plotFile = '';
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
validateattributes(opts.sourceCenter,{'numeric'}, ...
    {'vector','numel',2,'real','finite'});
validateattributes(opts.displayCenter,{'numeric'}, ...
    {'vector','numel',2,'real','finite'});
if ~isempty(opts.radialWindow)
    validateattributes(opts.radialWindow,{'numeric'}, ...
        {'vector','numel',2,'real','finite','positive','increasing'});
end
validateattributes(opts.maximumAngle,{'numeric'}, ...
    {'scalar','>',0,'<=',pi/2});
validateattributes(opts.fanThresholds,{'numeric'}, ...
    {'vector','>',0,'<',1,'increasing'});
validateattributes(opts.candidateExponents,{'numeric'}, ...
    {'vector','>',0,'<',1});
opts.coordinate = lower(char(opts.coordinate));
opts.sourceCenter = opts.sourceCenter(:)';
opts.displayCenter = opts.displayCenter(:)';
end

function [x,y,omega,fieldName] = select_field(profile,coordinate)
switch coordinate
    case 'rescaled'
        required = {'x','y','omega'};
        fieldName = 'omega=d_X1 R';
        xName = 'x';
        yName = 'y';
        omegaName = 'omega';
    case 'physical'
        required = {'physicalX','physicalY','physicalOmega'};
        fieldName = 'physicalOmega=d_x1 rho';
        xName = 'physicalX';
        yName = 'physicalY';
        omegaName = 'physicalOmega';
    otherwise
        error('ipm:PolarProfileCoordinate', ...
            'coordinate must be rescaled or physical.');
end
for index = 1:numel(required)
    if ~isfield(profile,required{index})
        error('ipm:PolarProfileField','Missing result field %s.', ...
            required{index});
    end
end
x = profile.(xName);
y = profile.(yName);
omega = profile.(omegaName);
end

function fit = fit_tail(r,omega,opts)
amplitude = abs(omega);
[peakAmplitude,peakIndex] = max(amplitude);
peakRadius = r(peakIndex);
valid = isfinite(amplitude) & amplitude > ...
    opts.minimumAmplitudeFraction*peakAmplitude;
if isempty(opts.radialWindow)
    lowerRadius = max(opts.tailStartFactor*peakRadius,r(peakIndex+1));
    upperRadius = opts.boundaryFraction*r(end);
else
    lowerRadius = opts.radialWindow(1);
    upperRadius = opts.radialWindow(2);
end
valid = valid & r >= lowerRadius & r <= upperRadius;
if nnz(valid) < opts.minimumFitPoints
    error('ipm:PolarProfileTail', ...
        'The proposed post-peak tail has too few usable nodes.');
end

rv = r(valid);
av = amplitude(valid);
logRadius = linspace(log(rv(1)),log(rv(end)),opts.radialSamples)';
logAmplitude = interp1(log(rv),log(av),logRadius,'linear');
if isempty(opts.radialWindow)
    selection = select_plateau(logRadius,logAmplitude,opts);
else
    selection = (1:numel(logRadius))';
end
fit = linear_power_fit(logRadius(selection),logAmplitude(selection));
fit.peakRadius = peakRadius;
fit.peakAmplitude = peakAmplitude;
fit.radius = exp(logRadius);
fit.amplitude = exp(logAmplitude);
fit.localAlpha = -gradient(logAmplitude,logRadius);
fit.selection = selection;
fit.radialWindow = exp([logRadius(selection(1)), ...
    logRadius(selection(end))]);
fit.logSpanDecades = diff(log10(fit.radialWindow));
fit.candidates = candidate_fits(logRadius(selection), ...
    logAmplitude(selection),opts.candidateExponents);
fit.acceptedPlateau = fit.rSquared >= 0.995 && ...
    fit.logSpanDecades >= opts.minimumTailDecades && ...
    fit.localSlopeStd <= 0.08;
end

function selection = select_plateau(logRadius,logAmplitude,opts)
n = numel(logRadius);
minimumCount = max(opts.minimumFitPoints,round(0.22*n));
minimumSpan = opts.minimumTailDecades*log(10);
bestScore = Inf;
selection = (1:n)';
for first = 1:4:max(1,n-minimumCount+1)
    for last = first+minimumCount-1:4:n
        indices = (first:last)';
        span = logRadius(last)-logRadius(first);
        if span < minimumSpan
            continue
        end
        candidate = linear_power_fit(logRadius(indices), ...
            logAmplitude(indices));
        score = candidate.logRmse+0.45*candidate.localSlopeStd+ ...
            0.025/span+0.01/sqrt(numel(indices));
        if score < bestScore
            bestScore = score;
            selection = indices;
        end
    end
end
end

function fit = linear_power_fit(logRadius,logAmplitude)
coefficients = polyfit(logRadius,logAmplitude,1);
prediction = polyval(coefficients,logRadius);
residual = logAmplitude-prediction;
centered = logAmplitude-mean(logAmplitude);
localSlope = -gradient(logAmplitude,logRadius);
fit.alpha = -coefficients(1);
fit.amplitudeCoefficient = exp(coefficients(2));
fit.logRmse = sqrt(mean(residual.^2));
fit.rSquared = 1-sum(residual.^2)/max(sum(centered.^2),eps);
fit.localSlopeStd = std(localSlope);
fit.predictedLogAmplitude = prediction;
end

function candidates = candidate_fits(logRadius,logAmplitude,exponents)
candidates = repmat(struct('alpha',0,'coefficient',0,'logRmse',0), ...
    numel(exponents),1);
for index = 1:numel(exponents)
    alpha = exponents(index);
    intercept = mean(logAmplitude+alpha*logRadius);
    residual = logAmplitude-(intercept-alpha*logRadius);
    candidates(index).alpha = alpha;
    candidates(index).coefficient = exp(intercept);
    candidates(index).logRmse = sqrt(mean(residual.^2));
end
end

function polar = sample_polar(x,y,omega,wallFit,opts)
rMin = wallFit.radialWindow(1);
rMax = wallFit.radialWindow(2);
rMax = min(rMax,x(end)-opts.sourceCenter(1));
if opts.maximumAngle > 0
    rMax = min(rMax,(y(end)-opts.sourceCenter(2))/ ...
        sin(opts.maximumAngle));
end
if rMax <= rMin
    error('ipm:PolarProfileDomain', ...
        'The polar fitting annulus does not fit inside the result domain.');
end
radii = logspace(log10(rMin),log10(rMax),opts.radialSamples)';
theta = linspace(0,opts.maximumAngle,opts.angularSamples);
[radiusGrid,thetaGrid] = ndgrid(radii,theta);
xQuery = opts.sourceCenter(1)+radiusGrid.*cos(thetaGrid);
yQuery = opts.sourceCenter(2)+radiusGrid.*sin(thetaGrid);
interpolant = griddedInterpolant({y,x},omega,'linear','none');
sample = interpolant(yQuery,xQuery);

angularNorm = sqrt(trapz(theta,sample.^2,2));
radialFit = linear_power_fit(log(radii),log(max(angularNorm,realmin)));
scaled = sample.*radii.^wallFit.alpha;
angularProfile = mean(scaled,1,'omitnan');
model = repmat(angularProfile,numel(radii),1)./ ...
    radii.^wallFit.alpha;
valid = isfinite(sample) & isfinite(model);
separabilityError = norm(sample(valid)-model(valid))/ ...
    max(norm(sample(valid)),eps);

fanAngles = nan(numel(radii),numel(opts.fanThresholds));
wallConnected = false(numel(radii),numel(opts.fanThresholds));
for radiusIndex = 1:numel(radii)
    row = abs(sample(radiusIndex,:));
    row = row/max(row,[],'omitnan');
    for thresholdIndex = 1:numel(opts.fanThresholds)
        [fanAngles(radiusIndex,thresholdIndex), ...
            wallConnected(radiusIndex,thresholdIndex)] = ...
            angular_crossing(theta,row,opts.fanThresholds(thresholdIndex));
    end
end
polar = struct();
polar.radii = radii;
polar.theta = theta;
polar.omega = sample;
polar.scaledOmega = scaled;
polar.angularProfile = angularProfile;
polar.angularL2Fit = radialFit;
polar.separabilityError = separabilityError;
polar.fanThresholds = opts.fanThresholds;
polar.fanAngles = fanAngles;
polar.medianFanAngleDegrees = median(fanAngles,1,'omitnan')*180/pi;
polar.wallConnectedFraction = mean(wallConnected,1);
polar.fanPower = fit_fan_power(radii,fanAngles);
end

function fits = fit_fan_power(radii,fanAngles)
fits = repmat(struct('sigma',NaN,'coefficient',NaN,'rSquared',NaN), ...
    1,size(fanAngles,2));
for index = 1:size(fanAngles,2)
    angle = fanAngles(:,index);
    valid = isfinite(angle) & angle > 0;
    if nnz(valid) < 4
        continue
    end
    fit = linear_power_fit(log(radii(valid)),log(angle(valid)));
    fits(index).sigma = -fit.alpha;
    fits(index).coefficient = fit.amplitudeCoefficient;
    fits(index).rSquared = fit.rSquared;
end
end

function [crossing,isWallConnected] = angular_crossing(theta,row,threshold)
mask = row >= threshold;
isWallConnected = mask(1);
if isWallConnected
    last = find(~mask,1,'first')-1;
    if isempty(last)
        crossing = theta(end);
        return
    end
else
    [~,maximumIndex] = max(row);
    last = maximumIndex;
    while last < numel(mask) && mask(last+1)
        last = last+1;
    end
end
if last >= numel(theta)
    crossing = theta(end);
else
    leftValue = row(last)-threshold;
    rightValue = row(last+1)-threshold;
    fraction = leftValue/max(leftValue-rightValue,eps);
    crossing = theta(last)+fraction*(theta(last+1)-theta(last));
end
end

function hypothesis = hypothesis_summary(alpha)
hypothesis = struct();
hypothesis.alpha = alpha;
hypothesis.beta = 1-alpha;
hypothesis.lambda = 1/alpha-1;
hypothesis.rateRatioCOmegaOverCL = 1-alpha;
hypothesis.wallOmegaForm = ...
    'd_x rho(T) ~ A*(x_1-a)_+^{-alpha}';
hypothesis.rescaledDensityFarField = ...
    'R(r,theta) ~ r^{1-alpha}*F(theta)';
end

function make_profile_plot(report,opts)
figureHandle = figure('Color','w','Position',[100,100,1200,820]);
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');
fit = report.wall;
indices = fit.selection;

nexttile(layout);
loglog(fit.radius,fit.amplitude,'k-','LineWidth',1.2);
hold on
fitRadius = fit.radius(indices);
fitAmplitude = fit.amplitudeCoefficient*fitRadius.^(-fit.alpha);
loglog(fitRadius,fitAmplitude,'r--','LineWidth',1.6);
xlabel('r from blow-up point');
ylabel('|\omega(r,0)|');
grid on
title(sprintf('wall tail: alpha=%.4f, R^2=%.5f', ...
    fit.alpha,fit.rSquared));

nexttile(layout);
semilogx(fit.radius,fit.localAlpha,'k-','LineWidth',1.2);
hold on
yline(fit.alpha,'r--','fit');
xline(fit.radialWindow(1),'Color',[0.5,0.5,0.5]);
xline(fit.radialWindow(2),'Color',[0.5,0.5,0.5]);
xlabel('r');
ylabel('-d log|\omega|/d log r');
grid on
title(sprintf('local exponent, std=%.3g',fit.localSlopeStd));

nexttile(layout);
[radiusGrid,thetaGrid] = ndgrid(report.polar.radii,report.polar.theta);
xPlot = opts.displayCenter(1)+radiusGrid.*cos(thetaGrid);
yPlot = opts.displayCenter(2)+radiusGrid.*sin(thetaGrid);
surface(xPlot,yPlot,zeros(size(xPlot)),report.polar.scaledOmega, ...
    'EdgeColor','none');
view(2)
axis equal tight
colorbar
xlabel('x_1');
ylabel('x_2');
title(sprintf('r^{alpha}|omega|, alpha=%.4f',fit.alpha));

nexttile(layout);
chosen = unique(round(linspace(1,numel(report.polar.radii),7)));
plot(report.polar.theta*180/pi, ...
    report.polar.scaledOmega(chosen,:),'LineWidth',1.0);
xlabel('theta (degrees)');
ylabel('r^{alpha} omega(r,theta)');
grid on
title(sprintf('angular collapse, relative error %.3g', ...
    report.polar.separabilityError));
title(layout,sprintf('IPM polar profile about (%g,%g)', ...
    opts.displayCenter(1),opts.displayCenter(2)));

if ~isempty(opts.plotFile)
    folder = fileparts(opts.plotFile);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    exportgraphics(figureHandle,opts.plotFile,'Resolution',180);
end
end
