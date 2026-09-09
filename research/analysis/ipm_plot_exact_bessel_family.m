function result = ipm_plot_exact_bessel_family(userOpts)
%IPM_PLOT_EXACT_BESSEL_FAMILY Plot exact far-decaying IPM steady modes.
%   RESULT = IPM_PLOT_EXACT_BESSEL_FAMILY(OPTS) evaluates
%
%     psi = A*exp(-lambda*X/2)*K_nu(|lambda|*r/2)*Theta(theta),
%     rho = lambda*psi,
%
%   on a polar sector about (1,0), together with analytic rho_x and rho_y.
%   The default gallery compares a wall-fed cosine mode with a sine mode that
%   vanishes at both sector rays. Set thetaStar=pi and modes [0,1] to compare
%   the global wall-fed and impermeable upper-half-plane members. Negative
%   lambda makes the slow tail point to the right. The center is omitted
%   because K_nu is singular there.

if nargin < 1 || isempty(userOpts)
    userOpts = struct();
end
opts = plot_options(userOpts);

radialCoordinate = logspace(log10(opts.innerRadius), ...
    log10(opts.outerRadius),opts.radialPoints);
angularCoordinate = linspace(0,opts.thetaStar,opts.angularPoints)';
[radius,theta] = meshgrid(radialCoordinate,angularCoordinate);
centeredX = radius.*cos(theta);
y = opts.center(2)+radius.*sin(theta);
x = opts.center(1)+centeredX;

numberOfModes = numel(opts.modeOrders);
modes = evaluate_mode(radius,theta,centeredX,opts.modeOrders(1), ...
    opts.modeCosine(1),opts.modeSine(1),opts);
modes = repmat(modes,numberOfModes,1);
for index = 2:numberOfModes
    modes(index) = evaluate_mode(radius,theta,centeredX, ...
        opts.modeOrders(index),opts.modeCosine(index), ...
        opts.modeSine(index),opts);
end

result = struct();
result.center = opts.center;
result.thetaStar = opts.thetaStar;
result.x = x;
result.y = y;
result.radius = radius;
result.theta = theta;
result.modes = modes;
result.options = opts;
result.maximumRelativePoissonResidual = ...
    max([modes.relativePoissonResidual]);
result.maximumRelativeTransportResidual = ...
    max([modes.relativeTransportResidual]);
result.maximumRelativeEdgeNormalVelocity = ...
    max([modes.relativeEdgeNormalVelocity]);

if opts.makePlot
    make_gallery(result,opts);
    make_detail(result,opts);
end
end

function mode = evaluate_mode(radius,theta,centeredX,order,cosineWeight, ...
        sineWeight,opts)
lambda = opts.lambda;
k = abs(lambda)/2;
argument = k*radius;
besselValue = besselk(order,argument);
besselDerivative = -0.5*(besselk(abs(order-1),argument)+ ...
    besselk(order+1,argument));
besselSecond = (1+order^2./argument.^2).*besselValue- ...
    besselDerivative./argument;

angularValue = cosineWeight*cos(order*theta)+ ...
    sineWeight*sin(order*theta);
angularDerivative = order*(-cosineWeight*sin(order*theta)+ ...
    sineWeight*cos(order*theta));
angularSecond = -order^2*angularValue;
exponential = exp(-0.5*lambda*centeredX);

psi = opts.amplitude*exponential.*besselValue.*angularValue;
rho = lambda*psi;
rhoX = lambda*opts.amplitude*exponential.*( ...
    -0.5*lambda*besselValue.*angularValue+ ...
    k*cos(theta).*besselDerivative.*angularValue- ...
    sin(theta)./radius.*besselValue.*angularDerivative);
rhoY = lambda*opts.amplitude*exponential.*( ...
    k*sin(theta).*besselDerivative.*angularValue+ ...
    cos(theta)./radius.*besselValue.*angularDerivative);
u1 = -rhoY/lambda;
u2 = rhoX/lambda;

phi = besselValue.*angularValue;
laplacianPhi = k^2*besselSecond.*angularValue+ ...
    k./radius.*besselDerivative.*angularValue+ ...
    besselValue./radius.^2.*angularSecond;
helmholtzResidual = laplacianPhi-k^2*phi;
poissonResidual = -opts.amplitude*exponential.*helmholtzResidual;
transportResidual = u1.*rhoX+u2.*rhoY;

edgeRows = [1,size(theta,1)];
edgeTheta = theta(edgeRows,:);
edgeNormalVelocity = -sin(edgeTheta).*u1(edgeRows,:)+ ...
    cos(edgeTheta).*u2(edgeRows,:);
edgePsi = psi(edgeRows,:);

poissonScale = max(abs(rhoX),[],'all');
transportScale = max(hypot(u1,u2).*hypot(rhoX,rhoY),[],'all');
velocityScale = max(hypot(u1,u2),[],'all');
psiScale = max(abs(psi),[],'all');
mode = struct();
mode.order = order;
mode.cosineWeight = cosineWeight;
mode.sineWeight = sineWeight;
mode.angularValue = angularValue;
mode.psi = psi;
mode.rho = rho;
mode.rhoX = rhoX;
mode.rhoY = rhoY;
mode.u1 = u1;
mode.u2 = u2;
mode.maximumPoissonResidual = max(abs(poissonResidual),[],'all');
mode.maximumTransportResidual = max(abs(transportResidual),[],'all');
mode.relativePoissonResidual = mode.maximumPoissonResidual/ ...
    max(poissonScale,eps);
mode.relativeTransportResidual = mode.maximumTransportResidual/ ...
    max(transportScale,eps);
mode.maximumEdgeNormalVelocity = max(abs(edgeNormalVelocity),[],'all');
mode.relativeEdgeNormalVelocity = mode.maximumEdgeNormalVelocity/ ...
    max(velocityScale,eps);
mode.maximumEdgePsi = max(abs(edgePsi),[],'all');
mode.relativeEdgePsi = mode.maximumEdgePsi/max(psiScale,eps);
end

function opts = plot_options(userOpts)
opts = struct();
opts.center = [1,0];
opts.lambda = -1;
opts.amplitude = -1;
opts.thetaStar = pi/3;
opts.innerRadius = 0.3;
opts.outerRadius = 14;
opts.radialPoints = 361;
opts.angularPoints = 281;
opts.modeOrders = [];
opts.modeCosine = [];
opts.modeSine = [];
opts.detailIndex = 1;
opts.radialCutAngles = [];
opts.colorQuantile = 0.995;
opts.logDisplayDecades = 9;
opts.makePlot = true;
opts.galleryFile = fullfile('result','verification', ...
    'exact_bessel_family_gallery.png');
opts.detailFile = fullfile('result','verification', ...
    'exact_bessel_rho_x_detail.png');
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
if isempty(opts.modeOrders)
    opts.modeOrders = [pi/(2*opts.thetaStar),pi/opts.thetaStar];
    opts.modeCosine = [1,0];
    opts.modeSine = [0,1];
end
if isempty(opts.radialCutAngles)
    opts.radialCutAngles = opts.thetaStar*[0,0.1,0.25,0.5,0.75,0.9,1];
end

validateattributes(opts.center,{'numeric'}, ...
    {'vector','numel',2,'real','finite'});
validateattributes(opts.lambda,{'numeric'}, ...
    {'scalar','real','finite','nonzero'});
validateattributes(opts.amplitude,{'numeric'}, ...
    {'scalar','real','finite','nonzero'});
validateattributes(opts.thetaStar,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<=',pi});
validateattributes(opts.innerRadius,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(opts.outerRadius,{'numeric'}, ...
    {'scalar','real','finite','>',opts.innerRadius});
validateattributes(opts.radialPoints,{'numeric'}, ...
    {'scalar','integer','>=',51});
validateattributes(opts.angularPoints,{'numeric'}, ...
    {'scalar','integer','>=',51});
validateattributes(opts.modeOrders,{'numeric'}, ...
    {'vector','real','finite','nonnegative','nonempty'});
validateattributes(opts.modeCosine,{'numeric'}, ...
    {'vector','real','finite','numel',numel(opts.modeOrders)});
validateattributes(opts.modeSine,{'numeric'}, ...
    {'vector','real','finite','numel',numel(opts.modeOrders)});
validateattributes(opts.detailIndex,{'numeric'}, ...
    {'scalar','integer','>=',1,'<=',numel(opts.modeOrders)});
validateattributes(opts.radialCutAngles,{'numeric'}, ...
    {'vector','real','finite','>=',0,'<=',opts.thetaStar});
validateattributes(opts.colorQuantile,{'numeric'}, ...
    {'scalar','real','finite','>',0.5,'<=',1});
validateattributes(opts.logDisplayDecades,{'numeric'}, ...
    {'scalar','real','finite','>=',2});
validateattributes(opts.makePlot,{'logical','numeric'},{'scalar'});
assert(all(hypot(opts.modeCosine,opts.modeSine)>0), ...
    'ipm:BesselMode','Every angular mode needs a nonzero coefficient.');
opts.center = opts.center(:)';
opts.modeOrders = opts.modeOrders(:)';
opts.modeCosine = opts.modeCosine(:)';
opts.modeSine = opts.modeSine(:)';
opts.radialCutAngles = opts.radialCutAngles(:)';
opts.makePlot = logical(opts.makePlot);
end

function make_gallery(result,opts)
numberOfModes = numel(result.modes);
figureHandle = figure('Color','w','Position',[60,80,520*numberOfModes,840]);
layout = tiledlayout(figureHandle,2,numberOfModes, ...
    'TileSpacing','compact','Padding','compact');
divergingMap = blue_white_red(256);
for index = 1:numberOfModes
    mode = result.modes(index);
    axisHandle = nexttile(layout,index);
    logDensity = normalized_log_magnitude(mode.rho,opts.logDisplayDecades);
    scalar_field(axisHandle,result.x,result.y,logDensity,opts,parula(256), ...
        [-opts.logDisplayDecades,0]);
    title(axisHandle,sprintf(['log_{10}(|\\rho|/max): \\nu=%.3g, ' ...
        '(a,b)=(%.3g,%.3g)'],mode.order,mode.cosineWeight, ...
        mode.sineWeight));

    axisHandle = nexttile(layout,numberOfModes+index);
    signedDerivative = signed_log_magnitude(mode.rhoX, ...
        opts.logDisplayDecades);
    scalar_field(axisHandle,result.x,result.y,signedDerivative,opts, ...
        divergingMap,opts.logDisplayDecades*[-1,1]);
    title(axisHandle,sprintf(['signed-log \\rho_x: PDE %.1e, ' ...
        'edge flux %.1e'],mode.relativePoissonResidual, ...
        mode.relativeEdgeNormalVelocity));
end
title(layout,sprintf(['exact decaying steady IPM modes: ' ...
    'polar sector \\theta_*=%.1f^\\circ, \\lambda=%.3g'], ...
    opts.thetaStar*180/pi,opts.lambda));
export_figure(figureHandle,opts.galleryFile);
end

function make_detail(result,opts)
mode = result.modes(opts.detailIndex);
figureHandle = figure('Color','w','Position',[80,70,1450,940]);
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');
divergingMap = blue_white_red(256);

axisHandle = nexttile(layout);
polar_field(axisHandle,result.x,result.y,mode.rho,opts,divergingMap);
title(axisHandle,'density \rho');

axisHandle = nexttile(layout);
polar_field(axisHandle,result.x,result.y,mode.rhoX,opts,divergingMap);
title(axisHandle,'analytic horizontal derivative \rho_x');

axisHandle = nexttile(layout);
signedDerivative = signed_log_magnitude(mode.rhoX, ...
    opts.logDisplayDecades);
scalar_field(axisHandle,result.x,result.y,signedDerivative,opts, ...
    divergingMap,opts.logDisplayDecades*[-1,1]);
title(axisHandle,'signed logarithmic structure of \rho_x');

axisHandle = nexttile(layout);
hold(axisHandle,'on');
radialCoordinate = result.radius(1,:);
for angle = opts.radialCutAngles
    [~,angularIndex] = min(abs(result.theta(:,1)-angle));
    cut = abs(mode.rhoX(angularIndex,:));
    cut = max(cut/max(abs(mode.rhoX),[],'all'),1e-16);
    loglog(axisHandle,radialCoordinate,cut,'LineWidth',1.25, ...
        'DisplayName',sprintf('\\theta=%g^\\circ',round(angle*180/pi)));
end
xlabel(axisHandle,'r');
ylabel(axisHandle,'|\rho_x|/max|\rho_x|');
grid(axisHandle,'on');
legend(axisHandle,'Location','southwest','NumColumns',2);
title(axisHandle,'radial cuts of |\rho_x|');

title(layout,sprintf(['mode \\nu=%.3g, (a,b)=(%.3g,%.3g), ' ...
    '\\theta_*=%.1f^\\circ, transport %.1e, edge flux %.1e'],mode.order, ...
    mode.cosineWeight,mode.sineWeight,opts.thetaStar*180/pi, ...
    mode.relativeTransportResidual,mode.relativeEdgeNormalVelocity));
export_figure(figureHandle,opts.detailFile);
end

function polar_field(axisHandle,x,y,value,opts,colorMap)
surface(axisHandle,x,y,zeros(size(value)),value,'EdgeColor','none');
view(axisHandle,2);
axis(axisHandle,'equal','tight');
xlabel(axisHandle,'x_1');
ylabel(axisHandle,'x_2');
colorbar(axisHandle);
limit = robust_limit(value,opts.colorQuantile);
clim(axisHandle,[-limit,limit]);
colormap(axisHandle,colorMap);
fan_edges(axisHandle,opts);
end

function scalar_field(axisHandle,x,y,value,opts,colorMap,colorLimits)
surface(axisHandle,x,y,zeros(size(value)),value,'EdgeColor','none');
view(axisHandle,2);
axis(axisHandle,'equal','tight');
xlabel(axisHandle,'x_1');
ylabel(axisHandle,'x_2');
colorbar(axisHandle);
clim(axisHandle,colorLimits);
colormap(axisHandle,colorMap);
fan_edges(axisHandle,opts);
end

function fan_edges(axisHandle,opts)
hold(axisHandle,'on');
radii = [opts.innerRadius,opts.outerRadius];
for angle = [0,opts.thetaStar]
    plot(axisHandle,opts.center(1)+radii*cos(angle), ...
        opts.center(2)+radii*sin(angle),'k-','LineWidth',0.7);
end
end

function transformed = normalized_log_magnitude(value,numberOfDecades)
relative = abs(value)/max(abs(value),[],'all');
transformed = log10(max(relative,10^(-numberOfDecades)));
end

function transformed = signed_log_magnitude(value,numberOfDecades)
relative = abs(value)/max(abs(value),[],'all');
visibleMagnitude = max(0,numberOfDecades+ ...
    log10(max(relative,10^(-numberOfDecades))));
transformed = sign(value).*visibleMagnitude;
end

function limit = robust_limit(value,quantileLevel)
samples = sort(abs(value(isfinite(value))));
index = max(1,min(numel(samples),ceil(quantileLevel*numel(samples))));
limit = samples(index);
if ~(isfinite(limit) && limit > 0)
    limit = 1;
end
end

function map = blue_white_red(numberOfColors)
half = ceil(numberOfColors/2);
blue = [linspace(0.08,1,half)',linspace(0.25,1,half)',ones(half,1)];
red = [ones(half,1),linspace(1,0.12,half)',linspace(1,0.08,half)'];
map = [blue;red(2:end,:)];
indices = round(linspace(1,size(map,1),numberOfColors));
map = map(indices,:);
end

function export_figure(figureHandle,fileName)
if isempty(fileName)
    return
end
folder = fileparts(fileName);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end
exportgraphics(figureHandle,fileName,'Resolution',210);
end
