function result = ipm_plot_compact_acute_tail(userOpts)
%IPM_PLOT_COMPACT_ACUTE_TAIL Audit a right-acute beta=0 far-field seed.
%   RESULT = IPM_PLOT_COMPACT_ACUTE_TAIL(OPTS) constructs smooth angular
%   functions supported inside a prescribed right acute sector such that
%
%       G0''+G0 = sin(theta)*F0',
%       integral sin(theta)^2*F0' dtheta = 0,
%       integral F0' dtheta = 0.
%
%   Consequently R=F0(theta)+O(1/r) and P=r*G0(theta)+O(1) satisfy the
%   leading beta=0 Dirichlet-Poisson/Fredholm conditions, while
%
%       R_X = -sin(theta)*F0'(theta)/r + O(1/r^2)
%
%   is active only in the requested acute sector and tends to zero in all
%   directions.  This function visualizes and audits only the leading far
%   tail.  It does not solve the nonlinear transport equation and must not
%   be classified as a fixed profile or a blow-up solution.
%
%   The equal-constant return condition is imposed with two disjoint Cinf
%   bumps.  If phi_1 and phi_2 are the bumps, then
%
%       G0 = phi_1-lambda*phi_2,
%       lambda = integral phi_1*csc(theta)^3 dtheta ...
%              / integral phi_2*csc(theta)^3 dtheta.
%
%   This uses the exact identity
%
%       integral (G0''+G0)/sin(theta) dtheta
%       = 2*integral G0*csc(theta)^3 dtheta.

if nargin < 1 || isempty(userOpts)
    userOpts = struct();
end
opts = tail_options(userOpts);

theta = linspace(0,pi,opts.numberOfAngles);
sectorWidth = opts.thetaInterval(2)-opts.thetaInterval(1);
firstInterval = opts.thetaInterval(1)+sectorWidth*[0.06,0.44];
secondInterval = opts.thetaInterval(1)+sectorWidth*[0.56,0.94];

[firstBump,firstBumpSecond] = compact_bump( ...
    theta,firstInterval(1),firstInterval(2));
[secondBump,secondBumpSecond] = compact_bump( ...
    theta,secondInterval(1),secondInterval(2));

safe = sin(theta) > sqrt(eps);
cubicCosecant = zeros(size(theta));
cubicCosecant(safe) = sin(theta(safe)).^(-3);
firstMoment = trapz(theta,firstBump.*cubicCosecant);
secondMoment = trapz(theta,secondBump.*cubicCosecant);
lambda = firstMoment/secondMoment;

rawG0 = firstBump-lambda*secondBump;
rawG0Second = firstBumpSecond-lambda*secondBumpSecond;
rawF0Prime = zeros(size(theta));
rawF0Prime(safe) = ...
    (rawG0Second(safe)+rawG0(safe))./sin(theta(safe));
rawRhoXCoefficient = -sin(theta).*rawF0Prime;
normalization = opts.amplitude/max(abs(rawRhoXCoefficient));
G0 = normalization*rawG0;
G0Second = normalization*rawG0Second;
F0Prime = normalization*rawF0Prime;
F0 = cumtrapz(theta,F0Prime);
F0 = F0-F0(1);
rhoXCoefficient = -sin(theta).*F0Prime;

poissonResidual = G0Second+G0-sin(theta).*F0Prime;
fredholmMoment = trapz(theta,sin(theta).^2.*F0Prime);
returnMoment = trapz(theta,F0Prime);
weightedReturnMoment = 2*trapz(theta,G0.*cubicCosecant);
outside = theta <= opts.thetaInterval(1) | ...
    theta >= opts.thetaInterval(2);
outsideActivity = max(abs(rhoXCoefficient(outside)));

diagnostics = struct();
diagnostics.lambda = lambda;
diagnostics.normalization = normalization;
diagnostics.poissonResidualInfinityNorm = max(abs(poissonResidual));
diagnostics.fredholmMoment = fredholmMoment;
diagnostics.returnMoment = returnMoment;
diagnostics.weightedReturnMoment = weightedReturnMoment;
diagnostics.outsideActivity = outsideActivity;
diagnostics.leftConstant = F0(1);
diagnostics.rightConstant = F0(end);
diagnostics.maximumRhoXCoefficient = max(abs(rhoXCoefficient));
diagnostics.isLeadingTailCompatible = ...
    diagnostics.poissonResidualInfinityNorm <= opts.auditTolerance && ...
    abs(fredholmMoment) <= opts.auditTolerance && ...
    abs(returnMoment) <= opts.auditTolerance && ...
    outsideActivity <= opts.auditTolerance;

result = struct();
result.theta = theta;
result.F0 = F0;
result.F0Prime = F0Prime;
result.G0 = G0;
result.G0Second = G0Second;
result.rhoXCoefficient = rhoXCoefficient;
result.firstBump = firstBump;
result.secondBump = secondBump;
result.diagnostics = diagnostics;
result.options = opts;

if opts.verbose
    fprintf(['compact acute tail: poisson %.3e, Fredholm %.3e, ', ...
        'return %.3e, outside %.3e, lambda %.6g, compatible=%d\n'], ...
        diagnostics.poissonResidualInfinityNorm,fredholmMoment, ...
        returnMoment,outsideActivity,lambda, ...
        diagnostics.isLeadingTailCompatible);
end
if opts.makePlot
    make_tail_figure(result);
end
end

function opts = tail_options(userOpts)
opts = struct();
opts.thetaInterval = pi/180*[8,38];
opts.numberOfAngles = 8001;
opts.numberOfRadii = 240;
opts.radialInterval = [0.75,20];
opts.amplitude = 1;
opts.auditTolerance = 2e-8;
opts.makePlot = true;
opts.plotFile = fullfile('result','verification', ...
    'compact_acute_beta0_tail.png');
opts.verbose = true;

names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end

validateattributes(opts.thetaInterval,{'numeric'}, ...
    {'vector','numel',2,'real','finite','increasing','positive'});
if opts.thetaInterval(2) >= pi/2
    error('ipm:AcuteTailAngle', ...
        'thetaInterval must lie strictly inside the right acute half-sector.');
end
validateattributes(opts.numberOfAngles,{'numeric'}, ...
    {'scalar','integer','>=',1001,'odd'});
validateattributes(opts.numberOfRadii,{'numeric'}, ...
    {'scalar','integer','>=',40});
validateattributes(opts.radialInterval,{'numeric'}, ...
    {'vector','numel',2,'real','finite','increasing','positive'});
validateattributes(opts.amplitude,{'numeric'}, ...
    {'scalar','real','finite','nonzero'});
validateattributes(opts.auditTolerance,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(opts.makePlot,{'logical','numeric'},{'scalar'});
validateattributes(opts.verbose,{'logical','numeric'},{'scalar'});
end

function [value,secondDerivative] = compact_bump(theta,left,right)
value = zeros(size(theta));
secondDerivative = zeros(size(theta));
inside = theta > left & theta < right;
scaled = (theta(inside)-left)/(right-left);
u = scaled.*(1-scaled);
logValue = 4-1./u;
localValue = exp(logValue);
logFirst = (1-2*scaled)./u.^2;
logSecond = -2./u.^2-2*(1-2*scaled).^2./u.^3;

value(inside) = localValue;
secondDerivative(inside) = localValue.* ...
    (logFirst.^2+logSecond)/(right-left)^2;
end

function make_tail_figure(result)
theta = result.theta;
opts = result.options;
diagnostics = result.diagnostics;
radius = linspace(opts.radialInterval(1), ...
    opts.radialInterval(2),opts.numberOfRadii).';
[angleGrid,radiusGrid] = meshgrid(theta,radius);
xGrid = radiusGrid.*cos(angleGrid);
yGrid = radiusGrid.*sin(angleGrid);
rhoX = repmat(result.rhoXCoefficient,opts.numberOfRadii,1)./radiusGrid;

figureHandle = figure('Color','w','Position',[70,70,1420,820]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile(layout);
plot(theta*180/pi,result.F0,'LineWidth',1.8);
hold on;
plot(theta*180/pi,result.G0,'LineWidth',1.8);
xline(opts.thetaInterval(1)*180/pi,'k--');
xline(opts.thetaInterval(2)*180/pi,'k--');
hold off;
xlim([0,90]);
xlabel('$\theta$ (degrees)','Interpreter','latex');
legend('$F_0$','$G_0$','Interpreter','latex','Location','best');
title('Direction-constant density and compact stream tail');
grid on;

nexttile(layout);
plot(theta*180/pi,result.F0Prime,'LineWidth',1.7);
hold on;
plot(theta*180/pi,result.rhoXCoefficient,'LineWidth',1.7);
hold off;
xlim([0,90]);
xlabel('$\theta$ (degrees)','Interpreter','latex');
legend('$F_0''$','$rR_X=-\sin\theta F_0''$', ...
    'Interpreter','latex','Location','best');
title('Signed return and horizontal-gradient coefficient');
grid on;

nexttile(layout);
surface(xGrid,yGrid,zeros(size(rhoX)),rhoX,'EdgeColor','none');
view(2);
axis equal tight;
xlim([0,opts.radialInterval(2)]);
ylim([0,opts.radialInterval(2)*sin(opts.thetaInterval(2))*1.2]);
xlabel('$X$','Interpreter','latex');
ylabel('$Y$','Interpreter','latex');
colorLimit = max(abs(rhoX),[],'all');
clim(colorLimit*[-1,1]);
colorbar;
title('$R_X=-\sin\theta F_0''/r$ (leading tail)', ...
    'Interpreter','latex');

nexttile(layout);
plot(theta*180/pi,result.G0Second+result.G0,'LineWidth',1.8);
hold on;
plot(theta*180/pi,sin(theta).*result.F0Prime,'--','LineWidth',1.6);
hold off;
xlim([0,90]);
xlabel('$\theta$ (degrees)','Interpreter','latex');
legend('$G_0''''+G_0$','$\sin\theta F_0''$', ...
    'Interpreter','latex','Location','best');
title(sprintf(['Poisson %.1e, Fredholm %.1e, return %.1e, ', ...
    'outside %.1e'],diagnostics.poissonResidualInfinityNorm, ...
    diagnostics.fredholmMoment,diagnostics.returnMoment, ...
    diagnostics.outsideActivity));
grid on;

title(layout,['Leading compact-acute beta=0 tail (not a nonlinear ', ...
    'fixed profile)'],'FontWeight','bold');

if ~isempty(opts.plotFile)
    folder = fileparts(opts.plotFile);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    exportgraphics(figureHandle,opts.plotFile,'Resolution',190);
end
end
