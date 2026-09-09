function result = ipm_mesh_exact_fan_profile(userOpts)
%IPM_MESH_EXACT_FAN_PROFILE Mesh the exact local IPM fan steady state.
%   The mesh is built in polar coordinates about (1,0). Upper and lower
%   branches are drawn separately so the accepted density jump at x_2=0 is
%   not interpolated away. The default angle is 60 degrees.

if nargin < 1 || isempty(userOpts)
    userOpts = struct();
end
opts = mesh_options(userOpts);
sineAngle = sin(opts.thetaStar);
cosineAngle = cos(opts.thetaStar);
a0 = opts.amplitude*sineAngle^2/(2*cosineAngle);

radialCoordinate = opts.radius*linspace(0,1,opts.radialPoints) ...
    .^opts.radialClustering;
angularCoordinate = linspace(0,opts.thetaStar,opts.angularPoints)';
[radius,theta] = meshgrid(radialCoordinate,angularCoordinate);
centeredX = radius.*cos(theta);
upperY = radius.*sin(theta);
q = sineAngle*centeredX-cosineAngle*upperY;

upperDensity = opts.amplitude*q;
upperPsi = opts.amplitude*sineAngle/(2*cosineAngle)*upperY.*q;
lowerDensity = -upperDensity;
lowerPsi = -upperPsi;
x = opts.center(1)+centeredX;
yUpper = opts.center(2)+upperY;
yLower = opts.center(2)-upperY;

upperU1 = -a0*centeredX+opts.amplitude*sineAngle*upperY;
upperU2 = a0*upperY;
lowerU1 = upperU1;
lowerU2 = -upperU2;
upperV1 = upperU1+a0*centeredX;
upperV2 = upperU2+a0*upperY;
lowerV1 = lowerU1+a0*centeredX;
lowerV2 = lowerU2-a0*upperY;

densityX = opts.amplitude*sineAngle;
densityY = -opts.amplitude*cosineAngle;
steadyResidual = upperV1*densityX+upperV2*densityY;
physicalResidual = a0*upperDensity+ ...
    upperU1*densityX+upperU2*densityY;
fanNormalResidual = sineAngle*upperU1(end,:)- ...
    cosineAngle*upperU2(end,:);

result = struct();
result.center = opts.center;
result.thetaStar = opts.thetaStar;
result.amplitude = opts.amplitude;
result.c_l = a0;
result.c_omega = 0;
result.c_rCentered = 0;
result.c_rOriginal = -a0;
result.x = x;
result.yUpper = yUpper;
result.yLower = yLower;
result.rhoUpper = upperDensity;
result.rhoLower = lowerDensity;
result.psiUpper = upperPsi;
result.psiLower = lowerPsi;
result.V1Upper = upperV1;
result.V2Upper = upperV2;
result.V1Lower = lowerV1;
result.V2Lower = lowerV2;
result.maximumSteadyResidual = max(abs(steadyResidual),[],'all');
result.maximumPhysicalResidual = max(abs(physicalResidual),[],'all');
result.maximumFanNormalVelocity = max(abs(fanNormalResidual),[],'all');
result.maximumWallNormalVelocity = max(abs(upperU2(1,:)),[],'all');
result.poissonResidual = 0;
result.divergenceResidual = 0;
result.options = opts;

if opts.makePlot
    make_mesh_figure(result,opts);
end
end

function opts = mesh_options(userOpts)
opts = struct();
opts.center = [1,0];
opts.thetaStar = pi/3;
opts.amplitude = 1;
opts.radius = 2.5;
opts.radialPoints = 181;
opts.angularPoints = 121;
opts.radialClustering = 1.35;
opts.makePlot = true;
opts.plotFile = fullfile('result','verification', ...
    'exact_fan_local_mesh.png');
names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
validateattributes(opts.center,{'numeric'}, ...
    {'vector','numel',2,'real','finite'});
validateattributes(opts.thetaStar,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<',pi/2});
validateattributes(opts.amplitude,{'numeric'}, ...
    {'scalar','real','finite','nonzero'});
validateattributes(opts.radius,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(opts.radialPoints,{'numeric'}, ...
    {'scalar','integer','>=',21});
validateattributes(opts.angularPoints,{'numeric'}, ...
    {'scalar','integer','>=',21});
validateattributes(opts.radialClustering,{'numeric'}, ...
    {'scalar','real','finite','>=',1});
opts.center = opts.center(:)';
end

function make_mesh_figure(result,opts)
figureHandle = figure('Color','w','Position',[80,120,1500,560]);
layout = tiledlayout(figureHandle,1,3,'TileSpacing','compact', ...
    'Padding','compact');
plotRadial = unique(round(linspace(1,opts.radialPoints,55)));
plotAngular = unique(round(linspace(1,opts.angularPoints,39)));

nexttile(layout);
hold on
mesh(result.x(plotAngular,plotRadial), ...
    result.yUpper(plotAngular,plotRadial), ...
    result.rhoUpper(plotAngular,plotRadial), ...
    'EdgeColor',[0.78,0.20,0.08],'FaceColor','none');
mesh(result.x(plotAngular,plotRadial), ...
    result.yLower(plotAngular,plotRadial), ...
    result.rhoLower(plotAngular,plotRadial), ...
    'EdgeColor',[0.05,0.30,0.75],'FaceColor','none');
plot3(result.x(1,:),result.yUpper(1,:),result.rhoUpper(1,:), ...
    'Color',[0.78,0.08,0.02],'LineWidth',2.2);
plot3(result.x(1,:),result.yLower(1,:),result.rhoLower(1,:), ...
    'Color',[0.02,0.20,0.78],'LineWidth',2.2);
plot3(opts.center(1),opts.center(2),0,'ko','MarkerFaceColor','k');
xlabel('x_1');
ylabel('x_2');
zlabel('R');
axis tight
view(42,28)
grid on
title('density mesh: wall jump retained');

nexttile(layout);
hold on
mesh(result.x(plotAngular,plotRadial), ...
    result.yUpper(plotAngular,plotRadial), ...
    result.psiUpper(plotAngular,plotRadial), ...
    'EdgeColor',[0.82,0.30,0.08],'FaceColor','none');
mesh(result.x(plotAngular,plotRadial), ...
    result.yLower(plotAngular,plotRadial), ...
    result.psiLower(plotAngular,plotRadial), ...
    'EdgeColor',[0.08,0.38,0.78],'FaceColor','none');
plot3(opts.center(1),opts.center(2),0,'ko','MarkerFaceColor','k');
xlabel('x_1');
ylabel('x_2');
zlabel('\Psi');
axis tight
view(42,28)
grid on
title('streamfunction mesh: continuous at wall');

nexttile(layout);
hold on
contourf(result.x,result.yUpper,result.rhoUpper,18, ...
    'LineStyle','none');
contourf(result.x,result.yLower,result.rhoLower,18, ...
    'LineStyle','none');
radialIndices = unique(round(linspace(8,opts.radialPoints,14)));
angularIndices = unique(round(linspace(1,opts.angularPoints,10)));
quiver(result.x(angularIndices,radialIndices), ...
    result.yUpper(angularIndices,radialIndices), ...
    result.V1Upper(angularIndices,radialIndices), ...
    result.V2Upper(angularIndices,radialIndices),0.6,'k');
quiver(result.x(angularIndices,radialIndices), ...
    result.yLower(angularIndices,radialIndices), ...
    result.V1Lower(angularIndices,radialIndices), ...
    result.V2Lower(angularIndices,radialIndices),0.6,'k');
plot(opts.center(1),opts.center(2),'ko','MarkerFaceColor','k');
axis equal tight
xlabel('x_1');
ylabel('x_2');
colorbar
title('R contours and U+c_l z');

title(layout,sprintf(['exact local fan steady state: theta_*=%.1f deg, ', ...
    'c_l=%.4g, c_omega=0'],opts.thetaStar*180/pi,result.c_l));
if ~isempty(opts.plotFile)
    folder = fileparts(opts.plotFile);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    exportgraphics(figureHandle,opts.plotFile,'Resolution',200);
end
end
