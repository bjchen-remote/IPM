function report = elliptic()
%IPMTESTS.BASELINE.ELLIPTIC Verify derivatives, Poisson, velocity, and Green data.

fprintf('Running IPM elliptic verification...\n');
nvals = [17,33,65,129];
errors = zeros(size(nvals));
residuals = zeros(size(nvals));
derivativeXErrors = zeros(size(nvals));
derivativeYErrors = zeros(size(nvals));
velocityErrors = zeros(size(nvals));
initialOmegaErrors = zeros(size(nvals));
ccfTailDerivativeErrors = zeros(size(nvals));
weightedSymmetryDefects = zeros(size(nvals));
for k = 1:numel(nvals)
    overrides = struct('nx',nvals(k),'ny',nvals(k),'xlim',[-1,1], ...
        'ymax',1,'rescalingMode','physical','saveResults',false, ...
        'makePlots',false,'verbose',false, ...
        'gridMode','stretched','gridStretch',3, ...
        'farBoundaryMode','dirichlet_zero');
    config = ipm.config.resolve(overrides);
    ops = ipm.mesh.build(config);
    psiExact = sin(pi*(ops.X+1)/2).*sin(pi*ops.Y);
    lambda = (pi/2)^2+pi^2;
    [psi,manufacturedRhs] = ipm.field.poisson(lambda*psiExact,ops);
    errors(k) = sqrt(sum((psi-psiExact).^2.*ops.weights,'all'));
    q = psi(2:end-1,2:end-1);
    residuals(k) = norm(ops.A*q(:)-manufacturedRhs(:),inf) / ...
        norm(manufacturedRhs(:),inf);
    weightedA = spdiags(ops.poissonWeights,0,numel(ops.poissonWeights), ...
        numel(ops.poissonWeights))*ops.A;
    weightedSymmetryDefects(k) = norm(weightedA-weightedA',inf) / ...
        norm(weightedA,inf);

    fieldXExact = (pi/2)*cos(pi*(ops.X+1)/2).*sin(pi*ops.Y);
    fieldYExact = pi*sin(pi*(ops.X+1)/2).*cos(pi*ops.Y);
    derivativeXErrors(k) = weighted_l2(psiExact*ops.Dx'-fieldXExact,ops);
    derivativeYErrors(k) = weighted_l2(ops.Dy*psiExact-fieldYExact,ops);

    rhoVelocityTest = -(2*lambda/pi)*cos(pi*(ops.X+1)/2).*sin(pi*ops.Y);
    [u1,u2] = ipm.field.velocity(rhoVelocityTest,ops);
    u1Exact = -pi*sin(pi*(ops.X+1)/2).*cos(pi*ops.Y);
    u2Exact = (pi/2)*cos(pi*(ops.X+1)/2).*sin(pi*ops.Y);
    velocityErrors(k) = weighted_l2(hypot(u1-u1Exact,u2-u2Exact),ops);

    absoluteInitialX = abs(ops.X);
    denominator = 1+absoluteInitialX.^4+ops.Y.^4;
    initialRho = absoluteInitialX.^4./denominator;
    initialOmega = 4*sign(ops.X).*absoluteInitialX.^3 .* ...
        (1+ops.Y.^4)./denominator.^2;
    initialOmegaErrors(k) = weighted_l2( ...
        initialRho*ops.Dx'-initialOmega,ops);

    ccfTailConfig = ipm.config.resolve(struct( ...
        'initialCondition','ccf_heavy_tail', ...
        'ccfTailOnset',0.35,'ccfTailRegularization',0.08, ...
        'ccfTailAmplitude',1.2,'ccfTailCutoffScale',2, ...
        'ccfTailVerticalScale',0.7, ...
        'ccfTailVerticalPower',4));
    ccfTailPhysics = ccfTailConfig.physics;
    ccfTailRho = ipm.field.initialDensity(ops,ccfTailPhysics);
    ccfTailZ = max(abs(ops.X)-ccfTailPhysics.ccfTailOnset,0);
    ccfTailEnvelope = 1./(1+ ...
        (ops.Y/ccfTailPhysics.ccfTailVerticalScale).^4);
    ccfTailBaseOmega = sign(ops.X).* ...
        ccfTailPhysics.ccfTailAmplitude.*ccfTailZ ./ ...
        (ccfTailZ.^2+ccfTailPhysics.ccfTailRegularization^2).^(3/4) .* ...
        ccfTailEnvelope;
    ccfTailBasePrimitive = 2*ccfTailPhysics.ccfTailAmplitude * ...
        ((ccfTailZ.^2+ccfTailPhysics.ccfTailRegularization^2).^(1/4) - ...
        sqrt(ccfTailPhysics.ccfTailRegularization));
    ccfTailSaturation = 2*ccfTailPhysics.ccfTailAmplitude * ...
        sqrt(ccfTailPhysics.ccfTailCutoffScale);
    ccfTailOmega = ccfTailBaseOmega ./ ...
        cosh(ccfTailBasePrimitive/ccfTailSaturation).^2;
    % The family is only C1 at z=0; use the same fixed smooth-tail mask.
    ccfTailSmoothMask = ccfTailZ > ...
        4*ccfTailPhysics.ccfTailRegularization & abs(ops.X) < 0.95;
    ccfTailDerivativeError = ccfTailRho*ops.Dx'-ccfTailOmega;
    ccfTailDerivativeErrors(k) = sqrt(sum( ...
        ccfTailDerivativeError(ccfTailSmoothMask).^2 .* ...
        ops.weights(ccfTailSmoothMask),'all'));
end
orders = log(errors(1:end-1)./errors(2:end))/log(2);
derivativeXOrders = log(derivativeXErrors(1:end-1)./ ...
    derivativeXErrors(2:end))/log(2);
derivativeYOrders = log(derivativeYErrors(1:end-1)./ ...
    derivativeYErrors(2:end))/log(2);
velocityOrders = log(velocityErrors(1:end-1)./velocityErrors(2:end))/log(2);
initialOmegaOrders = log(initialOmegaErrors(1:end-1)./ ...
    initialOmegaErrors(2:end))/log(2);
ccfTailDerivativeOrders = log(ccfTailDerivativeErrors(1:end-1)./ ...
    ccfTailDerivativeErrors(2:end))/log(2);

overrides = struct('nx',65,'ny',33,'xlim',[-4,4],'ymax',4, ...
    'rescalingMode','physical','saveResults',false, ...
    'makePlots',false,'verbose',false,'transportBoundaryMode','closed');
config = ipm.config.resolve(overrides);
ops = ipm.mesh.build(config);
rho = exp(-(ops.X.^2+(ops.Y-1.2).^2));
[rhs,flow] = ipm.evolve.rhs(rho,ops);
massRate = sum(rhs.*ops.weights,'all');
divU = flow.u1*ops.Dx'+ops.Dy*flow.u2;

boundaryChange = boundary_expansion_check('green');
zeroBoundaryChange = boundary_expansion_check('dirichlet_zero');
greenBoundaryError = green_boundary_check('double_odd_omega');
halfPlaneGreenBoundaryError = green_boundary_check('half_plane');
[greenCompressionBoundaryError,greenCompressionVelocityError] = ...
    green_compression_check();

report.poissonL2Errors = errors;
report.poissonOrders = orders;
report.poissonAlgebraicResiduals = residuals;
report.poissonWeightedSymmetryDefects = weightedSymmetryDefects;
report.derivativeXErrors = derivativeXErrors;
report.derivativeYErrors = derivativeYErrors;
report.derivativeXOrders = derivativeXOrders;
report.derivativeYOrders = derivativeYOrders;
report.velocityErrors = velocityErrors;
report.velocityOrders = velocityOrders;
report.initialOmegaErrors = initialOmegaErrors;
report.initialOmegaOrders = initialOmegaOrders;
report.ccfTailDerivativeErrors = ccfTailDerivativeErrors;
report.ccfTailDerivativeOrders = ccfTailDerivativeOrders;
report.divergenceInf = max(abs(divU),[],'all');
report.wallNormalInf = max(abs(flow.u2(1,:)),[],'all');
report.massRate = massRate;
report.boundaryExpansionRelativeChange = boundaryChange;
report.zeroBoundaryExpansionRelativeChange = zeroBoundaryChange;
report.greenBoundaryQuadratureError = greenBoundaryError;
report.halfPlaneGreenBoundaryQuadratureError = halfPlaneGreenBoundaryError;
report.greenCompressionBoundaryError = greenCompressionBoundaryError;
report.greenCompressionVelocityError = greenCompressionVelocityError;

assert(all(orders(end-1:end) > 1.9), ...
    'ipm:PoissonOrder','Poisson solver did not attain second-order accuracy.');
assert(max(residuals) < 1e-10,'ipm:PoissonResidual', ...
    'Poisson solve is inaccurate.');
assert(max(weightedSymmetryDefects) < 1e-12,'ipm:PoissonSymmetry', ...
    'Jacobian-weighted Poisson stiffness is not symmetric.');
assert(all(derivativeXOrders(end-1:end) > 1.9) && ...
    all(derivativeYOrders(end-1:end) > 1.9), ...
    'ipm:DerivativeOrder','First derivatives did not attain second-order accuracy.');
assert(all(velocityOrders(end-1:end) > 1.9), ...
    'ipm:VelocityOrder','IPM velocity operator did not attain second-order accuracy.');
assert(initialOmegaOrders(end) > 1.9, ...
    'ipm:InitialDerivativeOrder', ...
    'The rational initial density did not recover omega at second order.');
assert(ccfTailDerivativeOrders(end) > 1.9, ...
    'ipm:CcfTailDerivativeOrder', ...
    'The regularized CCF tail derivative was not second order away from its onset.');
assert(report.divergenceInf < 1e-11,'ipm:Divergence', ...
    'Velocity is not divergence free.');
assert(report.wallNormalInf < 1e-12,'ipm:WallNormalVelocity', ...
    'The wall-normal velocity must vanish on x_2=0.');
assert(max(greenBoundaryError,halfPlaneGreenBoundaryError) < 1e-12, ...
    'ipm:GreenBoundary', ...
    'Quarter-plane or half-plane Green boundary evaluation is inconsistent.');
assert(greenCompressionBoundaryError < 2e-3 && ...
    greenCompressionVelocityError < 1e-4, ...
    'ipm:GreenCompression', ...
    'Conservative Green-source aggregation is insufficiently accurate.');
assert(boundaryChange < 1e-2 && boundaryChange < 0.1*zeroBoundaryChange, ...
    'ipm:FarBoundaryConvergence', ...
    'Green far data did not sufficiently reduce doubled-box sensitivity.');
assert(abs(massRate) < 1e-11,'ipm:Mass', ...
    'Flux operator is not conservative.');

fprintf('  Poisson orders: %s\n',sprintf('%.3f ',orders));
fprintf('  D_x orders: %s\n',sprintf('%.3f ',derivativeXOrders));
fprintf('  D_y orders: %s\n',sprintf('%.3f ',derivativeYOrders));
fprintf('  Biot-Savart velocity orders: %s\n',sprintf('%.3f ',velocityOrders));
fprintf('  initial rho-to-omega orders: %s\n', ...
    sprintf('%.3f ',initialOmegaOrders));
fprintf('  regularized CCF-tail D_x orders away from onset: %s\n', ...
    sprintf('%.3f ',ccfTailDerivativeOrders));
fprintf('  max Poisson residual: %.3e\n',max(residuals));
fprintf('  max weighted-stiffness symmetry defect: %.3e\n', ...
    max(weightedSymmetryDefects));
fprintf('  ||div u||_inf: %.3e\n',report.divergenceInf);
fprintf('  wall ||u_2||_inf: %.3e\n',report.wallNormalInf);
fprintf('  conservative mass rate: %.3e\n',massRate);
fprintf('  core velocity change after doubling box: %.3e\n',boundaryChange);
fprintf('  same change with zero far boundary: %.3e\n',zeroBoundaryChange);
fprintf('  Green boundary quadrature error: %.3e\n',greenBoundaryError);
fprintf('  legacy half-plane Green boundary error: %.3e\n', ...
    halfPlaneGreenBoundaryError);
fprintf('  Green aggregation boundary/core-velocity errors: %.3e, %.3e\n', ...
    greenCompressionBoundaryError,greenCompressionVelocityError);
fprintf('IPM elliptic verification passed.\n');
end

function [boundaryError,velocityError] = green_compression_check()
overrides = struct('nx',129,'ny',65,'xlim',[-16,16],'ymax',16, ...
    'gridMode','stretched','targetCenterSpacing',0.05, ...
    'farBoundaryMode','green','greenMaxSources',3000, ...
    'rescalingMode','physical','saveResults',false,'makePlots',false, ...
    'verbose',false);
config = ipm.config.resolve(overrides);
compressedOps = ipm.mesh.build(config);
fullOps = compressedOps;
fullOps.greenMaxSources = numel(compressedOps.X);
denominator = 1+compressedOps.X.^4+compressedOps.Y.^4;
rho = abs(compressedOps.X).^4./denominator;
source = rho*compressedOps.Dx';
compressedBoundary = ipm.field.greenBoundary(source,compressedOps);
fullBoundary = ipm.field.greenBoundary(source,fullOps);
compressedVector = [compressedBoundary.left;compressedBoundary.right; ...
    compressedBoundary.top'];
fullVector = [fullBoundary.left;fullBoundary.right;fullBoundary.top'];
boundaryError = norm(compressedVector-fullVector,inf) / ...
    max(norm(fullVector,inf),eps);
[compressedU1,compressedU2] = ipm.field.biotSavart(source,compressedOps);
[fullU1,fullU2] = ipm.field.biotSavart(source,fullOps);
core = abs(compressedOps.X) <= 2 & compressedOps.Y <= 2;
velocityError = sqrt(sum(((compressedU1-fullU1).^2 + ...
    (compressedU2-fullU2).^2).*compressedOps.weights.*core,'all')) / ...
    max(sqrt(sum((fullU1.^2+fullU2.^2).*compressedOps.weights.*core, ...
    'all')),eps);
end

function relativeError = green_boundary_check(symmetryMode)
overrides = struct('nx',33,'ny',17,'xlim',[-6,6],'ymax',6, ...
    'gridMode','stretched','gridStretch',3,'farBoundaryMode','green', ...
    'symmetryMode',symmetryMode, ...
    'rescalingMode','physical','saveResults',false,'makePlots',false, ...
    'verbose',false);
config = ipm.config.resolve(overrides);
ops = ipm.mesh.build(config);
source = zeros(ops.ny,ops.nx);
ix = 20;
iy = 7;
source(iy,ix) = 1/ops.weights(iy,ix);
boundary = ipm.field.greenBoundary(source,ops);
sourceX = ops.x(ix);
sourceY = ops.y(iy);
if strcmp(symmetryMode,'double_odd_omega')
    kernel = @quarter_image_green;
else
    kernel = @image_green;
end
leftExact = kernel(ops.x(1)*ones(ops.ny,1),ops.y,sourceX,sourceY);
rightExact = kernel(ops.x(end)*ones(ops.ny,1),ops.y,sourceX,sourceY);
topExact = kernel(ops.x',ops.y(end)*ones(ops.nx,1),sourceX,sourceY)';
errorVector = [boundary.left-leftExact;boundary.right-rightExact; ...
    (boundary.top-topExact)'];
exactVector = [leftExact;rightExact;topExact'];
relativeError = norm(errorVector,inf)/max(norm(exactVector,inf),eps);
end

function value = quarter_image_green(targetX,targetY,sourceX,sourceY)
directSquared = (targetX-sourceX).^2+(targetY-sourceY).^2;
xImageSquared = (targetX+sourceX).^2+(targetY-sourceY).^2;
yImageSquared = (targetX-sourceX).^2+(targetY+sourceY).^2;
xyImageSquared = (targetX+sourceX).^2+(targetY+sourceY).^2;
value = log((xImageSquared.*yImageSquared) ./ ...
    (directSquared.*xyImageSquared))/(4*pi);
end

function value = image_green(targetX,targetY,sourceX,sourceY)
directSquared = (targetX-sourceX).^2+(targetY-sourceY).^2;
imageSquared = (targetX-sourceX).^2+(targetY+sourceY).^2;
value = log(imageSquared./directSquared)/(4*pi);
end

function value = weighted_l2(field,ops)
value = sqrt(sum(field.^2.*ops.weights,'all'));
end

function relativeChange = boundary_expansion_check(farBoundaryMode)
base = struct('rescalingMode','physical','saveResults',false, ...
    'makePlots',false,'verbose',false, ...
    'targetCenterSpacing',0.05, ...
    'symmetryMode','half_plane', ...
    'farBoundaryMode',farBoundaryMode);
small = base;
small.nx = 65; small.ny = 33; small.xlim = [-4,4]; small.ymax = 4;
large = base;
large.nx = 129; large.ny = 65; large.xlim = [-8,8]; large.ymax = 8;
smallConfig = ipm.config.resolve(small);
largeConfig = ipm.config.resolve(large);
os = ipm.mesh.build(smallConfig);
ol = ipm.mesh.build(largeConfig);
rhoS = exp(-(os.X.^2+((os.Y-1.2)/0.7).^2));
rhoL = exp(-(ol.X.^2+((ol.Y-1.2)/0.7).^2));
[us,vs] = ipm.field.velocity(rhoS,os);
[ul,vl] = ipm.field.velocity(rhoL,ol);
ul = interp2(ol.x,ol.y,ul,os.X,os.Y,'linear');
vl = interp2(ol.x,ol.y,vl,os.X,os.Y,'linear');
num = sqrt(sum((us-ul).^2+(vs-vl).^2,'all'));
den = sqrt(sum(ul.^2+vl.^2,'all'));
relativeChange = num/max(den,eps);
end
