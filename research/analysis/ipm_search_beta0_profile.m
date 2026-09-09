function result = ipm_search_beta0_profile(userOpts)
%IPM_SEARCH_BETA0_PROFILE Search for a bounded beta=0 IPM fixed profile.
%   RESULT = IPM_SEARCH_BETA0_PROFILE(OPTS) performs a finite-dimensional
%   least-squares search for
%
%       (grad^perp(P)+c*[X;Y]-b*e_1).grad(R) = 0,
%       -Delta(P) = R_X,             Y >= 0,
%       P(X,0) = 0.
%
%   The local density jet is fixed to the exact affine wall core
%
%       R = 2*c*k*(X-k*Y) + O(r^2),
%       P = c*X*Y-c*k*Y^2 + O(r^3).
%
%   The trial density vanishes only on the two far-field inflow sides.  No
%   density value is imposed on the top outflow.  Instead, R_X is matched
%   weakly to zero in an outer collar; the X-independent tail value remains
%   free, exactly as required by the weaker condition rho_x -> 0.
%   The total streamfunction is split as
%
%       P = a_infinity*X*Y+p,        p = 0 on the box boundary,
%
%   so the nondecaying harmonic strain needed by a beta=0 return geometry
%   is explicit rather than hidden in a Poisson boundary condition.
%
%   Poisson is eliminated exactly at the discrete level.  Only the
%   transport residual, the beta=0 wall-strain compatibility P_XY(0)=c,
%   the outer R_X matching, and a small coefficient regularizer enter the
%   Gauss-Newton solve.  No density values are prescribed on transport
%   outflow.
%
%   This is a falsifiable box/basis search, not an existence proof.  A
%   candidate must be repeated with larger boxes and richer bases while
%   keeping the bulk/wall residual and the outer R_X collar decreasing.
%
%   Reference noncandidate (2026-08-28): xlim=[-12,12], ymax=6,
%   nx-by-ny=129-by-49, modesX-by-modesY=32-by-14, a_infinity=1.001,
%   strainWeight=50, tailWeight=8, regularization=1e-10.  After 42 accepted
%   iterations the bulk/wall transport RMS values were 2.104e-1/1.044e-1,
%   the P_XY compatibility defect was 1.772e-5, and the outer-collar R_X
%   RMS/maximum values were 9.4e-3/2.424e-1.  The coefficient norm was about
%   1.12e3 and the wall trace remained oscillatory, so this run is explicitly
%   not classified as a profile candidate.  Its saved diagnostic artifacts
%   are result/verification/beta0_profile_search_wide.{mat,png}.

if nargin < 1 || isempty(userOpts)
    userOpts = struct();
end
opts = search_options(userOpts);
ops = search_operators(opts);

[~,originColumn] = min(abs(ops.x));
if abs(ops.x(originColumn)) > 100*eps(max(abs(ops.x)))
    error('ipm:BetaZeroGrid','The x grid must contain X=0.');
end
origin = sub2ind([ops.ny,ops.nx],1,originColumn);

xScale = max(abs(opts.xlim));
yScale = opts.ymax;
xHat = ops.X/xScale;
yHat = ops.Y/yScale;
envelope = max(1-xHat.^2,0).^opts.envelopePower;
affineDensity = 2*opts.c*opts.k*(ops.X-opts.k*ops.Y);
baseDensity = affineDensity.*envelope;
baseDensity = enforce_discrete_jet(baseDensity,envelope,ops,origin, ...
    [0,2*opts.c*opts.k,-2*opts.c*opts.k^2]);

[densityModes,modeLabels] = density_basis( ...
    ops,envelope,xHat,yHat,origin,opts);
numberOfModes = size(densityModes,3);

base = poisson_response(baseDensity,ops,origin);
modes = repmat(empty_response(ops),numberOfModes,1);
for modeIndex = 1:numberOfModes
    modes(modeIndex) = poisson_response( ...
        densityModes(:,:,modeIndex),ops,origin);
end

harmonicP = opts.outerStrain*ops.X.*ops.Y;
harmonicPx = opts.outerStrain*ops.Y;
harmonicPy = opts.outerStrain*ops.X;
harmonicPxy = opts.outerStrain;

weights = residual_weights(ops,opts);
coefficients = zeros(numberOfModes,1);
history = repmat(struct('iteration',0,'objective',0,'bulkRms',0, ...
    'wallRms',0,'tailRms',0,'strainDefect',0,'gradientInfinityNorm',0, ...
    'stepNorm',0,'lambda',0),opts.maxIterations+1,1);
lambda = opts.initialDamping;

[residual,jacobian,state,metrics] = profile_residual(coefficients, ...
    base,modes,harmonicP,harmonicPx,harmonicPy,harmonicPxy,ops, ...
    weights,opts,true);
objective = 0.5*(residual'*residual);
history(1) = history_entry(0,objective,metrics,0,lambda, ...
    norm(jacobian'*residual,inf));
acceptedIterations = 0;

for iteration = 1:opts.maxIterations
    gradient = jacobian'*residual;
    normalMatrix = jacobian'*jacobian;
    diagonalScale = max(full(diag(normalMatrix)),1e-10);
    accepted = false;
    step = zeros(numberOfModes,1);

    for dampingAttempt = 1:opts.maximumDampingAttempts
        step = -(normalMatrix+lambda*spdiags(diagonalScale,0, ...
            numberOfModes,numberOfModes))\gradient;
        trialCoefficients = coefficients+step;
        trialResidual = profile_residual( ...
            trialCoefficients,base,modes,harmonicP,harmonicPx, ...
            harmonicPy,harmonicPxy,ops,weights,opts,false);
        trialObjective = 0.5*(trialResidual'*trialResidual);
        if isfinite(trialObjective) && trialObjective < objective
            accepted = true;
            break;
        end
        lambda = min(lambda*opts.dampingIncrease,opts.maximumDamping);
    end

    if ~accepted
        if opts.verbose
            fprintf(['beta=0 search stopped: no decreasing LM step at ', ...
                'iteration %d (objective %.3e).\n'],iteration,objective);
        end
        break;
    end

    oldObjective = objective;
    coefficients = trialCoefficients;
    lambda = max(lambda/opts.dampingDecrease,opts.minimumDamping);
    [residual,jacobian,state,metrics] = profile_residual(coefficients, ...
        base,modes,harmonicP,harmonicPx,harmonicPy,harmonicPxy,ops, ...
        weights,opts,true);
    objective = 0.5*(residual'*residual);
    acceptedIterations = acceptedIterations+1;
    history(acceptedIterations+1) = history_entry(iteration,objective, ...
        metrics,norm(step),lambda,norm(jacobian'*residual,inf));

    if opts.verbose
        fprintf(['beta=0 iter %2d: J=%.3e, bulk=%.3e, wall=%.3e, ', ...
            'tail=%.3e, strain=%.3e, |dq|=%.3e\n'],iteration,objective, ...
            metrics.bulkRms,metrics.wallRms,metrics.tailRms, ...
            metrics.strainDefect, ...
            norm(step));
    end
    relativeDecrease = (oldObjective-objective)/max(oldObjective,eps);
    if norm(jacobian'*residual,inf) <= opts.gradientTolerance || ...
            (relativeDecrease <= opts.relativeObjectiveTolerance && ...
            norm(step) <= opts.stepTolerance*(1+norm(coefficients)))
        break;
    end
end

history = history(1:acceptedIterations+1);
diagnostics = final_diagnostics(state,ops,origin,opts,metrics);
diagnostics.acceptedIterations = acceptedIterations;
diagnostics.objective = objective;
diagnostics.numberOfModes = numberOfModes;
diagnostics.isBoxCandidate = diagnostics.bulkTransportRms <= ...
    opts.candidateBulkTolerance && diagnostics.wallTransportRms <= ...
    opts.candidateWallTolerance && abs(diagnostics.strainDefect) <= ...
    opts.candidateStrainTolerance && diagnostics.outerCollarRhoX <= ...
    opts.candidateTailTolerance;

result = struct();
result.x = ops.x;
result.y = ops.y;
result.X = ops.X;
result.Y = ops.Y;
result.R = state.R;
result.RX = state.RX;
result.RY = state.RY;
result.P = state.P;
result.pCorrection = state.pCorrection;
result.totalVelocityX = state.VX;
result.totalVelocityY = state.VY;
result.transportResidual = state.transport;
result.b = state.b;
result.coefficients = coefficients;
result.modeLabels = modeLabels;
result.history = history;
result.diagnostics = diagnostics;
result.options = opts;

if opts.makePlot
    make_search_figure(result,opts);
end
if ~isempty(opts.resultFile)
    folder = fileparts(opts.resultFile);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    save(opts.resultFile,'result','-v7.3');
end
end

function opts = search_options(userOpts)
opts = struct();
opts.nx = 65;
opts.ny = 41;
opts.xlim = [-6,6];
opts.ymax = 6;
opts.c = 1;
opts.k = 1;
opts.outerStrain = 1.01;
opts.modesX = 10;
opts.modesY = 8;
opts.envelopePower = 3;
opts.maxIterations = 35;
opts.maximumDampingAttempts = 9;
opts.initialDamping = 1e-2;
opts.minimumDamping = 1e-10;
opts.maximumDamping = 1e12;
opts.dampingIncrease = 10;
opts.dampingDecrease = 3;
opts.gradientTolerance = 1e-8;
opts.stepTolerance = 1e-8;
opts.relativeObjectiveTolerance = 1e-9;
opts.bulkWeight = 1;
opts.wallWeight = 2;
opts.tailWeight = 5;
opts.strainWeight = 40;
opts.regularization = 1e-8;
opts.residualRadialWeightPower = 0;
opts.outerCollarFraction = 0.78;
opts.candidateBulkTolerance = 2e-3;
opts.candidateWallTolerance = 2e-3;
opts.candidateStrainTolerance = 5e-3;
opts.candidateTailTolerance = 2e-3;
opts.makePlot = true;
opts.plotFile = fullfile('result','verification', ...
    'beta0_profile_search.png');
opts.resultFile = fullfile('result','verification', ...
    'beta0_profile_search.mat');
opts.verbose = true;

names = fieldnames(userOpts);
for index = 1:numel(names)
    opts.(names{index}) = userOpts.(names{index});
end
validateattributes(opts.nx,{'numeric'},{'scalar','integer','>=',17,'odd'});
validateattributes(opts.ny,{'numeric'},{'scalar','integer','>=',17});
validateattributes(opts.xlim,{'numeric'}, ...
    {'vector','numel',2,'real','finite','increasing'});
if abs(sum(opts.xlim)) > 100*eps(max(abs(opts.xlim)))
    error('ipm:BetaZeroDomain','xlim must be symmetric about zero.');
end
validateattributes(opts.ymax,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.c,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.k,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.outerStrain,{'numeric'}, ...
    {'scalar','real','finite','>',opts.c});
validateattributes(opts.modesX,{'numeric'},{'scalar','integer','>=',1});
validateattributes(opts.modesY,{'numeric'},{'scalar','integer','>=',1});
validateattributes(opts.envelopePower,{'numeric'}, ...
    {'scalar','integer','>=',2});
validateattributes(opts.maxIterations,{'numeric'}, ...
    {'scalar','integer','>=',0});
opts.makePlot = logical(opts.makePlot);
opts.verbose = logical(opts.verbose);
end

function ops = search_operators(opts)
solverOpts = struct();
solverOpts.nx = opts.nx;
solverOpts.ny = opts.ny;
solverOpts.xlim = opts.xlim;
solverOpts.ymax = opts.ymax;
solverOpts.gridMode = 'uniform';
solverOpts.gridStretchAutomatic = false;
solverOpts.gridStretch = [0,0];
solverOpts.farBoundaryMode = 'dirichlet_zero';
solverOpts.symmetryMode = 'half_plane';
solverOpts.rescalingMode = 'physical';
solverConfig = ipm.config.resolve(solverOpts);
ops = ipm.mesh.build(solverConfig);
end

function [basis,labels] = density_basis( ...
        ops,envelope,xHat,yHat,origin,opts)
radialSquare = xHat.^2+yHat.^2;
numberOfModes = opts.modesX*opts.modesY;
basis = zeros(ops.ny,ops.nx,numberOfModes);
labels = strings(numberOfModes,1);
modeIndex = 0;
for xOrder = 0:opts.modesX-1
    tx = chebyshev_value(xOrder,xHat);
    for yOrder = 0:opts.modesY-1
        ty = chebyshev_value(yOrder,2*yHat-1);
        value = envelope.*radialSquare.*tx.*ty;
        value = enforce_discrete_jet(value,envelope,ops,origin,[0,0,0]);
        scale = max(abs(value),[],'all');
        if scale == 0
            error('ipm:BetaZeroBasis','A trial basis function vanished.');
        end
        modeIndex = modeIndex+1;
        basis(:,:,modeIndex) = value/scale;
        labels(modeIndex) = sprintf('T_%d(X) T_%d(Y)',xOrder,yOrder);
    end
end
end

function value = enforce_discrete_jet(value,envelope,ops,origin,target)
jetValue = envelope;
jetX = ops.X.*envelope;
jetY = ops.Y.*envelope;
valueX = value*ops.Dx';
valueY = ops.Dy*value;
current = [value(origin),valueX(origin),valueY(origin)];

correctors = {jetValue,jetX,jetY};
jetMatrix = zeros(3,3);
for index = 1:3
    corrector = correctors{index};
    correctorX = corrector*ops.Dx';
    correctorY = ops.Dy*corrector;
    jetMatrix(:,index) = [corrector(origin);correctorX(origin); ...
        correctorY(origin)];
end
amplitudes = jetMatrix\(target(:)-current(:));
for index = 1:3
    value = value+amplitudes(index)*correctors{index};
end
end

function value = chebyshev_value(order,coordinate)
if order == 0
    value = ones(size(coordinate));
elseif order == 1
    value = coordinate;
else
    previous = ones(size(coordinate));
    value = coordinate;
    for degree = 2:order
        next = 2*coordinate.*value-previous;
        previous = value;
        value = next;
    end
end
end

function response = empty_response(ops)
zero = zeros(ops.ny,ops.nx);
response = struct('R',zero,'RX',zero,'RY',zero,'p',zero, ...
    'pX',zero,'pY',zero,'pXYOrigin',0,'pYOrigin',0);
end

function response = poisson_response(density,ops,origin)
response = empty_response(ops);
response.R = density;
response.RX = density*ops.Dx';
response.RY = ops.Dy*density;
response.p = ipm.field.poisson(response.RX,ops);
response.pX = response.p*ops.Dx';
response.pY = ops.Dy*response.p;
pXY = response.pY*ops.Dx';
response.pXYOrigin = pXY(origin);
response.pYOrigin = response.pY(origin);
end

function weights = residual_weights(ops,opts)
bulkMask = false(ops.ny,ops.nx);
bulkMask(2:end-1,2:end-1) = true;
bulkQuadrature = ops.weights(bulkMask);
radialScale = max(max(abs(ops.x)),max(ops.y));
radialWeight = (1+(hypot(ops.X(bulkMask),ops.Y(bulkMask))/ ...
    radialScale).^2).^(-opts.residualRadialWeightPower/2);
bulkQuadrature = bulkQuadrature.*radialWeight.^2;
weights.bulkMask = bulkMask;
weights.bulk = opts.bulkWeight*sqrt(bulkQuadrature/sum(bulkQuadrature));
weights.wallColumns = 2:ops.nx-1;
wallQuadrature = ops.hx(weights.wallColumns);
weights.wall = opts.wallWeight*sqrt(wallQuadrature/sum(wallQuadrature));
radius = hypot(ops.X,ops.Y);
outerRadius = opts.outerCollarFraction* ...
    min(max(abs(opts.xlim)),opts.ymax);
weights.tailMask = radius >= outerRadius & ops.Y > 0;
weights.tailMask(:,[1,end]) = false;
tailQuadrature = ops.weights(weights.tailMask);
weights.tail = opts.tailWeight*sqrt(tailQuadrature/sum(tailQuadrature));
end

function [residual,jacobian,state,metrics] = profile_residual( ...
        coefficients,base,modes,harmonicP,harmonicPx,harmonicPy, ...
        harmonicPxy,ops,weights,opts,needJacobian)
numberOfModes = numel(coefficients);
R = base.R;
RX = base.RX;
RY = base.RY;
p = base.p;
pX = base.pX;
pY = base.pY;
pXYOrigin = base.pXYOrigin;
pYOrigin = base.pYOrigin;
for modeIndex = 1:numberOfModes
    coefficient = coefficients(modeIndex);
    R = R+coefficient*modes(modeIndex).R;
    RX = RX+coefficient*modes(modeIndex).RX;
    RY = RY+coefficient*modes(modeIndex).RY;
    p = p+coefficient*modes(modeIndex).p;
    pX = pX+coefficient*modes(modeIndex).pX;
    pY = pY+coefficient*modes(modeIndex).pY;
    pXYOrigin = pXYOrigin+ ...
        coefficient*modes(modeIndex).pXYOrigin;
    pYOrigin = pYOrigin+coefficient*modes(modeIndex).pYOrigin;
end

P = harmonicP+p;
PX = harmonicPx+pX;
PY = harmonicPy+pY;
b = -pYOrigin;
VX = -PY+opts.c*ops.X-b;
VY = PX+opts.c*ops.Y;
transport = VX.*RX+VY.*RY;
strainDefect = harmonicPxy+pXYOrigin-opts.c;

bulkResidual = weights.bulk.*transport(weights.bulkMask);
wallResidual = weights.wall'.* ...
    transport(1,weights.wallColumns)';
tailResidual = weights.tail.*RX(weights.tailMask);
regularizationResidual = sqrt(opts.regularization)*coefficients;
residual = [bulkResidual;wallResidual;tailResidual; ...
    opts.strainWeight*strainDefect;regularizationResidual];

jacobian = [];
if needJacobian
    jacobian = zeros(numel(residual),numberOfModes);
    bulkCount = numel(bulkResidual);
    wallCount = numel(wallResidual);
    tailCount = numel(tailResidual);
    for modeIndex = 1:numberOfModes
        mode = modes(modeIndex);
        velocityModeX = -mode.pY+mode.pYOrigin;
        velocityModeY = mode.pX;
        transportMode = velocityModeX.*RX+velocityModeY.*RY+ ...
            VX.*mode.RX+VY.*mode.RY;
        jacobian(1:bulkCount,modeIndex) = ...
            weights.bulk.*transportMode(weights.bulkMask);
        jacobian(bulkCount+(1:wallCount),modeIndex) = ...
            weights.wall'.*transportMode(1,weights.wallColumns)';
        jacobian(bulkCount+wallCount+(1:tailCount),modeIndex) = ...
            weights.tail.*mode.RX(weights.tailMask);
        jacobian(bulkCount+wallCount+tailCount+1,modeIndex) = ...
            opts.strainWeight*mode.pXYOrigin;
        jacobian(bulkCount+wallCount+tailCount+1+modeIndex,modeIndex) = ...
            sqrt(opts.regularization);
    end
end

metrics = struct();
metrics.bulkRms = sqrt(sum((weights.bulk.* ...
    transport(weights.bulkMask)).^2))/max(opts.bulkWeight,eps);
metrics.wallRms = sqrt(sum((weights.wall'.* ...
    transport(1,weights.wallColumns)').^2))/max(opts.wallWeight,eps);
metrics.tailRms = sqrt(sum((weights.tail.*RX(weights.tailMask)).^2))/ ...
    max(opts.tailWeight,eps);
metrics.strainDefect = strainDefect;

state = struct('R',R,'RX',RX,'RY',RY,'P',P,'PX',PX,'PY',PY, ...
    'pCorrection',p,'VX',VX,'VY',VY,'transport',transport,'b',b);
end

function entry = history_entry(iteration,objective,metrics,stepNorm, ...
        lambda,gradientInfinityNorm)
entry = struct('iteration',iteration,'objective',objective, ...
    'bulkRms',metrics.bulkRms,'wallRms',metrics.wallRms, ...
    'tailRms',metrics.tailRms,'strainDefect',metrics.strainDefect, ...
    'gradientInfinityNorm',gradientInfinityNorm, ...
    'stepNorm',stepNorm,'lambda',lambda);
end

function diagnostics = final_diagnostics(state,ops,origin,opts,metrics)
negativeLaplacian = discrete_negative_laplacian(state.P,ops);
source = state.RX(2:end-1,2:end-1);
poissonResidual = negativeLaplacian-source;
radius = hypot(ops.X,ops.Y);
outerRadius = opts.outerCollarFraction* ...
    min(max(abs(opts.xlim)),opts.ymax);
outerMask = radius >= outerRadius;
outerMask(:,[1,end]) = false;

PYY = ops.Dy*(ops.Dy*state.P);
diagnostics = struct();
diagnostics.bulkTransportRms = metrics.bulkRms;
diagnostics.wallTransportRms = metrics.wallRms;
diagnostics.maximumTransportResidual = max(abs(state.transport),[],'all');
diagnostics.strainDefect = metrics.strainDefect;
diagnostics.maximumPoissonResidual = max(abs(poissonResidual),[],'all');
diagnostics.innerR = state.R(origin);
diagnostics.innerRX = state.RX(origin);
diagnostics.innerRY = state.RY(origin);
diagnostics.innerPXY = opts.c+metrics.strainDefect;
diagnostics.innerPYY = PYY(origin);
diagnostics.targetInnerRX = 2*opts.c*opts.k;
diagnostics.targetInnerRY = -2*opts.c*opts.k^2;
diagnostics.targetInnerPXY = opts.c;
diagnostics.targetInnerPYY = -2*opts.c*opts.k;
diagnostics.outerCollarRhoX = max(abs(state.RX(outerMask)),[],'all');
diagnostics.outerCollarRhoXRms = metrics.tailRms;
diagnostics.inflowSideDensity = max(abs([state.R(:,1)', ...
    state.R(:,end)']),[],'all');
diagnostics.artificialSideRhoX = max(abs([state.RX(:,1)', ...
    state.RX(:,end)']),[],'all');
diagnostics.topDensityVariation = max(state.R(end,:))-min(state.R(end,:));
diagnostics.topRhoX = max(abs(state.RX(end,:)),[],'all');
diagnostics.densityRange = [min(state.R,[],'all'),max(state.R,[],'all')];
diagnostics.maximumRhoX = max(abs(state.RX),[],'all');
diagnostics.pinnedTranslationRate = state.b;
end

function value = discrete_negative_laplacian(P,ops)
value = reshape(ops.A*reshape(P(2:end-1,2:end-1),[],1), ...
    ops.ny-2,ops.nx-2);
value(:,1) = value(:,1)+ ...
    ops.poissonBoundary.x(1)*P(2:end-1,1);
value(:,end) = value(:,end)+ ...
    ops.poissonBoundary.x(2)*P(2:end-1,end);
value(1,:) = value(1,:)+ ...
    ops.poissonBoundary.y(1)*P(1,2:end-1);
value(end,:) = value(end,:)+ ...
    ops.poissonBoundary.y(2)*P(end,2:end-1);
end

function make_search_figure(result,opts)
figureHandle = figure('Color','w','Position',[50,70,1520,850]);
layout = tiledlayout(figureHandle,2,3,'TileSpacing','compact', ...
    'Padding','compact');

nexttile(layout);
contourf(result.X,result.Y,result.R,31,'LineStyle','none');
hold on
contour(result.X,result.Y,result.R,[0,0],'k','LineWidth',1.3);
axis equal tight
colorbar
xlabel('X'); ylabel('Y'); title('R (zero contour in black)');

nexttile(layout);
contourf(result.X,result.Y,result.RX,31,'LineStyle','none');
axis equal tight
colorbar
xlabel('X'); ylabel('Y'); title('R_X');

nexttile(layout);
signedLog = sign(result.transportResidual).* ...
    log10(1+abs(result.transportResidual)/ ...
    max(result.diagnostics.maximumRhoX,eps));
contourf(result.X,result.Y,signedLog,31,'LineStyle','none');
axis equal tight
colorbar
xlabel('X'); ylabel('Y'); title('signed log transport residual');

nexttile(layout);
contourf(result.X,result.Y,result.pCorrection,31,'LineStyle','none');
axis equal tight
colorbar
xlabel('X'); ylabel('Y'); title('localized p=P-a_\inftyXY');

nexttile(layout);
contour(result.X,result.Y,result.R,17,'Color',[0.15,0.35,0.75]);
hold on
strideX = unique(round(linspace(1,numel(result.x),17)));
strideY = unique(round(linspace(1,numel(result.y),11)));
quiver(result.X(strideY,strideX),result.Y(strideY,strideX), ...
    result.totalVelocityX(strideY,strideX), ...
    result.totalVelocityY(strideY,strideX),0.8,'k');
axis equal tight
xlabel('X'); ylabel('Y'); title('R contours and total transport field');

nexttile(layout);
yyaxis left
plot(result.x,result.R(1,:),'LineWidth',1.7);
ylabel('R(X,0)');
yyaxis right
plot(result.x,result.RX(1,:),'LineWidth',1.3);
ylabel('R_X(X,0)');
xlabel('X'); grid on
title('wall trace');

diagnostics = result.diagnostics;
title(layout,sprintf(['beta=0 box search: bulk %.2e, wall %.2e, ', ...
    'strain %.2e, tail R_X %.2e, candidate=%d'], ...
    diagnostics.bulkTransportRms,diagnostics.wallTransportRms, ...
    diagnostics.strainDefect,diagnostics.outerCollarRhoX, ...
    diagnostics.isBoxCandidate));

if ~isempty(opts.plotFile)
    folder = fileparts(opts.plotFile);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    exportgraphics(figureHandle,opts.plotFile,'Resolution',180);
end
end
