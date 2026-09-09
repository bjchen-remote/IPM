function report = core()
%IPMTESTS.FOURTH.CORE Layered fourth-order-and-higher accuracy gates.
%   This laboratory verification is intentionally separate from IPMTESTS.BASELINE.SUITE
%   while the production solver is still evolving.  All reported spatial
%   norms include every node, including the one-sided boundary closures.

threshold = 3.8;
spatialGridSizes = [13,19,29,43];
transportGridSizes = [13,19,25,37];
manufacturedTime = 0.137;

fprintf('Running high-order IPM verification (required order %.1f)...\n', ...
    threshold);
report = struct('threshold',threshold,'manufacturedTime',manufacturedTime);
report.spatial.uniform = verify_spatial_family( ...
    'uniform',spatialGridSizes,manufacturedTime,threshold);
report.spatial.stretched = verify_spatial_family( ...
    'stretched',spatialGridSizes,manufacturedTime,threshold);
report.spatial.nonzeroDirichlet = verify_nonzero_dirichlet( ...
    spatialGridSizes,threshold);
report.ssprk54 = verify_ssprk54_conditions();
report.productionContracts = verify_production_contracts();

assert(~isempty(which('ipm.field.weno5FluxDerivative')), ...
    'ipm:HighOrderTransportUnavailable', ...
    ['The required WENO5-FD implementation is not on the MATLAB path; ' ...
    'the complete high-order certificate cannot skip transport.']);
report.transport.status = 'passed';
report.transport.uniform = verify_transport_family( ...
    'uniform',transportGridSizes,manufacturedTime,threshold);
report.transport.stretched = verify_transport_family( ...
    'stretched',transportGridSizes,manufacturedTime,threshold);
report.transport.mappedCriticalPoint = ...
    verify_mapped_critical_point(threshold);
report.transport.global = verify_forced_evolution( ...
    transportGridSizes,threshold);

function result = verify_mapped_critical_point(threshold)
% A varying metric can create a high-order critical zero in one LF split.
% This regression exposed a genuine second-order failure in the first
% WENO-Z prototype even though the coupled MMS happened to converge faster.
gridSizes = [33,65,129,257];
errors = zeros(2,numel(gridSizes));
h = 1./(gridSizes-1);
for level = 1:numel(gridSizes)
    options = struct('nx',gridSizes(level),'ny',9, ...
        'xlim',[-1,1],'ymax',1,'gridMode','stretched', ...
        'gridStretchAutomatic',false,'gridStretch',[1.2,0.5], ...
        'spatialDiscretization','high_order', ...
        'farBoundaryMode','dirichlet_zero','symmetryMode','half_plane', ...
        'transportBoundaryMode','open','transportScheme','weno5_fd', ...
        'adaptiveRemesh',false,'saveResults',false,'makePlots',false, ...
        'livePlot',false,'writeVideo',false,'verbose',false);
    ops = ipm.mesh.build(ipm.config.resolve(options));
    rho = sin(2*ops.X);
    numerical = ipm.field.transport(rho,ones(size(rho)), ...
        zeros(size(rho)),ops,'open');
    error = numerical+2*cos(2*ops.X);
    errors(1,level) = full_l2(error,ops);
    errors(2,level) = max(abs(error),[],'all');
end
orders = observed_orders(errors,h);
finalOrders = orders(:,end);
tailOrders = certify_convergence_tail(errors,orders,threshold, ...
    'ipm:HighOrderMappedCriticalPoint','mapped WENO critical point');
result = struct('gridSizes',gridSizes,'h',h,'errors',errors, ...
    'orders',orders,'tailOrders',tailOrders,'finalOrders',finalOrders);
fprintf('  mapped LF critical-point final orders: L2 %.3f, Linf %.3f\n', ...
    finalOrders(1),finalOrders(2));
end
report.passed = true;
fprintf('High-order IPM verification passed.\n');
end

function result = verify_nonzero_dirichlet(gridSizes,threshold)
% Exercise the exact production boundary-RHS orientation and metric factors.
modes = {'uniform','stretched'};
shapes = {'square','rectangular'};
kappas = [0.5,1,2];
metricNames = {'poissonL2','poissonLinf'};
result = struct('gridSizes',gridSizes,'kappas',kappas, ...
    'metricNames',{metricNames},'shapes',{shapes});
minimumFinalOrder = Inf;
for modeIndex = 1:numel(modes)
    mode = modes{modeIndex};
    modeMinimumOrder = Inf;
    for shapeIndex = 1:numel(shapes)
        shape = shapes{shapeIndex};
        errors = zeros(numel(metricNames),numel(gridSizes),numel(kappas));
        h = 1./(gridSizes-1);
        for kappaIndex = 1:numel(kappas)
            kappa = kappas(kappaIndex);
            for level = 1:numel(gridSizes)
                if strcmp(shape,'square')
                    ops = manufactured_operators(gridSizes(level),mode);
                else
                    ops = manufactured_rectangular_operators( ...
                        gridSizes(level),mode);
                end
                [psiExact,source] = nonzero_dirichlet_solution( ...
                    ops.X,ops.Y,kappa);
                boundary = struct('left',psiExact(:,1), ...
                    'right',psiExact(:,end),'bottom',psiExact(1,:), ...
                    'top',psiExact(end,:));
                psi = ipm.field.poisson(source,ops,kappa,boundary);
                error = psi-psiExact;
                errors(1,level,kappaIndex) = full_l2(error,ops);
                errors(2,level,kappaIndex) = max(abs(error),[],'all');
            end
        end
        orders = zeros(numel(metricNames),numel(gridSizes)-1,numel(kappas));
        for kappaIndex = 1:numel(kappas)
            orders(:,:,kappaIndex) = observed_orders( ...
                errors(:,:,kappaIndex),h);
        end
        finalOrders = squeeze(orders(:,end,:));
        tailOrders = certify_convergence_tail(errors,orders,threshold, ...
            'ipm:HighOrderNonzeroDirichlet', ...
            sprintf('%s %s nonzero-Dirichlet metric Poisson',mode,shape));
        shapeMinimumOrder = min(finalOrders,[],'all');
        modeMinimumOrder = min(modeMinimumOrder,shapeMinimumOrder);
        minimumFinalOrder = min(minimumFinalOrder,shapeMinimumOrder);
        result.(mode).(shape) = struct('h',h,'errors',errors, ...
            'orders',orders,'tailOrders',tailOrders, ...
            'finalOrders',finalOrders);
    end
    fprintf(['  %-9s nonzero-Dirichlet Poisson minimum final order: ' ...
        '%.3f (square/rectangular, kappa 0.5/1/2).\n'], ...
        mode,modeMinimumOrder);
end
result.minimumFinalOrder = minimumFinalOrder;
end

function [psi,source] = nonzero_dirichlet_solution(X,Y,kappa)
exponential = exp(0.3*X+0.2*Y);
psi = Y.*exponential+0.2*Y.^2.*sin(1.1*X);
psiXX = 0.09*Y.*exponential-0.2*1.1^2*Y.^2.*sin(1.1*X);
psiYY = (0.4+0.04*Y).*exponential+0.4*sin(1.1*X);
source = -(1/kappa)*psiXX-kappa*psiYY;
end

function result = verify_production_contracts()
% Cover the public assembler and the actual coupled production step, not
% only their isolated WENO/tableau ingredients.
ops = manufactured_operators(17,'stretched');
rhoConstant = ones(ops.ny,ops.nx);
zeroVelocity = zeros(size(rhoConstant));
[constantRhs,constantBase,transportU1,transportU2,sourceCoefficient] = ...
    ipm.evolve.assembleRhs(rhoConstant,zeroVelocity,zeroVelocity, ...
        0.3,-0.2,0.4,0.1,ops,'open');
freeStreamError = max([max(abs(constantBase),[],'all'), ...
    max(abs(constantRhs-0.4),[],'all'), ...
    max(abs(transportU1-(0.3*ops.X+0.1)),[],'all'), ...
    max(abs(transportU2+0.2*ops.Y),[],'all'), ...
    abs(sourceCoefficient-0.5)]);
assert(freeStreamError < 2e-12, ...
    'ipm:HighOrderRescaledFreeStream', ...
    'The anisotropic rescaled WENO free-stream/source contract failed.');

traceOps = ops;
traceOps.wallTransportMode = 'advective_upwind';
traceRho = 1+0.2*sin(1.1*ops.X).*exp(-ops.Y);
traceU1 = 0.4+0.1*cos(ops.X);
traceU2 = zeros(size(traceRho));
[openTraceRhs,traceTransportForm] = ipm.field.transport( ...
    traceRho,traceU1,traceU2,traceOps,'open');
closedTraceRhs = ipm.field.transport( ...
    traceRho,traceU1,traceU2,traceOps,'closed');
closedWallTraceError = max(abs( ...
    closedTraceRhs(1,:)-openTraceRhs(1,:)),[],'all');
assert(closedWallTraceError == 0 && ...
    ~traceTransportForm.bulkConservative && ...
    traceTransportForm.advectiveWall, ...
    'ipm:HighOrderAdvectiveWallProjection', ...
    ['Closed-boundary mass projection changed the explicitly selected ' ...
    'advective wall trace.']);

squareNodes = 17;
coordinate = linspace(-1,1,squareNodes);
squareState = repmat(sin(1.3*coordinate),squareNodes,1);
lineVelocity = linspace(0.4,1.2,squareNodes)';
lineMetric = linspace(0.8,1.1,squareNodes)';
columnBroadcast = ipm.field.weno5FluxDerivative(squareState, ...
    lineVelocity,lineMetric,coordinate(2)-coordinate(1),struct());
explicitBroadcast = ipm.field.weno5FluxDerivative(squareState, ...
    repmat(lineVelocity,1,squareNodes), ...
    repmat(lineMetric,1,squareNodes),coordinate(2)-coordinate(1),struct());
squareBroadcastError = max(abs(columnBroadcast-explicitBroadcast),[],'all');
assert(squareBroadcastError == 0, ...
    'ipm:HighOrderWenoSquareBroadcast', ...
    ['A square-grid nLines-by-1 velocity or metric did not retain its ' ...
    'documented line-constant broadcasting meaning.']);

userOptions = struct('nx',25,'ny',17,'xlim',[-2,2],'ymax',2, ...
    'gridMode','stretched','gridStretchAutomatic',false, ...
    'gridStretch',[0.7,0.5], ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','farBoundaryMode','dirichlet_zero', ...
    'transportBoundaryMode','open','symmetryMode','half_plane', ...
    'rescalingMode','dynamic','dynamicScaleGeometry','anisotropic', ...
    'anisotropicGaugeMode','fixed','anisotropicFixedCX',0.35, ...
    'anisotropicFixedCY',-0.2,'anisotropicFixedCOmega',0.1, ...
    'anisotropicFixedCR',0.05, ...
    'initialCondition',@(X,Y)2+exp(-((X-0.55).^2+0.7*Y.^2)), ...
    'adaptiveRemesh',false,'saveResults',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
config = ipm.config.resolve(userOptions);
stepOps = ipm.mesh.build(config);
rho = ipm.field.initialDensity(stepOps,config.physics);
stepOps = ipm.evolve.initializeScaling(rho,stepOps);
scale = struct('logC_l',0,'logC_x',0,'logC_y',0, ...
    'logC_omega',0,'physicalTime',0,'X_shift',0,'canonicalTime',0);
dt = 0.0125;
[rhoActual,flowActual,scaleActual] = ...
    ipm.evolve.stepSsprk54(rho,dt,stepOps,scale);
[rhoOracle,flowOracle,scaleOracle] = ...
    ssprk54_production_oracle(rho,dt,stepOps,scale);
stepError = max([max(abs(rhoActual-rhoOracle),[],'all'), ...
    max(abs(pack_anisotropic_scale(scaleActual)- ...
        pack_anisotropic_scale(scaleOracle))), ...
    abs(flowActual.aspect-flowOracle.aspect), ...
    abs(flowActual.poissonSolveInfo.kappa- ...
        flowOracle.poissonSolveInfo.kappa)]);
expectedAspect = exp(scaleActual.logC_y-scaleActual.logC_l);
aspectError = max([abs(flowActual.aspect-expectedAspect), ...
    abs(flowActual.poissonSolveInfo.kappa-expectedAspect)]);
assert(stepError < 5e-12 && aspectError < 5e-13, ...
    'ipm:HighOrderProductionSSPRK54', ...
    ['The production SSPRK54 packing, stage metric, or final flow ' ...
    'disagrees with the independent tableau oracle.']);

isotropicOptions = userOptions;
isotropicOptions.dynamicScaleGeometry = 'isotropic';
isotropicOptions.lengthGauge = 'local_strain';
isotropicOptions.cOmegaGauge = 'none';
isotropicConfig = ipm.config.resolve(isotropicOptions);
isotropicOps = ipm.mesh.build(isotropicConfig);
isotropicRho = ipm.field.initialDensity(isotropicOps,isotropicConfig.physics);
isotropicOps = ipm.evolve.initializeScaling(isotropicRho,isotropicOps);
isotropicScale = struct('logC_l',0,'logC_omega',0, ...
    'physicalTime',0,'X_shift',0,'canonicalTime',0);
[isotropicActual,isotropicFlowActual,isotropicScaleActual] = ...
    ipm.evolve.stepSsprk54(isotropicRho,dt,isotropicOps,isotropicScale);
[isotropicOracle,isotropicFlowOracle,isotropicScaleOracle] = ...
    ssprk54_isotropic_production_oracle( ...
        isotropicRho,dt,isotropicOps,isotropicScale);
isotropicStepError = max([ ...
    max(abs(isotropicActual-isotropicOracle),[],'all'), ...
    max(abs(pack_isotropic_scale(isotropicScaleActual)- ...
        pack_isotropic_scale(isotropicScaleOracle))), ...
    abs(isotropicFlowActual.aspect-isotropicFlowOracle.aspect)]);
isotropicRateMagnitude = max(abs([isotropicFlowActual.c_l, ...
    isotropicFlowActual.c_omega,isotropicFlowActual.c_r]));
assert(isotropicStepError < 5e-12 && isotropicRateMagnitude > 1e-6 && ...
    isotropicFlowActual.aspect == 1, ...
    'ipm:HighOrderProductionSSPRK54Isotropic', ...
    ['The nonzero-rate isotropic production SSPRK54 branch disagrees ' ...
    'with its independent tableau oracle.']);
result = struct('freeStreamError',freeStreamError, ...
    'closedWallTraceError',closedWallTraceError, ...
    'traceTransportForm',traceTransportForm, ...
    'squareBroadcastError',squareBroadcastError, ...
    'productionStepError',stepError,'stageAspectError',aspectError, ...
    'isotropicProductionStepError',isotropicStepError, ...
    'isotropicRateMagnitude',isotropicRateMagnitude, ...
    'dt',dt,'rates',[0.35,-0.2,0.1,0.05]);
fprintf(['  production free-stream / anisotropic SSPRK54 / aspect / ' ...
    'isotropic SSPRK54 errors: %.3e / %.3e / %.3e / %.3e.\n'], ...
    freeStreamError,stepError,aspectError,isotropicStepError);
end

function [rhoNew,flow,scaleNew] = ...
        ssprk54_production_oracle(rho,dt,ops,scale)
tableau = ipm.evolve.ssprk54Tableau();
z = pack_anisotropic_scale(scale);
stageRhoRates = cell(tableau.stages,1);
stageScaleRates = cell(tableau.stages,1);
for stage = 1:tableau.stages
    stageRho = rho;
    stageScaleVector = z;
    for previous = 1:stage-1
        stageRho = stageRho+dt*tableau.A(stage,previous) * ...
            stageRhoRates{previous};
        stageScaleVector = stageScaleVector + ...
            dt*tableau.A(stage,previous)*stageScaleRates{previous};
    end
    stageScale = unpack_anisotropic_scale(stageScaleVector);
    [stageRhoRates{stage},stageFlow] = ...
        ipm.evolve.flow(stageRho,ops,stageScale);
    stageScaleRates{stage} = anisotropic_scale_rate( ...
        stageScaleVector,stageFlow);
end
rhoNew = rho;
zNew = z;
for stage = 1:tableau.stages
    rhoNew = rhoNew+dt*tableau.b(stage)*stageRhoRates{stage};
    zNew = zNew+dt*tableau.b(stage)*stageScaleRates{stage};
end
scaleNew = unpack_anisotropic_scale(zNew);
[~,flow] = ipm.evolve.flow(rhoNew,ops,scaleNew);
end

function z = pack_anisotropic_scale(scale)
z = [scale.logC_l;scale.logC_y;scale.logC_omega; ...
    scale.physicalTime;scale.X_shift;scale.canonicalTime];
end

function scale = unpack_anisotropic_scale(z)
scale = struct('logC_l',z(1),'logC_x',z(1),'logC_y',z(2), ...
    'logC_omega',z(3),'physicalTime',z(4),'X_shift',z(5), ...
    'canonicalTime',z(6));
end

function rate = anisotropic_scale_rate(z,flow)
rate = [flow.c_x;flow.c_y;flow.c_omega; ...
    flow.timeSpeed*exp(z(3)-z(1)); ...
    flow.c_x*z(5)+flow.c_r;flow.timeSpeed];
end

function [rhoNew,flow,scaleNew] = ...
        ssprk54_isotropic_production_oracle(rho,dt,ops,scale)
tableau = ipm.evolve.ssprk54Tableau();
z = pack_isotropic_scale(scale);
stageRhoRates = cell(tableau.stages,1);
stageScaleRates = cell(tableau.stages,1);
for stage = 1:tableau.stages
    stageRho = rho;
    stageScaleVector = z;
    for previous = 1:stage-1
        stageRho = stageRho+dt*tableau.A(stage,previous) * ...
            stageRhoRates{previous};
        stageScaleVector = stageScaleVector + ...
            dt*tableau.A(stage,previous)*stageScaleRates{previous};
    end
    [stageRhoRates{stage},stageFlow] = ...
        ipm.evolve.flow(stageRho,ops, ...
            unpack_isotropic_scale(stageScaleVector));
    stageScaleRates{stage} = isotropic_scale_rate( ...
        stageScaleVector,stageFlow);
end
rhoNew = rho;
zNew = z;
for stage = 1:tableau.stages
    rhoNew = rhoNew+dt*tableau.b(stage)*stageRhoRates{stage};
    zNew = zNew+dt*tableau.b(stage)*stageScaleRates{stage};
end
scaleNew = unpack_isotropic_scale(zNew);
[~,flow] = ipm.evolve.flow(rhoNew,ops,scaleNew);
end

function z = pack_isotropic_scale(scale)
z = [scale.logC_l;scale.logC_omega;scale.physicalTime; ...
    scale.X_shift;scale.canonicalTime];
end

function scale = unpack_isotropic_scale(z)
scale = struct('logC_l',z(1),'logC_omega',z(2), ...
    'physicalTime',z(3),'X_shift',z(4),'canonicalTime',z(5));
end

function rate = isotropic_scale_rate(z,flow)
rate = [flow.c_l;flow.c_omega; ...
    flow.timeSpeed*exp(z(2)-z(1)); ...
    flow.c_l*z(4)+flow.c_r;flow.timeSpeed];
end

function family = verify_spatial_family(mode,gridSizes,time,threshold)
metricNames = {'DxL2','DxLinf','DyL2','DyLinf', ...
    'poissonL2','poissonLinf','velocityL2','velocityLinf'};
errors = zeros(numel(metricNames),numel(gridSizes));
h = 1./(gridSizes-1);

for level = 1:numel(gridSizes)
    ops = manufactured_operators(gridSizes(level),mode);
    exact = ipmtests.support.mms(ops.X,ops.Y,time);
    probeX = exact.derivativeProbe*ops.Dx';
    probeY = ops.Dy*exact.derivativeProbe;
    errors(1,level) = full_l2(probeX-exact.derivativeProbeX,ops);
    errors(2,level) = max(abs(probeX-exact.derivativeProbeX),[],'all');
    errors(3,level) = full_l2(probeY-exact.derivativeProbeY,ops);
    errors(4,level) = max(abs(probeY-exact.derivativeProbeY),[],'all');

    psi = ipm.field.poisson(exact.source,ops);
    errors(5,level) = full_l2(psi-exact.psi,ops);
    errors(6,level) = max(abs(psi-exact.psi),[],'all');

    [u1,u2] = ipm.field.velocity(exact.rho,ops);
    velocityError = hypot(u1-exact.u1,u2-exact.u2);
    errors(7,level) = full_l2(velocityError,ops);
    errors(8,level) = max(velocityError,[],'all');
end

orders = observed_orders(errors,h);
finalOrders = orders(:,end);
tailOrders = certify_convergence_tail(errors,orders,threshold, ...
    'ipm:HighOrderSpatialOrder',[mode,' all-node spatial family']);

family = struct('gridSizes',gridSizes,'h',h,'metricNames',{metricNames}, ...
    'errors',errors,'orders',orders,'tailOrders',tailOrders, ...
    'finalOrders',finalOrders);
fprintf('  %-9s D/Poisson/velocity final orders: %s\n',mode, ...
    format_named_values(metricNames,finalOrders));
end

function result = verify_ssprk54_conditions()
tableau = ipm.evolve.ssprk54Tableau();
A = tableau.A;
b = tableau.b(:);
c = tableau.c(:);
e = ones(tableau.stages,1);
values = [b'*e; b'*c; b'*(c.^2); b'*A*c; b'*(c.^3); ...
    b'*A*(c.^2); b'*(c.*(A*c)); b'*A*A*c];
targets = [1;1/2;1/3;1/6;1/4;1/12;1/8;1/24];
names = {'b''e','b''c','b''c^2','b''Ac','b''c^3', ...
    'b''A(c^2)','b''diag(c)Ac','b''A^2c'};
residuals = values-targets;
maximumResidual = max(abs(residuals));
assert(isequal(size(A),[5,5]) && all(abs(triu(A)) < eps,'all'), ...
    'ipm:SSPRK54Explicit','The SSPRK54 tableau is not explicit.');
assert(maximumResidual < 5e-13, ...
    'ipm:SSPRK54OrderConditions', ...
    'SSPRK54 fourth-order condition residual %.3e is too large.', ...
    maximumResidual);
result = struct('names',{names},'values',values,'targets',targets, ...
    'residuals',residuals,'maximumResidual',maximumResidual, ...
    'sspCoefficient',tableau.sspCoefficient);
fprintf('  SSPRK54 eight order-condition maximum residual: %.3e.\n', ...
    maximumResidual);
end

function family = verify_transport_family(mode,gridSizes,time,threshold)
metricNames = {'transportL2','transportLinf'};
errors = zeros(numel(metricNames),numel(gridSizes));
h = 1./(gridSizes-1);
for level = 1:numel(gridSizes)
    ops = manufactured_operators(gridSizes(level),mode);
    exact = ipmtests.support.mms(ops.X,ops.Y,time);
    numerical = weno_transport_rhs(exact.rho,exact.u1,exact.u2,ops);
    residual = numerical+exact.transport;
    errors(1,level) = full_l2(residual,ops);
    errors(2,level) = max(abs(residual),[],'all');
end
orders = observed_orders(errors,h);
finalOrders = orders(:,end);
tailOrders = certify_convergence_tail(errors,orders,threshold, ...
    'ipm:HighOrderTransportOrder',[mode,' WENO5-FD transport']);
family = struct('gridSizes',gridSizes,'h',h,'metricNames',{metricNames}, ...
    'errors',errors,'orders',orders,'tailOrders',tailOrders, ...
    'finalOrders',finalOrders);
fprintf('  %-9s WENO5-FD residual final orders: %s\n',mode, ...
    format_named_values(metricNames,finalOrders));
end

function result = verify_forced_evolution(gridSizes,threshold)
% Refine space and time together.  Velocity is recomputed from the evolving
% density at every stage, so this is a nonlinear forced IPM test rather than
% prescribed-velocity advection.
modes = {'uniform','stretched'};
finalTime = 0.15;
metricNames = {'rhoL2','rhoLinf'};
result = struct('gridSizes',gridSizes,'finalTime',finalTime);
for modeIndex = 1:numel(modes)
    mode = modes{modeIndex};
    h = 1./(gridSizes-1);
    errors = zeros(numel(metricNames),numel(gridSizes));
    stepCounts = (gridSizes-1)/2;
    for level = 1:numel(gridSizes)
        ops = manufactured_operators(gridSizes(level),mode);
        initial = ipmtests.support.mms(ops.X,ops.Y,0);
        rho = initial.rho;
        dt = finalTime/stepCounts(level);
        time = 0;
        for step = 1:stepCounts(level)
            rho = forced_ssprk54_step(rho,time,dt,ops);
            time = time+dt;
        end
        exact = ipmtests.support.mms(ops.X,ops.Y,finalTime);
        error = rho-exact.rho;
        errors(1,level) = full_l2(error,ops);
        errors(2,level) = max(abs(error),[],'all');
    end
    orders = observed_orders(errors,h);
    finalOrders = orders(:,end);
    tailOrders = certify_convergence_tail(errors,orders,threshold, ...
        'ipm:HighOrderForcedEvolutionOrder', ...
        [mode,' forced nonlinear IPM evolution']);
    result.(mode) = struct('h',h,'stepCounts',stepCounts, ...
        'dt',finalTime./stepCounts,'metricNames',{metricNames}, ...
        'errors',errors,'orders',orders,'tailOrders',tailOrders, ...
        'finalOrders',finalOrders);
    fprintf('  %-9s forced global final orders: %s\n',mode, ...
        format_named_values(metricNames,finalOrders));
end
end

function rhoNew = forced_ssprk54_step(rho,time,dt,ops)
tableau = ipm.evolve.ssprk54Tableau();
stageRates = cell(tableau.stages,1);
for stage = 1:tableau.stages
    stageRho = rho;
    for previous = 1:stage-1
        stageRho = stageRho+dt*tableau.A(stage,previous)* ...
            stageRates{previous};
    end
    stageTime = time+tableau.c(stage)*dt;
    stageRates{stage} = forced_rhs(stageRho,stageTime,ops);
end
rhoNew = rho;
for stage = 1:tableau.stages
    rhoNew = rhoNew+dt*tableau.b(stage)*stageRates{stage};
end
end

function rhs = forced_rhs(rho,time,ops)
[u1,u2] = ipm.field.velocity(rho,ops);
transportRhs = weno_transport_rhs(rho,u1,u2,ops);
exact = ipmtests.support.mms(ops.X,ops.Y,time);
rhs = transportRhs+exact.forcing;
end

function rhs = weno_transport_rhs(rho,u1,u2,ops)
rhs = ipm.field.transport(rho,u1,u2,ops,ops.transportBoundaryMode);
end

function ops = manufactured_operators(numberOfNodes,mode)
options = struct('nx',numberOfNodes,'ny',numberOfNodes, ...
    'xlim',[-1,1],'ymax',1,'gridMode',mode, ...
    'gridStretchAutomatic',false,'gridStretch',[0.9,0.7], ...
    'spatialDiscretization','high_order', ...
    'farBoundaryMode','dirichlet_zero','symmetryMode','double_odd_omega', ...
    'transportBoundaryMode','closed','transportScheme','weno5_fd', ...
    'adaptiveRemesh',false,'saveResults',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
ops = ipm.mesh.build(ipm.config.resolve(options));
end

function ops = manufactured_rectangular_operators(numberOfYNodes,mode)
numberOfXNodes = 2*numberOfYNodes-1;
options = struct('nx',numberOfXNodes,'ny',numberOfYNodes, ...
    'xlim',[-1,1],'ymax',1,'gridMode',mode, ...
    'gridStretchAutomatic',false,'gridStretch',[0.9,0.7], ...
    'spatialDiscretization','high_order', ...
    'farBoundaryMode','dirichlet_zero','symmetryMode','double_odd_omega', ...
    'transportBoundaryMode','closed','transportScheme','weno5_fd', ...
    'adaptiveRemesh',false,'saveResults',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
ops = ipm.mesh.build(ipm.config.resolve(options));
end

function value = full_l2(error,ops)
weights = ops.integrationWeights;
value = sqrt(sum(weights.*abs(error).^2,'all')/sum(weights,'all'));
end

function orders = observed_orders(errors,h)
orders = log(errors(:,1:end-1)./errors(:,2:end)) ./ ...
    log(h(1:end-1)./h(2:end));
end

function tailOrders = certify_convergence_tail( ...
        errors,orders,threshold,identifier,label)
tailSegmentCount = 2;
assert(size(orders,2) >= tailSegmentCount && ...
    all(isfinite(errors),'all') && all(errors > 0,'all'),identifier, ...
    '%s needs finite positive errors and at least two refinement segments.', ...
    label);
tailOrders = orders(:,end-tailSegmentCount+1:end,:);
assert(all(diff(errors,1,2) < 0,'all') && ...
    all(tailOrders >= threshold,'all'),identifier, ...
    ['%s must decrease monotonically and every one of its last two ' ...
    'observed orders must be at least %.1f; minimum tail order %.3f.'], ...
    label,threshold,min(tailOrders,[],'all'));
end

function text = format_named_values(names,values)
parts = cell(size(names));
for index = 1:numel(names)
    parts{index} = sprintf('%s %.3f',names{index},values(index));
end
text = strjoin(parts,', ');
end
