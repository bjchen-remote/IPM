function report = transport()
%IPMTESTS.BASELINE.TRANSPORT Accuracy and stability checks for transport schemes.

gridSizes = [33,65,129,257];
reconstructionErrors = zeros(size(gridSizes));
operatorErrors = zeros(size(gridSizes));
for level = 1:numel(gridSizes)
    x = test_axis(gridSizes(level),1.2);
    geometry = ipm.field.weno5Geometry(x);
    values = sin(2*x);
    target = 0.5*(x(geometry.faces)+x(geometry.faces+1));
    left = ipm.field.weno5Reconstruct(values,geometry.left,1e-12);
    right = ipm.field.weno5Reconstruct(values,geometry.right,1e-12);
    reconstructionErrors(level) = max([abs(left-sin(2*target)), ...
        abs(right-sin(2*target))]);

    ops = transport_operators(x,'weno5_nonuniform');
    rho = repmat(values,ops.ny,1);
    rhs = ipm.field.transport(rho,ones(size(rho)),zeros(size(rho)),ops);
    exact = -2*cos(2*x);
    interior = 4:ops.nx-3;
    operatorErrors(level) = max(abs(rhs(:,interior)-exact(interior)), ...
        [],'all');
end
reconstructionOrders = log2( ...
    reconstructionErrors(1:end-1)./reconstructionErrors(2:end));
operatorOrders = log2(operatorErrors(1:end-1)./operatorErrors(2:end));

x = test_axis(65,0.9);
ops = transport_operators(x,'weno5_nonuniform');
constantRho = 1.75*ones(ops.ny,ops.nx);
constantResidual = max(abs(ipm.field.transport(constantRho, ...
    ones(size(constantRho)),zeros(size(constantRho)),ops)),[],'all');
[X,Y] = meshgrid(x,ops.y);
rho = exp(-3*(X.^2+(Y-0.5).^2));
u1 = sin(pi*X).*cos(pi*Y);
u2 = -cos(pi*X).*sin(pi*Y);
closedRhs = ipm.field.transport(rho,u1,u2,ops,'closed');
closedMassRate = sum(closedRhs.*ops.weights,'all');

jumpX = test_axis(33,1.1);
jumpGeometry = ipm.field.weno5Geometry(jumpX);
jumpData = interval_values(jumpX,-0.55,0.13);
jumpValues = [ipm.field.weno5Reconstruct(jumpData, ...
    jumpGeometry.left,1e-12), ...
    ipm.field.weno5Reconstruct(jumpData,jumpGeometry.right,1e-12)];
jumpRange = [min(jumpValues),max(jumpValues)];

advectionX = test_axis(129,0.8);
initialInterval = [-0.55,-0.05];
finalTime = 0.3;
advection = struct();
schemes = {'muscl_minmod','weno5_nonuniform'};
for schemeNumber = 1:numel(schemes)
    scheme = schemes{schemeNumber};
    advectionOps = transport_operators(advectionX,scheme);
    initial = interval_values(advectionOps.x, ...
        initialInterval(1),initialInterval(2));
    evolved = advance_linear_advection(initial,advectionOps,finalTime,0.25);
    exact = interval_values(advectionOps.x, ...
        initialInterval(1)+finalTime,initialInterval(2)+finalTime);
    data = struct('l1Error',sum(abs(evolved-exact).*advectionOps.hx), ...
        'range',[min(evolved),max(evolved)], ...
        'totalVariation',sum(abs(diff(evolved))));
    advection.(scheme) = data;
end

assert(all(reconstructionOrders(end-1:end) > 4.2), ...
    'ipm:WenoReconstructionOrder', ...
    'Nonuniform WENO reconstruction did not attain fifth order.');
assert(all(operatorOrders(end-1:end) > 1.9), ...
    'ipm:WenoTransportOrder', ...
    'The nodal conservative WENO transport did not attain second order.');
assert(constantResidual < 1e-12, ...
    'ipm:WenoConstantState','WENO transport changed a constant state.');
assert(abs(closedMassRate) < 1e-12, ...
    'ipm:WenoMass','Closed WENO transport is not conservative.');
assert(jumpRange(1) > -1e-8 && jumpRange(2) < 1+1e-8, ...
    'ipm:WenoJump','WENO face reconstruction oscillated at a jump.');
assert(advection.weno5_nonuniform.range(1) > -2e-3 && ...
    advection.weno5_nonuniform.range(2) < 1+2e-3, ...
    'ipm:WenoAdvectionRange', ...
    'WENO-SSPRK3 generated a material overshoot in the step test.');
assert(advection.weno5_nonuniform.l1Error < ...
    advection.muscl_minmod.l1Error, ...
    'ipm:WenoAdvectionError', ...
    'WENO was not sharper than MUSCL in the step-advection test.');

report = struct('gridSizes',gridSizes, ...
    'reconstructionErrors',reconstructionErrors, ...
    'reconstructionOrders',reconstructionOrders, ...
    'operatorErrors',operatorErrors, ...
    'operatorOrders',operatorOrders, ...
    'constantResidual',constantResidual, ...
    'closedMassRate',closedMassRate, ...
    'jumpRange',jumpRange,'advection',advection);
fprintf(['Transport verification: WENO reconstruction/operator final ' ...
    'orders %.3f/%.3f, step L1 WENO/MUSCL %.3e/%.3e.\n'], ...
    reconstructionOrders(end),operatorOrders(end), ...
    advection.weno5_nonuniform.l1Error, ...
    advection.muscl_minmod.l1Error);
end

function x = test_axis(numberOfNodes,stretch)
reference = linspace(-1,1,numberOfNodes);
x = sinh(stretch*reference)/sinh(stretch);
end

function values = interval_values(axis,leftBound,rightBound)
values = double(axis >= leftBound & axis <= rightBound);
end

function ops = transport_operators(x,scheme)
ny = 9;
y = linspace(0,1,ny)';
edgesX = [x(1),0.5*(x(1:end-1)+x(2:end)),x(end)];
edgesY = [y(1);0.5*(y(1:end-1)+y(2:end));y(end)];
ops = struct('nx',numel(x),'ny',ny,'x',x,'y',y, ...
    'dxFaces',diff(x),'dyFaces',diff(y), ...
    'hx',diff(edgesX),'hy',diff(edgesY), ...
    'weights',diff(edgesY)*diff(edgesX), ...
    'transportBoundaryMode','open','transportScheme',scheme, ...
    'wallTransportMode','conservative_flux', ...
    'wenoEpsilon',1e-12,'edgesX',edgesX);
if strcmp(scheme,'weno5_nonuniform')
    ops.wenoX = ipm.field.weno5Geometry(x);
    ops.wenoY = ipm.field.weno5Geometry(y');
end
end

function evolved = advance_linear_advection(initial,ops,finalTime,cfl)
state = repmat(initial,ops.ny,1);
u1 = ones(size(state));
u2 = zeros(size(state));
time = 0;
while time < finalTime
    dt = min(cfl*min(ops.hx),finalTime-time);
    first = state+dt*ipm.field.transport(state,u1,u2,ops);
    second = 0.75*state+0.25*(first+dt* ...
        ipm.field.transport(first,u1,u2,ops));
    state = (state+2*(second+dt* ...
        ipm.field.transport(second,u1,u2,ops)))/3;
    time = time+dt;
end
evolved = state(ceil(ops.ny/2),:);
end
