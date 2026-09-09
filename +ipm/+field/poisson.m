function [psi,interiorRhs,boundary,solveInfo,poissonOperator] = ...
    poisson(source,ops,kappa,boundaryOverride)
%IPM.FIELD.POISSON Solve -Delta psi = source on a truncated upper half-plane.
%   The physical wall has psi=0. Artificial sides use the configured
%   quadrant/half-plane image Green formula or homogeneous Dirichlet data.
%   A fourth BOUNDARYOVERRIDE argument supplies explicit Dirichlet traces;
%   it is primarily an oracle hook for verifying the complete boundary RHS.

if nargin < 3
    kappa = 1;
end
validateattributes(kappa,{'numeric'},{'scalar','real','finite','positive'});

if ~isequal(size(source), [ops.ny, ops.nx])
    error('ipm:SizeMismatch', 'source must have size [ny,nx].');
end
if nargin >= 4 && ~isempty(boundaryOverride)
    boundary = validate_boundary_override(boundaryOverride,ops);
elseif kappa == 1
    boundary = ipm.field.greenBoundary(source,ops);
else
    boundary = ipm.field.greenBoundary(source,ops,kappa);
end
psi = zeros(ops.ny,ops.nx);
psi(:,1) = boundary.left;
psi(:,end) = boundary.right;
psi(1,:) = boundary.bottom;
psi(end,:) = boundary.top;
rhs = source(2:end-1,2:end-1);
highOrder = strcmp(ops.spatialDiscretization,'high_order') || ...
    strcmp(ops.spatialDiscretization,'sixth_order');
if highOrder
    xBoundaryValues = [boundary.left(2:end-1), ...
        boundary.right(2:end-1)];
    yBoundaryValues = [boundary.bottom(2:end-1); ...
        boundary.top(2:end-1)];
    rhs = rhs-(1/kappa)*xBoundaryValues*ops.poissonBoundary.x' - ...
        kappa*ops.poissonBoundary.y*yBoundaryValues;
    [poissonOperator,~] = ipm.field.poissonOperator(ops,kappa);
    if kappa == 1
        interior = ops.poisson\rhs(:);
        solver = 'cached_direct_high_order';
    else
        if ~strcmp(ops.anisotropicPoisson.solver,'direct')
            error('ipm:HighOrderPoissonSolver', ...
                ['The nonsymmetric high-order boundary closure currently ' ...
                'supports the direct anisotropic Poisson solver only.']);
        end
        factor = decomposition(poissonOperator,'lu');
        interior = factor\rhs(:);
        solver = 'direct_high_order';
    end
    relativeResidual = norm(poissonOperator*interior-rhs(:),inf) / ...
        max(norm(rhs(:),inf),eps);
    solveInfo = struct('solver',solver,'kappa',kappa,'flag',0, ...
        'relativeResidual',relativeResidual,'iterations',0);
elseif kappa == 1
    rhs(:,1) = rhs(:,1)-ops.poissonBoundary.x(1)*boundary.left(2:end-1);
    rhs(:,end) = rhs(:,end)- ...
        ops.poissonBoundary.x(2)*boundary.right(2:end-1);
    rhs(1,:) = rhs(1,:)-ops.poissonBoundary.y(1)*boundary.bottom(2:end-1);
    rhs(end,:) = rhs(end,:)-ops.poissonBoundary.y(2)*boundary.top(2:end-1);
    weightedRhs = ops.poissonWeights.*rhs(:);
    interior = ops.poisson \ weightedRhs;
    poissonOperator = ops.A;
    solveInfo = struct('solver','cached_direct','kappa',1,'flag',0, ...
        'relativeResidual',NaN,'iterations',0);
else
    inverseKappa = 1/kappa;
    rhs(:,1) = rhs(:,1)-inverseKappa*ops.poissonBoundary.x(1)* ...
        boundary.left(2:end-1);
    rhs(:,end) = rhs(:,end)-inverseKappa*ops.poissonBoundary.x(2)* ...
        boundary.right(2:end-1);
    rhs(1,:) = rhs(1,:)-kappa*ops.poissonBoundary.y(1)* ...
        boundary.bottom(2:end-1);
    rhs(end,:) = rhs(end,:)-kappa*ops.poissonBoundary.y(2)* ...
        boundary.top(2:end-1);
    weightedRhs = ops.poissonWeights.*rhs(:);
    [poissonOperator,weightedOperator] = ...
        ipm.field.poissonOperator(ops,kappa);
    solver = ops.anisotropicPoisson.solver;
    switch solver
        case 'direct'
            factor = decomposition(weightedOperator,'chol');
            interior = factor \ weightedRhs;
            flag = 0;
            iterations = 0;
            relativeResidual = norm(weightedOperator*interior-weightedRhs,inf) / ...
                max(norm(weightedRhs,inf),eps);
        case 'pcg'
            preconditioner = @(value) ops.poisson \ value;
            [interior,flag,relativeResidual,iterations] = pcg( ...
                weightedOperator,weightedRhs, ...
                ops.anisotropicPoisson.tolerance, ...
                ops.anisotropicPoisson.maxIterations,preconditioner);
            if flag ~= 0
                error('ipm:AnisotropicPoissonPcg', ...
                    ['Aspect-aware PCG failed with flag %d and relative ' ...
                    'residual %.3e after %d iterations.'], ...
                    flag,relativeResidual,iterations);
            end
        otherwise
            error('ipm:AnisotropicPoissonSolver', ...
                'Unknown aspect-aware Poisson solver ''%s''.',solver);
    end
    solveInfo = struct('solver',solver,'kappa',kappa,'flag',flag, ...
        'relativeResidual',relativeResidual,'iterations',iterations);
end
psi(2:end-1, 2:end-1) = reshape(interior,ops.ny-2,ops.nx-2);
interiorRhs = rhs;
end

function boundary = validate_boundary_override(boundary,ops)
if ~isstruct(boundary) || ~isscalar(boundary)
    error('ipm:PoissonBoundaryOverride', ...
        'The Poisson boundary override must be a scalar structure.');
end
names = {'left','right','bottom','top'};
sizes = {[ops.ny,1],[ops.ny,1],[1,ops.nx],[1,ops.nx]};
for index = 1:numel(names)
    name = names{index};
    if ~isfield(boundary,name) || ...
            ~isequal(size(boundary.(name)),sizes{index}) || ...
            ~isnumeric(boundary.(name)) || ...
            ~isreal(boundary.(name)) || ...
            any(~isfinite(boundary.(name)),'all')
        error('ipm:PoissonBoundaryOverride', ...
            'Boundary field "%s" must be a finite real array of size %s.', ...
            name,mat2str(sizes{index}));
    end
end
end
