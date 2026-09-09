function report = ipm_perflab_verify_sylvester()
%IPM_PERFLAB_VERIFY_SYLVESTER Bounded same-operator Schur/eigen screening.
%   Maintained admitted grids test MMS and unchanged image-Green boundary.
%   A separate raw-axis ladder tests maintained derivative matrices without
%   claiming that highly stretched tiny axes pass production grid gates.
methods = {'direct_cached','real_schur_cached', ...
    'fast_diagonalization','matlab_sylvester'};
cases = {};
rejectedGrids = {};
specifications = [4,17,0;4,33,0;4,65,0;4,65,2;4,65,4; ...
    6,33,0;6,65,0;6,65,2;6,65,4];
for specification = specifications.'
    order = specification(1);
    n = specification(2);
    stretch = specification(3);
    opts = small_options(order,n,stretch);
    try
        config = ipm.config.resolve(opts);
        ops = ipm.mesh.build(config);
    catch errorInfo
        rejectedGrids{end+1} = struct('order',order,'n',n,'stretch',stretch, ...
            'identifier',errorInfo.identifier,'message',errorInfo.message); %#ok<AGROW>
        fprintf('ORIGINAL GRID REJECT order%d n%d stretch%.3g %s\n', ...
            order,n,stretch,errorInfo.identifier);
        continue;
    end
    exact = mms(ops.X,ops.Y);
    boundary = struct('left',exact.psi(:,1),'right',exact.psi(:,end), ...
        'bottom',exact.psi(1,:),'top',exact.psi(end,:));
    [directMms,rhsMms] = ipm.field.poisson(exact.source,ops,1,boundary);
    rho = -exp(-ops.X.^2-ops.Y.^2);
    source = rho*ops.Dx.';
    [directGreen,rhsGreen] = ipm.field.poisson(source,ops);
    for methodIndex = 1:numel(methods)
        [cache,info] = ipm_perflab_sylvester_factor(ops.Ty,ops.Tx,methods{methodIndex});
        for problem = 1:2
            if problem == 1
                rhs = rhsMms;
                directPsi = directMms;
                problemName = 'admitted_nonzero_boundary_mms';
            else
                rhs = rhsGreen;
                directPsi = directGreen;
                problemName = 'admitted_image_green_frozen';
            end
            item = evaluate(cache,info,ops.A,rhs,directPsi(2:end-1,2:end-1));
            psi = directPsi;
            psi(2:end-1,2:end-1) = item.solution;
            item = rmfield(item,'solution');
            item.problem = problemName;
            item.order = order;
            item.n = n;
            item.stretch = stretch;
            item.rowScaleOrders = row_orders(ops.A);
            item.velocityRelativeDifference = velocity_difference(psi,directPsi,ops);
            item.boundaryExact = isequaln(psi([1,end],:),directPsi([1,end],:)) && ...
                isequaln(psi(:,[1,end]),directPsi(:,[1,end]));
            item.psiMmsRelativeError = NaN;
            item.velocityMmsRelativeError = NaN;
            if problem == 1
                item.psiMmsRelativeError = relative(psi,exact.psi);
                [u1,u2] = velocity(psi,ops);
                item.velocityMmsRelativeError = max(relative(u1,exact.u1),relative(u2,exact.u2));
            end
            item.screeningPass = item.finiteSolution && item.boundaryExact && ...
                item.originalRelativeResidual < 1e-10 && ...
                item.componentwiseBackwardError < 1e-10 && ...
                item.solutionRelativeDifference < 1e-9 && ...
                item.velocityRelativeDifference < 1e-8;
            cases{end+1} = item; %#ok<AGROW>
            print_item(item);
        end
    end
end
% Each raw case is the exact maintained fourth-order derivative operator,
% but deliberately has no claim of mesh.build admissibility or PDE accuracy.
rawCases = {};
n = 65;
H = 1e6;
for stretch = [0,4,8,12,16,20,24]
    sx = linspace(-1,1,n);
    sy = linspace(0,1,n).';
    if stretch == 0
        x = H*sx;
        y = H*sy;
    else
        x = H*sinh(stretch*sx)/sinh(stretch);
        y = H*sinh(stretch*sy)/sinh(stretch);
    end
    TxFull = -ipm.mesh.fdMatrix(x,2,7);
    TyFull = -ipm.mesh.fdMatrix(y,2,7);
    Tx = TxFull(2:end-1,2:end-1);
    Ty = TyFull(2:end-1,2:end-1);
    A = kron(speye(n-2),Ty)+kron(Tx,speye(n-2));
    [Xi,Yi] = meshgrid(sx(2:end-1),sy(2:end-1));
    exactDiscrete = sin(pi*(Xi+1)/2).*sin(pi*Yi) + ...
        0.1*sin(7*pi*(Xi+1)/2).*sin(5*pi*Yi);
    rhs = reshape(A*exactDiscrete(:),n-2,n-2);
    direct = reshape(decomposition(A,'lu')\rhs(:),n-2,n-2);
    for methodIndex = 1:numel(methods)
        [cache,info] = ipm_perflab_sylvester_factor(Ty,Tx,methods{methodIndex});
        item = evaluate(cache,info,A,rhs,direct);
        item.solutionRelativeError = relative(item.solution,exactDiscrete);
        item = rmfield(item,'solution');
        item.problem = 'raw_maintained_stencil_known_discrete_solution';
        item.order = 4;
        item.n = n;
        item.H = H;
        item.stretch = stretch;
        item.rowScaleOrders = row_orders(A);
        item.minimumSpacing = min([diff(x),diff(y).']);
        item.maximumSpacing = max([diff(x),diff(y).']);
        item.productionGridAdmissibilityClaimed = false;
        item.screeningPass = item.finiteSolution && ...
            item.componentwiseBackwardError < 1e-10 && ...
            item.solutionRelativeError < 1e-8 && ...
            item.solutionRelativeDifference < 1e-8;
        rawCases{end+1} = item; %#ok<AGROW>
        print_item(item);
    end
end
report = struct('schemaVersion',1,'kind','ipm_kronecker_sylvester_screening', ...
    'matlabVersion',version,'computer',computer,'computationalThreads',maxNumCompThreads, ...
    'cases',{cases},'originalRejectedGrids',{rejectedGrids},'rawCases',{rawCases}, ...
    'largeGridTestPerformed',false,'pdeEvolutionPerformed',false, ...
    'interpretation','Screening thresholds are research comparisons, not modified production gates.');
end

function item = evaluate(cache,info,A,rhs,direct)
lastwarn('');
errorIdentifier = '';
samples = nan(1,3);
try
    timer = tic;
    solution = ipm_perflab_sylvester_solve(cache,rhs);
    firstSolve = toc(timer);
    [~,warningIdentifier] = lastwarn;
    for index = 1:3
        timer = tic;
        ipm_perflab_sylvester_solve(cache,rhs);
        samples(index) = toc(timer);
    end
catch errorInfo
    solution = nan(size(rhs));
    firstSolve = NaN;
    warningIdentifier = '';
    errorIdentifier = errorInfo.identifier;
end
residual = A*solution(:)-rhs(:);
denominator = abs(A)*abs(solution(:))+abs(rhs(:));
rowMaximum = full(max(abs(A),[],2));
item = struct('method',cache.method,'factorInfo',info,'solution',solution, ...
    'firstSolveSeconds',firstSolve,'solveSamplesSeconds',samples, ...
    'solveMedianSeconds',median(samples),'warningIdentifier',warningIdentifier, ...
    'errorIdentifier',errorIdentifier, ...
    'originalRelativeResidual',norm(residual,inf)/max(norm(rhs(:),inf),eps), ...
    'componentwiseBackwardError',max(abs(residual)./max(denominator,realmin)), ...
    'rowNormalizedResidual',norm(residual./rowMaximum,inf)/ ...
        max(norm(rhs(:)./rowMaximum,inf),realmin), ...
    'solutionRelativeDifference',relative(solution,direct), ...
    'imaginaryRelativeMagnitude',max(abs(imag(solution)),[],'all')/ ...
        max(max(abs(solution),[],'all'),realmin), ...
    'finiteSolution',all(isfinite(solution),'all'));
end

function opts = small_options(order,n,stretch)
opts = struct('nx',n,'ny',n,'xlim',[-1,1],'ymax',1, ...
    'gridMode','stretched','gridStretchAutomatic',false,'gridStretch',[stretch,stretch], ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'rescalingMode','physical','farBoundaryMode','green', ...
    'transportBoundaryMode','closed','wallTransportMode','conservative_flux', ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'saveResults',false,'makePlots',false,'livePlot',false, ...
    'writeVideo',false,'verbose',false);
if order == 6
    opts = ipm.config.sixthOrder(opts);
    opts.spatialDiscretization = 'sixth_order';
    opts.transportScheme = 'weno7_fd';
    opts.timeIntegrator = 'rk6';
    opts.remeshTransferScheme = 'sixth_order';
    opts.farBoundaryMode = 'green';
end
end

function exact = mms(X,Y)
g = exp(0.3*X).*cos(1.2*X);
gx = exp(0.3*X).*(0.3*cos(1.2*X)-1.2*sin(1.2*X));
gxx = exp(0.3*X).*((0.3^2-1.2^2)*cos(1.2*X)-0.72*sin(1.2*X));
h = Y.*(1-Y);
exact = struct('psi',g.*h,'source',-gxx.*h+2*g, ...
    'u1',-g.*(1-2*Y),'u2',gx.*h);
end

function [u1,u2] = velocity(psi,ops)
u1 = -(ops.Dy*psi);
u2 = psi*ops.Dx.';
u2(1,:) = 0;
end

function difference = velocity_difference(actual,expected,ops)
[a1,a2] = velocity(actual,ops);
[e1,e2] = velocity(expected,ops);
difference = max(relative(a1,e1),relative(a2,e2));
end

function difference = relative(actual,expected)
difference = max(abs(actual-expected),[],'all')/ ...
    max(max(abs(expected),[],'all'),realmin);
end

function orders = row_orders(A)
maximum = full(max(abs(A),[],2));
orders = log10(max(maximum))-log10(min(maximum));
end

function print_item(item)
fprintf('%s order%d n%d stretch%.3g %s rowOrders%.3g dx%.3g back%.3g pass%d\n', ...
    item.problem,item.order,item.n,item.stretch,item.method,item.rowScaleOrders, ...
    item.solutionRelativeDifference,item.componentwiseBackwardError,item.screeningPass);
end
