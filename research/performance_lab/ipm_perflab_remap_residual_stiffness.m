function report = ipm_perflab_remap_residual_stiffness(outputDirectory)
%IPM_PERFLAB_REMAP_RESIDUAL_STIFFNESS Exact remainder need not have zero CFL.
% Frozen velocity V=alpha*X has identically zero residual velocity. The
% maintained mapped WENO F and the polynomial remap generator B still differ.
% Small one-row Jacobians expose the discrete remainder. No nonlinear IPM
% stability is inferred from this necessary linearized screening condition.
assert(~exist(outputDirectory,'dir'),'Pass a new directory.');
mkdir(outputDirectory);
items = {};
for nx = [17,33,65]
    opts = struct('nx',nx,'ny',9,'xlim',[-4,4],'ymax',4, ...
        'gridMode','uniform','gridStretchAutomatic',false,'gridStretch',[0,0], ...
        'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
        'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
        'transportBoundaryMode','open','wallTransportMode','advective_upwind', ...
        'farBoundaryMode','dirichlet_zero','symmetryMode','double_odd_omega', ...
        'saveResults',false,'makePlots',false,'verbose',false,'adaptiveRemesh',false);
    ops = ipm.mesh.build(ipm.config.resolve(opts));
    g = ipm_perflab_remap_generator(ops.x,6);
    alpha = 0.2;
    velocity = repmat(alpha*ops.x,ops.ny,1);
    row = 5;
    density = repmat(exp(-ops.x.^2),ops.ny,1);
    q = density(row,:).';
    B = -alpha*diag(ops.x)*full(g.D1);
    fullRhs = line_rhs(q,density,row,velocity,ops);
    fullRate = max(abs(velocity)./ops.hx,[],'all');
    jacobians = cell(1,2);
    for probe = 1:2
        delta = 1e-5/2^(probe-1);
        J = zeros(nx,nx);
        for column = 1:nx
            direction = zeros(nx,1); direction(column) = delta;
            plus = line_rhs(q+direction,density,row,velocity,ops);
            minus = line_rhs(q-direction,density,row,velocity,ops);
            J(:,column) = (plus-minus)/(2*delta);
        end
        jacobians{probe} = J;
    end
    J = jacobians{2}; N = J-B;
    factors = [1,2,4,8,10];
    amplifications = cell(1,numel(factors));
    for k = 1:numel(factors)
        h = factors(k)*0.5/fullRate;
        Ehalf = expm(0.5*h*B);
        lawson = expm(h*B)+h*Ehalf*N*Ehalf*(eye(nx)+0.5*h*N);
        explicitN = eye(nx)+h*N+0.5*h^2*N*N;
        exact = expm(h*J);
        amplifications{k} = struct('dt',h,'fullCfl',h*fullRate, ...
            'residualVelocityCfl',0, ...
            'lawsonSpectralRadius',max(abs(eig(lawson))), ...
            'lawsonNorm2',norm(lawson,2),'exactNorm2',norm(exact,2), ...
            'lawsonVersusExactNorm2',norm(lawson-exact,2), ...
            'explicitRemainderSpectralRadius',max(abs(eig(explicitN))));
    end
    item = struct('nx',nx,'alpha',alpha,'fullRate',fullRate, ...
        'residualVelocityRate',0,'fullRhsInf',norm(fullRhs,Inf), ...
        'remainderRhsInf',norm(fullRhs-B*q,Inf), ...
        'fullJacobianSpectralRadius',max(abs(eig(J))), ...
        'remainderJacobianSpectralRadius',max(abs(eig(N))), ...
        'fullJacobianNorm2',norm(J,2),'remainderJacobianNorm2',norm(N,2), ...
        'jacobianPerturbationDifferenceNorm2',norm(jacobians{2}-jacobians{1},2), ...
        'relativeJacobianPerturbationDifference',norm(jacobians{2}-jacobians{1},2)/norm(J,2), ...
        'amplifications',{amplifications});
    save(fullfile(outputDirectory,sprintf('linearized_row_%d.mat',nx)), ...
        'item','q','B','J','N','jacobians');
    items{end+1} = item; %#ok<AGROW>
    fprintf('Pure-background %d nodes: residual velocity 0, N RHS %.3e, spectral J/N %.3g/%.3g.\n', ...
        nx,item.remainderRhsInf,item.fullJacobianSpectralRadius,item.remainderJacobianSpectralRadius);
end
report = struct('schemaVersion',1,'kind','frozen_weno_remap_remainder_stiffness', ...
    'items',{items},'fullNonlinearStabilityClaimed',false, ...
    'interpretation','Zero residual velocity does not remove discrete full-WENO remainder.');
save(fullfile(outputDirectory,'residual_stiffness_report.mat'),'report');
fid = fopen(fullfile(outputDirectory,'residual_stiffness_report.json'),'w');
cleanup = onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
end

function value = line_rhs(q,density,row,velocity,ops)
density(row,:) = q.';
rhs = ipm.field.transport(density,velocity,zeros(size(velocity)),ops);
value = rhs(row,:).';
end
