function report = ipm_perflab_verify_background_remap(outputDirectory,frozenResultFile)
%IPM_PERFLAB_VERIFY_BACKGROUND_REMAP Frozen-background interpolation study.
%   Synthetic analytic data plus optional saved v2 rho. No elliptic solve,
%   nonlinear PDE stage, updated alpha, gauge, or native checkpoint writes.

assert(nargin >= 1 && ~exist(outputDirectory,'dir'),'Pass a new output directory.');
if nargin < 2
    frozenResultFile = '';
end
mkdir(outputDirectory);
methods = {'linear','pchip','poly6','poly8'};
dts = [0.004,0.008,0.016,0.032];
oracle = maintained_oracle();
synthetic = {};
invariants = {};
for nx = [65,129]
    ny = (nx+1)/2;
    for stretch = [0,2,4]
        x = axis_map(linspace(-1,1,nx),stretch,8);
        y = axis_map(linspace(0,1,ny).',stretch/2,4);
        [X,Y] = meshgrid(x,y);
        alpha = 0.2*(1-exp(-(y/0.2).^2));
        rho = analytic_density(X,Y);
        Dx = ipm.mesh.fdMatrix(x,1,7);
        weights = control_width(y)*control_width(x);
        regions = struct('lab',abs(X)<=2 & Y<=2, ...
            'peak',abs(abs(X)-1.1/sqrt(2))<=0.4 & Y<=0.5);
        for methodIndex = 1:numel(methods)
            method = methods{methodIndex};
            invariant = invariant_checks(rho,x,alpha,method);
            invariant.grid = [ny,nx];
            invariant.stretch = stretch;
            invariants{end+1} = invariant; %#ok<AGROW>
            for dt = dts
                [actual,fullInfo,samples] = timed_remap(rho,x,alpha,dt,method);
                exact = analytic_density(exp(-alpha*dt).*X,Y);
                exactRhoX = analytic_derivative(exp(-alpha*dt).*X,Y).*exp(-alpha*dt);
                half = ipm_perflab_background_remap(rho,x,alpha,dt/2,method);
                twice = ipm_perflab_background_remap(half,x,alpha,dt/2,method);
                [reverse,reverseInfo] = ipm_perflab_background_remap(actual,x,alpha,-dt,method);
                item = struct('kind','analytic_frozen_background','method',method, ...
                    'grid',[ny,nx],'stretch',stretch,'dt',dt, ...
                    'forward',pair_metrics(actual,exact,Dx,weights,regions), ...
                    'rhoXAnalytic',scalar_metrics(actual*Dx.',exactRhoX,weights,regions), ...
                    'semigroup',pair_metrics(twice,actual,Dx,weights,regions), ...
                    'reversibility',pair_metrics(reverse,rho,Dx,weights,regions), ...
                    'forwardInfo',fullInfo,'reverseInfo',reverseInfo, ...
                    'solveSamplesSeconds',samples,'solveMedianSeconds',median(samples), ...
                    'evenSymmetryRelativeDefect',even_defect(actual));
                synthetic{end+1} = item; %#ok<AGROW>
            end
        end
    end
end
frozen = {};
frozenFields = {};
frozenMetadata = struct();
if ~isempty(frozenResultFile)
    r = ipm.output.validate(frozenResultFile);
    assert(all(ipm.output.trustedMask(r)) && ...
        strcmp(r.config.transport.spatialDiscretization,'high_order') && ...
        strcmp(r.config.scaling.lengthGauge,'transport_anchor') && r.scale.rates.cr==0, ...
        'Frozen result must be paired, trusted, fourth-order and zero-translation anchor.');
    x = r.grid.x(:).';
    y = r.grid.y(:);
    [X,Y] = meshgrid(x,y);
    alpha = r.scale.rates.cx+interp1(x,r.state.velocity.x.',1,'linear').';
    assert(all(isfinite(alpha)) && all(alpha>=0), ...
        'Frozen forward tests require the screened nonnegative alpha.');
    rho = r.state.rho;
    Dx = ipm.mesh.fdMatrix(x,1,7);
    derivativePairingDefect = max(abs(rho*Dx.'-r.state.omega),[],'all');
    assert(derivativePairingDefect == 0,'Frozen rho and maintained rho_x must pair exactly.');
    weights = control_width(y)*control_width(x);
    peakX = r.history.gauge.omegaGaugeQuadraticPeakX(end);
    wx = r.history.mesh.trackedWallCoreWidth(end);
    wy = r.history.mesh.trackedVerticalCoreWidth(end);
    regions = struct('lab',abs((X-r.scale.Xshift)/r.scale.Cx)<=2 & Y/r.scale.Cy<=2, ...
        'peak',abs(abs(X)-peakX)<=3*wx & Y<=3*wy);
    frozenMetadata = struct('resultFile',frozenResultFile,'tau',r.state.canonicalTime, ...
        'physicalTime',r.state.physicalTime,'grid',size(rho),'x',x,'y',y,'alpha',alpha, ...
        'physicalLabWindow',[-2,2,2],'rescaledPeakX',peakX,'strictCoreWidths',[wx,wy], ...
        'labPoints',nnz(regions.lab),'peakPoints',nnz(regions.peak), ...
        'fieldsInRescaledUnits',true,'densityDerivativePairingDefect',derivativePairingDefect, ...
        'forwardReference','poly8 reconstruction of saved nodal data; not an exact solution', ...
        'semigroupReference','one frozen full remap versus two frozen half remaps', ...
        'reverseReference','saved nodal rho on covered interior; unknown inflow is NaN');
    for methodIndex = 1:numel(methods)
        invariant = invariant_checks(rho,x,alpha,methods{methodIndex});
        invariant.grid = size(rho);
        invariant.stretch = NaN;
        invariants{end+1} = invariant; %#ok<AGROW>
    end
    for dt = dts
        reference = ipm_perflab_background_remap(rho,x,alpha,dt,'poly8');
        for methodIndex = 1:numel(methods)
            method = methods{methodIndex};
            [actual,fullInfo,samples] = timed_remap(rho,x,alpha,dt,method);
            half = ipm_perflab_background_remap(rho,x,alpha,dt/2,method);
            twice = ipm_perflab_background_remap(half,x,alpha,dt/2,method);
            [reverse,reverseInfo] = ipm_perflab_background_remap(actual,x,alpha,-dt,method);
            item = struct('kind','frozen_saved_profile','method',method,'dt',dt, ...
                'referenceDifference',pair_metrics(actual,reference,Dx,weights,regions), ...
                'semigroup',pair_metrics(twice,actual,Dx,weights,regions), ...
                'reversibility',pair_metrics(reverse,rho,Dx,weights,regions), ...
                'forwardInfo',fullInfo,'reverseInfo',reverseInfo, ...
                'solveSamplesSeconds',samples,'solveMedianSeconds',median(samples), ...
                'evenSymmetryRelativeDefect',even_defect(actual));
            frozen{end+1} = item; %#ok<AGROW>
            % Preserve localized signed defects, rather than merely saying
            % a scalar norm passed. Coordinates/indices are included once.
            fieldItem = struct('method',method,'dt',dt, ...
                'referenceDifference',actual(regions.lab)-reference(regions.lab), ...
                'semigroupDifference',twice(regions.lab)-actual(regions.lab), ...
                'reverseDifference',reverse(regions.lab)-rho(regions.lab));
            dxActual = actual*Dx.';
            dxReference = reference*Dx.';
            dxTwice = twice*Dx.';
            dxReverse = reverse*Dx.';
            fieldItem.rhoXReferenceDifference = dxActual(regions.lab)-dxReference(regions.lab);
            fieldItem.rhoXSemigroupDifference = dxTwice(regions.lab)-dxActual(regions.lab);
            fieldItem.rhoXReverseDifference = dxReverse(regions.lab)-r.state.omega(regions.lab);
            frozenFields{end+1} = fieldItem; %#ok<AGROW>
            fprintf(['BACKGROUND_REMAP %s dt%.3g median%.6gs labRho/RhoX=%.3g/%.3g ' ...
                'peakRhoX=%.3g semiPeakRhoX=%.3g reversePeakRhoX=%.3g\n'], ...
                method,dt,item.solveMedianSeconds, ...
                item.referenceDifference.rho.lab.relativeL2, ...
                item.referenceDifference.rhoX.lab.relativeL2, ...
                item.referenceDifference.rhoX.peak.relativeL2, ...
                item.semigroup.rhoX.peak.relativeL2,item.reversibility.rhoX.peak.relativeL2);
        end
    end
    frozenMetadata.labLinearIndices = find(regions.lab);
    save(fullfile(outputDirectory,'frozen_remap_deviations.mat'), ...
        'frozenMetadata','frozenFields','-v7.3');
end
report = struct('schemaVersion',1,'kind','ipm_frozen_background_remap_screening', ...
    'matlabVersion',version,'computationalThreads',maxNumCompThreads, ...
    'dts',dts,'methods',{methods},'maintainedInterpolationOracle',oracle, ...
    'invariants',{invariants},'synthetic',{synthetic},'frozen',{frozen}, ...
    'frozenMetadata',rmfield_if_present(frozenMetadata,{'x','y','alpha','labLinearIndices'}), ...
    'poissonOperatorBuilt',false,'nonlinearPdeAdvanced',false,'nativeCheckpointWritten',false, ...
    'interpretation',['Frozen advective subproblem interpolation only. Poly8 is a ' ...
    'reconstruction comparison, not exact truth. Semigroup/reversal tests do not prove ' ...
    'nonlinear temporal order or authorize a larger full-PDE timestep.']);
save(fullfile(outputDirectory,'background_remap_report.mat'),'report');
file = fopen(fullfile(outputDirectory,'background_remap_report.json'),'w');
assert(file>=0,'Cannot write remap report.');
cleanup = onCleanup(@()fclose(file));
fprintf(file,'%s\n',jsonencode(report,'PrettyPrint',true));
end

function oracle = maintained_oracle()
x = axis_map(linspace(-1,1,65),4,8);
y = linspace(0,2,17).';
alpha = 0.2*(1-exp(-y));
rho = analytic_density(x,y);
oracle = struct('width',{},'dt',{},'maximumAbsDifference',{});
for width = [6,8]
    for dt = [0.004,0.032]
        actual = ipm_perflab_background_remap(rho,x,alpha,dt,['poly',num2str(width)]);
        expected = zeros(size(rho));
        for row = 1:numel(y)
            expected(row,:) = ipm.remesh.interpolate(x,rho(row,:), ...
                x*exp(-alpha(row)*dt),struct('stencilWidth',width));
        end
        defect = max(abs(actual-expected),[],'all');
        assert(defect<1e-12,'Vectorized remap disagreed with maintained interpolation.');
        oracle(end+1) = struct('width',width,'dt',dt,'maximumAbsDifference',defect); %#ok<AGROW>
    end
end
end

function item = invariant_checks(rho,x,alpha,method)
identity = ipm_perflab_background_remap(rho,x,alpha,0,method);
[constant,positiveInfo] = ipm_perflab_background_remap(ones(size(rho)),x,alpha,0.032,method);
[inflow,negativeInfo] = ipm_perflab_background_remap(ones(size(rho)),x,-alpha,0.032,method);
constantError = max(abs(constant-1),[],'all');
coveredInflowError = max(abs(inflow(isfinite(inflow))-1));
item = struct('method',method,'identityExact',isequaln(identity,rho), ...
    'constantMaximumAbsError',constantError, ...
    'outflowForwardMissingNodes',positiveInfo.outsideDomainNodes, ...
    'inflowReverseMissingNodes',negativeInfo.outsideDomainNodes, ...
    'missingNodesAreNaN',nnz(isnan(inflow))==negativeInfo.outsideDomainNodes, ...
    'inflowCoveredConstantMaximumAbsError',coveredInflowError, ...
    'wallRowsUnchanged',positiveInfo.wallRowsUnchanged && negativeInfo.wallRowsUnchanged);
assert(item.identityExact && constantError<1e-12 && coveredInflowError<1e-12 && ...
    item.missingNodesAreNaN && item.wallRowsUnchanged,'Remap invariants failed.');
end

function [actual,info,samples] = timed_remap(rho,x,alpha,dt,method)
[actual,info] = ipm_perflab_background_remap(rho,x,alpha,dt,method);
samples = zeros(1,3);
samples(1) = info.seconds;
for repeat = 2:3
    [~,timing] = ipm_perflab_background_remap(rho,x,alpha,dt,method);
    samples(repeat) = timing.seconds;
end
end

function metrics = pair_metrics(actual,reference,Dx,weights,regions)
metrics = struct('rho',scalar_metrics(actual,reference,weights,regions), ...
    'rhoX',scalar_metrics(actual*Dx.',reference*Dx.',weights,regions));
end

function metrics = scalar_metrics(actual,reference,weights,regions)
metrics = struct();
names = fieldnames(regions);
for index = 1:numel(names)
    name = names{index};
    requested = regions.(name);
    mask = requested & isfinite(actual) & isfinite(reference);
    assert(nnz(mask)>=9,'Too few covered comparison nodes.');
    delta = actual(mask)-reference(mask);
    ref = reference(mask);
    metrics.(name) = struct('points',nnz(mask),'uncoveredRequestedPoints',nnz(requested&~mask), ...
        'relativeL2',sqrt(sum(delta.^2.*weights(mask))/max(sum(ref.^2.*weights(mask)),eps)), ...
        'relativeInf',max(abs(delta))/max(max(abs(ref)),eps), ...
        'absoluteInf',max(abs(delta)));
end
end

function defect = even_defect(rho)
reflected = fliplr(rho);
mask = isfinite(rho) & isfinite(reflected);
defect = max(abs(rho(mask)-reflected(mask)))/max(max(abs(rho(mask))),eps);
end

function axis = axis_map(s,stretch,H)
if stretch == 0
    axis = H*s;
else
    axis = H*sinh(stretch*s)/sinh(stretch);
end
end

function rho = analytic_density(X,Y)
rho = -exp(-(X/1.1).^2-(Y/0.7).^2);
end

function rhoX = analytic_derivative(X,Y)
rhoX = 2*X/1.1^2.*exp(-(X/1.1).^2-(Y/0.7).^2);
end

function width = control_width(axis)
width = zeros(size(axis));
width([1,end]) = [axis(2)-axis(1),axis(end)-axis(end-1)]/2;
width(2:end-1) = (axis(3:end)-axis(1:end-2))/2;
end

function value = rmfield_if_present(value,names)
present = names(isfield(value,names));
if ~isempty(present)
    value = rmfield(value,present);
end
end
