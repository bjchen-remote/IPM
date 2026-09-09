function [report,details] = ipm_accellab_observer_sampling(fields,x,y,Dx,Dy,user)
%IPM_ACCELLAB_OBSERVER_SAMPLING No-LU observer quadrature/restriction study.
%   Fixed Omega, F_X, P/PPrime and coordinate rates. This is not PDE-grid
%   refinement. Every registered level and phase is reported without choice.
if nargin < 6, user = struct(); end
opts = struct('observerFactors',[1,2,4],'phases',[0,0.5],'nativeStrides',[1,2,4]);
names = fieldnames(user);
assert(all(ismember(names,fieldnames(opts))),'ipm:SamplingOptions','Unknown option.');
for k = 1:numel(names), opts.(names{k}) = user.(names{k}); end
validateattributes(opts.observerFactors,{'numeric'},{'vector','integer','positive','finite'});
validateattributes(opts.nativeStrides,{'numeric'},{'vector','integer','positive','finite'});
assert(all(ismember(opts.phases,[0,0.5])) && iscell(fields) && numel(fields) == 2, ...
    'ipm:SamplingProtocol','Two fixed fields and registered zero/half-step phases are required.');
branches = cell(1,2); details = cell(1,2);
for branch = 1:2
    data = fields{branch}; exact = data.exactCoordinates;
    [qx,wx] = gauss_axis((x-exact.peak.x)/exact.wallCoreWidth,[-2,-1,0,1,2]);
    [qy,wy] = gauss_axis(y/exact.verticalCoreWidth,[0,1.5,3]);
    [QX,QY] = meshgrid(qx,qy); qw = wy(:)*wx(:)';
    qcore = abs(QX) < 1 & QY < 1.5;
    assert(abs(sum(qw(qcore))-3) < 1e-11 && abs(sum(qw(~qcore))-9) < 1e-11, ...
        'ipm:SamplingQuadratureArea','Split-cell quadrature does not cover the registered geometric regions.');
    refs = pair(data,x,y,Dx,Dy,qx,qy);
    refMetrics = struct('bilinear',l2_metrics(refs.bilinear.G,qcore,qw), ...
        'hermite',l2_metrics(refs.hermite.G,qcore,qw));
    samples = cell(1,numel(opts.observerFactors)*numel(opts.phases)); counter = 0;
    for factor = opts.observerFactors
        for phase = opts.phases
            sx = sample_axis(-2,2,64*factor+1,phase,[-1,0,1]);
            sy = sample_axis(0,3,48*factor+1,phase,1.5);
            [SX,SY] = meshgrid(sx,sy); sw = trap_weights(sy(:))*trap_weights(sx(:))';
            core = abs(SX) <= 1 & SY <= 1.5;
            charts = pair(data,x,y,Dx,Dy,sx,sy);
            entry = struct('observerFactor',factor,'phase',phase,'nodeCount',[numel(sx),numel(sy)]);
            for method = {'bilinear','hermite'}
                name = method{1}; chart = charts.(name);
                masked = sample_metrics(chart.G,core,sw,SY);
                geometric = geometric_l2(chart.G,sx,sy);
                reference = refMetrics.(name);
                entry.(name) = struct('originalMaskedNodeMetrics',masked, ...
                    'geometricTrapezoidL2',geometric, ...
                    'maskedTrainingL2RelativeError',masked.trainingL2/reference.trainingL2-1, ...
                    'maskedHoldoutL2RelativeError',masked.holdoutL2/reference.holdoutL2-1, ...
                    'geometricTrainingL2RelativeError',geometric.trainingL2/reference.trainingL2-1, ...
                    'geometricHoldoutL2RelativeError',geometric.holdoutL2/reference.holdoutL2-1, ...
                    'minimumMovingCellMargin',chart.minimumMovingCellMargin);
            end
            counter = counter+1; samples{counter} = entry;
        end
    end
    restrictions = cell(1,numel(opts.nativeStrides));
    for k = 1:numel(opts.nativeStrides)
        stride = opts.nativeStrides(k);
        ix = unique([1:stride:numel(x),numel(x)]); iy = unique([1:stride:numel(y),numel(y)]);
        assert(numel(ix) >= 7 && numel(iy) >= 7,'ipm:SamplingRestriction','Restricted derivative stencil is too small.');
        limited = data; limited.omega = data.omega(iy,ix); limited.forcing = data.forcing(iy,ix);
        if stride == 1
            dx = Dx; dy = Dy;
        else
            dx = ipm.mesh.fdMatrix(x(ix),1,7); dy = ipm.mesh.fdMatrix(y(iy),1,7);
        end
        charts = pair(limited,x(ix),y(iy),dx,dy,qx,qy);
        entry = struct('nativeStride',stride,'restrictedNodeCount',[numel(ix),numel(iy)], ...
            'coordinateFunctionalsRecomputed',false,'rhsRecomputed',false);
        for method = {'bilinear','hermite'}
            name = method{1}; current = charts.(name); reference = refs.(name);
            entry.(name) = struct('splitCellL2',l2_metrics(current.G,qcore,qw), ...
                'valueDifferenceFromNative',l2_metrics(current.U-reference.U,qcore,qw), ...
                'derivativeDifferenceFromNative',l2_metrics(current.G-reference.G,qcore,qw));
        end
        restrictions{k} = entry;
    end
    branches{branch} = struct('splitCellGaussReference',refMetrics, ...
        'referenceQuadratureNodeCount',[numel(qx),numel(qy)], ...
        'sampling',{samples},'nativeRestriction',{restrictions}, ...
        'actualPeakPrime',exact.peakPrime,'fixedCoordinates',exact);
    details{branch} = struct('quadratureX',qx,'quadratureY',qy,'quadratureWeights',qw, ...
        'referenceObservers',refs);
end
ratios = cell(size(branches{1}.sampling));
for k = 1:numel(ratios)
    entry = struct('observerFactor',branches{1}.sampling{k}.observerFactor, ...
        'phase',branches{1}.sampling{k}.phase);
    for method = {'bilinear','hermite'}
        name = method{1}; a = branches{1}.sampling{k}.(name); b = branches{2}.sampling{k}.(name);
        r = struct();
        for metric = {'trainingL2','trainingInf','holdoutL2','holdoutInf','wallInf'}
            key = metric{1}; r.(key) = b.originalMaskedNodeMetrics.(key)/ ...
                max(a.originalMaskedNodeMetrics.(key),realmin);
        end
        entry.(name) = r;
    end
    ratios{k} = entry;
end
referenceRatios = struct();
for method = {'bilinear','hermite'}
    name = method{1}; a = branches{1}.splitCellGaussReference.(name); b = branches{2}.splitCellGaussReference.(name);
    referenceRatios.(name) = struct('trainingL2',b.trainingL2/max(a.trainingL2,realmin), ...
        'holdoutL2',b.holdoutL2/max(a.holdoutL2,realmin));
end
report = struct('kind','fixed_field_no_LU_observer_sampling_and_restriction', ...
    'options',opts,'baseline',branches{1},'rejectedCandidate',branches{2}, ...
    'candidateToBaselineSamplingRatios',{ratios},'candidateToBaselineGaussRatios',referenceRatios, ...
    'originalDecision','rejected_unchanged','accepted',false,'newRhsEvaluations',0, ...
    'newPdeSteps',0,'poissonOperatorBuilds',0,'parameterSelection',false, ...
    'quadratureMeaning','Within each native cell G is at most bicubic, so 4-point tensor Gauss integrates G squared through degree six exactly up to roundoff; cell and geometric holdout boundaries are split.', ...
    'supremumMeaning','Infinity metrics are finite sampled maxima; phase/refinement spread is diagnostic, not a certified continuum supremum bound.', ...
    'restrictionMeaning','Native Omega and F_X are restricted while all exact coordinates/rates remain fixed. This tests observer restriction sensitivity, not a new PDE mesh or spatial convergence of the PDE.');
end

function charts = pair(data,x,y,Dx,Dy,xi,eta)
charts = struct('bilinear',ipm_accellab_pullback(data.omega,data.forcing, ...
    x,y,data.exactCoordinates,xi,eta), ...
    'hermite',ipm_accellab_pullback_hermite(data.omega,data.forcing, ...
    x,y,Dx,Dy,data.exactCoordinates,xi,eta));
end

function result = l2_metrics(field,core,weights)
result = struct('trainingL2',sqrt(sum(weights(core).*field(core).^2)/sum(weights(core))), ...
    'holdoutL2',sqrt(sum(weights(~core).*field(~core).^2)/sum(weights(~core))));
end

function result = sample_metrics(field,core,weights,Y)
result = l2_metrics(field,core,weights);
result.trainingInf = max(abs(field(core))); result.holdoutInf = max(abs(field(~core)));
result.wallInf = max(abs(field(Y == 0)));
end

function result = geometric_l2(field,x,y)
wx = trap_weights(x(:)); wy = trap_weights(y(:));
total = sum((wy*wx').*field.^2,'all');
ix = abs(x) <= 1; iy = y <= 1.5;
core = sum((trap_weights(y(iy)')*trap_weights(x(ix))').*field(iy,ix).^2,'all');
result = struct('trainingL2',sqrt(core/3),'holdoutL2',sqrt(max(total-core,0)/9));
end

function axis = sample_axis(a,b,count,phase,breaks)
h = (b-a)/(count-1); values = a+((0:count-1)+phase)*h;
axis = unique([a,b,breaks,values(values > a & values < b)]);
end

function [nodes,weights] = gauss_axis(native,breaks)
a = min(breaks); b = max(breaks);
native = native(:)';
edges = unique([breaks(:)',native(native > a & native < b)]);
edges = edges(:)';
z = [-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7), ...
    sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
w = [(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
allNodes = (edges(1:end-1)'+edges(2:end)')/2+diff(edges)'*z/2;
allWeights = diff(edges)'*w/2;
[nodes,~,map] = unique(allNodes(:)); weights = accumarray(map,allWeights(:));
nodes = nodes(:)'; weights = weights(:)';
end

function weights = trap_weights(axis)
axis = axis(:);
weights = [diff(axis(1:2));axis(3:end)-axis(1:end-2);diff(axis(end-1:end))]/2;
end
