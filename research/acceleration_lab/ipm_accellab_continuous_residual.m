function [report,details] = ipm_accellab_continuous_residual(omega,forcing,x,y,Dx,Dy,rates)
%IPM_ACCELLAB_CONTINUOUS_RESIDUAL No-fit residual of the consistent pullback.
% Split-cell Gauss norms avoid choosing a favorable observation-node phase.
assert(rates.valid);
[qx,wx] = gauss_axis((x-rates.peak.x)/rates.wallCoreWidth,[-2,-1,0,1,2]);
[qy,wy] = gauss_axis(y/rates.verticalCoreWidth,[0,1.5,3]);
quadrature = evaluate(omega,forcing,x,y,Dx,Dy,rates,qx,qy);
[XI,ETA] = meshgrid(qx,qy); core = abs(XI) < 1 & ETA < 1.5;
weights = wy(:)*wx(:)';
assert(abs(sum(weights(core))-3) < 1e-10 && abs(sum(weights(~core))-9) < 1e-10);
gaussMetrics = metrics(quadrature,core,weights,false);
wall = evaluate(omega,forcing,x,y,Dx,Dy,rates,qx,0);
wallCore = abs(qx) < 1;
wallMetrics = metrics(wall,wallCore,wx,true);
samples = cell(1,2); sampleFields = cell(1,2);
for k = 1:2
    factor = 2^(k-1); sx = linspace(-2,2,64*factor+1); sy = linspace(0,3,48*factor+1);
    values = evaluate(omega,forcing,x,y,Dx,Dy,rates,sx,sy);
    [SX,SY] = meshgrid(sx,sy); mask = abs(SX) <= 1 & SY <= 1.5;
    sw = trap_weights(sy)'*trap_weights(sx);
    sample = metrics(values,mask,sw,false);
    sample.nodeCount = [numel(sx),numel(sy)];
    sample.wallInfinity = max(abs(values.G(1,:)));
    sample.geometricRootAndPeakRatesRecomputed = false;
    samples{k} = sample; sampleFields{k} = values;
end
report = struct('kind','consistent_continuous_geometry_pullback_residual', ...
    'coordinates',rates,'gauss',gaussMetrics,'wallGauss',wallMetrics, ...
    'sampledMetrics',{samples},'quadratureNodeCount',[numel(qx),numel(qy)], ...
    'fixedRectangles',struct('fullLocal',[-2,2,0,3],'core',[-1,1,0,1.5]), ...
    'fullLocalMeaning','The registered inner observation rectangle, not the entire global PDE domain.', ...
    'quadratureMeaning','G is bicubic within each mapped native cell; four-point tensor Gauss integrates its square exactly up to roundoff after native and rectangle boundaries are split.', ...
    'infinityMeaning','Finite sampled maxima only; the two node densities are both reported without selection.', ...
    'interpolantLinearWithFixedQueries',true,'completeGeometryPullbackLinearInField',false, ...
    'pdeSteps',0,'rhsEvaluations',0,'poissonOperatorBuilds',0, ...
    'interpretation','No-fit instantaneous separation of actual normalized field evolution, translation, and two-axis width motion. A small residual at one state is not a convergence rate, fixed point or singularity proof.');
details = struct('quadratureX',qx,'quadratureY',qy,'quadratureWeights',weights, ...
    'quadratureFields',quadrature,'wallFields',wall,'sampleFields',{sampleFields});
end

function values = evaluate(omega,forcing,x,y,Dx,Dy,r,xi,eta)
[XI,ETA] = meshgrid(xi,eta); QX = r.peak.x+r.wallCoreWidth*XI; QY = r.verticalCoreWidth*ETA;
o = ipm_accellab_tensor_hermite(omega,x,y,Dx,Dy,QX,QY);
f = ipm_accellab_tensor_hermite(forcing,x,y,Dx,Dy,QX,QY);
P = r.peak.value; U = o.value/P;
raw = f.value/P; amplitudeNormalized = raw-r.peakPrime/P*U;
translation = r.translationRate*o.derivativeX/P;
horizontal = r.wallCoreWidthPrime*XI.*o.derivativeX/P;
vertical = r.verticalCoreWidthPrime*ETA.*o.derivativeY/P;
values = struct('U',U,'G',amplitudeNormalized+translation+horizontal+vertical, ...
    'rawEulerian',raw,'amplitudeNormalizedEulerian',amplitudeNormalized, ...
    'translationOnly',amplitudeNormalized+translation, ...
    'translationContribution',translation,'horizontalContribution',horizontal,'verticalContribution',vertical, ...
    'xi',xi,'eta',eta,'queryX',QX,'queryY',QY);
end

function report = metrics(v,core,w,isWall)
names = {'G','U','rawEulerian','amplitudeNormalizedEulerian','translationOnly', ...
    'translationContribution','horizontalContribution','verticalContribution'};
report = struct('isWall',isWall);
for k = 1:numel(names)
    key = names{k}; f = v.(key);
    report.(key) = struct('coreL2',sqrt(sum(w(core).*f(core).^2)/sum(w(core))), ...
        'holdoutL2',sqrt(sum(w(~core).*f(~core).^2)/sum(w(~core))), ...
        'fullLocalL2',sqrt(sum(w(:).*f(:).^2)/sum(w(:))), ...
        'coreSampledInfinity',max(abs(f(core))),'holdoutSampledInfinity',max(abs(f(~core))), ...
        'fullLocalSampledInfinity',max(abs(f(:))));
end
for key = {'coreL2','holdoutL2','fullLocalL2'}
    name = key{1}; report.GOverAmplitudeNormalizedRhs.(name) = report.G.(name)/max(report.amplitudeNormalizedEulerian.(name),realmin);
    report.GOverProfile.(name) = report.G.(name)/max(report.U.(name),realmin);
end
end

function [nodes,weights] = gauss_axis(native,breaks)
native = native(:)'; edges = unique([breaks,native(native > min(breaks) & native < max(breaks))]);
z = [-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7), ...
    sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
w = [(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
q = (edges(1:end-1)'+edges(2:end)')/2+diff(edges)'*z/2;
weights0 = diff(edges)'*w/2; [nodes,~,map] = unique(q(:)); weights = accumarray(map,weights0(:));
nodes = nodes(:)'; weights = weights(:)';
end

function w = trap_weights(x)
x = x(:)'; h = diff(x); w = [h(1)/2,(h(1:end-1)+h(2:end))/2,h(end)/2];
end
