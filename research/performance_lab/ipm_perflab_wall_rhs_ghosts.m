function [rhs,baseRhs,transportU1,cache] = ipm_perflab_wall_rhs_ghosts(rho,u1,cX,cOmega,cR,ops,cache)
%IPM_PERFLAB_WALL_RHS_GHOSTS Wall reconstruction retaining full-size ghost GEMM.
% Independent research only. CACHE is valid for this identical rho/grid.
assert(isequal(size(rho),size(u1))&&strcmp(ops.transportScheme,'weno5_fd')&& ...
    strcmp(ops.wallTransportMode,'advective_upwind')&& ...
    any(strcmp(ops.transportBoundaryMode,{'open','closed'})));
if nargin<7
    metric=ops.metricX;[ny,nx]=size(rho);
    if isscalar(metric),metric=repmat(metric,ny,nx);
    elseif isequal(size(metric),[ny,1]),metric=repmat(metric,1,nx);
    elseif isvector(metric)&&numel(metric)==nx,metric=repmat(reshape(metric,1,nx),ny,1);
    end
    metric=cast(metric,'like',rho);assert(isequal(size(metric),size(rho)));
    cache=struct('rhoGhost',ghosts(rho),'unitGhost',ghosts(ones(size(rho),'like',rho)), ...
        'metricGhost',ghosts(metric),'metricWall',metric(1,:));
end
if cX==0&&cR==0
    a=u1;
else
    a=u1+cX*ops.X+cR;
end
transportU1=a(1,:);rhoWall=rho(1,:);
options=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',ops.wenoEpsilon,'wallQGhost',cache.rhoGhost, ...
    'wallAGhost',ghosts(a),'wallJGhost',cache.metricGhost);
drhoX=ipm_perflab_weno5_wall_ghosts(rhoWall,transportU1,cache.metricWall, ...
    ops.computationalSpacingX,options);
options.wallQGhost=cache.unitGhost;
divergenceX=ipm_perflab_weno5_wall_ghosts(ones(size(rhoWall),'like',rhoWall), ...
    transportU1,cache.metricWall,ops.computationalSpacingX,options);
baseRhs=-(drhoX-rhoWall.*divergenceX);
if cOmega==0,rhs=baseRhs;else,rhs=baseRhs+cOmega*rhoWall;end
end
function row=ghosts(values)
% Match both matrix dimensions and the original multiplication expression.
weights=lagrange_weights(0:5,-3:-1);
left=values(:,1:6)*weights.';
weights=lagrange_weights(-5:0,1:3);
right=values(:,end-5:end)*weights.';
row=[left(1,:),right(1,:)];
end
function weights=lagrange_weights(nodes,targets)
numberOfNodes=numel(nodes);numberOfTargets=numel(targets);
weights=ones(numberOfTargets,numberOfNodes);
for node=1:numberOfNodes
    others=[1:node-1,node+1:numberOfNodes];
    weights(:,node)=prod((targets(:)-nodes(others))./(nodes(node)-nodes(others)),2);
end
end
