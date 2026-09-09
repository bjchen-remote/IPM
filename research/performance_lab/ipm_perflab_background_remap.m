function [remapped,info] = ipm_perflab_background_remap(rho,x,alpha,dt,method)
%IPM_PERFLAB_BACKGROUND_REMAP Frozen advective background subproblem only.
%   rho_t + alpha(y)*x*rho_x = 0 gives q=x*exp(-alpha*dt).
%   Missing inflow data is returned as NaN, never silently extrapolated.
%   poly6/poly8 vectorize the maintained remesh local interpolation formula;
%   no global conservation correction is appropriate for this subproblem.

if nargin < 5
    method = 'poly6';
end
x = x(:).';
alpha = alpha(:);
assert(isequal(size(rho),[numel(alpha),numel(x)]) && ...
    all(isfinite(rho),'all') && all(isfinite(alpha)) && all(diff(x)>0), ...
    'Pass finite paired row signals, increasing x, and one alpha per row.');
validateattributes(dt,{'double'},{'scalar','real','finite'});
timer = tic;
q = exp(-alpha*dt).*x;
outside = q < x(1) | q > x(end);
if dt == 0 || all(alpha == 0)
    remapped = rho;
elseif ismember(method,{'linear','pchip','spline'})
    remapped = zeros(size(rho));
    for row = 1:size(rho,1)
        if alpha(row) == 0
            remapped(row,:) = rho(row,:);
        else
            remapped(row,:) = interp1(x,rho(row,:),q(row,:),method,NaN);
        end
    end
elseif ismember(method,{'poly6','poly8'})
    width = str2double(method(end));
    remapped = local_polynomial(rho,x,q,width);
    remapped(outside) = NaN;
    remapped(alpha==0,:) = rho(alpha==0,:);
else
    error('ipm:BackgroundRemapMethod','Unknown interpolation method.');
end
info = struct('method',method,'dt',dt,'wallRowsUnchanged', ...
    isequaln(remapped(alpha==0,:),rho(alpha==0,:)), ...
    'outsideDomainNodes',nnz(outside),'coveredNodes',nnz(~outside), ...
    'finiteOnCoveredDomain',all(isfinite(remapped(~outside))), ...
    'minimumBackwardScale',min(exp(-alpha*dt)), ...
    'maximumBackwardScale',max(exp(-alpha*dt)), ...
    'seconds',toc(timer),'boundaryPolicy','NaN_for_unprovided_inflow', ...
    'frozenAlpha',true,'pdeOrGaugeRecomputed',false);
end

function values = local_polynomial(rho,x,q,width)
ny = size(rho,1);
nx = numel(x);
assert(nx >= width,'Too few interpolation nodes.');
numberOfStencils = nx-width+1;
indices = (1:numberOfStencils).'+(0:width-1);
centers = 0.5*(x(indices(:,1))+x(indices(:,end))).';
scales = 0.5*(x(indices(:,end))-x(indices(:,1))).';
nodes = (x(indices)-centers)./scales;
weights = zeros(size(nodes));
for local = 1:width
    others = [1:local-1,local+1:width];
    weights(:,local) = 1./prod(nodes(:,local)-nodes(:,others),2);
end
cellIndex = discretize(q(:),[-Inf,x(2:end-1),Inf]);
first = min(max(cellIndex-(width/2-1),1),numberOfStencils);
normalizedQuery = (q(:)-centers(first))./scales(first);
distance = normalizedQuery-nodes(first,:);
[nearestDistance,nearestLocal] = min(abs(distance),[],2);
hit = nearestDistance <= 8*eps(max(1,abs(normalizedQuery)));
% Values on a node use exact source samples; harmless denominators avoid
% manufacturing intermediate Inf/NaN before these rows are replaced.
distance(hit,:) = 1;
factors = weights(first,:)./distance;
rows = repmat((1:ny).',nx,1);
samples = rho(rows+(first+(0:width-1)-1)*ny);
output = sum(factors.*samples,2)./sum(factors,2);
nearestGlobal = rows+(first+nearestLocal-2)*ny;
output(hit) = rho(nearestGlobal(hit));
values = reshape(output,size(q));
end
