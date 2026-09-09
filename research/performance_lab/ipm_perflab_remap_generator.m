function generator = ipm_perflab_remap_generator(x,width)
%IPM_PERFLAB_REMAP_GENERATOR Right-time derivatives of inward poly remapping.
% Research only. Alpha must be nonnegative and dt positive. At an original
% node the infinitesimal query uses the inward cell, not a symmetric stencil.
x = x(:).';
assert(width == 6 || width == 8,'Use the existing poly6/poly8 interpolant.');
assert(numel(x) >= width && all(diff(x)>0),'Invalid interpolation axis.');
nx = numel(x);
D1 = zeros(nx,nx);
D2 = zeros(nx,nx);
firstNodes = zeros(1,nx);
for node = 1:nx
    cellIndex = node-double(x(node)>0);
    first = min(max(cellIndex-(width/2-1),1),nx-width+1);
    stencil = first:first+width-1;
    localNode = node-first+1;
    localD1 = ipm.mesh.fdMatrix(x(stencil),1,width);
    localD2 = ipm.mesh.fdMatrix(x(stencil),2,width);
    D1(node,stencil) = localD1(localNode,:);
    D2(node,stencil) = localD2(localNode,:);
    firstNodes(node) = first;
end
generator = struct('x',x,'width',width,'D1',sparse(D1),'D2',sparse(D2), ...
    'firstNodes',firstNodes,'alphaContract','nonnegative', ...
    'interpretation','Right-time first/second derivative of inward interpolation.');
end
