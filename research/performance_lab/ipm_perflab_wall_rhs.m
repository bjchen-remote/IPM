function [rhs,baseRhs,transportU1] = ipm_perflab_wall_rhs(rhoWall,u1Wall,cX,cOmega,cR,ops)
%IPM_PERFLAB_WALL_RHS Research candidate for the original advective wall.
% Uses native WENO kernels. No Poisson solve and no geometry/projection.
assert(isrow(rhoWall)&&isequal(size(rhoWall),size(u1Wall)));
assert(strcmp(ops.transportScheme,'weno5_fd')&& ...
    strcmp(ops.wallTransportMode,'advective_upwind')&& ...
    any(strcmp(ops.transportBoundaryMode,{'open','closed'})));
if cX==0&&cR==0
    transportU1=u1Wall;
else
    transportU1=u1Wall+cX*ops.X(1,:)+cR;
end
options=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',ops.wenoEpsilon);
metricX=ops.metricX(1,:);
drhoX=ipm.field.weno5FluxDerivative(rhoWall,transportU1,metricX, ...
    ops.computationalSpacingX,options);
unit=ones(size(rhoWall),'like',rhoWall);
divergenceX=ipm.field.weno5FluxDerivative(unit,transportU1,metricX, ...
    ops.computationalSpacingX,options);
baseRhs=-(drhoX-rhoWall.*divergenceX);
if cOmega==0
    rhs=baseRhs;
else
    rhs=baseRhs+cOmega*rhoWall;
end
end
