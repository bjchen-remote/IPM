function [F,details]=ipm_accellab_local_lf_rhs(rho,flow,rates,ops,alphaMode)
%IPM_ACCELLAB_LOCAL_LF_RHS Research full 2D assembly; original open boundary.
% No Poisson. Line mode must reproduce maintained assembleRhs bitwise before
% any face-mode result can be used. Conservative wall y transport is kept.
if nargin<5,alphaMode='face';end
assert(strcmp(ops.spatialDiscretization,'high_order') && strcmp(ops.transportScheme,'weno5_fd') && ...
    strcmp(ops.transportBoundaryMode,'open') && ismember(ops.wallTransportMode,{'conservative_flux','advective_upwind'}));
cl=rates(1);cw=rates(2);cr=rates(3);
if cl==0 && cr==0,u1=flow.u1;else,u1=flow.u1+cl*ops.X+cr;end
if cl==0,u2=flow.u2;else,u2=flow.u2+cl*ops.Y;end
options=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',ops.wenoEpsilon,'alphaMode',alphaMode,'normalization','conditional_global');
[dx,qx]=ipm_accellab_local_lf_derivative(rho,u1,ops.metricX,ops.computationalSpacingX,options);
[dy,qy]=ipm_accellab_local_lf_derivative(rho',u2',ops.metricY',ops.computationalSpacingY,options);dy=dy';
unit=ones(size(rho),'like',rho);
[divx,vx]=ipm_accellab_local_lf_derivative(unit,u1,ops.metricX,ops.computationalSpacingX,options);
[divy,vy]=ipm_accellab_local_lf_derivative(unit',u2',ops.metricY',ops.computationalSpacingY,options);divy=divy';
base=-(dx+dy-rho.*(divx+divy));
if strcmp(ops.wallTransportMode,'advective_upwind')
    base(1,:)=-(dx(1,:)-rho(1,:).*divx(1,:));
end
if cw==0,F=base;else,F=base+cw*rho;end
details=struct('kind','independent_full_2D_research_transport','alphaMode',alphaMode,'normalization','conditional_global', ...
    'sourceWallMode',ops.wallTransportMode,'boundaryMode',ops.transportBoundaryMode,'rates',rates, ...
    'transportU1',u1,'transportU2',u2,'xDensityFluxDerivative',dx,'yDensityFluxDerivative',dy, ...
    'xVelocityDivergence',divx,'yVelocityDivergence',divy,'xDensityFlux',qx,'yDensityFlux',qy, ...
    'xUnitFlux',vx,'yUnitFlux',vy,'baseRhs',base, ...
    'xAdvectiveContribution',-(dx-rho.*divx),'yAdvectiveContribution',-(dy-rho.*divy), ...
    'interpretation','The free-stream correction, native ghost extrapolation and selected original wall mode are preserved. No closed-boundary mass projection is silently substituted. This is a different spatial operator only in face mode, not a new physical trajectory.');
end
