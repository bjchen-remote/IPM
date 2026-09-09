function [rhoNew,scaleNew,info] = ipm_perflab_lawson_remap_step( ...
        rho,scale,ops,dt,generator,variant,alphaMode)
%IPM_PERFLAB_LAWSON_REMAP_STEP Research Lawson midpoint with full IPM RHS.
% No solver integration, checkpoint, gauge projection, or CFL acceptance.
% B is held fixed within a step; F is recomputed at the actual midpoint.
% N=F-B*rho retains full-speed WENO and its nonlinear flux splitting.
% corrected uses R(h)+h^2/2*(B^2-R''(0)), giving an O(h^3) exponential
% approximation at fixed space while raw retains the remap semigroup defect.
assert(ops.nx <= 65 && ops.ny <= 65,'This probe is restricted to tiny grids.');
assert(strcmp(ops.dynamicScaleGeometry,'isotropic'),'Isotropic probe only.');
assert(ismember(variant,{'raw','corrected','generator_expm'}),'Unknown variant.');
[f0,flow0] = ipm.evolve.flow(rho,ops,scale);
cache0 = ipm.evolve.makeRhsCache(rho,f0,flow0,ops,scale);
z = pack_scale(scale);
if strcmp(alphaMode,'anchor')
    anchor = ops.rescaling.transportAnchorX;
    alpha = flow0.c_l+interp1(ops.x,flow0.u1.',anchor,'linear').'/anchor;
elseif strcmp(alphaMode,'prescribed_positive')
    alpha = 0.2*(1-exp(-(ops.y(:)/0.5).^2));
else
    error('ipm:LawsonAlpha','Unknown background selection.');
end
assert(all(isfinite(alpha)) && all(alpha>=0), ...
    'Negative alpha requires supplied inflow data and a different generator.');
b0 = apply_b(rho,alpha,generator);
n0 = f0-b0;
rhoMid = apply_e(rho+0.5*dt*n0,0.5*dt,alpha,generator,variant);
scaleMid = unpack_scale(z+0.5*dt*cache0.scaleRate);
[fMid,flowMid] = ipm.evolve.flow(rhoMid,ops,scaleMid);
cacheMid = ipm.evolve.makeRhsCache(rhoMid,fMid,flowMid,ops,scaleMid);
nMid = fMid-apply_b(rhoMid,alpha,generator);
rhoNew = apply_e(rho,dt,alpha,generator,variant)+ ...
    dt*apply_e(nMid,0.5*dt,alpha,generator,variant);
scaleNew = unpack_scale(z+dt*cacheMid.scaleRate);
scaleNew.canonicalTime = scale.canonicalTime+dt;
info = struct('variant',variant,'alphaMode',alphaMode, ...
    'alphaMinimum',min(alpha),'alphaMaximum',max(alpha), ...
    'fullNonlinearRhsEvaluations',2,'midpointVelocityRecomputed',true, ...
    'midpointGaugeRecomputed',true,'midpointScaleRecomputed',true, ...
    'initialFlow',flow0,'midpointFlow',flowMid, ...
    'initialAlgebraicSplitDefect',max(abs(b0+n0-f0),[],'all'), ...
    'dt',dt,'cflStabilityCertified',false);
end

function value = apply_b(rho,alpha,g)
value = -alpha.*(rho*g.D1.').*g.x;
end

function value = apply_e(rho,h,alpha,g,variant)
if strcmp(variant,'generator_expm')
    value = zeros(size(rho));
    B0 = -diag(g.x)*full(g.D1);
    for row = 1:size(rho,1)
        value(row,:) = (expm(h*alpha(row)*B0)*rho(row,:).').';
    end
    return;
end
value = ipm_perflab_background_remap(rho,g.x,alpha,h, ...
    sprintf('poly%d',g.width));
if strcmp(variant,'corrected')
    first = rho*g.D1.';
    remapSecond = alpha.^2.*(first.*g.x+(rho*g.D2.').*g.x.^2);
    generatorSecond = apply_b(apply_b(rho,alpha,g),alpha,g);
    value = value+0.5*h^2*(generatorSecond-remapSecond);
end
end

function z = pack_scale(s)
z = [s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime];
end

function s = unpack_scale(z)
s = struct('logC_l',z(1),'logC_omega',z(2), ...
    'physicalTime',z(3),'X_shift',z(4),'canonicalTime',z(5));
end
