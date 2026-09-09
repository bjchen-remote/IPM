function [rhoNew,zNew,cache,trace] = ipm_accellab_peak_step(rho,z,dt,ops,mode,cache)
%IPM_ACCELLAB_PEAK_STEP Independent tiny SSPRK54 with identical stage algebra.
% Native-mode parity is tested against maintained stepSsprk54 before use.
if nargin < 6 || isempty(cache)
    cache = ipm_accellab_peak_rhs(rho,z,ops,mode);
else
    assert(isequaln(cache.rho,rho) && isequaln(cache.z,z) && strcmp(cache.mode,mode), ...
        'ipm:PeakResearchCache','First-stage cache must match the exact field, scale and research rule.');
end
tableau = ipm.evolve.ssprk54Tableau(); A = tableau.A; b = tableau.b;
k1 = cache.rhoRate; s1 = cache.scaleRate; trace = cell(6,1); trace{1} = cache.diagnostic;
q2 = rho+dt*A(2,1)*k1; z2 = z+dt*A(2,1)*s1;
e = ipm_accellab_peak_rhs(q2,z2,ops,mode); k2 = e.rhoRate; s2 = e.scaleRate; trace{2} = e.diagnostic;
q3 = rho+dt*(A(3,1)*k1+A(3,2)*k2); z3 = z+dt*(A(3,1)*s1+A(3,2)*s2);
e = ipm_accellab_peak_rhs(q3,z3,ops,mode); k3 = e.rhoRate; s3 = e.scaleRate; trace{3} = e.diagnostic;
q4 = rho+dt*(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3);
z4 = z+dt*(A(4,1)*s1+A(4,2)*s2+A(4,3)*s3);
e = ipm_accellab_peak_rhs(q4,z4,ops,mode); k4 = e.rhoRate; s4 = e.scaleRate; trace{4} = e.diagnostic;
q5 = rho+dt*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4);
z5 = z+dt*(A(5,1)*s1+A(5,2)*s2+A(5,3)*s3+A(5,4)*s4);
e = ipm_accellab_peak_rhs(q5,z5,ops,mode); k5 = e.rhoRate; s5 = e.scaleRate; trace{5} = e.diagnostic;
rhoNew = rho+dt*(b(1)*k1+b(2)*k2+b(3)*k3+b(4)*k4+b(5)*k5);
zNew = z+dt*(b(1)*s1+b(2)*s2+b(3)*s3+b(4)*s4+b(5)*s5);
zNew(5) = z(5)+dt;
cache = ipm_accellab_peak_rhs(rhoNew,zNew,ops,mode); trace{6} = cache.diagnostic;
end
