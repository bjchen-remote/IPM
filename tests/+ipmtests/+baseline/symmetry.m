function report = symmetry()
%IPMTESTS.BASELINE.SYMMETRY Verify double-odd symmetry and the zero-origin frame.

fprintf('Running IPM symmetry verification...\n');
config = ipm.config.resolve(struct('nx',65,'ny',33,'xlim',[-4,4], ...
    'ymax',4,'initialCondition','degenerate','degeneratePower',4, ...
    'symmetryMode','double_odd_omega', ...
    'rescalingMode','dynamic','saveResults',false, ...
    'makePlots',false,'verbose',false, ...
    'lengthGauge','omega_peak_location'));
ops = ipm.mesh.build(config);
rho = ipm.field.initialDensity(ops,config.physics);
ops = ipm.evolve.initializeScaling(rho,ops);
[rhs,flow] = ipm.evolve.rhs(rho,ops);
[~,~,boundary] = ipm.field.poisson(flow.source,ops);

rhoDefect = max(abs(rho-fliplr(rho)),[],'all');
omegaDefect = max(abs(flow.source+fliplr(flow.source)),[],'all');
psiDefect = max(abs(flow.psi+fliplr(flow.psi)),[],'all');
u1Defect = max(abs(flow.u1+fliplr(flow.u1)),[],'all');
u2Defect = max(abs(flow.u2-fliplr(flow.u2)),[],'all');
rhsDefect = max(abs(rhs-fliplr(rhs)),[],'all');
boundaryDefect = max([ ...
    max(abs(boundary.left+boundary.right)), ...
    max(abs(boundary.top+fliplr(boundary.top)))]);
originVelocity = hypot(flow.u1(1,ops.rescaling.originIndex), ...
    flow.u2(1,ops.rescaling.originIndex));
scale = struct('logC_l',0,'logC_omega',0,'physicalTime',0, ...
    'X_shift',0,'canonicalTime',0);
[~,~,scaleNew] = ipm.evolve.stepSsprk3(rho,1e-4,ops,scale);

report.symmetricRhoDefect = rhoDefect;
report.symmetricOmegaDefect = omegaDefect;
report.symmetricPsiDefect = psiDefect;
report.symmetricU1Defect = u1Defect;
report.symmetricU2Defect = u2Defect;
report.symmetricRhsDefect = rhsDefect;
report.symmetricBoundaryDefect = boundaryDefect;
report.symmetricOriginVelocity = originVelocity;
report.symmetricCR = flow.c_r;
report.symmetricLengthGaugeResidual = flow.lengthGaugeResidual;
report.symmetricXShiftIncrement = scaleNew.X_shift;

assert(max([rhoDefect,omegaDefect,psiDefect, ...
    u1Defect,u2Defect,rhsDefect,boundaryDefect,originVelocity, ...
    abs(flow.c_r),abs(scaleNew.X_shift), ...
    abs(flow.lengthGaugeResidual)]) < 1e-11, ...
    'ipm:DoubleOddSymmetry', ...
    'The double-odd omega symmetry or zero-origin frame was not preserved.');

fprintf(['  rho/omega/psi/u1/u2/RHS defects: ' ...
    '%.1e %.1e %.1e %.1e %.1e %.1e\n'], ...
    rhoDefect,omegaDefect,psiDefect,u1Defect,u2Defect,rhsDefect);
fprintf('  origin speed/c_r/X_shift: %.1e %.1e %.1e\n', ...
    originVelocity,flow.c_r,scaleNew.X_shift);
fprintf('IPM symmetry verification passed.\n');
end
