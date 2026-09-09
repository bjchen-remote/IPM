function report = suite()
%IPMTESTS.SIXTH.SUITE Run the integrated sixth-order certificate.
%   This suite is intentionally separate from IPMTESTS.BASELINE.SUITE and
%   IPMTESTS.FOURTH.SUITE so selecting the sixth-order feature cannot
%   weaken the independent second- and fourth-order regression gates. It
%   covers complete spatial/elliptic/velocity, WENO7, RK6, remesh-transaction,
%   and unforced physical-convergence contracts within the documented scope.

started = tic;
fprintf('Running integrated sixth-order certificate...\n');
report = struct();
report.core = ipmtests.sixth.core();
report.time = ipmtests.sixth.time();
report.physicalConvergence = ...
    ipmtests.sixth.physicalConvergence();
report.elapsedSeconds = toc(started);
report.passed = report.core.passed && report.time.passed && ...
    report.physicalConvergence.passed;
assert(report.passed,'ipm:SixthOrderCertificate', ...
    'At least one integrated sixth-order verification layer failed.');
fprintf('Integrated sixth-order certificate passed in %.3f seconds.\n', ...
    report.elapsedSeconds);
end
