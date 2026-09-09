function report = suite(stabilityMode)
%IPMTESTS.FOURTH.SUITE Run every opt-in high-order numerical gate.
%   REPORT = IPMTESTS.FOURTH.SUITE() runs the spatial/coupled MMS,
%   temporal, unforced physical-convergence, Green, remap, and quick
%   stability certificates. Pass 'heavy' to extend the stability sweep. The
%   maintained second-order IPMTESTS.BASELINE.SUITE suite is intentionally separate so
%   each certificate remains independently useful.

if nargin < 1 || isempty(stabilityMode)
    stabilityMode = 'quick';
end
stabilityMode = validatestring(lower(char(string(stabilityMode))), ...
    {'quick','heavy'},mfilename,'stabilityMode',1);

started = tic;
fprintf('Running complete high-order certificate (%s stability)...\n', ...
    stabilityMode);
report = struct('stabilityMode',stabilityMode);
report.core = ipmtests.fourth.core();
report.time = ipmtests.fourth.time();
report.physicalConvergence = ...
    ipmtests.fourth.physicalConvergence();
report.green = ipmtests.fourth.green();
report.remap = ipmtests.fourth.remap();
report.stability = ipmtests.fourth.stability(stabilityMode);
report.elapsedSeconds = toc(started);
report.passed = report.core.passed && report.time.passed && ...
    report.physicalConvergence.passed && report.green.passed && ...
    report.remap.passed && report.stability.passed;
assert(report.passed,'ipm:HighOrderCertificate', ...
    'At least one high-order verification layer did not pass.');
fprintf('Complete high-order certificate passed in %.3f seconds.\n', ...
    report.elapsedSeconds);
end
