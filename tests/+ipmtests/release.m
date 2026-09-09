function report = release(stabilityMode)
%IPMTESTS.RELEASE Run maintained second-order and fourth-order certificates.
%   REPORT = IPMTESTS.RELEASE() uses the quick stability campaign. Pass
%   'heavy' for the extended fourth-order CFL and long-evolution gates.

if nargin < 1 || isempty(stabilityMode)
    stabilityMode = 'quick';
end
stabilityMode = validatestring(lower(char(string(stabilityMode))), ...
    {'quick','heavy'},mfilename,'stabilityMode',1);

started = tic;
fprintf('Running combined IPM release certificate (%s)...\n',stabilityMode);
report = struct('stabilityMode',stabilityMode);
report.secondOrder = ipmtests.baseline.suite();
report.highOrder = ipmtests.fourth.suite(stabilityMode);
report.elapsedSeconds = toc(started);
report.passed = report.highOrder.passed;
fprintf('Combined IPM release certificate passed in %.3f seconds.\n', ...
    report.elapsedSeconds);
end
