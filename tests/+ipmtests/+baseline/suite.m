function report = suite()
%IPMTESTS.BASELINE.SUITE Run the complete layered verification suite.

fprintf('Running layered IPM verification...\n');
report = struct();
report.config = ipmtests.baseline.config();
report.serverInterface = ipmtests.baseline.serverInterface();
report.grid = ipmtests.baseline.grid();
report.elliptic = ipmtests.baseline.elliptic();
report.numericCore = ipmtests.baseline.numericCore();
report.transport = ipmtests.baseline.transport();
report.symmetry = ipmtests.baseline.symmetry();
report.scaling = ipmtests.baseline.scaling();
report.remesh = ipmtests.baseline.remesh();
report.smoke = ipmtests.baseline.smoke();
report.output = ipmtests.baseline.output();
report.checkpoint = ipmtests.baseline.checkpoint();
report.integration = ipmtests.baseline.integration();
fprintf('Layered IPM verification passed.\n');
end
