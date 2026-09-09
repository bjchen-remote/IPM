function report = ipm_accellab_test_quotient_checkpoints(outputRoot)
%IPM_ACCELLAB_TEST_QUOTIENT_CHECKPOINTS Native validation and one-restore smoke.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['quotient_native_tiny_',token]);
mkdir(destination);
opts = struct('nx',129,'ny',65,'xlim',[-4,4],'ymax',4, ...
    'initialCondition','degenerate_primitive','degeneratePower',8, ...
    'symmetryMode','double_odd_omega','rescalingMode','dynamic', ...
    'lengthGauge','transport_anchor','transportAnchorX',1, ...
    'cOmegaGauge','wall_omega_quadratic_peak', ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'adaptiveRemesh',false,'maxDt',0.002,'finalTime',1,'physicalFinalTime',1, ...
    'saveResults',false,'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false);
state = ipm.evolve.initialize(opts); files = cell(1,5);
for k = 1:5
    if k > 1
        for step = 1:8
            [state,stop] = ipm.evolve.advance(state);
            assert(isempty(stop),'ipm:QuotientCheckpointSmoke','Tiny original history stopped.');
        end
    end
    state.runMetadata.resumeCount = k-1;
    [log,~] = ipm.output.record(ipm.output.initializeLog(),state);
    cursor = struct('nextOutput',state.scale.canonicalTime+0.01,'nextCheckpoint',Inf, ...
        'lastCheckpointStep',state.step,'lastCheckpointCanonicalTime',state.scale.canonicalTime);
    checkpoint = ipm.output.makeCheckpoint(state,log,cursor);
    files{k} = fullfile(destination,sprintf('tiny_native_%02d.mat',k));
    save(files{k},'checkpoint','-v7.3');
end
probe = ipm_accellab_quotient_checkpoints(files,struct( ...
    'outputRoot',destination,'expectedRemeshCount',0));
assert(probe.restoreCalls == 1 && all(probe.signatureValidated) && ...
    ~probe.parameterScan && probe.damping == 0.0625 && probe.trial.accepted, ...
    'ipm:QuotientCheckpointSmoke','The fixed native probe did not reproduce the tiny candidate.');
checkpoint = ipm.output.readCheckpoint(files{1});
checkpoint.payload.state.config.transport.wenoEpsilon = 2e-12;
checkpoint.signature = ipm.output.checkpointSignature(checkpoint.payload);
bad = fullfile(destination,'tiny_changed_numerics.mat');
save(bad,'checkpoint','-v7.3');
badFiles = files; badFiles{1} = bad;
rejected = false;
try
    ipm_accellab_quotient_checkpoints(badFiles,struct( ...
        'outputRoot',destination,'expectedRemeshCount',0));
catch exception
    rejected = strcmp(exception.identifier,'ipm:QuotientCheckpointNumerics');
end
assert(rejected,'ipm:QuotientCheckpointSmoke','A changed numerical operator passed the native preflight.');
report = struct('status','passed','outputDirectory',destination, ...
    'validProbe',probe,'changedNumericsRejected',rejected,'q512Run',false);
fprintf('QUOTIENT NATIVE PASS: one restore, fixed damping, numerical mismatch rejected.\n');
end
