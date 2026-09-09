function experiment = ipm_accellab_run_q512_fixed()
%IPM_ACCELLAB_RUN_Q512_FIXED Parent-window-only preregistered five-frame probe.
directory = fileparts(mfilename('fullpath'));
project = fileparts(fileparts(directory)); addpath(project,directory);
assert(maxNumCompThreads == 10,'ipm:QuotientNativeThreads', ...
    'Use the original ten-thread MATLAB environment for native signature validation.');
campaign = fullfile(project,'result','longtime','20260908_campaign_v1');
prefix = '20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7';
steps = [2135,2182,2226,2248,2270]; files = cell(1,5);
for k = 1:5
    if k <= 3
        folder = fullfile(campaign,'adaptive_campaign','stage_003');
    else
        folder = fullfile(campaign,'bridge_segment_001');
    end
    tag = 'checkpoint_';
    if k == 1, tag = 'regrid_checkpoint_'; end
    files{k} = fullfile(folder,sprintf('%s%s_step%010d.mat',tag,prefix,steps(k)));
end
fprintf('Q512_FIXED_REGISTRATION memory=4 damping=0.0625 threads=%d no_PDE_no_scan\n',maxNumCompThreads);
experiment = ipm_accellab_quotient_checkpoints(files,struct( ...
    'outputRoot',fullfile(campaign,'acceleration_lab'),'expectedRemeshCount',3));
fprintf('Q512_FIXED_OUTPUT %s\n',experiment.outputDirectory);
end
