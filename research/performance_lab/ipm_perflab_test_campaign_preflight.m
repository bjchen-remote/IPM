function summary=ipm_perflab_test_campaign_preflight(sourceRoot,historicalRoot,outDir,freshRegistrationFile)
%IPM_PERFLAB_TEST_CAMPAIGN_PREFLIGHT Registered launch-only regression cases.
% Uses only configuration, analytic samples, quadrature and pure axis factories.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
oldDirectory=pwd;oldPath=path;restore=onCleanup(@()restore_environment(oldDirectory,oldPath));
ambient=fullfile(outDir,'ambient_cwd_conflict');mkdir(ambient);
conflictFile=fullfile(ambient,'mesh_general_rounded_axis.m');fid=fopen(conflictFile,'w');assert(fid>=0);
fprintf(fid,'function varargout=mesh_general_rounded_axis(varargin)\nerror(''ipm:PreflightAmbientExecuted'',''Ambient constructor must never execute.'');\nend\n');fclose(fid);
addpath(ambient);cd(ambient);assert(strcmp(which('mesh_general_rounded_axis'),conflictFile));
options=struct('nx',49,'ny',25,'xlim',[-2,2],'ymax',2, ...
    'customX',linspace(-2,2,49),'customY',linspace(0,2,25)', ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','initialCondition','degenerate', ...
    'rescalingMode','dynamic','lengthGauge','transport_anchor', ...
    'cOmegaGauge','wall_omega_quadratic_peak','initialAnalyticRemesh',false, ...
    'saveResults',false,'storeSnapshots',false,'makePlots',false, ...
    'livePlot',false,'writeVideo',false,'verbose',false);
base=struct('sourceRoot',sourceRoot,'outputDirectory',fullfile(outDir,'unused_campaign'), ...
    'startKind','time_zero','entryPoint','ipm.solve','options',options, ...
    'additionalEntries',{{'run_fresh_profile_campaign'}});
goodProbe=struct('label','uniform rounded factory','functionName','mesh_general_rounded_axis', ...
    'arguments',{{1,2,49,1/12,8,.5,2}},'geometry','x');
base.axisProbes=goodProbe;
plans=repmat(base,1,7);labels={'complete_current_source','historical_v1_missing_dependency', ...
    'all_axis_families_fail','axis_quality_failure','existing_output','relative_output','late_driver_is_not_original_t0'};
expected={'','ipm:PreflightDependencies','ipm:PreflightAxisFamiliesExhausted', ...
    'ipm:PreflightAxisQuality','ipm:PreflightOutputExists','ipm:PreflightOutputPath','ipm:PreflightTimeZeroEntry'};
plans(2).sourceRoot=historicalRoot;
plans(3).axisProbes=repmat(goodProbe,1,2);
plans(3).axisProbes(1).label='insufficient budget A';plans(3).axisProbes(1).arguments={1,2,49,1/120,80,.5,2};
plans(3).axisProbes(2).label='insufficient budget B';plans(3).axisProbes(2).arguments={1,2,49,1/120,96,.5,2};
plans(4).options.customY(2)=plans(4).options.customY(2)/2;
plans(5).outputDirectory=outDir;plans(6).outputDirectory='relative_campaign_is_ambiguous';
plans(7).entryPoint='run_fresh_profile_campaign';
if nargin>=4
    registered=jsondecode(fileread(freshRegistrationFile));plans(8)=base;plans(8).startKind='fresh_checkpoint';
    plans(8).entryPoint='run_fresh_profile_campaign';plans(8).axisProbes=[];
    plans(8).bootstrapResultFile=registered.bootstrapResultFile;plans(8).initialDataFile=registered.initialDataFile;
    plans(8).qualificationFile=registered.qualificationFile;
    labels{8}='actual_fresh_3668_native_admission';expected{8}='';
end
rows=struct('label',{},'expectedIdentifier',{},'actualIdentifier',{},'matchedExpectation',{},'reportDirectory',{});
reports=cell(1,numel(plans));profile clear;profile on -timer real
for k=1:numel(plans)
    reportDir=fullfile(outDir,labels{k});r=ipm_perflab_preflight_fresh_campaign(plans(k),reportDir);
    reports{k}=r;actual='';if isfield(r,'failure'),actual=r.failure.identifier;end
    okay=strcmp(actual,expected{k})&&(r.passed==isempty(expected{k}));
    if k==2
        assert(any(strcmp({r.dependencies.unavailable.functionName},'mesh_general_rounded_axis')), ...
            'ipm:PreflightMissingHistoricalEvidence','The actual v1 missing dependency must be named.');
        missing=r.dependencies.unavailable(strcmp({r.dependencies.unavailable.functionName},'mesh_general_rounded_axis'));
        assert(strcmp(missing.resolvedFile,conflictFile),'ipm:PreflightAmbientMask','Ambient fallback must be detected as outside the frozen source.');
    end
    if k==3,assert(numel(r.axisFactories)==2&&~any([r.axisFactories.passed]));end
    assert(strcmp(pwd,ambient)&&strcmp(which('mesh_general_rounded_axis'),conflictFile), ...
        'ipm:PreflightEnvironmentRestore','Preflight must restore the caller cwd and MATLAB path.');
    rows(end+1)=struct('label',labels{k},'expectedIdentifier',expected{k},'actualIdentifier',actual, ...
        'matchedExpectation',okay,'reportDirectory',reportDir); %#ok<AGROW>
end
profile off;executionProfile=profile('info');
files={executionProfile.FunctionTable.FileName};
forbidden={'/+ipm/solve.m','/+ipm/+mesh/build.m','/+ipm/+field/poisson.m', ...
    '/+ipm/+field/velocity.m','/+ipm/+evolve/initialize.m','/+ipm/+evolve/flow.m', ...
    '/+ipm/+output/restoreCheckpoint.m','/+ipm/+output/writeCheckpoint.m'};
called=false(size(forbidden));for k=1:numel(forbidden),called(k)=any(endsWith(files,forbidden{k}));end
summary=struct('kind','no_LU_launch_regression','rows',rows,'allPassed',all([rows.matchedExpectation])&&~any(called), ...
    'forbiddenEntryFiles',{forbidden},'forbiddenEntriesExecuted',called, ...
    'currentSourceFileCount',numel(reports{1}.dependencies.manifest), ...
    'ambientConflictFile',conflictFile,'ambientCwdAndPathIsolationPassed',true, ...
    'scope','Static initial sampling and registered pure-axis constructors, not initial flow, field-transfer admission or long-time stability.');
save(fullfile(outDir,'regression.mat'),'summary','executionProfile','-v7.3');
fid=fopen(fullfile(outDir,'regression.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));
assert(summary.allPassed,'ipm:PreflightRegression','One or more registered launch regressions failed.');
fprintf('CAMPAIGN_PREFLIGHT_REGRESSION all=%d cases=%d noLU=%d\n',summary.allPassed,numel(rows),~any(called));
end
function restore_environment(directory,oldPath)
path(oldPath);cd(directory);
end
