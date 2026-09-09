function report=ipm_accellab_test_t0_tail_preflight(registrationFile)
%IPM_ACCELLAB_TEST_T0_TAIL_PREFLIGHT No-LU preparation and damaged-array guards.
reg=jsondecode(fileread(registrationFile));out=reg.outputDirectory;
assert(maxNumCompThreads==10&&~isfolder(out));mkdir(out);profile clear;profile on;
report=struct('kind','t0_tail_operator_no_lu_preflight_v1','positiveCases',{{}}, ...
    'testOnlyDamagedContexts',{{}},'LUCount',0,'PDECount',0);
for k=1:numel(reg.contextFiles)
    if k==1
        r=ipm_accellab_run_t0_tail_operator_pair(reg.contextFiles{k},fullfile(out,sprintf('positive_%d',k)));
    else
        r=ipm_accellab_run_t0_tail_operator_pair(reg.contextFiles{k},fullfile(out,sprintf('positive_%d',k)),false);
    end
    assert(strcmp(r.status,'prepared_only_no_lu')&&r.LUCount==0&&r.PoissonSolveCount==0);
    report.positiveCases{k}=r;
end
fields={'Omega','mass0','initialNativeCL'};
for k=1:numel(fields)
    loaded=load(reg.contextFiles{1},'context');context=loaded.context;key=fields{k};
    context.(key)(1)=context.(key)(1)+1;context.testOnlyDamagedContext=true;
    path=fullfile(out,['test_only_damaged_',key,'.mat']);save(path,'context','-v7.3');
    rejected=false;exceptionRecord=struct();
    try
        ipm_accellab_run_t0_tail_operator_pair(path,fullfile(out,['rejected_',key]),false);
    catch exception
        rejected=true;exceptionRecord=struct('identifier',exception.identifier,'message',exception.message);
    end
    assert(rejected);report.testOnlyDamagedContexts{k}=struct('field',key,'rejected',rejected,'exception',exceptionRecord);
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for key={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.field.poisson', ...
        'ipm.output.restoreCheckpoint','ipm.evolve.initializeScaling','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,key{1})),['Forbidden call: ',key{1}]);
end
report.allPassed=true;report.noLUCallGraphChecked=true;
save(fullfile(out,'report.mat'),'report','profileInfo');
fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));disp(jsonencode(report));
end
