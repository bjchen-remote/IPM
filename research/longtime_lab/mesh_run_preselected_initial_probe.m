function report=mesh_run_preselected_initial_probe(registrationFile,label)
%MESH_RUN_PRESELECTED_INITIAL_PROBE One explicitly preselected t0 native probe.
% Dormant until called by root; at most four accepted steps, no new solver path.
assert(maxNumCompThreads==10);
q=load(registrationFile,'registration','plans');k=find(strcmp({q.plans.label},label));assert(isscalar(k));p=q.plans(k);
assert(~isfolder(p.runDirectory));mkdir(p.runDirectory);timer=tic;
save(fullfile(p.runDirectory,'registration.mat'),'p','-v7.3');
try
    result=ipm.solve(p.options);seconds=toc(timer);
    report=mesh_audit_preselected_initial_probe(registrationFile,label, ...
        char(result.metadata.resultFile),char(result.metadata.latestCheckpointFile));
    report.seconds=seconds;
    save(fullfile(p.runDirectory,'report.mat'),'report','-v7.3');write_json(fullfile(p.runDirectory,'report.json'),report);
    assert(report.passed,'ipm:PreselectedInitialProbe','Native qualification failed; original outputs preserved.');
    fprintf('PRESELECTED_INITIAL_NATIVE label=%s steps=%d core=%.8g/%.8g passed=1\n',label,result.state.steps,report.endpointCore);
catch e
    failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack,'seconds',toc(timer));
    save(fullfile(p.runDirectory,'failure.mat'),'failure');write_json(fullfile(p.runDirectory,'failure.json'),failure);rethrow(e)
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
