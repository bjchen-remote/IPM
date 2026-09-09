function report=mesh_run_box_space_case(registrationFile,label)
%MESH_RUN_BOX_SPACE_CASE Explicitly run one previously registered original-t0 case.
% This runner is dormant until called by the root queue; never dispatches peers.
q=load(registrationFile,'registration','plans');r=q.registration;
assert(maxNumCompThreads==10);
assert(strcmp(which('ipm.solve'),fullfile(r.sourceDirectory,'+ipm','solve.m')));
assert(strcmp(which(mfilename),fullfile(r.sourceDirectory,'research',[mfilename,'.m'])));
mesh_verify_factorial_frozen_source(r.sourceDirectory);
i=find(strcmp({q.plans.label},label));assert(isscalar(i));p=q.plans(i);
assert(~p.reuseExistingBaseline,'A is immutable reused evidence; do not rerun it here.');
assert(p.pureGeometryLaunchable,'Initial geometry failed; no solver call is permitted.');
assert(~isfolder(p.runDirectory),'A new, empty run directory is required.');mkdir(p.runDirectory);
save(fullfile(p.runDirectory,'registration.mat'),'p','r','-v7.3');timer=tic;
try
    result=ipm.solve(p.options);wallSeconds=toc(timer);
    resultFile=char(result.metadata.resultFile);checkpointFile=char(result.metadata.latestCheckpointFile);
    clear result
    report=mesh_audit_box_space_case(resultFile,checkpointFile,r.physicalFinalTime,p.config);
    report.wallSeconds=wallSeconds;report.factorialLabel=label;
    save(fullfile(p.runDirectory,'endpoint_report.mat'),'report','-v7.3');write_json(fullfile(p.runDirectory,'endpoint_report.json'),report);
    assert(report.passed,'ipm:BoxSpaceEndpointRejected','The actual endpoint or original gates failed; preserved without promotion.');
    mesh_verify_factorial_frozen_source(r.sourceDirectory);
    fprintf('BOX_SPACE_CASE_COMPLETE label=%s step=%d t=%.17g seconds=%.3f passed=1\n',label,report.step,report.actualPhysicalTime,wallSeconds);
catch e
    failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack,'wallSeconds',toc(timer));
    save(fullfile(p.runDirectory,'failure.mat'),'failure','-v7.3');write_json(fullfile(p.runDirectory,'failure.json'),failure);rethrow(e)
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
