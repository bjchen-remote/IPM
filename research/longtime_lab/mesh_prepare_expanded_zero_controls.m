function report=mesh_prepare_expanded_zero_controls(factorialFile,timeFile,sourceDirectory,outDir)
%MESH_PREPARE_EXPANDED_ZERO_CONTROLS Register five fresh t0 runs at cap210000.
% Preserve old failures. Only the registered node cap changes numerically.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
assert(strcmp(which('ipm.solve'),fullfile(sourceDirectory,'+ipm','solve.m')));
q=load(factorialFile,'registration','plans');t=load(timeFile,'registration','plans');
labels={'A_H8_T32','B_H16_T32','C_H8_T42','D_H16_T42'};
assert(isequal({q.plans.label},labels)&&isscalar(t.plans)&& ...
    strcmp(t.plans.label,'E_H8_T32_HALF_TIME'));
assert(q.registration.maximumTotalNodes==110000&&t.registration.maximumTotalNodes==110000&& ...
    q.registration.physicalFinalTime==t.registration.physicalFinalTime&& ...
    isequaln(q.registration.limits,t.registration.limits));
b=jsondecode(fileread(fullfile(q.plans(2).runDirectory,'endpoint_report.json')));
assert(~b.passed&&~b.physicalEndpointReached&& ...
    strcmp(b.stopReason,'autonomous_mesh_axis_capacity'));
profile clear;profile on;
parentPlans=[q.plans,t.plans];expanded=parentPlans;families=cell(1,5);
for k=1:5
    p=parentPlans(k);o=p.options;
    assert(p.pureGeometryLaunchable&&o.autonomousMesh.nodeFamily.maximumTotalNodes==110000);
    o.autonomousMesh.nodeFamily.maximumTotalNodes=210000;
    o.autonomousMesh=ipm.config.autonomousMeshPolicy(o.autonomousMesh);
    group='factorial';if k==5,group='time';end
    runDir=fullfile(outDir,group,'runs',p.label);
    o.resultFile=fullfile(runDir,'result.mat');o.checkpoint.file=fullfile(runDir,'checkpoint.mat');
    o.caseMetadata.maximumTotalNodes=210000;
    o.caseMetadata.expandedCapacityControl=struct('parentRegistration',factorialFile, ...
        'parentTimeRegistration',timeFile,'originalPhysicalTime',0, ...
        'allCasesRestartFromOriginalDatum',true,'oldFailurePreserved',true);
    config=ipm.config.resolve(o);old=p.config;
    patched=old;patched.remesh.autonomousMesh.nodeFamily.maximumTotalNodes=210000;
    patched.output=config.output;
    assert(isequaln(config,patched),'ipm:ExpandedControl','Unexpected numerical change.');
    assert(isequaln(rmfield(config.output,{'resultFile','checkpoint','caseMetadata'}), ...
        rmfield(old.output,{'resultFile','checkpoint','caseMetadata'}))&& ...
        isequaln(rmfield(config.output.checkpoint,'file'),rmfield(old.output.checkpoint,'file')));
    index=k;if k==5,index=1;end
    candidate=q.plans(index).pureInitialAudit.selectedPureCandidate;
    a=load(fullfile(fileparts(factorialFile),sprintf('%s_candidate_%d.mat',labels{index},candidate)), ...
        'c','family','audit');
    assert(a.audit.passed);
    family=ipm.remesh.referenceAxisFamily(a.c.x,a.c.y,1,o.autonomousMesh);
    expected=a.family;
    for j=1:numel(expected.members),expected.members(j).resourceAdmitted=true;end
    assert(isequaln(family,expected)&&all([family.members.resourceAdmitted]));
    families{k}=family;
    p.options=o;p.config=config;p.runDirectory=runDir;p.reuseExistingBaseline=false;
    p.pureInitialAudit.reuseExistingBaseline=false;
    p.pureInitialAudit.nativeInitializationQualified=false;
    expanded(k)=p;
end
% The time control has exactly the new A numerical domains except its time step.
for group={'grid','physics','elliptic','transport','scaling','remesh','diagnostics'}
    assert(isequaln(expanded(1).config.(group{1}),expanded(5).config.(group{1})));
end
assert(expanded(5).config.time.cfl==expanded(1).config.time.cfl/2&& ...
    expanded(5).config.time.maxDt==expanded(1).config.time.maxDt/2);
registration=q.registration;registration.kind='from_zero_box_space_factorial_cap210000_v2';
registration.sourceDirectory=sourceDirectory;registration.maximumTotalNodes=210000;
registration.maximumNewRuns=4;registration.originalFactorialRegistration=factorialFile;
registration.allCasesRestartFromOriginalDatum=true;registration.oldCapacityFailure=b;
plans=expanded(1:4);mkdir(fullfile(outDir,'factorial'));
factorialNew=fullfile(outDir,'factorial','registration.mat');
save(factorialNew,'registration','plans','-v7.3');write_json(fullfile(outDir,'factorial','registration.json'),registration);
registration=t.registration;registration.kind='from_zero_half_time_cap210000_v2';
registration.sourceDirectory=sourceDirectory;registration.maximumTotalNodes=210000;
registration.baselineFactorialRegistration=factorialNew;
registration.originalTimeRegistration=timeFile;registration.baselineCasePending=true;
plans=expanded(5);mkdir(fullfile(outDir,'time'));
save(fullfile(outDir,'time','registration.mat'),'registration','plans','-v7.3');
write_json(fullfile(outDir,'time','registration.json'),registration);
profile off;pi=profile('info');names={pi.FunctionTable.FunctionName};
assert(~any(ismember(names,{'ipm.mesh.build','ipm.solve','ipm.evolve.flow', ...
    'ipm.evolve.advance','ipm.remesh.transfer','ipm.output.restoreCheckpoint'})));
report=struct('kind','expanded_five_zero_controls_preflight_v1','passed',true, ...
    'registeredNodeCap',210000,'originalFailurePreserved',true,'allFiveAreNewT0Runs',true, ...
    'onlyNumericalConfigChange','autonomousMesh.nodeFamily.maximumTotalNodes:110000->210000', ...
    'initialCandidateAxesAndAllReferenceGeometryExact',true,'allFourMembersNowResourceAdmitted',true, ...
    'allNativeInitializationsPending',true,'physicalFinalTime',q.registration.physicalFinalTime, ...
    'noLU',true,'noPDE',true,'memoryQualification','641x321 previously measured; no universal LU memory bound');
save(fullfile(outDir,'report.mat'),'report','families','pi','-v7.3');write_json(fullfile(outDir,'report.json'),report);
end

function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
