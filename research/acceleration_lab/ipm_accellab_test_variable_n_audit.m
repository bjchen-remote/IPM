function report=ipm_accellab_test_variable_n_audit(configFixture,outputDirectory)
%IPM_ACCELLAB_TEST_VARIABLE_N_AUDIT Pure state-view and early-interface tests.
% Polynomial samples and any ledger below are synthetic unit fixtures, not
% PDE trajectories, native checkpoints, or accepted numerical transactions.
assert(~isfolder(outputDirectory));mkdir(outputDirectory);
q=load(configFixture,'result');config=q.result.config;clear q
policy=ipm.config.autonomousMeshPolicy(struct('version',2, ...
    'nodeFamily',struct('maximumTotalNodes',110000)));
config.remesh.autonomousMesh=policy;
config.grid.nx=321;config.grid.ny=161;config.grid.xlim=[-4,4];config.grid.ymax=4;
config.scaling.transportAnchorX=1;
x=linspace(-4,4,321);y=linspace(0,4,161)';
family=ipm.remesh.referenceAxisFamily(x,y,1,policy);
registration=struct('kind','synthetic_variable_node_audit_and_transfer_guard_v1', ...
    'configFixture',configFixture,'originalNodeCount',[321,161], ...
    'maximumTotalNodes',110000,'syntheticPolynomialField',true,'syntheticLedger',true, ...
    'nativeTransactionClaim',false,'noPDE',true,'noLU',true, ...
    'transferImplementation',which('ipm.remesh.transfer'), ...
    'auditImplementation',which('ipm.remesh.auditCandidate'));
write_json(fullfile(outputDirectory,'registration.json'),registration);
base=make_view(config,family,1);
original=base;cases={};
profile clear;profile on
for target=1:4
    c=target_view(base,family.members(target));
    a=ipm.remesh.auditCandidate(base,c,policy,0);
    if target<4
        passed=a.passed && a.version==2 && a.sourceLevelId==1 && a.targetLevelId==target && ...
            a.sameBox && a.registeredNodeCounts && a.matchingReferenceMembers && a.admittedLevelTransition && ...
            a.sameBoxAndNodeCount==(target==1);
    else
        passed=~a.passed && any(strcmp(a.reasons,'unadmitted_node_family_transition'));
    end
    record(sprintf('registered_level_1_to_%d',target),passed,a);
    if target==2,a12=a;c2=c;end
end
assert(a12.passed,'ipm:AccellabFixture','The polynomial growth fixture must pass original numerical gates.');
% Test the next source membership without pretending to have advanced PDE.
source2=c2;source2.ops.remeshCount=1;
source2.runMetadata.autonomousMesh.transactions=a12;
source2.runMetadata.autonomousMesh.currentLevelId=2;
source2.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump=a12.relativePeakJump;
same2=ipm.remesh.auditCandidate(source2,source2,policy,a12.relativePeakJump);
record('registered_level_2_same_N',same2.passed && same2.sameBoxAndNodeCount && ...
    same2.sourceLevelId==2 && same2.targetLevelId==2,same2);
for target=[1,3,4]
    c=target_view(source2,family.members(target));
    a=ipm.remesh.auditCandidate(source2,c,policy,a12.relativePeakJump);
    record(sprintf('level_2_to_%d_rejected',target),~a.passed && ...
        any(strcmp(a.reasons,'unadmitted_node_family_transition')),a);
end
keys={'candidate_metadata','candidate_config','candidate_clock','candidate_mass0', ...
    'candidate_range0','premature_count','target_unregistered_base','source_base_mismatch', ...
    'source_level_mismatch','source_family_forged','source_initial_root_changed', ...
    'source_config_N_changed','bad_prior_budget','ledger_missing','ledger_wrong_factors', ...
    'ledger_budget_mismatch','target_weights_nan','target_Dx_inf','target_box', ...
    'target_anchor','target_field_complex','source_actual_box_mismatch'};
for k=1:numel(keys)
    b=base;c=target_view(base,family.members(2));prior=0;key=keys{k};
    switch key
        case 'candidate_metadata',c.runMetadata.testOnlyMutation=true;
        case 'candidate_config',c.config.grid.nx=c.config.grid.nx+2;
        case 'candidate_clock',c.scale.physicalTime=1;
        case 'candidate_mass0',c.mass0=c.mass0+1;
        case 'candidate_range0',c.rhoRange0(2)=c.rhoRange0(2)+1;
        case 'premature_count',c.ops.remeshCount=1;
        case 'target_unregistered_base',c.ops.baseX(2)=c.ops.baseX(2)+eps(c.ops.baseX(2));
        case 'source_base_mismatch',b.ops.baseX(2)=b.ops.baseX(2)+eps(b.ops.baseX(2));
        case 'source_level_mismatch'
            b.runMetadata.autonomousMesh.currentLevelId=2;c.runMetadata=b.runMetadata;
        case 'source_family_forged'
            b.runMetadata.autonomousMesh.referenceFamily.members(4).resourceAdmitted=true;c.runMetadata=b.runMetadata;
        case 'source_initial_root_changed'
            b.runMetadata.autonomousMesh.initialization.selectedBaseX(2)=0;c.runMetadata=b.runMetadata;
        case 'source_config_N_changed',b.config.grid.nx=641;c.config=b.config;
        case 'bad_prior_budget',prior=1e-4;
        case 'ledger_missing'
            b=source2;c=b;b.runMetadata.autonomousMesh.transactions=struct([]);c.runMetadata=b.runMetadata;
            prior=b.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump;
        case 'ledger_wrong_factors'
            b=source2;b.runMetadata.autonomousMesh.transactions.targetCellFactors=[1,1];c=b;
            prior=b.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump;
        case 'ledger_budget_mismatch'
            b=source2;b.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump=1e-4;c=b;prior=1e-4;
        case 'target_weights_nan',c.ops.integrationWeights(1,1)=NaN;
        case 'target_Dx_inf',c.ops.Dx(1,1)=Inf;
        case 'target_box',c.ops.x([1,end])=[-5,5];
        case 'target_anchor'
            j=find(c.ops.x==1);c.ops.x(j)=1+eps(1);c.ops.x(end+1-j)=-c.ops.x(j);
        case 'target_field_complex',c.rho(1,1)=c.rho(1,1)+1i;
        case 'source_actual_box_mismatch',b.ops.x([1,end])=[-5,5];c.ops.x([1,end])=[-5,5];
    end
    o=try_audit(b,c,policy,prior);record(['audit_reject_',key],~o.passed,o);
end
% V1 result structure and every arithmetic output remain exact on this
% independent same-N fixture, including nonzero defect/range tests.
legacy=base;legacy.config.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct());
legacy.runMetadata=struct('syntheticV1',true);p1=legacy.config.remesh.autonomousMesh;
for amplitude=[1,1+1e-4]
    c=legacy;c.rho=amplitude*c.rho;
    actual=ipm.remesh.auditCandidate(legacy,c,p1,0);
    expected=ipm_accellab_audit_candidate_v1_frozen(legacy,c,p1,0);
    record(sprintf('v1_frozen_exact_%.6g',amplitude),isequaln(actual,expected), ...
        struct('exact',isequaln(actual,expected),'actual',actual));
end
% Every transfer call is an intentionally invalid contract. No accepted
% transfer is invoked in this no-LU suite; real transfer is separately gated.
keys={'missing_registration','forged_family','source_level','target_level','id_type', ...
    'config_initial_N','resource_cap','target_count','source_base','box', ...
    'anchor','nonmonotone_axis','coarsening','cross_level','rho_shape','source_box'};
for k=1:numel(keys)
    b=base;cfg=config;proposal=struct('x',family.members(2).baseX,'y',family.members(2).baseY, ...
        'referenceFamily',family,'sourceLevelId',1,'targetLevelId',2);key=keys{k};
    switch key
        case 'missing_registration',proposal=rmfield(proposal,'referenceFamily');
        case 'forged_family',proposal.referenceFamily.members(4).resourceAdmitted=true;
        case 'source_level',proposal.sourceLevelId=2;
        case 'target_level',proposal.targetLevelId=5;
        case 'id_type',proposal.sourceLevelId=true;
        case 'config_initial_N',cfg.grid.nx=641;
        case 'resource_cap'
            proposal.targetLevelId=4;proposal.x=family.members(4).baseX;proposal.y=family.members(4).baseY;
        case 'target_count',proposal.targetLevelId=3;
        case 'source_base',b.ops.baseX(2)=0;
        case 'box',proposal.x([1,end])=[-5,5];
        case 'anchor'
            j=find(proposal.x==1);proposal.x(j)=1+eps(1);proposal.x(end+1-j)=-proposal.x(j);
        case 'nonmonotone_axis',proposal.x(2)=proposal.x(1);
        case 'coarsening'
            b=source2;proposal.sourceLevelId=2;proposal.targetLevelId=1;
            proposal.x=family.members(1).baseX;proposal.y=family.members(1).baseY;
        case 'cross_level'
            b=source2;proposal.sourceLevelId=2;proposal.targetLevelId=3;
            proposal.x=family.members(3).baseX;proposal.y=family.members(3).baseY;
        case 'rho_shape',b.rho=b.rho(2:end,:);
        case 'source_box',b.ops.x([1,end])=[-5,5];proposal.x([1,end])=[-5,5];
    end
    o=reject_transfer(b,cfg,proposal);record(['transfer_early_reject_',key], ...
        strcmp(o.identifier,'ipm:AutonomousMeshFamilyTransfer'),o);
end
profile off;prof=profile('info');calls={prof.FunctionTable.FileName};
suffixes={'+mesh/build.m','+evolve/flow.m','+evolve/advance.m', ...
    '+output/restoreCheckpoint.m','+remesh/interpolate.m'};found={};
for k=1:numel(suffixes),if any(endsWith(calls,suffixes{k})),found{end+1}=suffixes{k};end,end %#ok<AGROW>
record('no_build_flow_PDE_restore_or_interpolation',isempty(found),struct('forbiddenCalls',{found}));
record('source_fixture_unchanged',isequaln(original,base),struct('unchanged',isequaln(original,base)));
report=struct('registration',registration,'cases',{cases},'caseCount',numel(cases), ...
    'passedCount',sum(cellfun(@(v)v.passed,cases)), ...
    'allPassed',all(cellfun(@(v)v.passed,cases)), ...
    'actualTransferTested',false,'nativeTransactionQualified',false);
save(fullfile(outputDirectory,'report.mat'),'report','base','source2','family','policy','-v7.3');
write_json(fullfile(outputDirectory,'report.json'),report);
fprintf('VARIABLE_N_AUDIT all=%d passed=%d/%d noLU=1 actualTransfer=0\n',report.allPassed,report.passedCount,report.caseCount);
    function record(name,passed,outcome)
        cases{end+1}=struct('name',name,'passed',logical(passed),'outcome',outcome); %#ok<AGROW>
        write_json(fullfile(outputDirectory,'partial_cases.json'),struct('cases',{cases}));
        fprintf('VARIABLE_N_CASE %s passed=%d\n',name,passed);
    end
end
function state=make_view(config,family,id)
m=family.members(id);ops=sample_ops(m,1);[X,Y]=meshgrid(ops.x,ops.y);
rho=(X.^2-X.^4/32).*(1-Y.^2/16);
memory=struct('version',2,'policy',config.remesh.autonomousMesh, ...
    'referenceFamily',family,'currentLevelId',id,'transactions',struct([]), ...
    'cumulativeAbsolutePeakJump',0,'initialization',struct('selectedBaseX',family.rootX,'selectedBaseY',family.rootY));
state=struct('rho',rho,'ops',ops,'config',config,'step',0,'normalizedTime',0, ...
    'scale',struct('logC_l',0,'logC_omega',0,'physicalTime',0,'canonicalTime',0,'X_shift',0), ...
    'mass0',sum(rho.*ops.integrationWeights,'all'),'rhoRange0',[min(rho,[],'all'),max(rho,[],'all')], ...
    'runMetadata',struct('syntheticUnitFixture',true,'autonomousMesh',memory));
end
function candidate=target_view(source,member)
candidate=source;candidate.ops=sample_ops(member,1);candidate.ops.remeshCount=source.ops.remeshCount;
candidate.ops.rescaling=source.ops.rescaling;
[~,candidate.ops.rescaling.originIndex]=min(abs(candidate.ops.x));
[~,candidate.ops.rescaling.pinIndex]=min(abs(candidate.ops.x-1));
[X,Y]=meshgrid(candidate.ops.x,candidate.ops.y);
candidate.rho=(X.^2-X.^4/32).*(1-Y.^2/16);
end
function ops=sample_ops(member,anchor)
x=member.baseX;y=member.baseY;[~,i]=min(abs(x));[~,j]=min(abs(x-anchor));
ops=struct('x',x,'y',y,'nx',numel(x),'ny',numel(y), ...
    'baseX',x,'baseY',y,'Dx',ipm.mesh.fdMatrix(x,1,7), ...
    'integrationWeights',ipm.mesh.quadrature(y)*ipm.mesh.quadrature(x), ...
    'remeshCount',0,'rescaling',struct('pinX',anchor,'originIndex',i,'pinIndex',j, ...
        'transportAnchorX',anchor,'adaptiveLevels',[.1,.5,.9]));
end
function o=try_audit(a,b,p,prior)
try,v=ipm.remesh.auditCandidate(a,b,p,prior);o=struct('passed',v.passed,'audit',v);
catch e,o=struct('passed',false,'identifier',e.identifier,'message',e.message);end
end
function o=reject_transfer(b,c,p)
o=struct('identifier','','message','');
try,ipm.remesh.transfer(b.rho,b.ops,c,p);catch e,o.identifier=e.identifier;o.message=e.message;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
