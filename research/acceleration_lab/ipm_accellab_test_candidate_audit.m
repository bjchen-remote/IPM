function report=ipm_accellab_test_candidate_audit(preFile,postFile,nativeAuditFile,outputRoot)
%IPM_ACCELLAB_TEST_CANDIDATE_AUDIT No-LU native-data precommit regressions.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['candidate_audit_regression_',token]);mkdir(destination);
registration=struct('testRevision',2,'preFile',preFile,'postFile',postFile,'nativeAuditFile',nativeAuditFile, ...
    'outputDirectory',destination,'strictNativeRead',true,'testOnlyPrecommitView',true, ...
    'checkpointModified',false,'poissonBuilds',0,'flowEvaluations',0,'pdeSteps',0, ...
    'auditImplementation',which('ipm.remesh.auditCandidate'), ...
    'groups',{{'native_parity','signed_mass','invariants','precommit_bookkeeping', ...
    'off_wall_peak','actual_resolution','adaptive_level_scope','invalid_numeric'}});
write_json(fullfile(destination,'registration.json'),registration);
try
    pre=ipm.output.readCheckpoint(preFile);post=ipm.output.readCheckpoint(postFile);
    preOriginal=pre;postOriginal=post;
    o=pre.payload.state;n=post.payload.state;legacy=jsondecode(fileread(nativeAuditFile));
    assert(o.step==3668&&n.step==3668&&n.remeshCount==o.remeshCount+1);
    assert(isequaln(o.scale,n.scale)&&o.normalizedTime==n.normalizedTime&& ...
        isequal(o.mass0,n.mass0)&&isequal(o.rhoRange0,n.rhoRange0));
    old=state_view(o);candidate=old;
    nativeNew=state_view(n);candidate.rho=n.rho;candidate.ops=nativeNew.ops;
    candidate.ops.remeshCount=old.ops.remeshCount;
    candidate.testOnlyPrecommitView=true;old.testOnlyPrecommitView=true;
    % CP payloads are immutable. Only this explicitly synthetic state view
    % combines the new numerical grid with the unchanged caller references.
    policy=ipm.config.autonomousMeshPolicy(struct());
    prior=legacy.cumulativeAbsolutePeakJump-legacy.rhoXMaximumRelativeChange;
    actual=ipm.remesh.auditCandidate(old,candidate,policy,prior);
    expected=[legacy.xCoreCells,legacy.yCoreCells,legacy.rhoXMaximumRelativeChange, ...
        legacy.massRelativeDefect,legacy.relativeRangeViolation,legacy.actualLeftFrontCells];
    measured=[actual.coreCells,actual.relativePeakJump,actual.massRelativeDefect, ...
        actual.relativeRangeViolation,actual.leftFrontCells];
    parity=struct('fields',{{'coreX','coreY','wholeBoxPeakJump','massDefect','rangeViolation','frontCells'}}, ...
        'native',expected,'measured',measured,'difference',measured-expected, ...
        'bitwiseEqual',isequal(measured,expected),'withinRoundoff',all(abs(measured-expected)<=64*eps(max(1,abs(expected)))), ...
        'actualAuditPassed',actual.passed,'actualAudit',actual);
    identity=ipm.remesh.auditCandidate(candidate,candidate,policy,0);
    cases={};
    record('native_parity','actual_pre_post_measurement',parity.withinRoundoff&&actual.passed,parity);
    record('native_parity','identity_passing_fixture',identity.passed,identity);
    % A signed-cancellation field is a test fixture, not a PDE state.
    signed=candidate;w=signed.ops.integrationWeights;
    signed.rho=signed.rho-sum(w.*signed.rho,'all')/sum(w,'all');
    absMass=sum(w.*abs(signed.rho),'all');
    signed.rho=signed.rho+(.001*absMass-sum(w.*signed.rho,'all'))/sum(w,'all');
    signedTrial=signed;mass=sum(w.*signed.rho,'all');
    signedTrial.rho=signed.rho+1e-10*abs(mass)/sum(w,'all');
    outcome=call_audit(signed,signedTrial,policy,0);
    alternate=abs(sum(w.*signedTrial.rho,'all')-mass)/max(sum(w.*abs(signed.rho),'all'),realmin);
    outcome.l1NormalizedMassDefect=alternate;
    record('signed_mass','native_denominator_rejects_L1_would_pass', ...
        rejected_for(outcome,'mass_defect')&&alternate<policy.maximumMassRelativeDefect,outcome);
    fields={'mass0','rhoRange0','config','runMetadata','scale','step','normalizedTime','baseX','runtime_reference'};
    for k=1:numel(fields)
        c=candidate;key=fields{k};
        switch key
            case 'mass0',c.mass0=c.mass0+max(1,abs(c.mass0));
            case 'rhoRange0',c.rhoRange0(1)=c.rhoRange0(1)-1;
            case 'config',c.config.testOnlyMutation=true;
            case 'runMetadata',c.runMetadata.testOnlyMutation=true;
            case 'scale',c.scale.physicalTime=c.scale.physicalTime+1;
            case 'step',c.step=c.step+1;
            case 'normalizedTime',c.normalizedTime=c.normalizedTime+1;
            case 'baseX',c.ops.baseX(2)=c.ops.baseX(2)+eps(c.ops.baseX(2));
            case 'runtime_reference',c.ops.rescaling.pinX=c.ops.rescaling.pinX+1;
        end
        outcome=call_audit(candidate,c,policy,0);record('invariants',key,~outcome.returnedPassed,outcome);
    end
    fields={'origin_index','pin_index','missing_anchor','endpoint','premature_count'};
    for k=1:numel(fields)
        c=candidate;key=fields{k};
        switch key
            case 'origin_index',c.ops.rescaling.originIndex=c.ops.rescaling.originIndex+1;
            case 'pin_index',c.ops.rescaling.pinIndex=c.ops.rescaling.pinIndex+1;
            case 'missing_anchor'
                a=c.ops.rescaling.transportAnchorX;j=find(c.ops.x==a);jj=find(c.ops.x==-a);
                assert(isscalar(j)&&isscalar(jj));c.ops.x(j)=a+eps(a);c.ops.x(jj)=-c.ops.x(j);
            case 'endpoint',c.ops.x(end)=c.ops.x(end)+1;c.ops.x(1)=-c.ops.x(end);
            case 'premature_count',c.ops.remeshCount=c.ops.remeshCount+1;
        end
        outcome=call_audit(candidate,c,policy,0);record('precommit_bookkeeping',key,~outcome.returnedPassed,outcome);
    end
    outcome=call_audit(old,candidate,policy,policy.maximumCumulativeAbsolutePeakJump);
    record('precommit_bookkeeping','existing_budget_is_not_reset',rejected_for(outcome,'cumulative_peak_budget'),outcome);
    c=candidate;c.rho(2:end,:)=2*c.rho(2:end,:);
    outcome=call_audit(candidate,c,policy,0);outcome.wallTraceUnchanged=isequal(c.rho(1,:),candidate.rho(1,:));
    record('off_wall_peak','full_box_peak_increase_same_wall', ...
        outcome.wallTraceUnchanged&&rejected_for(outcome,'single_peak_jump'),outcome);
    % Deliberately under-resolved analytic mutation with fictitious padded
    % design scores: those scores must not replace actual field diagnostics.
    c=candidate;[X,Y]=meshgrid(c.ops.x,c.ops.y);a=c.ops.rescaling.transportAnchorX;
    [~,j]=min(abs(c.ops.x-a));dx=.5*(c.ops.x(j+1)-c.ops.x(j-1));dy=c.ops.y(2)-c.ops.y(1);
    c.rho=-exp(-((abs(X)-a)/(2*dx)).^2-(Y/(2*dy)).^2);
    c.proposedPaddedCoreCells=[999,999];c.proposedPaddedFrontCells=999;
    outcome=call_audit(candidate,c,policy,0);
    record('actual_resolution','padded_design_cannot_hide_actual_core',rejected_for(outcome,'actual_core_floor'),outcome);
    record('actual_resolution','padded_design_cannot_hide_actual_front',rejected_for(outcome,'actual_front_floor'),outcome);
    choices=struct('rescalingMode','dynamic','dynamicScaleGeometry','isotropic', ...
        'symmetryMode','double_odd_omega','lengthGauge','transport_anchor', ...
        'cOmegaGauge','wall_omega_quadratic_peak','spatialDiscretization','high_order', ...
        'transportScheme','weno5_fd','timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
        'adaptiveRemesh',true,'initialAnalyticRemesh',true,'adaptiveLevels',[.1,.5,.9]);
    valid=ipm.config.autonomousMeshPolicy(struct(),choices);assert(valid.enabled);
    choices.adaptiveLevels=[.1,.5,.8];scope=struct('rejected',false,'identifier','');
    try,ipm.config.autonomousMeshPolicy(struct(),choices);catch e,scope.rejected=true;scope.identifier=e.identifier;end
    record('adaptive_level_scope','nonstandard_levels',scope.rejected&&strcmp(scope.identifier,'ipm:AutonomousMeshCombination'),scope);
    fields={'rho_nan','candidate_weight_nan','source_weight_nan','derivative_nan','axis_nan', ...
        'candidate_weight_shape','candidate_weight_zero','candidate_weight_negative','candidate_weight_complex','candidate_weight_inf', ...
        'source_weight_shape','source_weight_zero','source_weight_negative','source_weight_complex','source_weight_inf', ...
        'derivative_inf','source_derivative_nan','source_derivative_inf','rho_complex','mass_overflow'};
    for k=1:numel(fields)
        b=candidate;c=candidate;key=fields{k};
        switch key
            case 'rho_nan',c.rho(1,1)=NaN;
            case 'candidate_weight_nan',c.ops.integrationWeights(1,1)=NaN;
            case 'source_weight_nan',b.ops.integrationWeights(1,1)=NaN;
            case 'derivative_nan',c.ops.Dx(1,1)=NaN;
            case 'axis_nan',c.ops.y(2)=NaN;
            case 'candidate_weight_shape',c.ops.integrationWeights=c.ops.integrationWeights(2:end,:);
            case 'candidate_weight_zero',c.ops.integrationWeights(1,1)=0;
            case 'candidate_weight_negative',c.ops.integrationWeights(1,1)=-1;
            case 'candidate_weight_complex',c.ops.integrationWeights(1,1)=1+1i;
            case 'candidate_weight_inf',c.ops.integrationWeights(1,1)=Inf;
            case 'source_weight_shape',b.ops.integrationWeights=b.ops.integrationWeights(2:end,:);
            case 'source_weight_zero',b.ops.integrationWeights(1,1)=0;
            case 'source_weight_negative',b.ops.integrationWeights(1,1)=-1;
            case 'source_weight_complex',b.ops.integrationWeights(1,1)=1+1i;
            case 'source_weight_inf',b.ops.integrationWeights(1,1)=Inf;
            case 'derivative_inf',c.ops.Dx(1,1)=Inf;
            case 'source_derivative_nan',b.ops.Dx(1,1)=NaN;
            case 'source_derivative_inf',b.ops.Dx(1,1)=Inf;
            case 'rho_complex',c.rho(1,1)=c.rho(1,1)+1i;
            case 'mass_overflow',c.ops.integrationWeights(:)=realmax;
        end
        outcome=call_audit(b,c,policy,0);record('invalid_numeric',key,~outcome.returnedPassed,outcome);
    end
    unchanged=isequaln(pre,preOriginal)&&isequaln(post,postOriginal);
    report=struct('registration',registration,'strictNativeReadsPassed',true, ...
        'nativeAuditParity',parity,'cases',{cases},'sourceCheckpointObjectsUnchanged',unchanged, ...
        'caseCount',numel(cases),'passedCount',sum(cellfun(@(c)c.passed,cases)), ...
        'allPassed',unchanged&&all(cellfun(@(c)c.passed,cases)), ...
        'meaning','Synthetic precommit views and deliberate negative mutations only; no native checkpoint or physical solution was modified.');
    save(fullfile(destination,'report.mat'),'report');write_json(fullfile(destination,'report.json'),report);
    save(fullfile(destination,'test_only_precommit_fixture.mat'),'old','candidate','policy','prior','-v7.3');
catch exception
    failure=struct('registration',registration,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('CANDIDATE_AUDIT_REGRESSION %s all=%d cases=%d passed=%d parity=%d\n', ...
    destination,report.allPassed,report.caseCount,report.passedCount,parity.bitwiseEqual);
    function record(group,name,passed,outcome)
        cases{end+1}=struct('group',group,'name',name,'passed',logical(passed),'outcome',outcome); %#ok<AGROW>
        fprintf('AUDIT_CASE %s %s passed=%d\n',group,name,passed);
    end
end
function state=state_view(payload)
assert(strcmp(payload.config.transport.spatialDiscretization,'high_order'));
x=payload.x(:)';y=payload.y(:);Dx=ipm.mesh.fdMatrix(x,1,7);
wx=ipm.mesh.quadrature(x);wy=ipm.mesh.quadrature(y);
ops=struct('x',x,'y',y,'Dx',Dx,'integrationWeights',wy*wx, ...
    'baseX',payload.baseX,'baseY',payload.baseY,'rescaling',payload.rescaling,'remeshCount',payload.remeshCount);
state=struct('rho',payload.rho,'ops',ops,'scale',payload.scale,'step',payload.step, ...
    'normalizedTime',payload.normalizedTime,'config',payload.config,'runMetadata',payload.runMetadata, ...
    'mass0',payload.mass0,'rhoRange0',payload.rhoRange0);
end
function outcome=call_audit(b,c,p,prior)
try
    audit=ipm.remesh.auditCandidate(b,c,p,prior);
    outcome=struct('returnedPassed',audit.passed,'threw',false,'identifier','','audit',audit);
catch exception
    outcome=struct('returnedPassed',false,'threw',true,'identifier',exception.identifier,'message',exception.message);
end
end
function passed=rejected_for(outcome,reason)
passed=~outcome.returnedPassed&&isfield(outcome,'audit')&&any(strcmp(outcome.audit.reasons,reason));
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
