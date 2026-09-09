function report=ipm_accellab_test_controller_plan(fixtureFile,outputRoot)
%IPM_ACCELLAB_TEST_CONTROLLER_PLAN Pure geometry plus synthetic telemetry.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);out=fullfile(outputRoot,['controller_plan_regression_',token]);mkdir(out);
registration=struct('kind','synthetic_telemetry_pure_controller_regression', ...
    'fixtureFile',fixtureFile,'outputDirectory',out,'telemetryIsSynthetic',true, ...
    'trajectoryClaim',false,'nativeCheckpointClaim',false,'noLU',true,'noPDE',true, ...
    'timeIncrement',1e-4,'highSequence',[42,41,40,39],'lowSequence',[34,33,32,31], ...
    'planImplementation',which('ipm.evolve.planAutonomousMesh'),'applyImplementation',which('ipm.evolve.applyAutonomousMesh'));
write_json(fullfile(out,'registration.json'),registration);
try
    d=load(fixtureFile,'result');r=d.result;
    assert(r.snapshots.normalizedTime(1)==0);
    x=r.snapshots.x{1}(:)';y=r.snapshots.y{1}(:);rho=r.snapshots.rho{1};
    [X,Y]=meshgrid(x,y);geometry=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y), ...
        'symmetryMode',r.config.physics.symmetryMode);
    assert(isequal(rho,ipm.field.initialDensity(geometry,r.config.physics)), ...
        'ipm:AccellabFixture','Saved initial samples must equal the actual configured datum.');
    policy=ipm.config.autonomousMeshPolicy(struct());
    init=r.metadata.autonomousMesh.initialization;
    assert(isequal(init.selectedBaseX,x)&&isequal(init.selectedBaseY,y));
    ops=struct('x',x,'y',y,'Dx',ipm.mesh.fdMatrix(x,1,7), ...
        'baseX',init.selectedBaseX,'baseY',init.selectedBaseY,'remeshCount',0);
    base=struct('rho',rho,'ops',ops,'config',r.config,'step',0,'normalizedTime',0, ...
        'scale',struct('canonicalTime',0,'physicalTime',0), ...
        'flow',struct('coreGridPoints',42,'verticalCoreGridPoints',42,'safetyFactor',.2), ...
        'runMetadata',struct('testOnlySyntheticTelemetry',true));
    base.config.remesh.autonomousMesh=policy;
    oldRho=rho;oldAxes={x,y};cases={};
    profile clear;profile on
    [high,highPlan]=sequence(base,[42,41,40,39],[42,41,40,39],policy);
    record('above_target_forecast_does_not_request',~highPlan.requested&&isempty(highPlan.candidates)&& ...
        isempty(highPlan.stopReason)&&any(highPlan.predictedCoreCells<policy.predictedCoreBuffer),plan_summary(highPlan));
    [low,lowPlan]=sequence(base,[34,33,32,31],[34,33,32,31],policy);
    changed=true;
    for k=1:numel(lowPlan.candidates)
        a=lowPlan.candidates(k);changed=changed&&~a.unchanged&&(~isequal(a.x,x)||~isequal(a.y,y));
    end
    record('below_target_forecast_excludes_unchanged',lowPlan.requested&&strcmp(lowPlan.reason,'forecast_core_trigger')&& ...
        isempty(lowPlan.stopReason)&&~isempty(lowPlan.candidates)&&changed,plan_summary(lowPlan));
    [~,crossPlan]=sequence(base,[42,41,40,39],[31,31,31,31],policy);
    record('forecast_and_deficit_must_share_axis',~crossPlan.requested&&isempty(crossPlan.candidates)&& ...
        crossPlan.predictedCoreCells(1)<policy.predictedCoreBuffer(1)&&crossPlan.coreCells(2)<policy.targetCoreCells(2),plan_summary(crossPlan));
    % Deliberately minimal apply input has no transfer-scheme or grid fields.
    % If the unchanged guard regressed, transfer would error before any build.
    applyState=base;applyState.config=struct('remesh',struct('autonomousMesh',policy));
    for flag=[true,false]
        plan=struct('initial',false,'candidates',struct('x',x,'y',y,'unchanged',flag));
        o=reject(@()ipm.evolve.applyAutonomousMesh(applyState,plan));
        o.claimedUnchangedFlag=flag;
        record(sprintf('apply_same_axes_flag_%d_rejected_before_build',flag), ...
            strcmp(o.identifier,'ipm:AutonomousMeshUnchangedTransaction'),o);
    end
    memory=high.runMetadata.autonomousMesh;
    repeated=ipm.remesh.controllerTelemetry(memory,high,policy);
    record('identical_time_step_replaces_single_row',isequaln(repeated,memory),struct('windowLength',numel(repeated.window.time)));
    b=high;b.step=b.step+1;
    o=reject(@()ipm.remesh.controllerTelemetry(memory,b,policy));record('same_time_different_step_rejected',strcmp(o.identifier,'ipm:AutonomousMeshMemory'),o);
    b=high;b.scale.canonicalTime=b.scale.canonicalTime-1e-4;b.normalizedTime=b.scale.canonicalTime;
    o=reject(@()ipm.remesh.controllerTelemetry(memory,b,policy));record('backward_time_rejected',strcmp(o.identifier,'ipm:AutonomousMeshMemory'),o);
    b=high;b.ops.remeshCount=1;
    o=reject(@()ipm.remesh.controllerTelemetry(memory,b,policy));record('missing_transaction_rejected',strcmp(o.identifier,'ipm:AutonomousMeshMemory'),o);
    altered=policy;altered.maximumReviewInterval=policy.maximumReviewInterval/2;
    o=reject(@()ipm.remesh.controllerTelemetry(memory,high,altered));record('changed_policy_rejected',strcmp(o.identifier,'ipm:AutonomousMeshMemory'),o);
    broken=memory;broken.cumulativeAbsolutePeakJump=1e-4;
    o=reject(@()ipm.remesh.controllerTelemetry(broken,high,policy));record('nonreconstructible_budget_rejected',strcmp(o.identifier,'ipm:AutonomousMeshMemory'),o);
    % This minimal ledger entry is solely a controller unit fixture, not a
    % persisted/accepted transaction. validateController is not bypassed for CP.
    nextMemory=memory;nextMemory.transactions=struct('relativePeakJump',1e-4);nextMemory.cumulativeAbsolutePeakJump=1e-4;
    b=high;b.ops.remeshCount=1;b.flow.coreGridPoints=35;b.flow.verticalCoreGridPoints=36;
    reset=ipm.remesh.controllerTelemetry(nextMemory,b,policy);
    record('same_time_new_epoch_resets_trend',numel(reset.window.time)==1&&reset.window.time==high.scale.canonicalTime&& ...
        reset.window.step==high.step&&reset.window.remeshCount==1&&isequal(reset.window.coreCells,[35,36]), ...
        struct('window',reset.window,'ledgerIsSynthetic',true));
    profile off;profiling=profile('info');forbidden={};
    suffixes={fullfile('+mesh','build.m'),fullfile('+evolve','flow.m'),fullfile('+remesh','transfer.m'), ...
        fullfile('+evolve','advance.m'),fullfile('+output','restoreCheckpoint.m')};
    for k=1:numel(profiling.FunctionTable)
        name=profiling.FunctionTable(k).FileName;
        for j=1:numel(suffixes)
            if endsWith(name,suffixes{j}),forbidden{end+1}=name;end %#ok<AGROW>
        end
    end
    record('profile_confirms_no_build_flow_transfer_or_PDE',isempty(forbidden),struct('forbiddenCalls',{forbidden}));
    unchanged=isequal(base.rho,oldRho)&&isequal(base.ops.x,oldAxes{1})&&isequal(base.ops.y,oldAxes{2})&& ...
        isequal(high.rho,oldRho)&&isequal(low.rho,oldRho)&&base.ops.remeshCount==0&& ...
        high.ops.remeshCount==0&&low.ops.remeshCount==0;
    record('original_field_axes_and_remesh_count_unchanged',unchanged,struct('passed',unchanged));
    report=struct('registration',registration,'cases',{cases},'allPassed',all(cellfun(@(c)c.passed,cases)), ...
        'caseCount',numel(cases),'passedCount',sum(cellfun(@(c)c.passed,cases)), ...
        'interpretation','Actual saved analytic initial samples with deliberately synthetic core/time telemetry. Pure controller/axis tests; no physical trajectory or native transaction qualification.');
    save(fullfile(out,'report.mat'),'report','highPlan','lowPlan','crossPlan');write_json(fullfile(out,'report.json'),report);
catch e
    profile off;failure=struct('registration',registration,'identifier',e.identifier,'message',e.message,'stack',e.stack);
    write_json(fullfile(out,'failure.json'),failure);rethrow(e);
end
fprintf('CONTROLLER_PLAN_REGRESSION %s all=%d passed=%d/%d\n',out,report.allPassed,report.passedCount,report.caseCount);
    function record(name,passed,outcome)
        cases{end+1}=struct('name',name,'passed',logical(passed),'outcome',outcome);
        fprintf('PLAN_CASE %s passed=%d\n',name,passed);
    end
end
function [state,plan]=sequence(base,xs,ys,policy)
state=base;memory=[];
for k=1:4
    t=(k-1)*1e-4;state.step=k-1;state.normalizedTime=t;
    state.scale.canonicalTime=t;state.scale.physicalTime=t;
    state.flow.coreGridPoints=xs(k);state.flow.verticalCoreGridPoints=ys(k);
    if k<4,memory=ipm.remesh.controllerTelemetry(memory,state,policy);
    else,state.runMetadata.autonomousMesh=memory;[state,plan]=ipm.evolve.planAutonomousMesh(state,false);end
end
end
function o=reject(f)
o=struct('identifier','','message','');
try,f();catch e,o.identifier=e.identifier;o.message=e.message;end
end
function s=plan_summary(p)
s=rmfield(p,{'candidates','axisReport'});s.candidateCount=numel(p.candidates);
if isfield(p.axisReport,'unchangedPairsExcludedFromEvolution'),s.excludedUnchanged=p.axisReport.unchangedPairsExcludedFromEvolution;else,s.excludedUnchanged=0;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
