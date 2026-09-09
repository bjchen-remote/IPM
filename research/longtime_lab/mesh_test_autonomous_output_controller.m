function report=mesh_test_autonomous_output_controller(oldResultFile,oldCheckpointFile,enabledFiles,outDir)
%MESH_TEST_AUTONOMOUS_OUTPUT_CONTROLLER Pure result-to-controller hook replay.
% Existing solved data are immutable inputs. Mutations below are rejected
% contract fixtures, never asserted to be native trajectories or checkpoints.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
report=struct('kind','autonomous_output_controller_hook_nolu_v1', ...
    'noLU',true,'noPDE',true,'positiveChecks',struct([]),'negativeChecks',struct([]));
for name={'ipm.output.validateV2','ipm.remesh.validateController',mfilename}
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
loaded=load(oldResultFile,'result');old=loaded.result;clear loaded;
assert(~isfield(old.config.remesh,'autonomousMesh') && ~isfield(old.metadata,'autonomousMesh'));
report.positiveChecks=append_row(report.positiveChecks,positive('old_result_exact',isequaln(ipm.output.validate(old),old)));
loaded=load(oldCheckpointFile,'checkpoint');cp=loaded.checkpoint;clear loaded;
report.positiveChecks=append_row(report.positiveChecks,positive('old_native_checkpoint_exact',isequaln(ipm.output.readCheckpoint(cp),cp)));clear cp;
fixture=old;fixture.config.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct('enabled',false));
report.positiveChecks=append_row(report.positiveChecks,positive('disabled_no_ledger_exact',isequaln(ipm.output.validate(fixture),fixture)));
for k=1:numel(enabledFiles)
    loaded=load(enabledFiles{k},'result');enabled=loaded.result;clear loaded;
    [~,label]=fileparts(enabledFiles{k});
    report.positiveChecks=append_row(report.positiveChecks,positive([label '_actual_exact'],isequaln(ipm.output.validate(enabled),enabled)));
    fixture=enabled;fixture.metadata=rmfield(fixture.metadata,'autonomousMesh');
    report.negativeChecks=append_row(report.negativeChecks,reject([label '_missing_ledger'],fixture,'ipm:AutonomousMeshController'));
end
fixture=enabled;fixture.config.remesh.autonomousMesh.enabled=false;
report.negativeChecks=append_row(report.negativeChecks,reject('disabled_with_ledger',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.config.remesh=rmfield(fixture.config.remesh,'autonomousMesh');
report.negativeChecks=append_row(report.negativeChecks,reject('absent_policy_with_ledger',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.window.coreCells(end,1)=fixture.metadata.autonomousMesh.window.coreCells(end,1)+.01;
report.negativeChecks=append_row(report.negativeChecks,reject('window_core_mismatch',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.window.safety(end)=fixture.metadata.autonomousMesh.window.safety(end)+.01;
report.negativeChecks=append_row(report.negativeChecks,reject('window_safety_mismatch',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.window.step(end)=fixture.metadata.autonomousMesh.window.step(end)+1;
report.negativeChecks=append_row(report.negativeChecks,reject('window_step_mismatch',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.initialization.originalPhysicalTime=.1;
report.negativeChecks=append_row(report.negativeChecks,reject('fresh_epoch_impostor',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.initialization=rmfield(fixture.metadata.autonomousMesh.initialization,'operatorGridRepresentation');
report.negativeChecks=append_row(report.negativeChecks,reject('missing_operator_representation',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.cumulativeAbsolutePeakJump=.001;
report.negativeChecks=append_row(report.negativeChecks,reject('unearned_cumulative_budget',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.policy.enabled=1;
report.negativeChecks=append_row(report.negativeChecks,reject('ledger_noncanonical_type',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.config.remesh.autonomousMesh.enabled=1;
report.negativeChecks=append_row(report.negativeChecks,reject('config_noncanonical_type',fixture,'ipm:ResultConfig'));
fixture=enabled;fixture.metadata.autonomousMesh.initialization.selectedBaseX(end)=fixture.metadata.autonomousMesh.initialization.selectedBaseX(end)+1;
report.negativeChecks=append_row(report.negativeChecks,reject('selected_box_changed',fixture,'ipm:AutonomousMeshController'));
fixture=enabled;fixture.metadata.autonomousMesh.initialization.audit.xQuality.maximumAdjacentCellRatio=1.09;
fixture.metadata.autonomousMesh.initialization.attempts(end).audit=fixture.metadata.autonomousMesh.initialization.audit;
report.negativeChecks=append_row(report.negativeChecks,reject('initial_quality_relaxed',fixture,'ipm:AutonomousMeshController'));
report.allPassed=all([report.positiveChecks.passed])&&all([report.negativeChecks.passed]);
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report');
assert(report.allPassed,'ipm:ControllerHookRegression','At least one output/controller boundary case failed.');
fprintf('AUTONOMOUS_OUTPUT_CONTROLLER_COMPLETE positive=%d negative=%d noLU=1 noPDE=1\n',numel(report.positiveChecks),numel(report.negativeChecks));
end
function row=positive(label,passed)
row=struct('label',label,'passed',passed);
end
function row=reject(label,fixture,expected)
identifier='';message='';
try
    ipm.output.validate(fixture);
catch exception
    identifier=exception.identifier;message=exception.message;
end
row=struct('label',label,'expected',expected,'identifier',identifier,'message',message,'passed',strcmp(identifier,expected));
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
