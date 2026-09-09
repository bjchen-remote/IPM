function report=mesh_test_autonomous_result_policy(resultFile,checkpointFile,outDir)
%MESH_TEST_AUTONOMOUS_RESULT_POLICY Persisted configuration tests, no LU/PDE.
% Enabled synthetic fixtures exercise configuration only; they are not runs
% or native checkpoints and do not assert an autonomous-controller history.
% With the later controller hook, reaching its missing-ledger rejection is
% the expected evidence that configuration validation accepted the fixture.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
report=struct('kind','autonomous_result_policy_contract_nolu_v1', ...
    'resultFile',resultFile,'checkpointFile',checkpointFile, ...
    'syntheticFixturesAreConfigurationOnly',true,'noLU',true,'noPDE',true, ...
    'positiveChecks',struct([]),'negativeChecks',struct([]));
for name={'ipm.output.validateV2','ipm.config.autonomousMeshPolicy',mfilename}
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
loaded=load(resultFile,'result');old=loaded.result;clear loaded;
assert(~isfield(old.config.remesh,'autonomousMesh'));
validated=ipm.output.validate(old);
report.positiveChecks=append_row(report.positiveChecks,positive('actual_old_result_exact',isequaln(validated,old)));
loaded=load(checkpointFile,'checkpoint');originalCheckpoint=loaded.checkpoint;clear loaded;
validatedCheckpoint=ipm.output.readCheckpoint(originalCheckpoint);
report.positiveChecks=append_row(report.positiveChecks,positive('actual_old_checkpoint_exact',isequaln(validatedCheckpoint,originalCheckpoint)));
clear validatedCheckpoint originalCheckpoint;
legal=old;legal.config.remesh.adaptiveRemesh=true;legal.config.remesh.initialAnalyticRemesh=true;
legal.config.scaling.adaptiveLevels=[.1,.5,.9];
policy=ipm.config.autonomousMeshPolicy(struct());legal.config.remesh.autonomousMesh=policy;
for target={[32,32],[42,56],[256,256]}
    fixture=legal;fixture.config.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct('targetCoreCells',target{1}));
    passed=configuration_accepted(fixture);
    report.positiveChecks=append_row(report.positiveChecks,positive(sprintf('legal_target_%d_%d',target{1}),passed));
end
fixture=old;fixture.config.remesh.autonomousMesh=ipm.config.autonomousMeshPolicy(struct('enabled',false));
value=ipm.output.validate(fixture);
report.positiveChecks=append_row(report.positiveChecks,positive('disabled_unchanged_choices',isequaln(value,fixture)));
for version=1:3
    fixture=legacy_fixture(old,version);value=ipm.output.validate(fixture);
    report.positiveChecks=append_row(report.positiveChecks,positive(sprintf('synthetic_legacy_%d_exact',version),isequaln(value,fixture)));
    fixture.config.remesh.autonomousMesh=policy;
    report.negativeChecks=append_row(report.negativeChecks,reject(sprintf('legacy_%d_policy_forbidden',version),fixture));
end
mutations={ ...
    'partial_policy',{'autonomousMesh'},struct(); ...
    'numeric_enabled',{'autonomousMesh','enabled'},1; ...
    'integer_target_type',{'autonomousMesh','targetCoreCells'},uint32([32,32]); ...
    'column_target',{'autonomousMesh','targetCoreCells'},[32;32]; ...
    'string_time_unit',{'autonomousMesh','timeUnit'},"native_canonical"; ...
    'uppercase_time_unit',{'autonomousMesh','timeUnit'},'NATIVE_CANONICAL'; ...
    'partial_quality',{'autonomousMesh','qualityLimits'},struct(); ...
    'partial_search',{'autonomousMesh','search'},struct(); ...
    'column_search_vector',{'autonomousMesh','search','fineCells'},[32;48;64;80;96;112;128]; ...
    'relaxed_axis_ratio',{'autonomousMesh','qualityLimits','maxAdjacentCellRatio'},1.09; ...
    'relaxed_peak_gate',{'autonomousMesh','maximumSinglePeakJump'},.003; ...
    'wrong_derived_floor',{'autonomousMesh','transactionMinimumCoreCells'},[30,31]; ...
    'wrong_target_bounds',{'autonomousMesh','targetCoreCells'},[31,32]; ...
    'unknown_field',{'autonomousMesh','unknown'},true; ...
    'disabled_partial_still_forbidden',{'autonomousMesh'},struct('enabled',false); ...
    'adaptive_disabled',{'adaptiveRemesh'},false; ...
    'analytic_initial_disabled',{'initialAnalyticRemesh'},false; ...
    'wrong_transfer',{'remeshTransferScheme'},'pchip'};
for k=1:size(mutations,1)
    fixture=legal;fixture.config.remesh=set_nested(fixture.config.remesh,mutations{k,2},mutations{k,3});
    report.negativeChecks=append_row(report.negativeChecks,reject(mutations{k,1},fixture));
end
choices={ ...
    'physical_mode','scaling','rescalingMode','physical'; ...
    'anisotropic','scaling','dynamicScaleGeometry','anisotropic'; ...
    'symmetry','physics','symmetryMode','half_plane'; ...
    'length_gauge','scaling','lengthGauge','density_width'; ...
    'omega_gauge','scaling','cOmegaGauge','none'; ...
    'wrong_levels','scaling','adaptiveLevels',[.1,.5,.8]; ...
    'spatial','transport','spatialDiscretization','legacy_second_order'; ...
    'transport','transport','transportScheme','muscl_minmod'; ...
    'integrator','time','timeIntegrator','ssprk3'};
for k=1:size(choices,1)
    fixture=legal;fixture.config.(choices{k,2}).(choices{k,3})=choices{k,4};
    report.negativeChecks=append_row(report.negativeChecks,reject(choices{k,1},fixture));
end
report.allPassed=all([report.positiveChecks.passed])&&all([report.negativeChecks.passed]);
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report');
assert(report.allPassed,'ipm:PolicyRegression','At least one persisted-policy regression failed.');
fprintf('AUTONOMOUS_RESULT_POLICY_COMPLETE positive=%d negative=%d noLU=1 noPDE=1\n',numel(report.positiveChecks),numel(report.negativeChecks));
end
function value=set_nested(value,path,replacement)
if isscalar(path),value.(path{1})=replacement;else,value.(path{1})=set_nested(value.(path{1}),path(2:end),replacement);end
end
function passed=configuration_accepted(fixture)
try
    passed=isequaln(ipm.output.validate(fixture),fixture);
catch exception
    passed=~isfield(fixture.metadata,'autonomousMesh') && ...
        strcmp(exception.identifier,'ipm:AutonomousMeshController') && ...
        strcmp(exception.message,'An enabled policy requires its original controller ledger.');
    if ~passed,rethrow(exception);end
end
end
function fixture=legacy_fixture(value,version)
fixture=value;fixture.config.schemaVersion=version;
s=rmfield(fixture.config.scaling,'scalingContract');s.cOmegaGauge='none';
s.adaptiveLengthScaling=false;s.widthExpansionStrength=1;s.widthContractionOnset=1.35;
s.widthContractionStrength=28;s.maxWidthRateCorrection=16;
if version<=2,s=rmfield(s,'omegaGaugeWindowRadius');s.maxDynamicRate=2;end
fixture.config.scaling=s;
if version==1
    fixture.config.time=rmfield(fixture.config.time,'timeIntegrator');
    fixture.config.transport=rmfield(fixture.config.transport,'spatialDiscretization');
    fixture.config.transport.transportScheme='muscl_minmod';
    fixture.config.remesh=rmfield(fixture.config.remesh,'remeshTransferScheme');
end
end
function row=positive(label,passed)
row=struct('label',label,'passed',passed);
end
function row=reject(label,fixture)
identifier='';message='';
try
    ipm.output.validate(fixture);
catch exception
    identifier=exception.identifier;message=exception.message;
end
row=struct('label',label,'identifier',identifier,'message',message,'passed',strcmp(identifier,'ipm:ResultConfig'));
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
