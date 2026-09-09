function report=ipm_perflab_test_budget_chain(campaignRoot,outDir)
% Actual no-LU chain plus explicitly synthetic local scheduling negatives.
assert(~isfolder(outDir));mkdir(outDir);
directory=fullfile(campaignRoot,'fresh_adaptive_campaign_v4');
stageDir=fullfile(directory,'stage_006');
endpoint=jsondecode(fileread(fullfile(stageDir,'endpoint_audit.json')));
sourceReview=jsondecode(fileread(fullfile(stageDir,'source_review.json')));
r=jsondecode(fileread(fullfile(directory,'registration.json')));
q=jsondecode(fileread(fullfile(campaignRoot,'fresh_adaptive_launch_v4','caller_qualification.json')));
q.bootstrapResultFile=endpoint.resultFile;q.bootstrapCheckpointFile=endpoint.checkpointFile;
q.campaignDirectories{end+1}=directory;
q.budgetLimitedCheckpointFiles={endpoint.checkpointFile};
q.newChainEvidenceDirectory=fullfile(outDir,'actual_v1_to_v4');
q.purpose='Research regression: explicit budget-limited step7510 linkage, no numerical advancement.';
qualificationFile=fullfile(outDir,'test_caller_qualification.json');write_json(qualificationFile,q);
profile clear;profile on
actual=fresh_profile_campaign_admit(qualificationFile,endpoint.resultFile,r.initialDataFile,endpoint.checkpointFile);
profile off;profileInfo=profile('info');save(fullfile(outDir,'profile.mat'),'profileInfo');
files={profileInfo.FunctionTable.FileName};
assert(~any(contains(files,{'/+mesh/build.m','/+field/poisson.m','/+evolve/flow.m','/+ipm/solve.m'})), ...
    'ipm:ResearchBudgetChainNoLU','A forbidden numerical constructor or PDE entry was invoked.');
loaded=load(fullfile(q.newChainEvidenceDirectory,'chain_report.mat'),'chain','contract');
chain=loaded.chain;contract=loaded.contract;
assert(chain.allPassed && chain.links(end).numericalStateQualified && ...
    ~chain.links(end).reviewHorizonReached && ~chain.links(end).savedStageAccepted && ...
    strcmp(chain.links(end).stopReason,'maximum_steps') && chain.targetReview.step==7510);
assert(all([chain.links(1:end-1).reviewHorizonReached]) && all([chain.links(1:end-1).savedStageAccepted]));
aCP=ipm.output.readCheckpoint(endpoint.sourceCheckpoint);bCP=ipm.output.readCheckpoint(endpoint.checkpointFile);
a=aCP.payload.state;b=bCP.payload.state;next=chain.targetReview;
lines=splitlines(fileread(fullfile(directory,'campaign_manifest.jsonl')));started=[];
for k=1:numel(lines)
    if strlength(strtrim(lines(k)))==0,continue,end
    row=jsondecode(lines{k});
    if isfield(row,'event') && strcmp(row.event,'solve_started') && row.stage==6,started=row;end
end
args={r,6,sourceReview,endpoint,a,b,next,contract,'maximum_steps', ...
    chain.cumulativeAbsolutePeakJump,{endpoint.checkpointFile},started};
v=fresh_chain_segment_test_entry(args{:});assert(~v.reviewHorizonReached && v.numericalStateQualified);
tests=struct('name',{},'args',{});
z=args;z{11}={};tests(end+1)=test_case('default_rejects_actual_partial',z);
z=args;z{11}={[endpoint.checkpointFile '_wrong']};tests(end+1)=test_case('wrong_allowlist',z);
z=args;z{9}='minimum_timestep';z{4}.stopReason=z{9};tests(end+1)=test_case('nonbudget_stop_reason',z);
z=args;z{6}.step=b.step-1;z{4}.step=z{6}.step;z{4}.segmentSteps=z{6}.step-a.step;z{7}.step=z{6}.step;
tests(end+1)=test_case('budget_not_exhausted',z);
z=args;z{4}.targetParentEquivalentTau=z{4}.targetParentEquivalentTau+eps(z{4}.targetParentEquivalentTau);
tests(end+1)=test_case('tampered_saved_target',z);
z=args;z{12}.targetLocalCanonicalTime=z{12}.targetLocalCanonicalTime+eps(z{12}.targetLocalCanonicalTime);
tests(end+1)=test_case('tampered_manifest_target',z);
z=args;z{6}.config.time.finalTime=b.config.time.finalTime+eps(b.config.time.finalTime);
tests(end+1)=test_case('tampered_native_target',z);
z=args;z{1}.maximumAdditionalStepsPerStage=513;tests(end+1)=test_case('tampered_old_registration_budget',z);
z=args;z{6}.config.time.maxSteps=b.config.time.maxSteps+1;tests(end+1)=test_case('tampered_native_budget',z);
z=args;z{3}.localCanonicalTime=z{3}.localCanonicalTime+eps(z{3}.localCanonicalTime);
tests(end+1)=test_case('tampered_source_time',z);
z=args;z{4}.stageAccepted=true;tests(end+1)=test_case('partial_claims_saved_stage_accepted',z);
z=args;z{1}.firstEquivalentTauStep=.1;tests(end+1)=test_case('old_kind_short_first_interval',z);
z=args;z{7}.coreCells(1)=19.99;tests(end+1)=test_case('unchanged_core20_floor',z);
z=args;z{7}.safety=.70;tests(end+1)=test_case('unchanged_safety_gate',z);
% New-kind scheduling fixtures are synthetic projections only, never native
% CPs and never passed through the full chain/native state qualification.
full=full_projection(args,.2,512,6);
fullVerdict=fresh_chain_segment_test_entry(full{:});assert(fullVerdict.reviewHorizonReached);
short=full_projection(args,.013032182165655,4096,1);
short{1}.kind='independent_fresh_profile_campaign_v3';short{1}.maximumAdditionalStepsPerStage=4096;
short{1}.firstEquivalentTauStep=.013032182165655;
shortVerdict=fresh_chain_segment_test_entry(short{:});assert(shortVerdict.reviewHorizonReached);
later=full_projection(short,.2,4096,2);laterVerdict=fresh_chain_segment_test_entry(later{:});
assert(laterVerdict.reviewHorizonReached && laterVerdict.registeredEquivalentTauIncrement==.2);
z=short;z{1}.maximumAdditionalStepsPerStage=512;tests(end+1)=test_case('new_kind_requires4096',z);
z=short;z{1}.firstEquivalentTauStep=.2001;tests(end+1)=test_case('new_first_interval_too_long',z);
z=short;z{1}.firstEquivalentTauStep=0;tests(end+1)=test_case('new_first_interval_zero',z);
negatives=struct([]);
for k=1:numel(tests)
    caught=false;identifier='';
    try
        fresh_chain_segment_test_entry(tests(k).args{:});
    catch exception
        caught=true;identifier=exception.identifier;
    end
    assert(caught && startsWith(identifier,'ipm:FreshChain'), ...
        'ipm:ResearchBudgetChainNegative','Invalid case %s did not reach a chain contract rejection.',tests(k).name);
    item=struct('name',tests(k).name,'rejected',true,'identifier',identifier);
    if isempty(negatives),negatives=item;else,negatives(end+1)=item;end %#ok<AGROW>
end
report=struct('kind','actual_budget_limited_chain_and_synthetic_boundary_tests_v1','allPassed',true, ...
    'actualAdmission',actual,'actualLinks',numel(chain.links),'actualStep',chain.targetReview.step, ...
    'actualEquivalentTau',chain.targetReview.parentEquivalentTau,'actualPartialLink',chain.links(end), ...
    'actualNumericalStateQualified',true,'actualReviewHorizonReached',false,'originalStageAcceptedPreserved',false, ...
    'syntheticSchedulingPositiveCount',3,'negativeTests',negatives,'noLU',true,'noPDE',true, ...
    'syntheticTestsAreNativeEvidence',false, ...
    'partialBoundaryTestUsesActualSourceText',true);
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report','-v7.3');
fprintf('BUDGET_CHAIN_PASS actualLinks=%d actualStep=%d negatives=%d\n',numel(chain.links),chain.targetReview.step,numel(negatives));
end

function z=full_projection(args,increment,budget,stage)
z=args;lambda=z{8}.lineage.canonicalCovarianceFactor;
target=z{5}.scale.canonicalTime+increment/lambda;
equivalent=z{8}.lineage.parentCanonicalTime+lambda*target;
z{2}=stage;z{6}.step=z{5}.step+100;z{6}.scale.canonicalTime=target;
z{6}.config.time.finalTime=target;z{6}.config.time.maxSteps=z{5}.step+budget;
z{4}.stage=stage;z{4}.segmentSteps=100;z{4}.step=z{6}.step;z{4}.localCanonicalTime=target;
z{4}.parentEquivalentTau=equivalent;z{4}.targetParentEquivalentTau=equivalent;
z{4}.parentEquivalentEndpointError=0;z{4}.stageAccepted=true;z{4}.stopReason='final_time';
z{7}.step=z{6}.step;z{7}.localCanonicalTime=target;z{7}.parentEquivalentTau=equivalent;
z{9}='final_time';z{11}={};z{12}.stage=stage;
z{12}.targetLocalCanonicalTime=target;z{12}.targetParentEquivalentTau=equivalent;
end
function t=test_case(name,args)
t=struct('name',name,'args',{args});
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
