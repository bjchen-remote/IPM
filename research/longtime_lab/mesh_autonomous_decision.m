function decision=mesh_autonomous_decision(samples,policy)
%MESH_AUTONOMOUS_DECISION Research scheduling from accepted mesh telemetry.
% This is not a numerical step, native checkpoint validator, mesh constructor,
% or guarantee against unresolved changes between observations. A future
% integration must inspect each accepted step and validate every transaction.
% time and nextStableStep must have the explicitly named SAME time unit.
% nextStableStep is the solver's current unconstrained stability limit, not
% the last endpoint-clipped accepted step. Both returned limits must be used.
required={'time','step','remeshCount','coreCells','safety','trusted', ...
    'nextStableStep','cumulativeAbsolutePeakJump','timeUnit'};
assert(isstruct(samples) && isscalar(samples) && all(isfield(samples,required)), ...
    'ipm:AutonomousMeshTelemetry','Missing accepted-state telemetry.');
defaults=struct('maximumReviewInterval',.2,'trendWindow',.35, ...
    'maximumStepsPerReview',512,'stepBudgetFraction',.75, ...
    'historyCoreFloor',14,'endpointCoreFloor',20,'regridCoreTrigger',26, ...
    'predictedCoreBuffer',22,'maximumSafety',.70, ...
    'maximumCumulativeAbsolutePeakJump',.02);
assert(isstruct(policy) && isscalar(policy) && isfield(policy,'timeUnit') && ...
    isempty(setdiff(fieldnames(policy),[fieldnames(defaults);{'timeUnit'}])), ...
    'ipm:AutonomousMeshPolicy','The caller must explicitly register its time unit.');
assert(strcmp(samples.timeUnit,policy.timeUnit) && ...
    any(strcmp(policy.timeUnit,{'native_canonical','parent_equivalent_canonical'})), ...
    'ipm:AutonomousMeshClock','The telemetry and policy time units must agree.');
for name=fieldnames(defaults)'
    if ~isfield(policy,name{1}),policy.(name{1})=defaults.(name{1});end
    validateattributes(policy.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
assert(policy.historyCoreFloor<policy.endpointCoreFloor && ...
    policy.endpointCoreFloor<policy.predictedCoreBuffer && ...
    policy.predictedCoreBuffer<policy.regridCoreTrigger && ...
    policy.stepBudgetFraction<1 && mod(policy.maximumStepsPerReview,1)==0);
t=samples.time(:);step=samples.step(:);epoch=samples.remeshCount(:);
safety=samples.safety(:);trusted=samples.trusted(:);cores=samples.coreCells;n=numel(t);
assert(n>=1 && numel(step)==n && numel(epoch)==n && numel(safety)==n && ...
    numel(trusted)==n && isequal(size(cores),[n,2]) && ...
    all(isfinite([t;step;epoch;safety;cores(:)])) && ...
    all(diff(t)>0) && all(diff(step)>0) && all(diff(epoch)>=0) && ...
    all(step>=0 & mod(step,1)==0) && all(epoch>=0 & mod(epoch,1)==0) && ...
    all(cores(:)>0) && all(safety>=0) && islogical(trusted), ...
    'ipm:AutonomousMeshTelemetry','Malformed or nonchronological accepted samples.');
validateattributes(samples.nextStableStep,{'numeric'},{'scalar','finite','positive'});
validateattributes(samples.cumulativeAbsolutePeakJump,{'numeric'},{'scalar','finite','nonnegative'});
decision=struct('kind','research_autonomous_mesh_decision_v1','timeUnit',policy.timeUnit, ...
    'action','advance','reason','resolved_with_current_budget','policy',policy, ...
    'currentStep',step(end),'currentRemeshCount',epoch(end),'currentCoreCells',cores(end,:), ...
    'sameMeshRecords',0,'trendAvailable',false,'decayUpperEstimate',[NaN,NaN], ...
    'predictedCoreCells',[NaN,NaN],'maximumReviewInterval',0, ...
    'maximumAdvanceSteps',0,'historicalPrefixPassed',false, ...
    'controllerIsIntegrated',false,'nativeSignaturesCheckedHere',false, ...
    'pdeAdvanced',false,'errorOrConvergenceGuarantee',false);
if ~all(trusted) || any(cores(:)<policy.historyCoreFloor) || any(safety>=policy.maximumSafety)
    decision.action='preserve_and_reject';decision.reason='historical_resolution_or_trust_failure';return
end
decision.historicalPrefixPassed=true;
if samples.cumulativeAbsolutePeakJump>policy.maximumCumulativeAbsolutePeakJump
    decision.action='preserve_and_reject';decision.reason='cumulative_transfer_budget_exceeded';return
end
if any(cores(end,:)<policy.endpointCoreFloor)
    decision.action='regrid';decision.reason='endpoint_buffer_already_exhausted';return
end
% Use only the current remesh epoch. Cross-transfer count changes are not
% physical contraction rates. A fresh epoch initially gets one-step review.
use=find(epoch==epoch(end) & t>=t(end)-policy.trendWindow);
decision.sameMeshRecords=numel(use);
decision.trendAvailable=numel(use)>=4;
if any(cores(end,:)<policy.regridCoreTrigger)
    decision.action='regrid';decision.reason='current_core_trigger';return
end
if ~decision.trendAvailable
    decision.reason='one_step_observation_for_new_mesh_epoch';
    decision.maximumAdvanceSteps=1;
    decision.maximumReviewInterval=min(policy.maximumReviewInterval,samples.nextStableStep);return
end
u=t(use);v=cores(use,:);decay=zeros(1,2);
for k=1:2
    slope=polyfit(u-u(end),log(v(:,k)),1);
    secants=-diff(log(v(:,k)))./diff(u);
    % An empirical upper estimate, not an a priori bound. Keep the steepest
    % observed decrease; do not smooth away a recently accelerating collapse.
    decay(k)=max([0;-slope(1);secants]);
end
decision.decayUpperEstimate=decay;
% A fixed review interval can exceed a fixed step budget late in the run.
% Shorten the scheduling horizon, never enlarge the PDE stable step/CFL.
% Old endpoint clipping can make an accepted dt arbitrarily small; it must
% not create a persistent artificial restriction on later review intervals.
stepEstimate=samples.nextStableStep;
budgetSteps=max(1,floor(policy.stepBudgetFraction*policy.maximumStepsPerReview));
decision.maximumAdvanceSteps=budgetSteps;
decision.maximumReviewInterval=min(policy.maximumReviewInterval,budgetSteps*stepEstimate);
decision.predictedCoreCells=cores(end,:).*exp(-decision.maximumReviewInterval*decay);
if any(decision.predictedCoreCells<policy.predictedCoreBuffer)
    decision.action='regrid';decision.reason='forecast_core_trigger';
    decision.maximumAdvanceSteps=0;decision.maximumReviewInterval=0;return
end
if decision.maximumReviewInterval<policy.maximumReviewInterval
    decision.reason='shorten_horizon_for_step_budget';
end
end
