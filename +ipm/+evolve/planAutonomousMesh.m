function [state,plan]=planAutonomousMesh(state,initial)
%IPM.EVOLVE.PLANAUTONOMOUSMESH Decide and prepare axes without a new LU.
policy=state.config.remesh.autonomousMesh;
assert(policy.enabled,'ipm:AutonomousMeshPlan','An enabled policy is required.');
if initial
    memory=[];
else
    assert(isfield(state.runMetadata,'autonomousMesh'), ...
        'ipm:AutonomousMeshMemory','An evolved state must retain its original controller ledger.');
    memory=state.runMetadata.autonomousMesh;
end
memory=ipm.remesh.controllerTelemetry(memory,state,policy);
w=memory.window;core=w.coreCells(end,:);decay=[0,0];
trend=numel(w.time)>=4;
if trend
    for k=1:2
        fit=polyfit(w.time-w.time(end),log(w.coreCells(:,k)),1);
        secants=-diff(log(w.coreCells(:,k)))./diff(w.time);
        decay(k)=max([0;-fit(1);secants]);
    end
end
predicted=core.*exp(-policy.maximumReviewInterval*decay);
reason='resolved';requested=false;
if initial
    requested=true;reason='analytic_zero_time_mesh_selection';
elseif any(core<policy.regridCoreTrigger)
    requested=true;reason='current_core_trigger';
elseif trend && any(predicted<policy.predictedCoreBuffer & core<policy.targetCoreCells)
    requested=true;reason='forecast_core_trigger';
elseif state.flow.safetyFactor>=policy.maximumSafety
    requested=true;reason='safety_buffer';
end
plan=struct('initial',initial,'requested',requested,'reason',reason,'candidates',struct([]), ...
    'axisReport',struct(),'trendAvailable',trend,'coreCells',core,'predictedCoreCells',predicted, ...
    'decayEstimate',decay,'stopReason','');
if policy.version==4
    plan.axisSearchEvidence=ipm.remesh.searchEvidence('empty',policy,initial);
end
if any(policy.version == [2,3,4])
    plan.sourceStep=state.step;plan.sourceCanonicalTime=state.scale.canonicalTime;
    plan.sourcePhysicalTime=state.scale.physicalTime;plan.sourceNormalizedTime=state.normalizedTime;
    plan.sourceRemeshCount=state.ops.remeshCount;plan.sourceLevelId=memory.currentLevelId;
    plan.sourceNodeCount=[state.ops.nx,state.ops.ny];
    plan.decayEvidence=decay_evidence(w,trend);
end
if ~initial && (any(core<policy.historyCoreFloor) || state.flow.safetyFactor>=policy.maximumSafety)
    plan.stopReason='autonomous_mesh_resolution_failure';plan.requested=false;
elseif requested && state.ops.remeshCount>=state.config.remesh.maxRemeshes && ~initial
    plan.stopReason='autonomous_mesh_count_budget';plan.requested=false;
elseif requested
    ops=state.ops;
    view=struct('rho',state.rho,'x',ops.x,'y',ops.y,'Dx',ops.Dx, ...
        'source',state.rho*ops.Dx','trusted',true);
    if any(policy.version == [2,3,4]) && ~initial
        if policy.version==4
            [plan.candidates,plan.axisReport,plan.axisSearchEvidence]=family_candidates( ...
                view,memory,state.config.scaling.transportAnchorX,policy);
        else
        [plan.candidates,plan.axisReport]=family_candidates( ...
            view,memory,state.config.scaling.transportAnchorX,policy);
        end
    else
        reference=struct('x',ops.baseX,'y',ops.baseY);
        [plan.candidates,plan.axisReport]=ipm.remesh.plannedAxisPairs( ...
            view,reference,state.config.scaling.transportAnchorX,policy);
        if any(policy.version == [2,3,4])
            for k=1:numel(plan.candidates),plan.candidates(k).nodeFamilyIndex=1;end
        end
    end
    if ~initial && ~isempty(plan.candidates)
        unchanged=[plan.candidates.unchanged];
        plan.axisReport.unchangedPairsExcludedFromEvolution=sum(unchanged);
        plan.candidates=plan.candidates(~unchanged);
    end
    if isempty(plan.candidates),plan.stopReason='autonomous_mesh_axis_capacity';end
end
if policy.version==4
    ipm.remesh.searchEvidence('validate',plan.axisSearchEvidence,policy, ...
        memory.referenceFamily,memory.currentLevelId,initial,plan.requested);
end
memory.lastDecision=rmfield(plan,{'candidates','axisReport'});
state.runMetadata.autonomousMesh=memory;
end

function rows=decay_evidence(w,trend)
% Record the actual endpoints behind the most restrictive observed secant.
rows=repmat(struct('axis',0,'fitDecay',0,'maximumSecant',0, ...
    'secantFromStep',w.step(end),'secantToStep',w.step(end), ...
    'secantFromTime',w.time(end),'secantToTime',w.time(end)),1,2);
for k=1:2
    rows(k).axis=k;
    if trend
        fit=polyfit(w.time-w.time(end),log(w.coreCells(:,k)),1);
        [rate,index]=max(-diff(log(w.coreCells(:,k)))./diff(w.time));
        rows(k).fitDecay=-fit(1);rows(k).maximumSecant=rate;
        rows(k).secantFromStep=w.step(index);rows(k).secantToStep=w.step(index+1);
        rows(k).secantFromTime=w.time(index);rows(k).secantToTime=w.time(index+1);
    end
end
end

function [candidates,report,evidence]=family_candidates(view,memory,anchor,policy)
% Each level receives its own registered pair budget. An exhausted smaller
% level cannot consume the attempts reserved for a larger admitted member.
family=memory.referenceFamily;source=family.members(memory.currentLevelId);
factors=vertcat(family.members.cellFactors);counts=vertcat(family.members.nodeCount);
eligible=find(all(factors>=source.cellFactors,2) & [family.members.resourceAdmitted]' & ...
    [family.members.qualityPassed]');
[~,order]=sortrows([prod(counts(eligible,:),2),eligible],[1,2]);eligible=eligible(order);
assert(~isempty(eligible) && eligible(1)==memory.currentLevelId, ...
    'ipm:AutonomousMeshFamily','The admitted current level must be planned first.');
candidates=struct([]);reports=struct([]);evidence=[];
if policy.version==4
    evidence=ipm.remesh.searchEvidence('empty',policy,false);evidence.phase='evolved_requested';
end
for index=eligible(:)'
    member=family.members(index);
    if policy.version==4
        [pairs,detail,searchAudit]=ipm.remesh.hierarchicalAxisPairs(view, ...
            struct('x',member.baseX,'y',member.baseY),anchor,policy,policy.nodeFamily.maximumTotalNodes);
        phase='primary';if searchAudit.fallbackUsed,phase='refined';end
        summary=struct('levelId',index,'nodeCount',member.nodeCount,'phase',phase, ...
            'evaluatedAxisTrials',searchAudit.axisTrialsEvaluated, ...
            'primaryPairCount',numel(searchAudit.primaryCandidates), ...
            'primaryYIndices',find([searchAudit.primaryReport.yTrials.admissible]), ...
            'completedRefinementLevel',searchAudit.lastCompletedLevel, ...
            'returnedPairCount',numel(pairs),'status',searchAudit.status);
        if isempty(evidence.members),evidence.members=summary;else,evidence.members(end+1)=summary;end
        for k=1:numel(pairs)
            descriptor=struct('targetLevelId',index,'localPairIndex',k,'unchanged',pairs(k).unchanged, ...
                'searchPhase',phase,'refinementLevel',searchAudit.lastCompletedLevel, ...
                'xScheduleIndex',pairs(k).xIndex,'yTrialIndex',pairs(k).yIndex);
            pairs(k).searchDescriptor=descriptor;
            if isempty(evidence.candidates),evidence.candidates=descriptor;else,evidence.candidates(end+1)=descriptor;end
        end
    else
    [pairs,detail]=ipm.remesh.plannedAxisPairs(view, ...
        struct('x',member.baseX,'y',member.baseY),anchor,policy);
    end
    for k=1:numel(pairs),pairs(k).nodeFamilyIndex=index;end
    if any(policy.version == [3,4])
        for k=1:numel(pairs)
            pairs(k).targetLevelId=index;pairs(k).localPairIndex=k;
        end
    end
    if ~isempty(pairs)
        if isempty(candidates),candidates=pairs;else,candidates=[candidates,pairs];end %#ok<AGROW>
    end
    row=struct('levelId',index,'nodeCount',member.nodeCount, ...
        'status',detail.status,'xAdmitted',sum([detail.xTrials.admissible]), ...
        'yAdmitted',sum([detail.yTrials.admissible]),'proposedPairs',numel(pairs), ...
        'xRejectionReasons',{{detail.xTrials.reasons}}, ...
        'yRejectionReasons',{{detail.yTrials.reasons}});
    if isempty(reports),reports=row;else,reports(end+1)=row;end %#ok<AGROW>
end
report=struct('kind','registered_node_family_planning_v2','sourceLevelId',memory.currentLevelId, ...
    'eligibleLevelOrder',eligible(:)','maximumPairsPerLevel',policy.search.maximumPairCandidates, ...
    'levels',reports,'noLU',true,'noTransfer',true);
if policy.version == 3,report.kind='registered_node_family_planning_v3';end
if policy.version==4
    report.kind='registered_node_family_planning_v4';
    if ~isempty(evidence.candidates),evidence.filteredCandidateIndices=find(~[evidence.candidates.unchanged]);end
    report.fullTrialArtifactWritten=false;
    report.evidenceBoundary='Compact identities persisted; full trial construction arrays are transient request-time diagnostics.';
end
end
