function [state,success,attempts]=applyAutonomousMesh(state,plan)
%IPM.EVOLVE.APPLYAUTONOMOUSMESH Try planned axes on the original source only.
% Caller releases every reference to the old LU before entering this function.
% Failure returns the unchanged source data without a factor; the caller
% rebuilds the source factor from its exact stored matrix when rolling back.
assert(~isfield(state.ops,'poisson'),'ipm:AutonomousMeshFactorLifetime', ...
    'Release the old factor in the owning caller before building a candidate.');
policy=state.config.remesh.autonomousMesh;success=false;
if policy.version==2
    assert(plan.requested && isempty(plan.stopReason) && ...
        plan.sourceStep==state.step && plan.sourceCanonicalTime==state.scale.canonicalTime && ...
        plan.sourcePhysicalTime==state.scale.physicalTime && ...
        plan.sourceNormalizedTime==state.normalizedTime && ...
        plan.sourceRemeshCount==state.ops.remeshCount && ...
        plan.sourceLevelId==state.runMetadata.autonomousMesh.currentLevelId && ...
        isequal(plan.sourceNodeCount,[state.ops.nx,state.ops.ny]) && ...
        isequal(plan.coreCells,[state.flow.coreGridPoints,state.flow.verticalCoreGridPoints]), ...
        'ipm:AutonomousMeshStalePlan','A candidate plan must belong to the exact source state.');
end
attempts=struct('candidateIndex',{},'passed',{},'audit',{}, ...
    'errorIdentifier',{},'errorMessage',{});
for k=1:numel(plan.candidates)
    axes=plan.candidates(k);
    assert(plan.initial || (~axes.unchanged && ...
        (~isequal(axes.x,state.ops.x) || ~isequal(axes.y,state.ops.y))), ...
        'ipm:AutonomousMeshUnchangedTransaction', ...
        'An evolved mesh transaction must actually change at least one axis.');
    proposal=struct('x',axes.x,'y',axes.y);
    if policy.version==2 && ~plan.initial
        proposal.referenceFamily=state.runMetadata.autonomousMesh.referenceFamily;
        proposal.sourceLevelId=state.runMetadata.autonomousMesh.currentLevelId;
        proposal.targetLevelId=axes.nodeFamilyIndex;
    end
    row=struct('candidateIndex',k,'passed',false,'audit',struct(), ...
        'errorIdentifier','','errorMessage','');
    trial=state;
    try
        if plan.initial
            grid=state.config.grid;grid.customX=axes.x;grid.customY=axes.y;
            trial.ops=ipm.mesh.build(state.config,grid);
            trial.ops.remeshCount=0;
            trial.ops.baseX=trial.ops.x;trial.ops.baseY=trial.ops.y;
            trial.rho=ipm.field.initialDensity(trial.ops,state.config.physics);
            trial.ops=ipm.evolve.initializeScaling(trial.rho,trial.ops);
            [rhoRate,trial.flow]=ipm.evolve.flow(trial.rho,trial.ops,trial.scale);
            trial.mass0=sum(trial.rho.*trial.ops.integrationWeights,'all');
            trial.rhoRange0=[min(trial.rho,[],'all'),max(trial.rho,[],'all')];
            row.audit=initial_audit(trial,policy);
        else
            [trial.rho,trial.ops]=ipm.remesh.transfer(state.rho,state.ops,state.config,proposal);
            row.audit=ipm.remesh.auditCandidate(state,trial,policy, ...
                state.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump);
            if row.audit.passed
                [rhoRate,trial.flow]=ipm.evolve.flow(trial.rho,trial.ops,trial.scale);
            end
        end
        if row.audit.passed
            finite=ipm.evolve.isFinite(trial);
            safe=trial.flow.safetyFactor<policy.maximumSafety && ...
                all([trial.flow.coreGridPoints,trial.flow.verticalCoreGridPoints]>=policy.transactionMinimumCoreCells);
            if ~finite || ~safe
                row.audit.passed=false;
                row.audit.reasons{end+1}='post_transfer_flow_resolution_or_finiteness';
            end
        end
        row.passed=row.audit.passed;
        attempts(k)=row;
        if row.passed
            memory=state.runMetadata.autonomousMesh;
            if plan.initial
                memory.initialization=struct('sourceBaseX',state.ops.baseX,'sourceBaseY',state.ops.baseY, ...
                    'selectedCandidateIndex',k,'audit',row.audit,'attempts',attempts, ...
                    'analyticDatumResampled',true,'originalPhysicalTime',0, ...
                    'originalStep',0,'originalCanonicalTime',0,'originalRemeshCount',0, ...
                    'selectedBaseX',trial.ops.baseX,'selectedBaseY',trial.ops.baseY, ...
                    'operatorGridRepresentation','custom_axes', ...
                    'newPhysicalEpochCreated',false);
                if isfield(state.config.remesh,'initialMeshObservationFallback') && ...
                        state.config.remesh.initialMeshObservationFallback.enabled
                    if isfield(plan,'observationFallback')
                        evidence=plan.observationFallback;
                        evidence.selectedLocalCandidateIndex=k;
                        evidence.selectedOriginalPairIndex=evidence.admittedOriginalPairIndices(k);
                        evidence.status='native_initialization_accepted';
                    else
                        evidence=struct('version',1, ...
                            'policy',state.config.remesh.initialMeshObservationFallback, ...
                            'used',false,'attemptCount',0,'selectedPhase','original');
                    end
                    memory.initialization.observationFallback=evidence;
                end
                if policy.version==2
                    memory.referenceFamily=ipm.remesh.referenceAxisFamily( ...
                        trial.ops.baseX,trial.ops.baseY, ...
                        state.config.scaling.transportAnchorX,policy);
                    memory.currentLevelId=1;
                end
            else
                if policy.version==2
                    assert(row.audit.targetLevelId==axes.nodeFamilyIndex, ...
                        'ipm:AutonomousMeshFamily','The audited candidate must match its planned level.');
                    row.audit.controllerDecision=rmfield(plan,{'candidates','axisReport'});
                    row.audit.candidateIndex=k;
                    row.audit.attemptSummary=attempt_summary(attempts);
                    memory.currentLevelId=row.audit.targetLevelId;
                    attempts(k).audit=row.audit;
                end
                if isempty(memory.transactions),memory.transactions=row.audit;
                else,memory.transactions(end+1)=row.audit;end
                memory.cumulativeAbsolutePeakJump=sum([memory.transactions.relativePeakJump]);
                trial.ops.remeshCount=state.ops.remeshCount+1;
            end
            if plan.initial
                trial.runMetadata.gaugeContract=ipm.evolve.gaugeContract(trial.config,trial.ops);
            end
            trial.runMetadata.autonomousMesh=ipm.remesh.controllerTelemetry(memory,trial,policy);
            trial.rhsCache=ipm.evolve.makeRhsCache(trial.rho,rhoRate,trial.flow,trial.ops,trial.scale);
            state=trial;success=true;return
        end
    catch exception
        row.passed=false;
        if isfield(row.audit,'passed')
            row.audit.passed=false;
            row.audit.reasons{end+1}='candidate_exception';
        end
        row.errorIdentifier=exception.identifier;row.errorMessage=exception.message;
        attempts(k)=row;
        clear trial rhoRate
        if ~recoverable_geometry(exception.identifier)
            % Report a programming/configuration failure separately. Do not
            % disguise it as exhaustion of all geometrically valid candidates.
            return
        end
    end
    clear trial rhoRate
end
end

function summary=attempt_summary(attempts)
summary=repmat(struct('candidateIndex',0,'passed',false,'errorIdentifier','', ...
    'auditReasons',{{}}),1,numel(attempts));
for k=1:numel(attempts)
    row=attempts(k);summary(k).candidateIndex=row.candidateIndex;
    summary(k).passed=row.passed;summary(k).errorIdentifier=row.errorIdentifier;
    if isfield(row.audit,'reasons'),summary(k).auditReasons=row.audit.reasons;end
end
end

function audit=initial_audit(state,policy)
% Selecting the t=0 mesh resamples the analytic datum, not a migrated field.
% There is no inherited mass or evolution history to correct at this point.
ops=state.ops;
view=struct('rho',state.rho,'x',ops.x,'y',ops.y,'Dx',ops.Dx, ...
    'source',state.rho*ops.Dx','trusted',true);
features=ipm.diagnostics.meshFeatureIntervals(view);
q={ipm.mesh.quality(ops.x,7,ipm.mesh.quadrature(ops.x)), ...
   ipm.mesh.quality(ops.y,7,ipm.mesh.quadrature(ops.y))};
limits=policy.qualityLimits;qualityPassed=true;
for k=1:2
    a=q{k};qualityPassed=qualityPassed && a.quadratureWeightsStrictlyPositive && ...
        a.maximumAdjacentCellRatio<=limits.maxAdjacentCellRatio && ...
        a.maximumLogSpacingCurvature<=limits.maxLogSpacingCurvature && ...
        a.minimumStencilRcond>=limits.minStencilRcond && ...
        a.minimumQuadratureWeightRatio>=limits.minQuadratureWeightRatio && ...
        a.minimumQuadratureWeightToControlWidthRatio>=limits.minWeightToControlWidth && ...
        a.maximumQuadratureWeightToControlWidthRatio<=limits.maxWeightToControlWidth;
end
anchor=state.config.scaling.transportAnchorX;
zero=state.step==0 && state.scale.canonicalTime==0 && state.scale.physicalTime==0 && state.normalizedTime==0;
counts=all(features.actualCoreCells>=policy.transactionMinimumCoreCells) && ...
    features.leftFrontCells>=policy.minimumFrontCells;
passed=zero && counts && qualityPassed && any(ops.x==anchor) && any(ops.x==-anchor);
audit=struct('kind','analytic_zero_time_axis_selection_audit','passed',passed,'reasons',{{}}, ...
    'coreCells',features.actualCoreCells,'leftFrontCells',features.leftFrontCells, ...
    'qualityPassed',qualityPassed,'exactZeroTime',zero,'xQuality',q{1},'yQuality',q{2}, ...
    'massMigrationClaim',false,'initialDensityExactlyResampled',true);
if ~passed,audit.reasons={'initial_axis_or_analytic_resolution_gate'};end
end
function yes=recoverable_geometry(identifier)
yes=any(strcmp(identifier,{'ipm:HighOrderQuadratureWeights','ipm:HighOrderGridMetric', ...
    'ipm:FiniteDifferenceScale','ipm:HighOrderRemeshBubble'}));
end
