function [checkpoint,audit,Dx,Dy] = ipm_accellab_verify_rejected_source(source,saved)
%IPM_ACCELLAB_VERIFY_REJECTED_SOURCE Reproduce the complete rejected array.
%   Five strict native reads and one-dimensional derivative matrices only.
assert(~source.trial.accepted && isscalar(source.trial.attempts) && ...
    numel(source.checkpointFiles) == 5 && all(source.signatureValidated), ...
    'ipm:RejectedSourceProtocol','The original fixed five-frame rejected experiment is required.');
t = source.trial;
assert(isequaln(saved.report.sourceResidualAudit,t.attempts.audit) && ...
    isequaln(saved.report.guardFailures,t.attempts.guardFailures), ...
    'ipm:RejectedSourceAudit','Materialized field audit does not match the original evaluated rejection.');
assert(t.resolvedSecantOptions.memory == 4 && t.attempts.damping == 0.0625 && ...
    numel(t.coefficients) == 4 && source.options.expectedRemeshCount == 3, ...
    'ipm:RejectedSourceProtocol','Recorded memory, damping or remesh registration changed.');
states = cell(1,5);
for k = 1:5
    checkpoint = ipm.output.readCheckpoint(source.checkpointFiles{k});
    states{k} = checkpoint.payload.state;
end
reference = states{end};
assert(strcmp(reference.config.transport.spatialDiscretization,'high_order') && ...
    strcmp(reference.config.scaling.dynamicScaleGeometry,'isotropic'), ...
    'ipm:RejectedSourceDerivative','This fixed observer comparison requires maintained isotropic seven-node derivatives.');
Dx = ipm.mesh.fdMatrix(reference.x,1,7); Dy = ipm.mesh.fdMatrix(reference.y,1,7);
history = zeros(numel(reference.rho),5); measuredPeaks = zeros(1,5);
for k = 1:5
    state = states{k};
    for name = {'x','y','baseX','baseY','rescaling','remeshCount','mass0','rhoRange0'}
        assert(isequaln(state.(name{1}),reference.(name{1})), ...
            'ipm:RejectedSourceGrid','Native frame %d changed %s.',k,name{1});
    end
    assert(state.remeshCount == 3 && state.step == source.sourceSteps(k) && ...
        state.scale.canonicalTime == source.sourceCanonicalTimes(k), ...
        'ipm:RejectedSourceIdentity','Frame step/time/remesh identity differs from the original experiment.');
    a = rmfield(state.config,'output'); b = rmfield(reference.config,'output');
    a.time = rmfield(a.time,{'finalTime','physicalFinalTime','maxSteps'});
    b.time = rmfield(b.time,{'finalTime','physicalFinalTime','maxSteps'});
    assert(isequaln(a,b),'ipm:RejectedSourceNumerics','Native numerical configuration changed.');
    allowed = {'latestCheckpointFile','resumeCount','resumedFromStep', ...
        'resumedFromCanonicalTime','resumedFromPhysicalTime','resumedFromCheckpoint', ...
        'resultFile','videoFile','caseMetadata'};
    a = rmfield(state.runMetadata,intersect(fieldnames(state.runMetadata),allowed));
    b = rmfield(reference.runMetadata,intersect(fieldnames(reference.runMetadata),allowed));
    assert(isequaln(a,b),'ipm:RejectedSourceLineage','Native case/gauge lineage changed.');
    if k > 1
        assert(state.scale.canonicalTime > states{k-1}.scale.canonicalTime, ...
            'ipm:RejectedSourceTime','Native times are not strictly ordered.');
    end
    measuredPeaks(k) = peak_value(state.rho,reference,Dx);
    history(:,k) = t.historyAmplitudeFactors(k)*state.rho(:);
end
assert(abs(measuredPeaks(end)-t.baselineAudit.peak) < 1e-13 && ...
    max(abs(measuredPeaks(end)./measuredPeaks-t.historyAmplitudeFactors)) < 1e-13, ...
    'ipm:RejectedSourceAmplitude','Recorded baseline/history amplitude normalization is not reproduced.');
direction = -diff(history,1,2)*t.coefficients;
vector = history(:,end)+t.attempts.damping*direction;
rho = reshape(vector,size(reference.rho));
vector = (t.baselineAudit.peak/peak_value(rho,reference,Dx))*vector;
assert(isequal(size(saved.rhoTrial),size(reference.rho)) && ...
    isequal(saved.x,reference.x) && isequal(saved.y,reference.y), ...
    'ipm:RejectedSourceSavedGrid','Materialized rejected field axes or shape differ.');
maximumDifference = max(abs(vector-saved.rhoTrial(:)));
tolerance = 64*eps(max(1,max(abs(vector))));
assert(maximumDifference <= tolerance,'ipm:RejectedSourceField', ...
    'The saved rejected field does not reproduce the recorded full-array formula.');
jump = norm(vector-history(:,end))/max(norm(history(:,end)),realmin);
assert(abs(jump-t.attempts.relativeJump) < 1e-13, ...
    'ipm:RejectedSourceJump','The recorded attempt jump is not reproduced.');
audit = struct('nativeSignaturesPassed',true(1,5),'nativeSteps',source.sourceSteps, ...
    'nativeCanonicalTimes',source.sourceCanonicalTimes,'completeArrayExact',isequal(vector,saved.rhoTrial(:)), ...
    'completeArrayMaximumDifference',maximumDifference,'completeArrayTolerance',tolerance, ...
    'reproducedJump',jump,'measuredQuadraticPeaks',measuredPeaks, ...
    'restoreCalls',0,'newRhsEvaluations',0,'originalDecision','rejected_unchanged');
end

function peak = peak_value(rho,state,Dx)
omega = rho*Dx';
tracking = abs(state.x-state.rescaling.pinX) <= state.rescaling.peakTrackingHalfWidth;
indices = find(tracking); [~,index] = max(omega(1,tracking));
functional = ipm.evolve.quadraticPeakFunctional(omega(1,:),state.x,indices(index));
peak = functional.value;
end
