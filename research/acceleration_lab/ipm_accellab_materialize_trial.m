function report = ipm_accellab_materialize_trial(experimentFile)
%IPM_ACCELLAB_MATERIALIZE_TRIAL Save a tested field from its recorded formula.
%   Native signature reads and a one-dimensional sparse Dx only. No restore,
%   Poisson construction/solve, fresh residual call, or PDE step is performed.
directory = fileparts(mfilename('fullpath'));
addpath(directory,fileparts(fileparts(directory)));
loaded = load(experimentFile,'experiment'); experiment = loaded.experiment;
t = experiment.trial;
assert(numel(t.attempts) == 1 && ~t.accepted && numel(experiment.checkpointFiles) == 5, ...
    'ipm:MaterializeProtocol','This helper preserves the fixed five-frame rejected attempt only.');
assert(maxNumCompThreads == 10,'ipm:MaterializeThreads','Native frames require their original ten-thread environment.');
history = [];
for k = 1:5
    checkpoint = ipm.output.readCheckpoint(experiment.checkpointFiles{k});
    saved = checkpoint.payload.state;
    assert(strcmp(saved.config.transport.spatialDiscretization,'high_order'), ...
        'ipm:MaterializeDerivative','Only the recorded seven-node high-order Dx is supported.');
    if k == 1
        history = zeros(numel(saved.rho),5);
        x = saved.x; y = saved.y;
    else
        assert(isequal(saved.x,x) && isequal(saved.y,y),'ipm:MaterializeGrid','Saved axes changed.');
    end
    history(:,k) = t.historyAmplitudeFactors(k)*saved.rho(:);
end
direction = -diff(history,1,2)*t.coefficients;
vector = history(:,end)+t.attempts.damping*direction;
rhoTrial = reshape(vector,size(saved.rho));
Dx = ipm.mesh.fdMatrix(x,1,7);
omega = rhoTrial*Dx';
tracking = abs(x-saved.rescaling.pinX) <= saved.rescaling.peakTrackingHalfWidth;
indices = find(tracking); [~,index] = max(omega(1,tracking));
peak = ipm.evolve.quadraticPeakFunctional(omega(1,:),x,indices(index));
vector = (t.baselineAudit.peak/peak.value)*vector;
rhoTrial = reshape(vector,size(saved.rho));
relativeJump = norm(vector-history(:,end))/max(norm(history(:,end)),realmin);
assert(abs(relativeJump-t.attempts.relativeJump) < 1e-13, ...
    'ipm:MaterializeReproduction','Reconstructed field jump differs from the evaluated attempt.');
[~,token] = fileparts(tempname);
destination = fullfile(experiment.outputDirectory,['materialized_rejection_',token]);
mkdir(destination);
report = struct('kind','recorded_rejected_full_field_reconstruction', ...
    'experimentFile',experimentFile,'outputDirectory',destination, ...
    'relativeJump',relativeJump,'recordedRelativeJump',t.attempts.relativeJump, ...
    'freshResidualComputed',false,'sourceResidualAudit',t.attempts.audit, ...
    'guardFailures',{t.attempts.guardFailures},'restoreCalls',0,'newRhsEvaluations',0, ...
    'pdeAcceptedSteps',0,'accepted',false,'isNativeCheckpoint',false);
save(fullfile(destination,'rejected_trial.mat'),'rhoTrial','x','y','report','-v7.3');
fprintf('REJECTED_FIELD_SAVED %s jump %.12g\n',destination,relativeJump);
end
