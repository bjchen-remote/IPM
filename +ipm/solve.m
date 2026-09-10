function result = solve(userOpts,restartCheckpoint)
%IPM.SOLVE Solve the physical 2-D incompressible porous media equation.
%   RESULT = IPM.SOLVE(OPTIONS) is the only solver entry point. It resolves
%   configuration once, advances a shared numerical path, and returns the
%   grouped version-2 result contract.
%   RESULT = IPM.SOLVE(OVERRIDES,CHECKPOINT) resumes a validated same-grid
%   accepted-state checkpoint. Only terminal horizons and output controls
%   may be overridden; every numerical and physical choice remains frozen.
%
%   The maintained equation is
%       rho_t + u dot grad(rho) = 0,
%       u = grad^perp (-Delta)^(-1) d_x rho
%   on a nonperiodic truncation of R x R_+. Isotropic scaling is the exact
%   C_x=C_y specialization of the two-scale geometry.

if nargin < 1
    userOpts = [];
end
if nargin < 2 || isempty(restartCheckpoint)
    state = ipm.evolve.initialize(userOpts);
    log = ipm.output.initializeLog();
    [log,~] = ipm.output.record(log,state);
    cursor = initial_cursor(state.config);
else
    [state,log,cursor] = ipm.output.restoreCheckpoint( ...
        restartCheckpoint,userOpts);
end
time = state.config.time;
diagnostics = state.config.diagnostics;
remesh = state.config.remesh;
output = state.config.output;

viz = [];
if output.livePlot || output.writeVideo
    viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
        state.scale,log.history,output,state.runMetadata);
    if output.writeVideo && isfield(viz,'videoFile')
        state.runMetadata.videoFile = string(viz.videoFile);
    end
end
videoCleanup = onCleanup(@()ipm.output.closeLive(viz));

stopReason = 'final_time';
while ipm.evolve.isActive(state)
    [state,stepStopReason,meshPlan] = ipm.evolve.advance(state);
    if ~isempty(meshPlan)
        % Own the factor lifetime here: advance has returned, so neither its
        % pre-step rollback state nor this state can retain a second old LU.
        state.ops = rmfield(state.ops,'poisson');
        [state,applied,attempts] = ipm.evolve.applyAutonomousMesh(state,meshPlan);
        if ~applied
            state.ops.poisson = decomposition(state.ops.A,'lu');
            state.runMetadata.autonomousMesh.lastFailure = struct( ...
                'kind','candidate_transaction_failure','attempts',attempts, ...
                'step',state.step,'canonicalTime',state.scale.canonicalTime);
            if any(state.config.remesh.autonomousMesh.version == [2,3,4])
                state.runMetadata.autonomousMesh.lastFailure.controllerDecision= ...
                    rmfield(meshPlan,{'candidates','axisReport'});
                state.runMetadata.autonomousMesh.lastFailure.axisReport=meshPlan.axisReport;
            end
            stepStopReason = 'autonomous_mesh_candidates_rejected';
            if ~isempty(attempts) && ~isempty(attempts(end).errorIdentifier)
                stepStopReason = 'autonomous_mesh_error';
            end
        end
    end
    if ~isempty(stepStopReason)
        stopReason = stepStopReason;
        [log,recorded] = ipm.output.recordIfNew(log,state);
        [state,cursor] = ipm.output.maybeCheckpoint( ...
            state,log,cursor,true);
        if recorded && ~isempty(viz)
            viz = ipm.output.visualize(viz,state.rho,state.flow, ...
                state.ops,state.scale,log.history,output,state.runMetadata);
        end
        break;
    end

    if ~ipm.output.shouldRecord(state,cursor.nextOutput)
        continue;
    end
    [log,~] = ipm.output.record(log,state);
    cursor.nextOutput = cursor.nextOutput+time.outputEvery;
    stopReasonAtOutput = ipm.diagnostics.stopPolicy(log.history, ...
        log.history.common.physicalGradInf(end),state.flow,state.ops, ...
        diagnostics,remesh);
    [state,cursor] = ipm.output.maybeCheckpoint( ...
        state,log,cursor,~isempty(stopReasonAtOutput));
    if ~isempty(viz)
        viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
            state.scale,log.history,output,state.runMetadata);
    end
    if ~isempty(stopReasonAtOutput)
        stopReason = stopReasonAtOutput;
        break;
    end
end

if strcmp(stopReason,'final_time')
    if state.scale.physicalTime >= time.physicalFinalTime
        stopReason = 'physical_final_time';
    elseif state.step >= time.maxSteps && ...
            state.scale.canonicalTime < time.finalTime
        stopReason = 'maximum_steps';
    end
end

% Always preserve the last finite field, even when it falls between output
% intervals. Snapshot grids remain paired with their corresponding fields.
[log,recorded] = ipm.output.recordIfNew(log,state);
[state,~] = ipm.output.maybeCheckpoint(state,log,cursor,true);
if recorded && ~isempty(viz)
    viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
        state.scale,log.history,output,state.runMetadata);
end

result = ipm.output.finalize(state,log,stopReason);
result = ipm.output.write(result,output);
ipm.output.report(result);
if output.makePlots && isempty(viz)
    ipm.output.plotResult(result);
end
end

function cursor = initial_cursor(config)
nextCheckpoint = Inf;
if isfield(config.output,'checkpoint') && ...
        config.output.checkpoint.enabled
    nextCheckpoint = config.output.checkpoint.every;
end
cursor = struct('nextOutput',config.time.outputEvery, ...
    'nextCheckpoint',nextCheckpoint,'lastCheckpointStep',-1, ...
    'lastCheckpointCanonicalTime',NaN);
end
