function [state,stopReason,meshPlan] = advance(state)
%IPM.EVOLVE.ADVANCE Select dt, advance the configured RK scheme, and remesh.

meshPlan = [];
[dt,stopReason,timeStep] = ipm.evolve.selectTimestep( ...
    state.flow,state.scale,state.ops,state.config.time);
if ~isempty(stopReason)
    return;
end

% Commit a step only after its complete state is known to be finite.  This
% keeps the last accepted field available for diagnostics and persistence
% when an attempted Runge--Kutta stage overflows.
acceptedState = state;
initialRhsCache = [];
if isfield(state,'rhsCache')
    initialRhsCache = state.rhsCache;
end
switch state.config.time.timeIntegrator
    case 'ssprk3'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepSsprk3( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    case 'ssprk54'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepSsprk54( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    case 'rk6'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepRk6( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    otherwise
        error('ipm:TimeIntegrator','Unknown time integrator ''%s''.', ...
            state.config.time.timeIntegrator);
end
% normalizedTime is retained as a result/checkpoint compatibility alias.
state.normalizedTime = state.scale.canonicalTime;
state.step = state.step+1;

if ~ipm.evolve.isFinite(state)
    state = acceptedState;
    stopReason = 'non_finite_solution';
    return;
end

if isfield(state.config.remesh,'autonomousMesh') && ...
        state.config.remesh.autonomousMesh.enabled
    % Return the pure plan to solve, which owns the last external reference
    % to the old factor. No new LU may be constructed in this call scope.
    try
        [state,meshPlan] = ipm.evolve.planAutonomousMesh(state,false);
    catch exception
        failedStep = state.step;
        failureClock = [state.scale.canonicalTime,state.scale.physicalTime,state.normalizedTime];
        state = acceptedState;
        state.runMetadata.autonomousMesh.lastFailure = struct( ...
            'kind','planning_exception','attemptedStep',failedStep, ...
            'identifier',exception.identifier,'message',exception.message);
        if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
            state.runMetadata.autonomousMesh.lastFailure.sourceCanonicalTime=failureClock(1);
            state.runMetadata.autonomousMesh.lastFailure.sourcePhysicalTime=failureClock(2);
            state.runMetadata.autonomousMesh.lastFailure.sourceNormalizedTime=failureClock(3);
        end
        stopReason = 'autonomous_mesh_error';
        meshPlan = [];
        return;
    end
    if ~isempty(meshPlan.stopReason)
        failedPlan = rmfield(meshPlan,{'candidates','axisReport'});
        if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
            failedPlan.axisReport=meshPlan.axisReport;
        end
        state = acceptedState;
        state.runMetadata.autonomousMesh.lastFailure = failedPlan;
        stopReason = meshPlan.stopReason;
        meshPlan = [];
        return;
    end
    policy = state.config.remesh.autonomousMesh;
    if any([state.flow.coreGridPoints,state.flow.verticalCoreGridPoints] < ...
            policy.endpointCoreFloor)
        failureDecision=rmfield(meshPlan,{'candidates','axisReport'});
        state = acceptedState;
        if any(policy.version == [2,3,4,5])
            state.runMetadata.autonomousMesh.lastFailure=failureDecision;
            state.runMetadata.autonomousMesh.lastFailure.stopReason='autonomous_mesh_resolution_failure';
        end
        stopReason = 'autonomous_mesh_resolution_failure';
        meshPlan = [];
        return;
    end
    state.timeStep = timeStep;
    if ~meshPlan.requested
        meshPlan = [];
    end
    return;
end

[state.rho,state.ops,state.flow,remeshInfo,remeshedRhoRate] = ...
    ipm.evolve.remeshIfNeeded( ...
    state.rho,state.ops,state.flow,state.scale,state.config);
if remeshInfo.applied
    state.rhsCache = ipm.evolve.makeRhsCache( ...
        state.rho,remeshedRhoRate,state.flow,state.ops,state.scale);
end
if ~ipm.evolve.isFinite(state)
    state = acceptedState;
    stopReason = 'non_finite_solution';
    return;
end
if remeshInfo.applied
    state = ipm.evolve.publishRemeshWarnings(state);
end
state.timeStep = timeStep;
end
