function [rho,ops,flow,info,rhoRate] = ...
    remeshIfNeeded(rho,ops,flow,scale,config)
%IPM.EVOLVE.REMESHIFNEEDED Apply one accepted adaptive-grid proposal when due.

remesh = config.remesh;
output = config.output;
info = struct('applied',false);
rhoRate = [];
requested = remesh.adaptiveRemesh && isfinite(flow.safetyFactor) && ...
    flow.safetyFactor > remesh.remeshSafetyTrigger && ...
    ops.remeshCount < remesh.maxRemeshes;
if ~requested
    return;
end

[rho,ops,info] = ipm.remesh.adapt( ...
    rho,ops,flow,config,exp(scale.logC_omega));
if ~info.applied
    return;
end
[rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
if output.verbose
    fprintf(['  adaptive remesh %d at physical t=%.6g: ' ...
        'peak dx %.3e -> %.3e, wall dy %.3e -> %.3e\n'], ...
        info.count,scale.physicalTime,info.oldPeakSpacing, ...
        info.newPeakSpacing,info.oldWallSpacing,info.newWallSpacing);
end
end
