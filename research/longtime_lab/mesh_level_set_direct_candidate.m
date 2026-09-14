function [candidate,report] = mesh_level_set_direct_candidate(view,reference,anchor,policy)
%MESH_LEVEL_SET_DIRECT_CANDIDATE One geometry proposal from measured widths.
% Research-only, no LU, transfer, time advance, or controller mutation.
% The 90% wall core and the left-front interval set the fine platform;
% the 90% vertical width sets one rounded-geometric y axis. No tuple sweep.
feature=ipm.diagnostics.meshFeatureIntervals(view);
limits=policy.qualityLimits;
targetX=policy.targetCoreCells(1);
targetY=policy.targetCoreCells(2);
targetFront=policy.minimumFrontCells;
core=feature.coreInterval;
front=feature.frontInterval;
padding=policy.search.xPadding;
h=min(diff(core)/targetX,diff(front)/targetFront)/padding;
center=feature.coreCenter;
left=min(core(1),front(1));
right=max(core(2),front(2));
fineCells=max(8,ceil(diff(core)/h));
fineFraction=(center-left)/(right-left);
leftCells=max(1,min(fineCells-1,round(fineCells*fineFraction)));
rightCells=fineCells-leftCells;
fineFraction=leftCells/fineCells;
% If the outer slope reaches the registered ratio cap, this many rounded
% cells keep the initial log-spacing curvature below its registered cap.
rounding=max(8,ceil(log(limits.maxAdjacentCellRatio)/ ...
    limits.maxLogSpacingCurvature));
halfCells=(numel(reference.x)-1)/2;
report=struct('kind','level_set_direct_geometry_research_v1', ...
    'feature',feature,'fineSpacing',h,'fineCells',fineCells, ...
    'leftFineCells',leftCells,'rightFineCells',rightCells, ...
    'fineFraction',fineFraction,'roundingCells',rounding, ...
    'targetNodeCount',[numel(reference.x),numel(reference.y)], ...
    'predictedCells',[NaN,NaN,NaN], ...
    'quality',struct(),'reasons',{{}},'qualified',false, ...
    'noLU',true,'noTransfer',true);
candidate=struct([]);
if ~(isfinite(h)&&h>0&&left>0&&right<reference.x(end)&& ...
        fineCells+2*(rounding+2)<=halfCells)
    report.reasons={'direct_platform_capacity_or_geometry'};
    return
end
try
    [x,xInfo]=ipm.remesh.corePatchAxis(anchor,center,reference.x(end), ...
        numel(reference.x),h,fineCells,fineFraction,rounding, ...
        policy.search.maximumWarpFraction);
    [y,yInfo]=ipm.remesh.roundedGeometricAxis(reference.y(:), ...
        feature.yCoreWidth,targetY*policy.search.yPadding,rounding);
catch exception
    known=strcmp(exception.identifier,'ipm:RoundedGeometricInfeasible') || ...
        (strcmp(exception.identifier,'MATLAB:assertion:failed') && ...
        any(strcmp(exception.message, ...
        {'No feasible rounded-log cell split.','Insufficient rounded-branch cell budget.'})));
    if ~known,rethrow(exception);end
    report.reasons={exception.identifier};
    return
end
qx=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));
qy=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
counts=[interval_count(x(x>=0),core), ...
    interval_count(y,[0,feature.yCoreWidth]), ...
    interval_count(x(x>=0),front)];
reasons=[quality_reasons(qx,limits,'x'),quality_reasons(qy,limits,'y')];
if counts(1)<targetX-1e-6,reasons{end+1}='x_core';end
if counts(2)<targetY-1e-6,reasons{end+1}='y_core';end
if counts(3)<targetFront-1e-6,reasons{end+1}='x_front';end
if ~any(x==anchor)||~any(x==-anchor),reasons{end+1}='exact_anchor';end
report.predictedCells=counts;
report.quality=struct('x',qx,'y',qy);
report.construction=struct('x',xInfo,'y',yInfo);
report.reasons=reasons;
report.qualified=isempty(reasons);
candidate=struct('x',x,'y',y,'predictedCells',counts, ...
    'quality',report.quality,'qualified',report.qualified);
end

function count=interval_count(axis,interval)
v=axis(:)';
count=sum(max(0,min(v(2:end),interval(2))- ...
    max(v(1:end-1),interval(1)))./diff(v));
end

function reasons=quality_reasons(q,limits,prefix)
reasons={};
if q.maximumAdjacentCellRatio>limits.maxAdjacentCellRatio
    reasons{end+1}=[prefix '_adjacent_ratio'];
end
if q.maximumLogSpacingCurvature>limits.maxLogSpacingCurvature
    reasons{end+1}=[prefix '_log_curvature'];
end
if q.minimumStencilRcond<limits.minStencilRcond
    reasons{end+1}=[prefix '_stencil'];
end
if ~q.quadratureWeightsStrictlyPositive || ...
        q.minimumQuadratureWeightRatio<limits.minQuadratureWeightRatio
    reasons{end+1}=[prefix '_quadrature'];
end
if q.minimumQuadratureWeightToControlWidthRatio<limits.minWeightToControlWidth || ...
        q.maximumQuadratureWeightToControlWidthRatio>limits.maxWeightToControlWidth
    reasons{end+1}=[prefix '_local_weight'];
end
end
