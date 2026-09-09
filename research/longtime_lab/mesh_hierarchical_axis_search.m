function [candidates,report]=mesh_hierarchical_axis_search(view,referenceAxes,anchor,policy,maximumTotalNodes)
%MESH_HIERARCHICAL_AXIS_SEARCH Research-only bounded pure axis refinement.
% No clock, step number, field transfer, flow, checkpoint or initial sampler.
validateattributes(maximumTotalNodes,{'numeric'},{'scalar','integer','positive'});
assert(numel(referenceAxes.x)*numel(referenceAxes.y)<=maximumTotalNodes, ...
    'ipm:HierarchicalResourceCap','Reference exceeds the separately registered node cap.');
[schedule,registration]=mesh_hierarchical_axis_schedule(policy.search);
registration.maximumTotalNodes=maximumTotalNodes;
[candidates,primary]=ipm.remesh.plannedAxisPairs(view,referenceAxes,anchor,policy);
assert(numel(primary.xTrials)==nnz(schedule(:,1)==0));
report=struct('kind','research_hierarchical_axis_search_v1','registration',registration, ...
    'primaryReport',primary,'primaryCandidates',candidates,'schedule',schedule, ...
    'extraTrials',struct([]),'stages',struct([]),'fallbackUsed',false, ...
    'axisTrialsEvaluated',numel(primary.xTrials),'lastCompletedLevel',0, ...
    'status','primary_qualified','noLU',true,'noFlow',true,'noTransfer',true, ...
    'resourceNodeCount',[numel(referenceAxes.x),numel(referenceAxes.y)], ...
    'globalFeasibilityClaim',false);
if ~isempty(candidates),return;end
yIndices=find([primary.yTrials.admissible]);
if isempty(yIndices),report.status='primary_y_family_exhausted_no_x_refinement';return;end
report.fallbackUsed=true;f=primary.feature;limits=policy.qualityLimits;
spacing=min(diff(f.coreInterval)/policy.targetCoreCells(1), ...
    diff(f.frontInterval)/policy.minimumFrontCells)/policy.search.xPadding;
n=numel(referenceAxes.x);passedAxes={};passedRows=struct([]);passedIndices=[];
for level=1:registration.refinementDepth
    indices=find(schedule(:,1)==level);before=numel(report.extraTrials);
    for i=indices(:)'
        tuple=schedule(i,:);fine=tuple(2);rounding=tuple(3);fraction=tuple(4);
        row=struct('positiveFineCells',fine,'roundingCells',rounding, ...
            'coreFineCellFraction',fraction,'admissible',false,'reasons',{{}}, ...
            'quality',struct(),'coreCells',NaN,'frontCells',NaN,'qualityMargin',-Inf, ...
            'construction',struct(),'failure',struct());axis=[];
        if tuple(5)>(n-1)/2
            row.reasons={'rounded_branch_node_budget'};
        else
            try
                [axis,info]=ipm.remesh.corePatchAxis(anchor,f.coreCenter,view.x(end),n, ...
                    spacing,fine,fraction,rounding,policy.search.maximumWarpFraction);
                [q,reasons,margin]=axis_quality(axis,limits);
                core=interval_count(axis(axis>=0),f.coreInterval);
                front=interval_count(axis(axis>=0),f.frontInterval);
                if core<policy.targetCoreCells(1)-1e-6,reasons{end+1}='x_core_target';end %#ok<AGROW>
                if front<policy.minimumFrontCells-1e-6,reasons{end+1}='x_front_target';end %#ok<AGROW>
                if ~any(axis==anchor)||~any(axis==-anchor),reasons{end+1}='exact_anchor';end %#ok<AGROW>
                row.quality=q;row.reasons=reasons;row.qualityMargin=margin;
                row.coreCells=core;row.frontCells=front;row.construction=info;
                row.admissible=isempty(reasons);
            catch e
                known_geometry_rejection(e);
                row.reasons={'geometry_construction_rejected'};
                row.failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
            end
        end
        report.extraTrials=append(report.extraTrials,struct('scheduleIndex',i,'level',level,'row',row));
        if row.admissible
            passedAxes{end+1}=axis;passedRows=append(passedRows,row);passedIndices(end+1)=i; %#ok<AGROW>
        end
    end
    report.axisTrialsEvaluated=numel(primary.xTrials)+numel(report.extraTrials);
    report.lastCompletedLevel=level;
    report.stages=append(report.stages,struct('level',level, ...
        'newTrials',numel(report.extraTrials)-before,'xAdmitted',numel(passedRows), ...
        'completeScheduledLevel',true));
    if ~isempty(passedRows)
        pairs=zeros(numel(passedRows)*numel(yIndices),4);k=0;
        for i=1:numel(passedRows)
            for j=yIndices
                k=k+1;pairs(k,:)=[i,j,min(passedRows(i).qualityMargin,primary.yTrials(j).qualityMargin),passedIndices(i)];
            end
        end
        pairs=sortrows(pairs,[-3,4,2]);pairs=pairs(1:min(3,size(pairs,1)),:);
        candidates=struct([]);meshLimits=rmfield(limits,{'minWeightToControlWidth','maxWeightToControlWidth'});
        for k=1:size(pairs,1)
            i=pairs(k,1);j=pairs(k,2);sigma=primary.yTrials(j).sigma;
            [y,~]=ipm.remesh.equalizedAxis(referenceAxes.y(:)/anchor, ...
                [0,f.yCoreWidth/anchor],policy.targetCoreCells(2)*policy.search.yPadding, ...
                struct('geometry','positive','focusCenters',0,'sigma',sigma,'power',2,'meshLimits',meshLimits));
            y=y(:)*anchor;y(1)=view.y(1);y(end)=view.y(end);
            [qy,reasons]=axis_quality(y,limits);assert(isempty(reasons)&&isequaln(qy,primary.yTrials(j).quality));
            row=struct('x',passedAxes{i},'y',y,'unchanged',false,'xIndex',passedIndices(i), ...
                'yIndex',j,'predictedCells',[passedRows(i).coreCells,primary.yTrials(j).coreCells,passedRows(i).frontCells], ...
                'quality',struct('x',passedRows(i).quality,'y',qy));
            candidates=append(candidates,row);
        end
        report.status='refined_axes_proposed_require_actual_transfer';return;
    end
end
report.status='registered_search_budget_exhausted_not_global_infeasibility';
end

function [q,reasons,margin]=axis_quality(axis,limits)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
v=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio,q.maximumQuadratureWeightToControlWidthRatio];
assert(all(isfinite(v)),'ipm:HierarchicalNonfiniteQuality','Nonfinite quality cannot enter the ranking.');
bad=[v(1)>limits.maxAdjacentCellRatio,v(2)>limits.maxLogSpacingCurvature,v(3)<limits.minStencilRcond, ...
    v(4)<limits.minQuadratureWeightRatio,v(5)<limits.minWeightToControlWidth,v(6)>limits.maxWeightToControlWidth];
names={'adjacent_ratio','log_spacing_curvature','stencil_rcond','global_quadrature','local_quadrature_min','local_quadrature_max'};
reasons=names(bad);if ~q.quadratureWeightsStrictlyPositive,reasons{end+1}='nonpositive_quadrature';end
margin=min([limits.maxAdjacentCellRatio/v(1),limits.maxLogSpacingCurvature/max(v(2),realmin), ...
    v(3)/limits.minStencilRcond,v(4)/limits.minQuadratureWeightRatio, ...
    v(5)/limits.minWeightToControlWidth,limits.maxWeightToControlWidth/v(6)]);
end
function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end
function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function known_geometry_rejection(e)
known={'No feasible rounded-log cell split.','Insufficient rounded-branch cell budget.'};
fromFactory=~isempty(e.stack)&&any(strcmp(e.stack(1).name, ...
    {'ipm.remesh.corePatchAxis','ipm.remesh.roundedAxis','ipm.remesh.roundedAxis>growth', ...
    'corePatchAxis','roundedAxis','roundedAxis>growth'}));
if ~(fromFactory&&strcmp(e.identifier,'MATLAB:assertion:failed'))&&~any(strcmp(e.message,known)),rethrow(e);end
end
