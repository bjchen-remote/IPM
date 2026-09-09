function [candidate,report]=mesh_design_paired_snapshot(snapshot,referenceAxes,actualAnchor,controls)
%MESH_DESIGN_PAIRED_SNAPSHOT Bounded no-LU same-box/same-N mesh design.
% Inputs are actual paired samples and immutable axes in the SAME units.
% trusted is a caller eligibility assertion, not native signature validation.
% No checkpoint lineage, gauge reference, physical epoch or history is made.
if nargin<4,controls=struct();end
c=resolve_controls(controls);validate_inputs(snapshot,referenceAxes,actualAnchor);
x=snapshot.x(:)';y=snapshot.y(:);snapshot.x=x;snapshot.y=y;
snapshot.datasetIndex=1;snapshot.sourceLabel='caller supplied paired snapshot';
dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','general_paired_design', ...
    'snapshotCount',1,'snapshots',snapshot,'xLimits',x([1,end]),'yLimits',y([1,end])');
profiles=ipm_gridlab_feature_profiles(dataset);p=profiles.records;
limits=struct('maxAdjacentCellRatio',1.08,'maxLogSpacingCurvature',.01, ...
    'minStencilRcond',1e-9,'minQuadratureWeightRatio',1e-8);
core=[max(0,p.safetyCoreCenter-p.safetyCoreWidthLeft),min(x(end),p.safetyCoreCenter+p.safetyCoreWidthRight)];
front=[max(0,p.leftFrontCenter-p.leftFrontHalfWidth),min(x(end),p.leftFrontCenter+p.leftFrontHalfWidth)];
spacing=min(diff(core)/c.targetXCoreCells,diff(front)/c.minimumXFrontCells)/c.resolutionPadding;
report=struct('kind','bounded_general_paired_mesh_design_v1','controls',c,'meshLimits',limits, ...
    'snapshotProvenance',snapshot.provenance,'referenceProvenance',referenceAxes.provenance, ...
    'actualAnchor',actualAnchor,'nodeCount',[numel(x),numel(y)],'xLimits',x([1,end]),'yLimits',y([1,end])', ...
    'inputAnchorExact',any(x==actualAnchor)&&any(x==-actualAnchor), ...
    'sourceCanonicalTime',snapshot.canonicalTime,'sourcePhysicalTime',snapshot.physicalTime, ...
    'feature',struct('coreCenter',p.safetyCoreCenter,'coreInterval',core,'frontInterval',front, ...
    'yCoreWidth',p.yCoreWidth,'anchorGapOverCoreWidth',abs(actualAnchor-p.safetyCoreCenter)/diff(core)), ...
    'xTrials',struct([]),'yTrials',struct([]),'pairTrials',struct([]),'status','registered', ...
    'nativeSignatureValidated',false,'noLU',true,'noPDE',true,'boxChanged',false,'nodeCountChanged',false, ...
    'historyOrRuntimeReferencesCreated',false,'globalFeasibilityClaim',false);
candidate=struct('transactionReady',false,'status','not_selected','candidateX',[],'candidateY',[]);
xRows=struct([]);xAxes={};yRows=struct([]);yAxes={};
% Enumerate all registered y sigmas, including their limiting-gate reports.
for sigma=c.yEqualizationSigmas
    [v,info]=ipm_gridlab_equalize_axis(referenceAxes.y(:)/actualAnchor, ...
        [0,p.yCoreWidth/actualAnchor],c.targetYCoreCells*c.verticalResolutionPadding, ...
        struct('geometry','positive','focusCenters',0,'sigma',sigma,'power',2,'meshLimits',limits));
    v=v(:)*actualAnchor;v(1)=y(1);v(end)=y(end);
    [q,reasons,margin]=quality(v,limits);count=interval_count(v,[0,min(v(end),p.yCoreWidth)]);
    reasons=with_reasons(reasons,count<c.targetYCoreCells-1e-6,{'y_core_target'});
    row=struct('sigma',sigma,'quality',q,'coreCells',count,'admissible',isempty(reasons), ...
        'reasons',{reasons},'qualityMargin',margin,'equalizerInfo',info);
    yRows=append_row(yRows,row);yAxes{end+1}=v; %#ok<AGROW>
end
% A numerical construction rejection is separate from an input/program error.
for fine=c.positiveFineCells
 for rounding=c.roundingCells
  for fraction=c.coreFineCellFractions
    row=struct('positiveFineCells',fine,'roundingCells',rounding,'coreFineCellFraction',fraction, ...
        'admissible',false,'reasons',{{}},'quality',struct(),'coreCells',NaN,'frontCells',NaN, ...
        'qualityMargin',-Inf,'construction',struct(),'failure',struct());v=[];
    if fine+2*(rounding+2)>(numel(x)-1)/2
        row.reasons={'rounded_branch_node_budget'};
    elseif p.safetyCoreCenter<=0 || core(2)<=core(1) || front(2)<=front(1)
        row.reasons={'unsupported_or_degenerate_core_geometry'};
    else
        try
            [v,construction]=mesh_core_patch_axis(actualAnchor,p.safetyCoreCenter,x(end),numel(x), ...
                spacing,fine,fraction,rounding,c.anchorWarpRadius);
            [q,reasons,margin]=quality(v,limits);
            nx=interval_count(v(v>=0),core);nf=interval_count(v(v>=0),front);
            reasons=with_reasons(reasons,[nx<c.targetXCoreCells-1e-6,nf<c.minimumXFrontCells-1e-6, ...
                ~any(v==actualAnchor)||~any(v==-actualAnchor)],{'x_core_target','x_front_target','exact_anchor'});
            row.quality=q;row.coreCells=nx;row.frontCells=nf;row.qualityMargin=margin;
            row.construction=construction;row.reasons=reasons;row.admissible=isempty(reasons);
        catch exception
            assert_geometry_rejection(exception);
            row.reasons={'geometry_construction_rejected'};
            row.failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
        end
    end
    xRows=append_row(xRows,row);xAxes{end+1}=v; %#ok<AGROW>
  end
 end
end
report.xTrials=xRows;report.yTrials=yRows;
% Rank by the weakest dimensionless quality margin, then enumeration order.
pairIndex=[];ranks=[];
if c.includeUnchangedPair
    [~,rx]=quality(x,limits);[~,ry]=quality(y,limits);
    if any(x==actualAnchor)&&any(x==-actualAnchor)&&isempty(rx)&&isempty(ry)&&interval_count(x(x>=0),core)>=c.targetXCoreCells-1e-6 && ...
            interval_count(x(x>=0),front)>=c.minimumXFrontCells-1e-6 && ...
            interval_count(y,[0,min(y(end),p.yCoreWidth)])>=c.targetYCoreCells-1e-6
        pairIndex=[0,0];ranks=Inf;
    end
end
for i=find([xRows.admissible])
 for j=find([yRows.admissible])
    pairIndex(end+1,:)=[i,j];ranks(end+1,1)=min(xRows(i).qualityMargin,yRows(j).qualityMargin); %#ok<AGROW>
 end
end
if isempty(pairIndex)
    if ~any([yRows.admissible]),report.status='registered_y_family_infeasible';
    else,report.status='registered_x_family_infeasible';end
    candidate.status=report.status;return
end
[~,order]=sort(ranks,'descend');report.admittedAxisPairCount=numel(order);
for k=1:min(c.maximumFullPairScores,numel(order))
    ij=pairIndex(order(k),:);unchanged=all(ij==0);
    if unchanged,vx=x;vy=y;else,vx=xAxes{ij(1)};vy=yAxes{ij(2)};end
    score=ipm_gridlab_score_frozen_pair(dataset,vx,vy,struct('trustedOnly',true,'profiles',profiles, ...
        'meshLimits',limits,'minimumXCoreCells',c.targetXCoreCells,'minimumYCoreCells',c.targetYCoreCells, ...
        'minimumXFrontCells',c.minimumXFrontCells));
    a=score.aggregate;ready=false;reasons=score.rejectionReasons;
    if score.admissible
        reasons=with_reasons(reasons,[a.minimumXCoreCells<c.targetXCoreCells-1e-6, ...
            a.minimumYCoreCells<c.targetYCoreCells-1e-6,a.minimumXLeftFrontCells<c.minimumXFrontCells-1e-6, ...
            a.worst.rhoXMaximumRelativeChange>.002,a.worst.conservationRelativeDefect>5e-12, ...
            a.worst.relativeRangeViolation>2e-4], ...
            {'paired_x_core','paired_y_core','paired_front','paired_peak_jump','paired_mass','paired_range'});
        ready=isempty(reasons);
    end
    item=struct('xTrial',ij(1),'yTrial',ij(2),'unchanged',unchanged,'passed',ready,'reasons',{reasons},'score',score);
    report.pairTrials=append_row(report.pairTrials,item);
    if ready
        candidate=struct('kind','general_paired_core_patch_candidate','transactionReady',true, ...
            'status','frozen_pair_passed','candidateX',vx,'candidateY',vy,'unchanged',unchanged, ...
            'anchorPosition',actualAnchor,'anchorExact',any(vx==actualAnchor)&&any(vx==-actualAnchor), ...
            'referenceAxes',referenceAxes,'snapshotProvenance',snapshot.provenance,'controls',c, ...
            'selectedXTrial',ij(1),'selectedYTrial',ij(2),'pairScore',score, ...
            'completePairScoredInOriginalUnits',true,'currentTransactionNodeCountSupported',true, ...
            'nativeTransactionPerformed',false,'pdeAdvanced',false);
        report.status='frozen_pair_passed';return
    end
end
report.status='registered_pair_budget_exhausted';candidate.status=report.status;
end

function c=resolve_controls(c)
d=struct('targetXCoreCells',32,'targetYCoreCells',32,'minimumXFrontCells',20, ...
    'positiveFineCells',[32,48,64,80,96,112,128],'roundingCells',[16,24,32,40,48], ...
    'coreFineCellFractions',[.5,.65],'yEqualizationSigmas',[.18,.25,.35,.5], ...
    'resolutionPadding',1.10,'verticalResolutionPadding',1.15,'anchorWarpRadius',.25, ...
    'maximumAxisCandidates',70,'maximumFullPairScores',3,'includeUnchangedPair',true);
assert(isstruct(c)&&isscalar(c)&&isempty(setdiff(fieldnames(c),fieldnames(d))), ...
    'ipm:PairedDesignControls','Unknown design controls; hard mesh gates cannot be overridden.');
for n=fieldnames(d)',if ~isfield(c,n{1}),c.(n{1})=d.(n{1});end,end
for n={'targetXCoreCells','targetYCoreCells','minimumXFrontCells','maximumAxisCandidates','maximumFullPairScores'}
    validateattributes(c.(n{1}),{'numeric'},{'scalar','integer','positive'});
end
assert(c.targetXCoreCells>=21&&c.targetYCoreCells>=21&&c.minimumXFrontCells>=20);
for n={'positiveFineCells','roundingCells'}
    validateattributes(c.(n{1}),{'numeric'},{'row','integer','finite','positive','nonempty'});
end
assert(all(c.positiveFineCells>=8)&&all(c.roundingCells>=2));
validateattributes(c.coreFineCellFractions,{'numeric'},{'row','finite','>',0,'<',1});
validateattributes(c.yEqualizationSigmas,{'numeric'},{'row','finite','positive'});
validateattributes(c.resolutionPadding,{'numeric'},{'scalar','finite','>=',1});
validateattributes(c.verticalResolutionPadding,{'numeric'},{'scalar','finite','>=',1});
validateattributes(c.anchorWarpRadius,{'numeric'},{'scalar','finite','>',0,'<',.9});
assert(islogical(c.includeUnchangedPair)&&isscalar(c.includeUnchangedPair));
assert(numel(c.positiveFineCells)*numel(c.roundingCells)*numel(c.coreFineCellFractions)<=c.maximumAxisCandidates && ...
    c.maximumAxisCandidates<=500 && c.maximumFullPairScores<=20 && numel(c.yEqualizationSigmas)<=16, ...
    'ipm:PairedDesignBudget','The entire registered search must fit its finite axis/pair budgets.');
end

function validate_inputs(s,r,a)
assert(isstruct(s)&&isscalar(s)&&all(isfield(s,{'rho','x','y','scale','canonicalTime','physicalTime','trusted','provenance'})), ...
    'ipm:PairedDesignInput','A complete caller-validated paired snapshot is required.');
assert(isstruct(r)&&isscalar(r)&&all(isfield(r,{'x','y','provenance'})), ...
    'ipm:PairedDesignInput','Immutable actual reference axes and provenance are required.');
validateattributes(a,{'numeric'},{'scalar','finite','positive'});
for ax={s.x,s.y,r.x,r.y},validateattributes(ax{1},{'numeric'},{'vector','finite','real','increasing'});end
x=s.x(:)';y=s.y(:);rx=r.x(:)';ry=r.y(:);
assert(isequal(size(s.rho),[numel(y),numel(x)])&&all(isfinite(s.rho),'all')&&isreal(s.rho)&& ...
    islogical(s.trusted)&&isscalar(s.trusted)&&s.trusted,'ipm:PairedDesignInput','Field, axes and eligibility must pair.');
assert(mod(numel(x),2)==1 && numel(x)>=17 && numel(y)>=8 && isequal(x,-fliplr(x)) && ...
    isequal(rx,-fliplr(rx)) && y(1)==0 && ry(1)==0 && ...
    numel(x)==numel(rx)&&numel(y)==numel(ry)&&isequal(x([1,end]),rx([1,end]))&& ...
    isequal(y([1,end]),ry([1,end]))&&a<x(end), ...
    'ipm:PairedDesignInput','Fixed N, fixed box, symmetry and an interior positive actual anchor are mandatory.');
assert(isstruct(s.scale)&&all(isfield(s.scale,{'Cx','Cy','Comega'})));
validateattributes([s.scale.Cx,s.scale.Cy,s.scale.Comega],{'numeric'},{'finite','positive','numel',3});
validateattributes([s.canonicalTime,s.physicalTime],{'numeric'},{'finite','nonnegative','numel',2});
end

function [q,reasons,margin]=quality(axis,limits)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
values=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio,q.maximumQuadratureWeightToControlWidthRatio];
bad=[values(1)>limits.maxAdjacentCellRatio,values(2)>limits.maxLogSpacingCurvature, ...
    values(3)<limits.minStencilRcond,values(4)<limits.minQuadratureWeightRatio,values(5)<.35,values(6)>1.65];
names={'adjacent_ratio','log_spacing_curvature','stencil_rcond','global_quadrature','local_quadrature_min','local_quadrature_max'};
reasons=names(bad);if ~q.quadratureWeightsStrictlyPositive,reasons{end+1}='nonpositive_quadrature';end
margin=min([limits.maxAdjacentCellRatio/values(1),limits.maxLogSpacingCurvature/max(values(2),realmin), ...
    values(3)/limits.minStencilRcond,values(4)/limits.minQuadratureWeightRatio,values(5)/.35,1.65/values(6)]);
end

function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end

function reasons=with_reasons(reasons,failed,names)
reasons=[reasons,names(failed)];
end

function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end

function assert_geometry_rejection(e)
% Only known geometric infeasibility may become a rejected candidate.
known={'No feasible rounded-log cell split.','Insufficient rounded-branch cell budget.'};
fromFactory=~isempty(e.stack)&&any(strcmp(e.stack(1).name, ...
    {'mesh_core_patch_axis','mesh_general_rounded_axis','mesh_general_rounded_axis>growth'}));
if ~(fromFactory&&strcmp(e.identifier,'MATLAB:assertion:failed')) && ~any(strcmp(e.message,known))
    rethrow(e)
end
end
