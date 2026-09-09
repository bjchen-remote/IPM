function report=ipm_perflab_fresh_stage5_geometry(planFile,checkpointFile,outDir)
%IPM_PERFLAB_FRESH_STAGE5_GEOMETRY Two preregistered, bounded no-LU protocols.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
registration=struct('kind','fresh_stage5_gap_then_protected_warp_geometry_v1', ...
    'originalAxisPlanFile',planFile,'sourceCheckpointFile',checkpointFile, ...
    'protocol1FineCells',[36,40,44,48],'protocol1RoundingCells',[8,10,12,14,16], ...
    'protocol1Fractions',[.5,.65],'protocol1TotalRows',40,'protocol1MaximumNewTrials',38, ...
    'protocol1WarpUnchanged',.25,'protocol2OnlyIfProtocol1AllRejected',true, ...
    'protocol2OriginalIndices',[1,2,3],'protocol2GapFractions',[.5,.75,.9], ...
    'protocol2MaximumTrials',9,'protocol2GapDefinition', ...
    'base nearest-anchor node minus rightmost current core/front endpoint, in anchor units', ...
    'coreFrontCountsAndFullQualityGatesUnchanged',true, ...
    'sourceBaseOrRuntimeReferencesChanged',false,'currentAxisFallbackUsed',false, ...
    'sameBoxSameNodeCount',true,'noLU',true,'noPDE',true,'noTransfer',true, ...
    'nativeOrFrozenFieldTransferQualified',false,'globalFeasibilityClaim',false);
write_json(fullfile(outDir,'registration.json'),registration);
issues=checkcode(which(mfilename),'-id');write_json(fullfile(outDir,'checkcode.json'),issues);
assert(isempty(issues),jsonencode(issues));
timer=tic;profile clear;profile on;
loaded=load(planFile,'axisReport');original=loaded.axisReport;clear loaded;
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;
assert(s.step==9205&&numel(s.x)==895&&numel(s.y)==386&& ...
    all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
Dx=ipm.mesh.fdMatrix(s.x,1,7);
v=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(v);assert(isequaln(feature,original.feature));
assert(isequal(original.xLimits,s.x([1,end]))&&isequal(original.yLimits,s.y([1,end])')&& ...
    isequal(original.nodeCount,[numel(s.x),numel(s.y)])&& ...
    original.anchor==s.rescaling.transportAnchorX);
p=original.policy;anchor=original.anchor;
assert(p.search.maximumWarpFraction==.25&&~any([original.xTrials.admissible])&& ...
    nnz([original.yTrials.admissible])==2);
spacing=min(diff(feature.coreInterval)/p.targetCoreCells(1), ...
    diff(feature.frontInterval)/p.minimumFrontCells)/p.search.xPadding;
rows=struct([]);axes=cell(1,40);reused=0;
for fine=registration.protocol1FineCells
 for rounding=registration.protocol1RoundingCells
  for fraction=registration.protocol1Fractions
    oldIndex=find([original.xTrials.positiveFineCells]==fine & ...
        [original.xTrials.roundingCells]==rounding & [original.xTrials.coreFineCellFraction]==fraction);
    if ~isempty(oldIndex)
        assert(isscalar(oldIndex));old=original.xTrials(oldIndex);reused=reused+1;
        row=make_row(fine,rounding,fraction,.25,old.quality,old.coreCells,old.frontCells, ...
            old.reasons,old.construction,oldIndex,false);
    else
        [axis,info]=ipm.remesh.corePatchAxis(anchor,feature.coreCenter,s.x(end),numel(s.x), ...
            spacing,fine,fraction,rounding,.25);
        row=score_axis(axis,info,fine,rounding,fraction,.25,p,feature,0);
        axes{numel(rows)+1}=axis;
    end
    rows=append(rows,row);
  end
 end
end
assert(numel(rows)==40&&reused==2&&nnz([rows.recomputed])==38);
write_json(fullfile(outDir,'protocol1_rows.json'),rows);
secondRows=struct([]);secondAxes={};secondExecuted=~any([rows.admissible]);
if secondExecuted
    protected=[min(feature.coreInterval(1),feature.frontInterval(1)), ...
        max(feature.coreInterval(2),feature.frontInterval(2))]/anchor;
    for index=registration.protocol2OriginalIndices
        old=original.xTrials(index);assert(old.qualityMargin>=1&&~old.admissible);
        center=old.construction.normalizedAnchorNodeBeforeWarp;gap=center-protected(2);assert(gap>0);
        for fraction=registration.protocol2GapFractions
            radius=fraction*gap;
            [axis,info]=ipm.remesh.corePatchAxis(anchor,feature.coreCenter,s.x(end),numel(s.x), ...
                spacing,old.positiveFineCells,old.coreFineCellFraction,old.roundingCells,radius);
            assert(info.normalizedAnchorNodeBeforeWarp==center&&center-radius>protected(2));
            row=score_axis(axis,info,old.positiveFineCells,old.roundingCells, ...
                old.coreFineCellFraction,radius,p,feature,index);
            oldShift=warp_shift(feature.coreCenter/anchor,center,.25,1-center)*anchor;
            newShift=warp_shift(feature.coreCenter/anchor,center,radius,1-center)*anchor;
            assert(newShift==0);
            row.protection=struct('protectedInterval',protected*anchor, ...
                'supportInterval',(center+[-radius,radius])*anchor, ...
                'gapFraction',fraction,'disjoint',true,'coreCenterOriginalShift',oldShift, ...
                'coreCenterOriginalShiftInCoreWidths',oldShift/diff(feature.coreInterval), ...
                'coreCenterNewShift',newShift,'originalCoreCells',old.coreCells, ...
                'originalFrontCells',old.frontCells);
            secondRows=append(secondRows,row);secondAxes{end+1}=axis; %#ok<AGROW>
        end
    end
    assert(numel(secondRows)==9);
    write_json(fullfile(outDir,'protocol2_rows.json'),secondRows);
end
report=registration;report.sourceNativeSignatureValidated=true;report.actualFeaturesExactlyPaired=true;
report.originalBaseX=s.baseX;report.originalBaseY=s.baseY;
report.originalRuntimeReferences=s.rescaling;report.sourceRemeshCount=s.remeshCount;
report.originalParentLineage=s.runMetadata.caseMetadata.latePhysicalBoxBranch;
report.originalXRejections=rejection_counts(original.xTrials);
report.originalYRejections=rejection_counts(original.yTrials);
report.protocol1=summary(rows);report.protocol1.reusedOriginalRows=reused;
report.protocol2Executed=secondExecuted;report.protocol2=summary(secondRows);
report.sourceStateUnchanged=isequaln(s,cp.payload.state);assert(report.sourceStateUnchanged);
report.feature=feature;report.wallSeconds=toc(timer);report.completed=true;
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.output.restoreCheckpoint','ipm.evolve.flow', ...
        'ipm.field.velocity','ipm.field.poisson','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})));
end
save(fullfile(outDir,'report.mat'),'report','rows','axes','secondRows','secondAxes','original','profileInfo','-v7.3');
write_json(fullfile(outDir,'report.json'),report);
fprintf('FRESH_STAGE5_GEOMETRY_COMPLETE gapPassed=%d protectedExecuted=%d protectedPassed=%d noLU=1\n', ...
    report.protocol1.admitted,secondExecuted,report.protocol2.admitted);
end
function row=score_axis(axis,info,fine,rounding,fraction,radius,p,f,oldIndex)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));l=p.qualityLimits;
values=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio, ...
    q.maximumQuadratureWeightToControlWidthRatio];assert(all(isfinite(values)));
bad=[values(1)>l.maxAdjacentCellRatio,values(2)>l.maxLogSpacingCurvature, ...
    values(3)<l.minStencilRcond,values(4)<l.minQuadratureWeightRatio, ...
    values(5)<l.minWeightToControlWidth,values(6)>l.maxWeightToControlWidth];
labels={'adjacent_ratio','log_spacing_curvature','stencil_rcond','global_quadrature','local_quadrature_min','local_quadrature_max'};
reasons=labels(bad);if ~q.quadratureWeightsStrictlyPositive,reasons{end+1}='nonpositive_quadrature';end
core=interval_count(axis(axis>=0),f.coreInterval);front=interval_count(axis(axis>=0),f.frontInterval);
if core<p.targetCoreCells(1)-1e-6,reasons{end+1}='x_core_target';end
if front<p.minimumFrontCells-1e-6,reasons{end+1}='x_front_target';end
if ~any(axis==info.actualTransportAnchor)||~any(axis==-info.actualTransportAnchor),reasons{end+1}='exact_anchor';end
row=make_row(fine,rounding,fraction,radius,q,core,front,reasons,info,oldIndex,true);
end
function row=make_row(fine,rounding,fraction,radius,q,core,front,reasons,info,index,recomputed)
row=struct('positiveFineCells',fine,'roundingCells',rounding,'coreFineCellFraction',fraction, ...
    'warpRadius',radius,'quality',q,'coreCells',core,'frontCells',front, ...
    'admissible',isempty(reasons),'reasons',{reasons},'construction',info, ...
    'originalIndex',index,'recomputed',recomputed,'protection',struct());
end
function v=warp_shift(x,center,radius,delta)
u=(x-center)/radius;if abs(u)<1,v=delta*(1-u^2)^3;else,v=0;end
end
function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end
function row=summary(rows)
row=struct('trialCount',numel(rows),'admitted',0,'rejections',struct([]),'passingIndices',[]);
if ~isempty(rows)
    row.admitted=nnz([rows.admissible]);row.rejections=rejection_counts(rows);
    row.passingIndices=find([rows.admissible]);
end
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct('reason',{},'count',{});
for k=1:numel(labels)
    rows(k)=struct('reason',labels{k},'count',nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials)));
end
end
function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
