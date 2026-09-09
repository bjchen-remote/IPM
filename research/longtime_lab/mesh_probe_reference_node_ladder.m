function report=mesh_probe_reference_node_ladder(checkpointFile,resultFile,outDir)
%MESH_PROBE_REFERENCE_NODE_LADDER Strict source + four research-only N levels.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
for name={'mesh_plan_reference_node_count',mfilename}
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;
r=ipm.output.validate(resultFile);policy=s.config.remesh.autonomousMesh;
pairing=isequaln(s.rho,r.state.rho)&&isequal(s.x,r.grid.x)&&isequal(s.y,r.grid.y)&& ...
    s.step==r.state.steps&&s.scale.canonicalTime==r.state.canonicalTime&& ...
    s.scale.physicalTime==r.state.physicalTime&&s.normalizedTime==r.state.normalizedTime&& ...
    isequaln(cp.payload.log.history,r.history)&&isequaln(s.config,r.config)&& ...
    strcmp(s.runMetadata.caseId,r.metadata.caseId)&&s.remeshCount==r.grid.remeshCount&& ...
    r.scale.Cx==exp(s.scale.logC_l)&&r.scale.Cy==exp(s.scale.logC_l)&& ...
    r.scale.Comega==exp(s.scale.logC_omega)&&r.scale.Xshift==s.scale.X_shift;
assert(pairing&&all(ipm.output.trustedMask(r)),'ipm:ResearchNativePair','Source result/native history must match exactly and be wholly trusted.');
assert(isequal([numel(s.x),numel(s.y)],[321,161]));
anchor=s.config.scaling.transportAnchorX;
reference=struct('x',s.baseX,'y',s.baseY);
assert(isequal(reference.x,s.runMetadata.autonomousMesh.initialization.selectedBaseX)&& ...
    isequal(reference.y,s.runMetadata.autonomousMesh.initialization.selectedBaseY));
Dx=ipm.mesh.fdMatrix(s.x,1,7);
view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',s.rho*Dx','trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);
report=struct('kind','research_immutable_reference_node_ladder_v2','sourceCheckpoint',checkpointFile, ...
    'sourceResult',resultFile,'sourcePairingExact',pairing,'sourceStep',s.step, ...
    'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime, ...
    'originalNativeReferenceExact',true,'policy',policy,'nodeLadder',[321,161;321,321;641,161;641,321], ...
    'feature',feature,'noLU',true,'noPDE',true,'nativeTransactionPerformed',false, ...
    'nativeVariableNodeCountSupported',false,'sourceFieldResampledForGeometry',false, ...
    'sourceStopReason',r.state.stopReason,'savedLastFailure',s.runMetadata.autonomousMesh.lastFailure, ...
    'levels',struct([]),'selectedResearchLevel',0,'selectedResearchCandidate',0);
write_json(fullfile(outDir,'registration.json'),report);
[original,originalReport]=ipm.remesh.plannedAxisPairs(view,reference,anchor,policy);
[same,sameReport]=mesh_plan_reference_node_count(feature,reference_axes(view),reference,anchor,policy);
report.sameNodePlannerBitwise=isequaln(original,same)&&isequaln(originalReport.xTrials,sameReport.xTrials)&& ...
    isequaln(originalReport.yTrials,sameReport.yTrials);
assert(report.sameNodePlannerBitwise,'ipm:ResearchPlannerEquivalence','Same-N wrapper must reproduce every original candidate/axis score exactly.');
snapshot=struct('rho',s.rho,'x',s.x,'y',s.y,'scale', ...
    struct('Cx',r.scale.Cx,'Cy',r.scale.Cy,'Comega',r.scale.Comega), ...
    'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime,'trusted',true, ...
    'datasetIndex',1,'sourceLabel','strict native source before node refinement');
dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','actual_321x161_source', ...
    'snapshotCount',1,'snapshots',snapshot,'xLimits',s.x([1,end]),'yLimits',s.y([1,end])');
profiles=ipm_gridlab_feature_profiles(dataset);
for level=1:size(report.nodeLadder,1)
    if level==1
        targetReference=reference;axes=original;axisReport=originalReport;
        referenceAudit=struct('oldNodesPreservedExactly',true,'exactAnchor',true,'sameEndpoints',true,'method','identity');
    else
        targetReference=refine_reference(reference,anchor,report.nodeLadder(level,:));
        strides=(report.nodeLadder(level,:)-1)./([numel(reference.x),numel(reference.y)]-1);
        referenceAudit=struct('oldNodesPreservedExactly',isequal(targetReference.x(1:strides(1):end),reference.x)&& ...
            isequal(targetReference.y(1:strides(2):end),reference.y),'exactAnchor',any(targetReference.x==anchor)&&any(targetReference.x==-anchor), ...
            'sameEndpoints',isequal(targetReference.x([1,end]),reference.x([1,end]))&&isequal(targetReference.y([1,end]),reference.y([1,end])), ...
            'method','normalized_index_pchip_from_original_selected_base; old knots explicitly retained; positive x mirrored');
        assert(referenceAudit.oldNodesPreservedExactly&&referenceAudit.exactAnchor&&referenceAudit.sameEndpoints);
        [axes,axisReport]=mesh_plan_reference_node_count(feature,reference_axes(view),targetReference,anchor,policy);
    end
    current=struct('nodeCount',[numel(targetReference.x),numel(targetReference.y)], ...
        'referenceAudit',referenceAudit,'axisStatus',axisReport.status, ...
        'admittedX',nnz([axisReport.xTrials.admissible]),'admittedY',nnz([axisReport.yTrials.admissible]), ...
        'xRejectionCounts',rejection_counts(axisReport.xTrials),'yRejectionCounts',rejection_counts(axisReport.yTrials), ...
        'candidateCount',numel(axes),'pairTrials',struct([]),'researchTransferPassed',false);
    save(fullfile(outDir,sprintf('level_%d_axes.mat',level)),'axes','axisReport','targetReference','referenceAudit','-v7.3');
    fprintf('NODE_LADDER_AXES level=%d N=%dx%d x=%d y=%d candidates=%d\n',level,current.nodeCount,current.admittedX,current.admittedY,numel(axes));
    for k=1:numel(axes)
        candidate=axes(k);meshLimits=rmfield(policy.qualityLimits,{'minWeightToControlWidth','maxWeightToControlWidth'});
        score=ipm_gridlab_score_frozen_pair(dataset,candidate.x,candidate.y, ...
            struct('trustedOnly',true,'profiles',profiles,'meshLimits',meshLimits, ...
            'minimumXCoreCells',32,'minimumYCoreCells',32,'minimumXFrontCells',20));
        rhoAfterX=ipm.remesh.interpolate(s.x,s.rho,candidate.x, ...
            struct('sampleDimension',2,'conservation','constant','stencilWidth',6));
        rhoNew=ipm.remesh.interpolate(s.y,rhoAfterX,candidate.y, ...
            struct('sampleDimension',1,'conservation','constant','stencilWidth',6));
        newDx=ipm.mesh.fdMatrix(candidate.x,1,7);
        transferredView=struct('rho',rhoNew,'x',candidate.x,'y',candidate.y,'Dx',newDx,'source',rhoNew*newDx','trusted',true);
        transferredFeature=ipm.diagnostics.meshFeatureIntervals(transferredView);
        transferPassed=false;peak=NaN;mass=NaN;range=NaN;
        if score.admissible
            a=score.aggregate;peak=a.worst.rhoXMaximumRelativeChange;mass=a.worst.conservationRelativeDefect;range=a.worst.relativeRangeViolation;
            transferPassed=all(transferredFeature.actualCoreCells>=policy.transactionMinimumCoreCells)&& ...
                transferredFeature.leftFrontCells>=policy.minimumFrontCells&&peak<=policy.maximumSinglePeakJump&& ...
                mass<=policy.maximumMassRelativeDefect&&range<=policy.maximumRelativeRangeViolation&& ...
                s.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump+peak<=policy.maximumCumulativeAbsolutePeakJump;
        end
        observation=struct('index',k,'frozenPairAdmissible',score.admissible,'scoreReasons',{score.rejectionReasons}, ...
            'actualTransferredCoreCells',transferredFeature.actualCoreCells,'actualTransferredFrontCells',transferredFeature.leftFrontCells, ...
            'relativePeakJump',peak,'massRelativeDefect',mass,'relativeRangeViolation',range, ...
            'researchTransferPassed',transferPassed,'nativeTransactionSupported',false,'noLU',true,'noPDE',true);
        current.pairTrials=append_row(current.pairTrials,observation);
        candidate.referenceFamily=targetReference;candidate.sourceCheckpoint=checkpointFile;
        candidate.nodeCountChanged=level>1;candidate.nativeTransactionSupported=false;
        candidate.nativeTransactionPerformed=false;candidate.researchTransferPassed=transferPassed;
        save(fullfile(outDir,sprintf('level_%d_candidate_%d.mat',level,k)), ...
            'candidate','score','observation','transferredFeature','rhoNew','snapshot','-v7.3');
        fprintf('NODE_LADDER_PAIR level=%d k=%d actualCore=%g/%g front=%g peak=%g pass=%d\n',level,k,transferredFeature.actualCoreCells,transferredFeature.leftFrontCells,peak,transferPassed);
        if transferPassed
            current.researchTransferPassed=true;
            if level>1&&report.selectedResearchLevel==0
                report.selectedResearchLevel=level;report.selectedResearchCandidate=k;
            end
            break;
        end
    end
    report.levels=append_row(report.levels,current);
    write_json(fullfile(outDir,sprintf('level_%d.json',level)),current);
end
report.protocolCompleted=true;report.nodeLadderResearchPassed=report.selectedResearchLevel>0;
report.minimumQualifiedNodeCount=[];
if report.nodeLadderResearchPassed,report.minimumQualifiedNodeCount=report.levels(report.selectedResearchLevel).nodeCount;end
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report','reference','feature','snapshot','-v7.3');
fprintf('NODE_LADDER_COMPLETE ladderResearchPassed=%d selectedLevel=%d nativeTransaction=0 noLU=1\n',report.nodeLadderResearchPassed,report.selectedResearchLevel);
end
function axes=reference_axes(view)
axes=struct('x',view.x,'y',view.y);
end
function refined=refine_reference(reference,anchor,targetCount)
x=reference.x;y=reference.y;
if targetCount(1)>numel(x)
    p=x(x>=0);u=linspace(0,1,numel(p));v=linspace(0,1,2*numel(p)-1);
    positive=pchip(u,p,v);positive(1:2:end)=p;positive(1)=0;
    x=[-fliplr(positive(2:end)),positive];
end
if targetCount(2)>numel(y)
    y=pchip(linspace(0,1,numel(y)),y,linspace(0,1,2*numel(y)-1))';
    y(1:2:end)=reference.y;y(1)=0;
end
assert(isequal([numel(x),numel(y)],targetCount)&&any(x==anchor)&&any(x==-anchor)&&all(diff(x)>0)&&all(diff(y)>0));
refined=struct('x',x,'y',y);
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct([]);
for k=1:numel(labels)
    count=nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials));
    rows=append_row(rows,struct('reason',labels{k},'count',count));
end
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
