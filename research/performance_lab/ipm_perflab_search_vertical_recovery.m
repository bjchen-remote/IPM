function report=ipm_perflab_search_vertical_recovery(xCandidateFile,outputDirectory)
%IPM_PERFLAB_SEARCH_VERTICAL_RECOVERY Root-y, then audited single-hop y.
% Fixed accepted-field x candidate; no evolution, LU, or shared mutation.
assert(~exist(outputDirectory,'dir'),'Use a new output directory.');mkdir(outputDirectory);
loaded=load(xCandidateFile,'candidate');xParent=loaded.candidate;p=xParent.provenance;
current=ipm.output.validate(p.currentFile);assert(all(ipm.output.trustedMask(current)));
loaded=load(p.baseDesignFile,'design');d=loaded.design;
loaded=load(p.rootReferenceFile,'dataset');root=loaded.dataset;
loaded=load(p.parentPlatformFile,'candidate');parent=loaded.candidate;
loaded=load(parent.referenceFile,'dataset');parentRoot=loaded.dataset;
assert(parent.rootAxisFactoryReproducedExactly&&parent.transactionReady&& ...
    strcmp(parent.kind,'ipm_longtime_moving_platform_candidate'));
assert(isequal(parent.candidateX(:),current.grid.x(:))&&isequal(parent.candidateY(:),current.grid.y(:)));
assert(isequal(root.snapshots(1).x,parentRoot.snapshots(1).x)&&isequal(root.snapshots(1).y,parentRoot.snapshots(1).y));
assert(isequal(d.dataset.snapshots(1).rho,current.state.rho)&& ...
    isequal(d.dataset.snapshots(1).x(:),current.grid.x(:))&&isequal(d.dataset.snapshots(1).y(:),current.grid.y(:)));
assert(isscalar(d.profiles.records)&&isscalar(d.dataset.snapshots));
assert(xParent.incrementalGeneration==1&&xParent.currentIncrementalReferenceDepth==0);
x=xParent.candidateX;assert(any(x==1)&&any(x==-1));
assert(xParent.pairScore.aggregate.minimumXCoreCells>=21-1e-6&& ...
    xParent.pairScore.aggregate.minimumXLeftFrontCells>=20-1e-6);
protocol=struct('sigmas',[.18,.25,.35,.5,.75,1,1.5,2,3],'targetY',21*1.15,'minimumY',21, ...
    'targetX',21,'frontMinimum',20,'localQuadratureLimits',[.35,1.65],'meshLimits',d.options.meshLimits, ...
    'sourceOrder',{{'immutable_external_root','trusted_current_y_single_hop'}}, ...
    'nonmonotoneSamples',33,'maximumAmplitudeCap',16,'maximumFullPairScores',3, ...
    'maximumPeakJump',2e-3,'maximumMassDefect',5e-12,'maximumRangeViolation',2e-4, ...
    'xCandidateFile',xCandidateFile,'xReferenceDepth',1,'entryYReferenceDepth',0, ...
    'currentYDepthEvidence','Actual entry y equals the last accepted globally rooted factory platform y exactly', ...
    'runtimeReferencesModifiedOrReplaced',false,'sameRoundRecursiveReference',false);
write_json(fullfile(outputDirectory,'frozen_protocol.json'),protocol);
references={root.snapshots(1).y(:)',current.grid.y(:)'};
labels=protocol.sourceOrder;interval=[0,d.profiles.records.yCoreWidth];
items=cell(36,1);axes=cell(36,1);infos=cell(36,1);itemCount=0;referenceAudit=cell(2,1);
for source=1:2
    yref=references{source};q=ipm.mesh.quality(yref,7,ipm.mesh.quadrature(yref));
    referenceAudit{source}=struct('source',labels{source},'quality',q,'admitted',admit(q,protocol.meshLimits), ...
        'sourceAxis',yref,'rootAndFactoryParentVerified',true, ...
        'entryReferenceDepth',0,'candidateIncrementalDepth',source-1);
    assert(referenceAudit{source}.admitted,'Reference y itself fails unchanged axis gates.');
    for sigma=protocol.sigmas
        o=struct('geometry','positive','focusCenters',0,'sigma',sigma,'power',2, ...
            'meshLimits',protocol.meshLimits,'maximumAmplitude',protocol.maximumAmplitudeCap);
        [y,info]=ipm_gridlab_equalize_axis(yref,interval,protocol.targetY,o);
        item=summarize(y,info,source,sigma,'maintained_target_search',interval,protocol);
        itemCount=itemCount+1;items{itemCount}=item;axes{itemCount}=y;infos{itemCount}=info;
        fprintf('Y_ORIGINAL source=%s sigma=%.2g count=%.9g status=%s qmin=%.9g ratio=%.9g\n', ...
            labels{source},sigma,item.cells,info.status,item.minimumQuadratureRatio,item.ratio);
    end
    write_json(fullfile(outputDirectory,'partial_y_axes.json'),items(1:itemCount));
    if any(cellfun(@(v)v.admitted&&v.cells>=protocol.targetY-1e-6,items(1:itemCount))),break;end
end
originalCount=itemCount;nonmonotoneUsed=false;
% Only when no maintained construction meets the actual floor, inspect the
% entire bounded path. Every A starts from its original root/current axis.
if ~any(cellfun(@(v)v.minimumPassed,items(1:itemCount)))
    nonmonotoneUsed=true;
    for k=1:originalCount
        original=items{k};source=original.sourceIndex;yref=references{source};sigma=original.sigma;
        upper=original.amplitude;trialRecords=cell(protocol.nonmonotoneSamples,1);trialIndex=0;best=original;bestY=axes{k};bestInfo=infos{k};
        for amplitude=linspace(0,upper,protocol.nonmonotoneSamples)
            o=struct('geometry','positive','focusCenters',0,'sigma',sigma,'power',2,'meshLimits',protocol.meshLimits);
            if amplitude==0,target=1e-6;else,o.initialAmplitude=amplitude;o.maximumAmplitude=amplitude;target=numel(yref);end
            [y,info]=ipm_gridlab_equalize_axis(yref,interval,target,o);
            assert(info.amplitude==amplitude);
            item=summarize(y,info,source,sigma,'nonmonotone_33_sample',interval,protocol);
            trialIndex=trialIndex+1;trialRecords{trialIndex}=item;
            if item.admitted&&(~best.admitted||item.cells>best.cells),best=item;bestY=y;bestInfo=info;end
        end
        itemCount=itemCount+1;items{itemCount}=best;axes{itemCount}=bestY;infos{itemCount}=bestInfo;
        file=fullfile(outputDirectory,sprintf('y_path_source%d_sigma%04d.mat',source,round(100*sigma)));
        save(file,'trialRecords','best','bestY','bestInfo','-v7.3');
        fprintf('Y_PATH source=%s sigma=%.2g best=%.9g atA=%.9g passes=%d\n',labels{source},sigma,best.cells,best.amplitude,best.minimumPassed);
    end
end
items=items(1:itemCount);axes=axes(1:itemCount);infos=infos(1:itemCount);
save(fullfile(outputDirectory,'axis_candidates.mat'),'items','axes','infos','referenceAudit','-v7.3');
write_json(fullfile(outputDirectory,'axis_candidates.json'),items);
eligible=find(cellfun(@(v)v.minimumPassed,items));
selected=[];
if ~isempty(eligible)
    counts=cellfun(@(v)v.cells,items(eligible));depths=cellfun(@(v)v.sourceIndex-1,items(eligible));
    % Prefer an immutable-root solution; within that source, prefer padding
    % margin, then the lower adjacent ratio. No field scores enter selection.
    ratios=cellfun(@(v)v.ratio,items(eligible));
    [~,order]=sortrows([depths(:),-counts(:),ratios(:)],[1,2,3]);
    selected=eligible(order(1:min(numel(order),protocol.maximumFullPairScores)));
end
pairItems=cell(numel(selected),1);
for rank=1:numel(selected)
    k=selected(rank);y=axes{k};item=items{k};
    o=d.options.pairScoreOptions;o.profiles=d.profiles;o.minimumXCoreCells=21;o.minimumYCoreCells=21;
    o.minimumXFrontCells=20;o.meshLimits=protocol.meshLimits;
    score=ipm_gridlab_score_frozen_pair(d.dataset,x,y,o);
    localPassed=local_quality(score.quality.x)&&local_quality(score.quality.y);
    resolutionPassed=false;transferPassed=false;
    if score.admissible&&~isempty(fieldnames(score.aggregate))
        a=score.aggregate;
        assert(abs(a.minimumYCoreCells-item.cells)<1e-10);
        resolutionPassed=a.minimumXCoreCells>=21-1e-6&&a.minimumYCoreCells>=21-1e-6&&a.minimumXLeftFrontCells>=20-1e-6;
        transferPassed=a.worst.rhoXMaximumRelativeChange<=2e-3&&a.worst.conservationRelativeDefect<=5e-12&&a.worst.relativeRangeViolation<=2e-4;
    end
    ready=score.admissible&&localPassed&&resolutionPassed&&transferPassed;
    yProvenance=struct('source',labels{item.sourceIndex},'referenceAxis',references{item.sourceIndex}, ...
        'parentPlatformFile',p.parentPlatformFile,'rootReferenceFile',p.rootReferenceFile, ...
        'sourceReferenceDepth',0,'candidateIncrementalDepth',item.sourceIndex-1, ...
        'sameRoundRecursiveReference',false,'parentAxesMatchActualEntryExactly',true);
    candidate=struct('kind','ipm_vertical_recovery_single_hop_candidate','protocol',protocol,'provenance',p, ...
        'xProvenance',struct('file',xCandidateFile,'candidateIncrementalDepth',1),'yProvenance',yProvenance, ...
        'candidateX',x,'candidateY',y,'anchorPosition',1,'currentIncrementalReferenceDepth',0,'incrementalGeneration',1, ...
        'yEqualizationInfo',infos{k},'yAxisSelection',item,'pairScore',score,'localQuadraturePassed',localPassed, ...
        'resolutionPassed',resolutionPassed,'transferPassed',transferPassed,'transactionReady',ready, ...
        'pdeAdvanced',false,'checkpointRegridCommitted',false);
    file=fullfile(outputDirectory,sprintf('candidate_rank%d_source%d_sigma%04d.mat',rank,item.sourceIndex,round(100*item.sigma)));
    save(file,'candidate','-v7.3');
    summary=struct('file',file,'ready',ready,'source',labels{item.sourceIndex},'sigma',item.sigma,'coreY',item.cells, ...
        'localPassed',localPassed,'resolutionPassed',resolutionPassed,'transferPassed',transferPassed);
    if score.admissible,summary.aggregate=score.aggregate;end
    pairItems{rank}=summary;
    fprintf('Y_PAIR rank=%d source=%s sigma=%.2g ready=%d y=%.9g\n',rank,labels{item.sourceIndex},item.sigma,ready,item.cells);
end
report=struct('kind','vertical_recovery_frozen_search','protocol',protocol,'provenance',p, ...
    'originalAxisCount',originalCount,'nonmonotoneUsed',nonmonotoneUsed,'axisItems',{items}, ...
    'selectedAxisIndices',selected,'pairItems',{pairItems},'readyCount',nnz(cellfun(@(v)v.ready,pairItems)), ...
    'pdeAdvanced',false,'poissonBuilt',false);
save(fullfile(outputDirectory,'vertical_report.mat'),'report','-v7.3');write_json(fullfile(outputDirectory,'vertical_report.json'),report);
end

function item=summarize(y,info,source,sigma,method,interval,p)
q=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
overlap=max(0,min(y(2:end),interval(2))-max(y(1:end-1),interval(1)));cells=sum(overlap./diff(y));
admitted=admit(q,p.meshLimits)&&all(diff(y)>0)&&y(1)==0&&all(isfinite(y));
item=struct('sourceIndex',source,'sigma',sigma,'method',method,'amplitude',info.amplitude, ...
    'status',info.status,'cells',cells,'admitted',admitted,'minimumPassed',admitted&&cells>=p.minimumY-1e-6, ...
    'paddingPassed',admitted&&cells>=p.targetY-1e-6,'ratio',q.maximumAdjacentCellRatio, ...
    'curvature',q.maximumLogSpacingCurvature,'minimumStencilRcond',q.minimumStencilRcond, ...
    'minimumQuadratureRatio',q.minimumQuadratureWeightRatio,'minimumSpacing',q.minimumSpacing, ...
    'localQuadratureRatios',[q.minimumQuadratureWeightToControlWidthRatio,q.maximumQuadratureWeightToControlWidthRatio]);
end
function value=admit(q,l)
value=q.maximumAdjacentCellRatio<=l.maxAdjacentCellRatio&&q.maximumLogSpacingCurvature<=l.maxLogSpacingCurvature&& ...
    q.minimumStencilRcond>=l.minStencilRcond&&q.minimumQuadratureWeightRatio>=l.minQuadratureWeightRatio&&local_quality(q);
end
function passed=local_quality(q)
passed=q.quadratureWeightsStrictlyPositive&&q.minimumQuadratureWeightToControlWidthRatio>=.35&&q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,PrettyPrint=true),'char');
end
