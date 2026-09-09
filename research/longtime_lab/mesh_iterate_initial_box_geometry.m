function report=mesh_iterate_initial_box_geometry(previousScreenRun,outDir)
%MESH_ITERATE_INITIAL_BOX_GEOMETRY Bounded t=0 analytic observation iteration.
% Rejected candidate axes are observations only. The original uniform
% reference and physical datum remain fixed; no native history is created.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
issues=checkcode(which(mfilename),'-id');assert(isempty(issues),jsonencode(issues));
policy=ipm.config.autonomousMeshPolicy(struct('version',2, ...
    'nodeFamily',struct('maximumTotalNodes',110000)));
priorRegistration=jsondecode(fileread(fullfile(previousScreenRun,'registration.json')));
boxes=priorRegistration.halfWidthAndYmax;
validateattributes(boxes,{'double'},{'real','finite','positive','2d','ncols',2,'nonempty'});
assert(size(boxes,1)<=8 && all(boxes(:,1)>1) && size(unique(boxes,'rows'),1)==size(boxes,1));
registration=struct('kind','bounded_t0_analytic_geometry_iteration_v1', ...
    'previousScreenRun',previousScreenRun,'maximumRounds',3,'maximumProposalsPerRound',3, ...
    'maximumAnalyticCandidateMeasurementsPerBox',9,'halfWidthAndYmax',boxes, ...
    'policy',policy,'referenceRule','original configuration uniform axes fixed across rounds', ...
    'datumRule','same native analytic primitive k8 at physical and canonical time zero', ...
    'retryRule','all candidate rejection reasons must be actual core or front', ...
    'observationSelection','maximum min(actualCore/31,actualFront/20); original proposal order breaks ties', ...
    'stopRule','first candidate passing all pure analytic geometry and full family gates', ...
    'physicalTime',0,'canonicalTime',0,'pdeSteps',0,'noLU',true,'noFlow',true,'noMeshBuild',true, ...
    'noSourceInterpolation',true,'nativeInitializationQualification',false,'dynamicBoxErrorQualification',false);
write_json(fullfile(outDir,'registration.json'),registration);
report=registration;report.boxes=struct([]);profile clear;profile on;
for boxIndex=1:size(boxes,1)
    previous=load(fullfile(previousScreenRun,sprintf('box_%d_source_and_design.mat',boxIndex)));
    config=previous.config;assert(isequaln(config.remesh.autonomousMesh,policy) && ...
        isequal(config.grid.xlim,[-boxes(boxIndex,1),boxes(boxIndex,1)]) && config.grid.ymax==boxes(boxIndex,2));
    reference=struct('x',previous.source.x,'y',previous.source.y);
    assert(isequal([numel(reference.x),numel(reference.y)],[321,161])&& ...
        isequal(reference.x,linspace(config.grid.xlim(1),config.grid.xlim(2),321))&& ...
        isequal(reference.y,linspace(0,config.grid.ymax,161)')&& ...
        strcmp(config.physics.initialCondition,'degenerate_primitive')&&config.physics.degeneratePower==8);
    [source,sourceMeasure]=analytic_view(reference.x,reference.y,config);
    assert(isequaln(source,previous.source)&&isequaln(sourceMeasure,previous.sourceMeasure));
    box=struct('boxIndex',boxIndex,'xlim',config.grid.xlim,'ymax',config.grid.ymax, ...
        'rounds',struct([]),'acceptedRound',0,'acceptedCandidate',0,'status','pending', ...
        'firstRoundAxesBitwise',false,'pureGeometryQualified',false,'nativeInitializationQualification',false);
    sourceFromRound=0;sourceFromCandidate=0;
    for roundIndex=1:registration.maximumRounds
        [candidates,axisReport]=ipm.remesh.plannedAxisPairs(source,reference,1,policy);
        assert(isequal(reference.x,previous.source.x)&&isequal(reference.y,previous.source.y));
        if roundIndex==1
            box.firstRoundAxesBitwise=isequaln(candidates,previous.candidates);
            assert(box.firstRoundAxesBitwise,'ipm:InitialIterationRoundOne','The initial proposal family must remain bitwise identical.');
        end
        round=struct('roundIndex',roundIndex,'physicalTime',0,'canonicalTime',0,'acceptedStep',0, ...
            'analyticObservationOnly',true,'sourceFromRound',sourceFromRound, ...
            'sourceFromCandidate',sourceFromCandidate,'sourceMeasure',sourceMeasure, ...
            'sourcePadding',padding_record(axisReport.feature),'plannerStatus',axisReport.status, ...
            'xAdmitted',nnz([axisReport.xTrials.admissible]),'yAdmitted',nnz([axisReport.yTrials.admissible]), ...
            'xRejections',rejection_counts(axisReport.xTrials),'yRejections',rejection_counts(axisReport.yTrials), ...
            'proposalCount',numel(candidates),'candidateAudits',struct([]), ...
            'nextObservationCandidate',0,'status','evaluating');
        candidateViews=cell(1,numel(candidates));candidateMeasures=cell(1,numel(candidates));
        for candidateIndex=1:numel(candidates)
            candidate=candidates(candidateIndex);[view,measurement]=analytic_view(candidate.x,candidate.y,config);
            if roundIndex==1
                assert(isequaln(measurement,previous.row.candidateAudits(candidateIndex).analyticMeasurement));
            end
            feature=ipm.diagnostics.meshFeatureIntervals(view);
            [qx,rx]=quality(candidate.x,policy.qualityLimits,'x');[qy,ry]=quality(candidate.y,policy.qualityLimits,'y');
            qualityPassed=isempty(rx)&&isempty(ry);anchorExact=any(candidate.x==1)&&any(candidate.x==-1);
            reasons=[rx,ry];
            if any(measurement.actualCoreCells<policy.transactionMinimumCoreCells),reasons{end+1}='actual_analytic_core_floor';end %#ok<AGROW>
            if measurement.leftFrontCells<policy.minimumFrontCells,reasons{end+1}='actual_analytic_front_floor';end %#ok<AGROW>
            if ~anchorExact,reasons{end+1}='exact_anchor';end %#ok<AGROW>
            family=struct();familyPassed=false;failure=struct();
            try
                family=ipm.remesh.referenceAxisFamily(candidate.x,candidate.y,1,policy);
                familyPassed=all([family.members.qualityPassed])&&family.members(1).resourceAdmitted;
            catch e
                failure=struct('identifier',e.identifier,'message',e.message);
                reasons{end+1}='reference_family_rejected'; %#ok<AGROW>
            end
            passed=isempty(reasons)&&familyPassed;
            score=min([measurement.actualCoreCells./policy.transactionMinimumCoreCells, ...
                measurement.leftFrontCells/policy.minimumFrontCells]);
            audit=struct('candidateIndex',candidateIndex,'physicalTime',0,'canonicalTime',0,'acceptedStep',0, ...
                'analyticObservationOnly',true,'xIndex',candidate.xIndex,'yIndex',candidate.yIndex, ...
                'predictedCellsFromCurrentSource',candidate.predictedCells,'analyticMeasurement',measurement, ...
                'padding',padding_record(feature),'actualMinusPredictedCells', ...
                [measurement.actualCoreCells,measurement.leftFrontCells]-candidate.predictedCells, ...
                'qualityPassed',qualityPassed,'exactAnchor',anchorExact,'xQuality',qx,'yQuality',qy, ...
                'referenceFamilyPassed',familyPassed,'referenceFamilyFailure',failure, ...
                'observationSelectionScore',score,'passed',passed,'reasons',{reasons});
            round.candidateAudits=append_row(round.candidateAudits,audit);
            candidateViews{candidateIndex}=view;candidateMeasures{candidateIndex}=measurement;
            save(fullfile(outDir,sprintf('box_%d_round_%d_candidate_%d.mat',boxIndex,roundIndex,candidateIndex)), ...
                'candidate','view','measurement','feature','audit','family','config','reference','-v7.3');
            fprintf('INITIAL_ITERATION_CANDIDATE H=%g round=%d candidate=%d actual=%.8g/%.8g front=%.8g score=%.8g purePassed=%d\n', ...
                config.grid.xlim(2),roundIndex,candidateIndex,measurement.actualCoreCells,measurement.leftFrontCells,score,passed);
            if passed
                box.acceptedRound=roundIndex;box.acceptedCandidate=candidateIndex;box.pureGeometryQualified=true;
                box.status='pure_geometry_qualified_requires_native_initialization';round.status='accepted';break
            end
        end
        if ~box.pureGeometryQualified
            if isempty(candidates)
                round.status='planner_capacity';box.status='planner_capacity';
            elseif roundIndex==registration.maximumRounds
                round.status='registered_round_budget_exhausted';box.status='registered_round_budget_exhausted';
            else
                allowed={'actual_analytic_core_floor','actual_analytic_front_floor'};
                eligible=arrayfun(@(a)~a.passed&&a.qualityPassed&&a.exactAnchor&&a.referenceFamilyPassed&& ...
                    ~isempty(a.reasons)&&all(ismember(a.reasons,allowed)),round.candidateAudits);
                if all(eligible)
                    [~,next]=max([round.candidateAudits.observationSelectionScore]);
                    round.nextObservationCandidate=next;round.status='resolution_only_retry';
                else
                    round.status='nonresolution_rejection';box.status='nonresolution_rejection';
                end
            end
        end
        save(fullfile(outDir,sprintf('box_%d_round_%d_source_and_design.mat',boxIndex,roundIndex)), ...
            'source','sourceMeasure','reference','config','candidates','axisReport','round','-v7.3');
        box.rounds=append_row(box.rounds,round);
        write_json(fullfile(outDir,sprintf('box_%d_round_%d.json',boxIndex,roundIndex)),round);
        if box.pureGeometryQualified||round.nextObservationCandidate==0,break;end
        source=candidateViews{round.nextObservationCandidate};sourceMeasure=candidateMeasures{round.nextObservationCandidate};
        sourceFromRound=roundIndex;sourceFromCandidate=round.nextObservationCandidate;
    end
    report.boxes=append_row(report.boxes,box);write_json(fullfile(outDir,sprintf('box_%d.json',boxIndex)),box);
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})),'ipm:InitialIterationUnexpectedSolver','This is analytic observation only.');
end
report.protocolCompleted=true;report.allBoxesPureGeometryQualified=all([report.boxes.pureGeometryQualified]);
save(fullfile(outDir,'report.mat'),'report','profileInfo','-v7.3');write_json(fullfile(outDir,'report.json'),report);
fprintf('INITIAL_ITERATION_COMPLETE allPureQualified=%d nativeQualification=0 noLU=1\n',report.allBoxesPureGeometryQualified);
end

function record=padding_record(f)
raw=f.rawCoreBounds;center=f.coreCenter;
record=struct('coreCenter',center,'rawCoreBounds',raw,'paddedCoreInterval',f.coreInterval, ...
    'rawCoreWidth',diff(raw),'paddedCoreWidth',diff(f.coreInterval), ...
    'widthAmplification',diff(f.coreInterval)/diff(raw), ...
    'leftRawHalfWidth',center-raw(1),'leftPaddedHalfWidth',center-f.coreInterval(1), ...
    'rightRawHalfWidth',raw(2)-center,'rightPaddedHalfWidth',f.coreInterval(2)-center);
end

function [view,measure]=analytic_view(x,y,config)
[X,Y]=meshgrid(x,y);ops=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y), ...
    'symmetryMode',config.physics.symmetryMode);
rho=ipm.field.initialDensity(ops,config.physics);Dx=ipm.mesh.fdMatrix(x,1,7);source=rho*Dx';
assert(all(isfinite(rho),'all')&&isreal(rho));
view=struct('rho',rho,'x',x,'y',y,'Dx',Dx,'source',source,'trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);k=config.physics.degeneratePower;
exact=sign(X).*abs(X).^(k-1)./(1+abs(X).^k+Y.^k);
window=X>=.5&X<=1.8&Y<=1;
measure=struct('actualCoreCells',feature.actualCoreCells,'leftFrontCells',feature.leftFrontCells, ...
    'coreInterval',feature.coreInterval,'frontInterval',feature.frontInterval, ...
    'yCoreWidth',feature.yCoreWidth,'coreCenter',feature.coreCenter, ...
    'wallNodalPeak',max(source(1,:)),'exactContinuousWallPeakX',(k-1)^(1/k), ...
    'exactContinuousWallPeak',(k-1)^((k-1)/k)/k, ...
    'pairedDerivativeInnerRelativeInf',max(abs(source(window)-exact(window)))/max(abs(exact(window))), ...
    'minimumDx',min(diff(x)),'minimumDy',min(diff(y)), ...
    'trustedViewMeansFinitePairedAnalyticDatumOnly',true);
end
function [q,reasons]=quality(axis,limits,label)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
bad=[q.maximumAdjacentCellRatio>limits.maxAdjacentCellRatio,q.maximumLogSpacingCurvature>limits.maxLogSpacingCurvature, ...
    q.minimumStencilRcond<limits.minStencilRcond,q.minimumQuadratureWeightRatio<limits.minQuadratureWeightRatio, ...
    q.minimumQuadratureWeightToControlWidthRatio<limits.minWeightToControlWidth, ...
    q.maximumQuadratureWeightToControlWidthRatio>limits.maxWeightToControlWidth,~q.quadratureWeightsStrictlyPositive];
names={'adjacent_ratio','spacing_curvature','stencil_rcond','global_quadrature','local_weight_min','local_weight_max','positive_weights'};
reasons=cellfun(@(s)[label,':',s],names(bad),'UniformOutput',false);
end
function rows=rejection_counts(trials)
labels={};for k=1:numel(trials),labels=union(labels,trials(k).reasons,'stable');end
rows=struct([]);
for k=1:numel(labels)
    rows=append_row(rows,struct('reason',labels{k},'count',nnz(arrayfun(@(v)any(strcmp(v.reasons,labels{k})),trials))));
end
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
