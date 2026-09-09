function report=mesh_scan_fresh_resolution_family(checkpointFile,resultFile,initialDataFile,outputDirectory)
%MESH_SCAN_FRESH_RESOLUTION_FAMILY Finer same-count designs; no LU or PDE.
% Every original-unit candidate uses the immutable real fresh initial axes.
% Synthetic unit tests are separately labelled and never native candidates.
assert(maxNumCompThreads==10 && ~isfolder(outputDirectory));mkdir(outputDirectory);
sourceRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
dependencies={'ipm.solve','ipm_gridlab_collect_dataset','ipm_gridlab_feature_profiles', ...
    'ipm_gridlab_equalize_axis','ipm_gridlab_score_frozen_pair','mesh_design_fresh_core_patch', ...
    'mesh_fresh_core_patch_candidate','mesh_core_patch_axis','mesh_general_rounded_axis', ...
    'fresh_profile_campaign_review','fresh_profile_value_equal'};
resolved=cell(size(dependencies));
for k=1:numel(dependencies)
    resolved{k}=which(dependencies{k});assert(startsWith(resolved{k},[sourceRoot,filesep]), ...
        'ipm:FreshResolutionFreeze','A required dependency resolves outside the complete frozen source.');
end
sourceReview=fresh_profile_campaign_review(resultFile,initialDataFile);
assert(strcmp(checkpointFile,sourceReview.checkpointFile));
registration=struct('kind','same_fresh_epoch_finer_grid_frozen_screen_v1', ...
    'checkpointFile',checkpointFile,'resultFile',resultFile,'initialDataFile',initialDataFile, ...
    'targets',[42,56],'fineCells',[64,96,128],'fractions',[.65,.5], ...
    'immutableNodeCount',[895,386],'frontMinimum',20,'originalMeshGatesUnchanged',true, ...
    'unitFactors',[.25,.37,4],'normalizedAxisTolerance',5e-11,'invariantMetricTolerance',1e-8, ...
    'unitTestScope','spatial construction and paired frozen transfer; not PDE or gauge-rate covariance', ...
    'dependencies',{dependencies},'resolvedDependencies',{resolved},'sourceReview',sourceReview, ...
    'noEllipticOperatorsBuilt',true,'noPDEAdvanced',true,'candidateIsNativeTransaction',false);
write_json(fullfile(outputDirectory,'registration.json'),registration);
records=struct([]);selected=struct([]);covariance=struct([]);
for target=registration.targets
    seedDirectory=fullfile(outputDirectory,sprintf('target%d_fine64_fraction65',target));
    controls=struct('targetXCoreCells',target,'targetYCoreCells',target, ...
        'minimumXFrontCells',20,'positiveFineCells',64,'coreFineCellFraction',.65);
    firstFailure=struct();
    try
        mesh_design_fresh_core_patch(checkpointFile,resultFile,initialDataFile,seedDirectory,controls);
    catch exception
        firstFailure=failure_value(exception);write_json(fullfile(seedDirectory,'direct_design_failure.json'),firstFailure);
    end
    designFile=fullfile(seedDirectory,'design.mat');yFile=fullfile(seedDirectory,'y_axis_input.mat');
    assert(isfile(designFile) && isfile(yFile),'ipm:FreshResolutionReference','Source/reference/y design must exist before family fallback.');
    for fine=registration.fineCells
        for fraction=registration.fractions
            directory=fullfile(outputDirectory,sprintf('target%d_fine%d_fraction%02d',target,fine,round(100*fraction)));
            candidateFile=fullfile(directory,'candidate.mat');
            if ~isfolder(directory),mkdir(directory);end
            record=struct('target',target,'positiveFineCells',fine,'coreFineCellFraction',fraction, ...
                'candidateFile',candidateFile,'frozenPassed',false,'failure',struct(), ...
                'summary',struct(),'axisCovariancePassed',false,'axisCovariance',struct([]));
            try
                if fine==64 && fraction==.65
                    if ~isempty(fieldnames(firstFailure)),error('ipm:FreshResolutionSeedDesign','%s',firstFailure.message);end
                    loaded=load(candidateFile,'candidate');candidate=loaded.candidate;
                else
                    candidate=mesh_fresh_core_patch_candidate(designFile,yFile,candidateFile, ...
                        struct('positiveFineCells',fine,'coreFineCellFraction',fraction));
                end
                record.frozenPassed=candidate.transactionReady;record.summary=candidate_summary(candidate);
                if record.frozenPassed
                    record.axisCovariance=axis_covariance(candidate,registration.unitFactors);
                    record.axisCovariancePassed=all([record.axisCovariance.passed]);
                end
            catch exception
                record.failure=failure_value(exception);
            end
            write_json(fullfile(directory,'frozen_trial_report.json'),record);if isempty(records),records=record;else,records(end+1)=record;end %#ok<AGROW>
            fprintf('FRESH_RESOLUTION_SCREEN target=%d fine=%d fraction=%.2f frozen=%d axisCovariance=%d\n', ...
                target,fine,fraction,record.frozenPassed,record.axisCovariancePassed);
        end
    end
    eligible=find([records.target]==target & [records.frozenPassed] & [records.axisCovariancePassed]);
    if ~isempty(eligible)
        chosen=records(eligible(1));if isempty(selected),selected=chosen;else,selected(end+1)=chosen;end %#ok<AGROW>
        test=paired_unit_covariance(designFile,yFile,chosen,outputDirectory,registration.unitFactors);
        if isempty(covariance),covariance=test;else,covariance(end+1)=test;end %#ok<AGROW>
    end
end
report=struct('kind',registration.kind,'registrationFile',fullfile(outputDirectory,'registration.json'), ...
    'trials',records,'selectedSmallestFineBudget',selected,'selectedPairedUnitCovariance',covariance, ...
    'allRequestedTrialsRecorded',numel(records)==12,'successfulCandidateCount',nnz([records.frozenPassed]), ...
    'allAdmittedAxisCovariancePassed',all([records([records.frozenPassed]).axisCovariancePassed]), ...
    'selectedPairedUnitCovariancePassed',~isempty(covariance) && all([covariance.passed]), ...
    'noLU',true,'noPDE',true,'nativeTransactionNotRun',true,'dynamicQualificationNotGranted',true);
save(fullfile(outputDirectory,'report.mat'),'report','-v7.3');write_json(fullfile(outputDirectory,'report.json'),report);
fprintf('FRESH_RESOLUTION_SCREEN_COMPLETE admitted=%d/12 pairedUnitCovariance=%d\n', ...
    report.successfulCandidateCount,report.selectedPairedUnitCovariancePassed);
end

function summary=candidate_summary(c)
summary=struct('anchorExact',c.anchorExact,'nodeCount',[numel(c.candidateX),numel(c.candidateY)], ...
    'quality',c.pairScore.quality,'resolutionPassed',c.resolutionPassed,'transferPassed',c.transferPassed, ...
    'localQuadraturePassed',c.localQuadraturePassed,'rejectionReasons',{c.pairScore.rejectionReasons}, ...
    'originalUnitCompletePair',c.completePairScoredInOriginalUnits,'aggregate',c.pairScore.aggregate);
end

function tests=axis_covariance(c,factors)
p=c.profiles.records;controls=c.corePatchControls;anchor=c.anchorPosition;
spacing=min((p.safetyCoreWidthLeft+p.safetyCoreWidthRight)/c.controls.targetXCoreCells, ...
    2*p.leftFrontHalfWidth/c.controls.minimumXFrontCells)/controls.resolutionPadding;
baseline=quality_invariants(c.pairScore.quality.x);tests=struct([]);
for factor=factors
    x=mesh_core_patch_axis(anchor*factor,p.safetyCoreCenter*factor,c.candidateX(end)*factor, ...
        numel(c.candidateX),spacing*factor,controls.positiveFineCells,controls.coreFineCellFraction, ...
        controls.roundingCells,controls.anchorWarpRadius);
    q=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));
    axisError=max(abs(x/(anchor*factor)-c.candidateX/anchor)./max(1,abs(c.candidateX/anchor)));
    qualityError=max(abs(quality_invariants(q)-baseline));
    test=struct('factor',factor,'normalizedAxisError',axisError,'qualityInvariantError',qualityError, ...
        'exactAnchor',any(x==anchor*factor) && any(x==-anchor*factor), ...
        'passed',axisError<=5e-11 && qualityError<=1e-8 && hard_quality(q));
    if isempty(tests),tests=test;else,tests(end+1)=test;end %#ok<AGROW>
end
end

function report=paired_unit_covariance(designFile,yFile,chosen,outputDirectory,factors)
loaded=load(designFile,'design');original=loaded.design;loaded=load(yFile,'candidate');seed=loaded.candidate;
loaded=load(chosen.candidateFile,'candidate');baseline=loaded.candidate;
directory=fullfile(outputDirectory,sprintf('target%d_selected_unit_covariance',chosen.target));mkdir(directory);
tests=struct([]);
for factor=factors
    d=fullfile(directory,sprintf('factor_%g',factor));mkdir(d);design=original;
    design.anchor=original.anchor*factor;design.dataset.xLimits=original.dataset.xLimits*factor;
    design.dataset.yLimits=original.dataset.yLimits*factor;
    design.dataset.tag='synthetic_spatial_unit_covariance_only';
    if isfield(design.dataset,'sources'),design.dataset=rmfield(design.dataset,'sources');end
    for k=1:numel(design.dataset.snapshots)
        for name={'x','y','physicalX','physicalY'}
            design.dataset.snapshots(k).(name{1})=original.dataset.snapshots(k).(name{1})*factor;
        end
        design.dataset.snapshots(k).sourceLabel='synthetic coordinate-unit transform, no native state';
    end
    if isfield(design.dataset.snapshots,'signature'),design.dataset.snapshots=rmfield(design.dataset.snapshots,'signature');end
    design.reference.x=original.reference.x*factor;design.reference.y=original.reference.y*factor;
    design.reference.kind='synthetic_unit_transform_of_immutable_fresh_axes';
    design.reference.provenance=struct('unitCovarianceOnly',true,'factor',factor, ...
        'originalDesignFile',designFile,'nativeStateOrInitialDataModified',false);
    if isfield(design.reference,'runtimeReferences'),design.reference=rmfield(design.reference,'runtimeReferences');end
    candidate=seed;candidate.reference=design.reference;candidate.anchorPosition=design.anchor;
    candidate.candidateY=seed.candidateY*factor;
    transformedDesignFile=fullfile(d,'synthetic_design.mat');transformedYFile=fullfile(d,'synthetic_y.mat');
    save(transformedDesignFile,'design','-v7.3');save(transformedYFile,'candidate','-v7.3');
    candidate=mesh_fresh_core_patch_candidate(transformedDesignFile,transformedYFile, ...
        fullfile(d,'synthetic_candidate.mat'),baseline.corePatchControls);
    xError=max(abs(candidate.candidateX/design.anchor-baseline.candidateX/baseline.anchorPosition)./ ...
        max(1,abs(baseline.candidateX/baseline.anchorPosition)));
    metricError=Inf;
    if candidate.pairScore.admissible
        a=aggregate_invariants(candidate.pairScore.aggregate);b=aggregate_invariants(baseline.pairScore.aggregate);
        metricError=max(abs(a-b)./max(1,abs(b)));
    end
    test=struct('factor',factor,'normalizedAxisError',xError,'pairedInvariantError',metricError, ...
        'frozenGatePassed',candidate.transactionReady,'syntheticOnly',true, ...
        'passed',candidate.transactionReady && xError<=5e-11 && metricError<=1e-8);
    write_json(fullfile(d,'unit_test_report.json'),test);if isempty(tests),tests=test;else,tests(end+1)=test;end %#ok<AGROW>
end
report=struct('target',chosen.target,'actualCandidateFile',chosen.candidateFile,'tests',tests,'passed',all([tests.passed]));
end

function v=aggregate_invariants(a)
v=[a.minimumXCoreCells,a.minimumYCoreCells,a.minimumXLeftFrontCells, ...
    a.worst.rhoXMaximumRelativeChange,a.worst.conservationRelativeDefect,a.worst.relativeRangeViolation];
end
function v=quality_invariants(q)
v=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio,q.maximumQuadratureWeightToControlWidthRatio];
end
function passed=hard_quality(q)
passed=q.maximumAdjacentCellRatio<=1.08 && q.maximumLogSpacingCurvature<=.01 && ...
    q.minimumStencilRcond>=1e-9 && q.minimumQuadratureWeightRatio>=1e-8 && ...
    q.quadratureWeightsStrictlyPositive && q.minimumQuadratureWeightToControlWidthRatio>=.35 && ...
    q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function failure=failure_value(exception)
failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
