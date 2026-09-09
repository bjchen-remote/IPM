function report=mesh_test_paired_designer_actual(projectRoot,outDir)
%MESH_TEST_PAIRED_DESIGNER_ACTUAL Three actual input regimes; no LU/PDE.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
c=fullfile(projectRoot,'result','longtime','20260908_campaign_v1');
rootDir=fullfile(projectRoot,'result','verification','q512_fresh_exact_gauge_roots_v4', ...
    'root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f');
rootFile=fullfile(rootDir,'checkpoint_20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7_step0000001818.mat');
fullResultFile=fullfile(c,'adaptive_campaign_v4','stage_002','result.mat');
freshResultFile=fullfile(c,'fresh_adaptive_campaign_v2','stage_012','result.mat');
initialFile=fullfile(projectRoot,'result','verification','performance_lab_20260908', ...
    'late_box_v4_stage002_epoch_v1','dynamic','crop_natural','initial_physical_data.mat');
controls=struct('positiveFineCells',[32,48,64,80,96,112,128],'roundingCells',[16,24,32,40,48], ...
    'coreFineCellFractions',[.5,.65],'yEqualizationSigmas',[.18,.25,.35,.5], ...
    'targetXCoreCells',32,'targetYCoreCells',32,'minimumXFrontCells',20, ...
    'maximumAxisCandidates',70,'maximumFullPairScores',3);
registration=struct('kind','actual_three_state_general_designer_regression_v1', ...
    'controls',controls,'rootNative',rootFile,'fullResult',fullResultFile,'freshResult',freshResultFile, ...
    'freshInitialData',initialFile,'noLU',true,'noPDE',true,'threads',10, ...
    'unitFactors',[.37,4],'unitRelativeTolerance',1e-8,'candidateNotNativeQualified',true);
write_json(fullfile(outDir,'registration.json'),registration);
root=ipm.output.readCheckpoint(rootFile);s=root.payload.state;
assert(s.remeshCount==0 && ~s.config.remesh.initialAnalyticRemesh);
reference=struct('x',s.baseX,'y',s.baseY,'provenance',struct('kind','native_immutable_original_base_axes','checkpoint',rootFile));
assert(isequal(s.x,s.baseX)&&isequal(s.y,s.baseY));
[X,Y]=meshgrid(reference.x,reference.y);
geometry=struct('X',X,'Y',Y,'nx',numel(reference.x),'ny',numel(reference.y),'symmetryMode',s.config.physics.symmetryMode);
rho=ipm.field.initialDensity(geometry,s.config.physics);
provenance=struct('kind','actual_original_analytic_initial_density_resampled_without_LU', ...
    'sourceNative',rootFile,'physics',s.config.physics,'initialAnalyticRemeshDisabled',true, ...
    'zeroTimeNativeCheckpointClaim',false,'PDEInitialQualificationClaim',false);
snapshot=struct('rho',rho,'x',reference.x,'y',reference.y,'scale',struct('Cx',1,'Cy',1,'Comega',1), ...
    'canonicalTime',0,'physicalTime',0,'trusted',true,'provenance',provenance);
anchor=s.config.scaling.transportAnchorX;clear root s rho X Y geometry
report=struct('registration',registration,'cases',struct([]),'negativeInputs',struct([]),'unitChecks',struct([]));
[candidate,design]=mesh_design_paired_snapshot(snapshot,reference,anchor,controls);
report.cases=append_row(report.cases,save_case(outDir,'original_t0',candidate,design,snapshot,reference,anchor));
clear candidate design snapshot
fullResult=ipm.output.validate(fullResultFile);full=ipm.output.readCheckpoint(char(fullResult.metadata.latestCheckpointFile));
assert_native_pair(full,fullResult);s=full.payload.state;
assert(isequal(s.baseX,reference.x)&&isequal(s.baseY,reference.y)&&s.step==3326);
snapshot=native_snapshot(s,char(fullResult.metadata.latestCheckpointFile));anchor=s.config.scaling.transportAnchorX;
clear full fullResult s
[candidate,design]=mesh_design_paired_snapshot(snapshot,reference,anchor,controls);
report.cases=append_row(report.cases,save_case(outDir,'full_box_step3326',candidate,design,snapshot,reference,anchor));
clear candidate design snapshot reference
[freshReview,fresh,contract]=fresh_profile_campaign_review(freshResultFile,initialFile);
assert(freshReview.step==3668);s=fresh.payload.state;
reference=struct('x',contract.initialX,'y',contract.initialY, ...
    'provenance',struct('kind','actual_immutable_fresh_initial_axes','initialFile',initialFile,'lineage',contract.lineage));
snapshot=native_snapshot(s,freshReview.checkpointFile);anchor=s.config.scaling.transportAnchorX;clear fresh s
[candidate,design]=mesh_design_paired_snapshot(snapshot,reference,anchor,controls);
report.cases=append_row(report.cases,save_case(outDir,'fresh_step3668',candidate,design,snapshot,reference,anchor));
% Negative input checks are inexpensive and must fail before constructing axes.
bad=snapshot;bad.rho=bad.rho(:,1:end-1);
report.negativeInputs=append_row(report.negativeInputs,reject('mismatched_field',@()mesh_design_paired_snapshot(bad,reference,anchor,controls),'ipm:PairedDesignInput'));
bad=snapshot;bad.trusted=false;
report.negativeInputs=append_row(report.negativeInputs,reject('untrusted_input',@()mesh_design_paired_snapshot(bad,reference,anchor,controls),'ipm:PairedDesignInput'));
badRef=reference;badRef.y(end)=2*badRef.y(end);
report.negativeInputs=append_row(report.negativeInputs,reject('changed_box',@()mesh_design_paired_snapshot(snapshot,badRef,anchor,controls),'ipm:PairedDesignInput'));
badControls=controls;badControls.minQuadratureWeightRatio=1e-10;
report.negativeInputs=append_row(report.negativeInputs,reject('relaxed_hard_gate',@()mesh_design_paired_snapshot(snapshot,reference,anchor,badControls),'ipm:PairedDesignControls'));
badControls=controls;badControls.maximumAxisCandidates=1;
report.negativeInputs=append_row(report.negativeInputs,reject('unregistered_search_expansion',@()mesh_design_paired_snapshot(snapshot,reference,anchor,badControls),'ipm:PairedDesignBudget'));
assert(candidate.transactionReady && ~candidate.unchanged,'The latest actual fresh state must have an admitted new candidate.');
choice=design.xTrials(candidate.selectedXTrial);yc=design.yTrials(candidate.selectedYTrial);
single=controls;single.positiveFineCells=choice.positiveFineCells;single.roundingCells=choice.roundingCells;
single.coreFineCellFractions=choice.coreFineCellFraction;single.yEqualizationSigmas=yc.sigma;single.includeUnchangedPair=false;
for factor=registration.unitFactors
    u=snapshot;u.x=factor*u.x;u.y=factor*u.y;u.scale.Cx=factor*u.scale.Cx;u.scale.Cy=factor*u.scale.Cy;
    ref=reference;ref.x=factor*ref.x;ref.y=factor*ref.y;
    [unitCandidate,unitReport]=mesh_design_paired_snapshot(u,ref,factor*anchor,single);
    assert(unitCandidate.transactionReady,'Unit covariance changed candidate admission.');
    axisError=max([max(abs(unitCandidate.candidateX/factor-candidate.candidateX))/max(abs(candidate.candidateX)), ...
        max(abs(unitCandidate.candidateY/factor-candidate.candidateY))/max(abs(candidate.candidateY))]);
    a=candidate.pairScore.aggregate;b=unitCandidate.pairScore.aggregate;
    valuesA=[a.minimumXCoreCells,a.minimumYCoreCells,a.minimumXLeftFrontCells];
    valuesB=[b.minimumXCoreCells,b.minimumYCoreCells,b.minimumXLeftFrontCells];
    metricError=max(abs(valuesA-valuesB)./max(1,abs(valuesA)));
    transferA=cell2mat(struct2cell(a.worst));transferB=cell2mat(struct2cell(b.worst));
    transferError=max(abs(transferA-transferB)./max(1,abs(transferA)));
    passed=max([axisError,metricError,transferError])<=registration.unitRelativeTolerance;
    report.unitChecks=append_row(report.unitChecks,struct('factor',factor,'axisRelativeError',axisError, ...
        'coreFrontRelativeError',metricError,'allPairedTransferMetricError',transferError,'passed',passed));
    save(fullfile(outDir,sprintf('unit_factor_%g.mat',factor)),'unitCandidate','unitReport','-v7.3');
end
report.allChecksPassed=all([report.negativeInputs.passed])&&all([report.unitChecks.passed]);
report.interpretation='Actual input eligibility and frozen design tests; no new native transaction, time integration or box change.';
save(fullfile(outDir,'report.mat'),'report','-v7.3');write_json(fullfile(outDir,'report.json'),report);
assert(report.allChecksPassed);fprintf('PAIRED_DESIGNER_ACTUAL_COMPLETE cases=3 negative=%d units=%d noLU=1\n',numel(report.negativeInputs),numel(report.unitChecks));
end

function row=save_case(out,label,candidate,design,snapshot,reference,anchor)
file=fullfile(out,[label '.mat']);save(file,'candidate','design','snapshot','reference','anchor','-v7.3');
write_json(fullfile(out,[label '_design.json']),design);
row=struct('label',label,'file',file,'status',design.status,'candidateReady',candidate.transactionReady, ...
    'xAxisTrials',numel(design.xTrials),'xAdmitted',nnz([design.xTrials.admissible]), ...
    'yAxisTrials',numel(design.yTrials),'yAdmitted',nnz([design.yTrials.admissible]), ...
    'pairScores',numel(design.pairTrials),'selectedControls',struct(),'aggregate',struct());
if candidate.transactionReady
    row.aggregate=candidate.pairScore.aggregate;
    if ~candidate.unchanged
        row.selectedControls=struct('x',design.xTrials(candidate.selectedXTrial), ...
            'ySigma',design.yTrials(candidate.selectedYTrial).sigma);
    end
end
fprintf('PAIRED_DESIGNER_CASE %s ready=%d status=%s x=%d/%d y=%d/%d\n', ...
    label,row.candidateReady,row.status,row.xAdmitted,row.xAxisTrials,row.yAdmitted,row.yAxisTrials);
end

function s=native_snapshot(n,file)
s=struct('rho',n.rho,'x',n.x,'y',n.y,'scale',struct('Cx',exp(n.scale.logC_l), ...
    'Cy',exp(n.scale.logC_l),'Comega',exp(n.scale.logC_omega)), ...
    'canonicalTime',n.scale.canonicalTime,'physicalTime',n.scale.physicalTime,'trusted',true, ...
    'provenance',struct('kind','strict_native_result_paired_terminal','checkpoint',file,'caseId',n.runMetadata.caseId,'step',n.step));
end

function assert_native_pair(cp,r)
n=cp.payload.state;assert(isequal(n.rho,r.state.rho)&&isequal(n.x(:),r.grid.x(:))&&isequal(n.y(:),r.grid.y(:))&& ...
    n.step==r.state.steps&&n.scale.canonicalTime==r.state.canonicalTime&&n.scale.physicalTime==r.state.physicalTime&& ...
    isequaln(cp.payload.log.history,r.history)&&strcmp(n.runMetadata.caseId,r.metadata.caseId)&& ...
    all(ipm.output.trustedMask(r))&&all(ipm.output.trustedMask(cp.payload.log.history,n.config)));
end

function row=reject(label,f,expected)
identifier='';
try
    f();
catch e
    identifier=e.identifier;
end
row=struct('label',label,'expected',expected,'actual',identifier,'passed',strcmp(identifier,expected));assert(row.passed);
end

function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
