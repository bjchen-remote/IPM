function report=ipm_perflab_extend_box_to_reference(shortReportFile,referenceResultFile,referenceCheckpointFile,outputDirectory,controls)
%IPM_PERFLAB_EXTEND_BOX_TO_REFERENCE Three serial natural-epsilon cases.
% Continue full/crop native fresh histories, initialize only the middle IVP,
% and remeasure covariance against the actual native endpoint at this time.
defaults=struct('absolutePhysicalTarget',NaN,'middleFraction',.1,'maximumAdditionalSteps',128,'allowLarge',false);
assert(isstruct(controls)&&isscalar(controls)&&isempty(setdiff(fieldnames(controls),fieldnames(defaults))));
for name=fieldnames(defaults)'
 if ~isfield(controls,name{1}),controls.(name{1})=defaults.(name{1});end
end
validateattributes(controls.absolutePhysicalTarget,{'numeric'},{'scalar','finite','positive'});
validateattributes(controls.maximumAdditionalSteps,{'numeric'},{'scalar','integer','positive'});
assert(~isfolder(outputDirectory));assert(maxNumCompThreads==10);
old=jsondecode(fileread(shortReportFile));r=old.registration;
if ~iscell(old.branches),old.branches=num2cell(old.branches);end
assert(old.completed&&old.fullCovariance.passed&&all(cellfun(@(b)b.audit.passed,old.branches)));
assert(~r.large||controls.allowLarge);assert(~r.large||old.shortIntervalComparisonPassed);
targetElapsed=controls.absolutePhysicalTarget-r.absolutePhysicalEpoch;
assert(targetElapsed>r.elapsedPhysicalTime&&controls.middleFraction>r.fraction&&controls.middleFraction<1);
if r.large
 static=jsondecode(fileread(r.staticAssessment.file));assert(strcmp(static.registration.checkpointFile,r.sourceFile));
 for f=[controls.middleFraction,r.fraction]
  i=find([static.assessment.fraction]==f);assert(isscalar(i));
  assert(static.assessment(i).passed&&static.assessment(i).eligibleForDynamicTest);
 end
 clear static
end
source=ipm.output.readCheckpoint(r.sourceFile);
reference=ipm.output.validate(referenceResultFile);referenceCP=ipm.output.readCheckpoint(referenceCheckpointFile);
assert(strcmp(char(source.payload.state.runMetadata.caseId),r.parentCaseId)&&source.payload.state.step==r.parentStep);
referenceAudit=paired_audit(reference,referenceCP);
[referencePrefix,referencePrefixAudit]=ipm_gaugelab_history_prefix(source.payload.log.history,reference.history);
referenceAudit.sourceHistoryPrefixExact=referencePrefix;
referenceAudit.physicalTargetError=abs(reference.state.physicalTime-controls.absolutePhysicalTarget);
referenceAudit.sameGridAsEpoch=isequal(reference.grid.x(:),source.payload.state.x(:))&&isequal(reference.grid.y(:),source.payload.state.y(:));
referenceAudit.sameNativeCase=strcmp(char(reference.metadata.caseId),r.parentCaseId);
referenceAudit.completedStopReason=any(strcmp(reference.state.stopReason,{'final_time','physical_final_time'}));
referenceAudit.productionResolution=all([reference.history.mesh.coreGridPoints(end),reference.history.mesh.verticalCoreGridPoints(end)]>=r.limits.productionMinimumCoreCells)&& ...
 reference.history.mesh.safetyFactor(end)<r.limits.productionMaximumSafety;
referenceAudit.passed=referenceAudit.passed&&referencePrefix&&referenceAudit.sameGridAsEpoch&&referenceAudit.sameNativeCase&& ...
 referenceAudit.completedStopReason&&referenceAudit.physicalTargetError<=r.limits.physicalEndpointAbsolute&& ...
 (referenceAudit.productionResolution||~r.large);
assert(referenceAudit.passed,'ipm:BoxReferenceGate','Native endpoint reference failed; no new solve was started.');
referenceObservation=observe(reference,1);clear reference referenceCP source
labels={'full_natural','crop_natural','middle_natural'};fractions=[1,r.fraction,controls.middleFraction];seeds=[2,3,0];
registration=struct('kind','three_natural_boxes_with_actual_endpoint_covariance','shortReportFile',shortReportFile, ...
 'parentRegistration',r,'referenceResultFile',referenceResultFile,'referenceCheckpointFile',referenceCheckpointFile, ...
 'referenceAudit',referenceAudit,'referencePrefixAudit',referencePrefixAudit,'controls',controls, ...
 'absolutePhysicalTarget',controls.absolutePhysicalTarget,'targetElapsedPhysicalTime',targetElapsed,'caseOrder',{labels}, ...
 'fractions',fractions,'seeds',seeds,'threads',10,'pdeCases',3,'ellipticAssemblies',3,'referenceOperatorRestores',0, ...
 'epsilonChoice','Natural public epsilon only, preregistered; prior common-epsilon control is not extended to this longer time', ...
 'samePhysicalWindow',r.physicalWindow,'limits',r.limits,'nativeCovarianceCheckedAtTarget',true, ...
 'productionDomainModified',false,'dynamicBoxConvergenceEstablished',false);
mkdir(outputDirectory);write_json(fullfile(outputDirectory,'registration.json'),registration);
save(fullfile(outputDirectory,'registration.mat'),'registration','-v7.3');
save(fullfile(outputDirectory,'native_reference_observation.mat'),'referenceObservation','-v7.3');
branches=cell(1,3);observations=cell(1,3);comparisons=struct();covariance=struct();
for k=1:3
 directory=fullfile(outputDirectory,labels{k});timer=tic;
 seedHistory=[];seedReferences=[];seedLineage=[];seedConfig=[];
 if seeds(k)>0
  seed=old.branches{seeds(k)};cp=ipm.output.readCheckpoint(seed.checkpointFile);s=cp.payload.state;
  lineage=s.runMetadata.caseMetadata.latePhysicalBoxBranch;
  assert(strcmp(char(s.runMetadata.caseId),seed.newCaseId)&&strcmp(lineage.parentCheckpoint,r.sourceFile)&& ...
   lineage.absolutePhysicalEpoch==r.absolutePhysicalEpoch&&lineage.boxFraction==fractions(k)&& ...
   s.scale.physicalTime==seed.elapsedPhysicalTime&&~s.config.remesh.adaptiveRemesh&& ...
   s.config.transport.wenoEpsilon==r.naturalPublicWenoEpsilon&& ...
   all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
  initialStep=s.step;initialElapsed=s.scale.physicalTime;
  seedHistory=cp.payload.log.history;seedReferences=s.rescaling;seedLineage=lineage;seedConfig=s.config;
  clear cp s
  mkdir(directory);
  opts=struct('finalTime',Inf,'physicalFinalTime',targetElapsed,'maxSteps',initialStep+controls.maximumAdditionalSteps, ...
   'saveResults',true,'resultFile',fullfile(directory,'result.mat'),'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false, ...
   'checkpoint',struct('enabled',true,'file',fullfile(directory,'checkpoint.mat'),'every',Inf,'atExit',true));
  result=ipm.solve(opts,seed.checkpointFile);
  files=dir(fullfile(directory,'checkpoint_*_step*.mat'));assert(isscalar(files));
  branch=struct('label',labels{k},'newCaseId',char(result.metadata.caseId),'resultFile',char(result.metadata.resultFile), ...
   'checkpointFile',fullfile(files.folder,files.name),'lineage',lineage,'seedCheckpoint',seed.checkpointFile, ...
   'existingFreshHistoryContinued',true);
  assert(strcmp(branch.newCaseId,seed.newCaseId));
 else
  [branch,result]=mesh_fresh_box_branch(r.sourceFile,fractions(k),targetElapsed,directory, ...
   struct('maximumSteps',controls.maximumAdditionalSteps,'allowLarge',controls.allowLarge));
  initialStep=0;initialElapsed=0;branch.label=labels{k};branch.seedCheckpoint='';branch.existingFreshHistoryContinued=false;
 end
 terminal=ipm.output.readCheckpoint(branch.checkpointFile);
 branch.wallSeconds=toc(timer);branch.segmentAcceptedSteps=result.state.steps-initialStep;
 branch.elapsedPhysicalTime=result.state.physicalTime;branch.segmentElapsedPhysicalTime=result.state.physicalTime-initialElapsed;
 branch.absolutePhysicalTime=r.absolutePhysicalEpoch+result.state.physicalTime;
 a=paired_audit(result,terminal);
 a.physicalEndpointError=abs(branch.absolutePhysicalTime-controls.absolutePhysicalTarget);
 a.stopReason=result.state.stopReason;a.physicalEndpointStop=strcmp(a.stopReason,'physical_final_time');
 a.freshHistoryStartsAtZero=result.history.common.physicalTime(1)==0&&result.history.common.acceptedStep(1)==0;
 a.coreCells=[result.history.mesh.coreGridPoints(end),result.history.mesh.verticalCoreGridPoints(end)];
 a.safety=result.history.mesh.safetyFactor(end);a.productionResolutionPassed=all(a.coreCells>=r.limits.productionMinimumCoreCells)&&a.safety<r.limits.productionMaximumSafety;
 a.finitePhysicalFields=all(isfinite(result.physical.rho),'all')&&all(isfinite(result.physical.rhoX),'all')&&all(isfinite(result.physical.rhoY),'all');
 a.withinAdditionalStepCap=branch.segmentAcceptedSteps<=controls.maximumAdditionalSteps;
 a.historyPrefixExact=true;a.runtimeReferencesExact=true;a.lineageExact=true;a.numericalSpecificationExact=true;
 if seeds(k)>0
  [a.historyPrefixExact,a.historyPrefixAudit]=ipm_gaugelab_history_prefix(seedHistory,terminal.payload.log.history);
  a.runtimeReferencesExact=references_equal(seedReferences,terminal.payload.state.rescaling);
  a.lineageExact=isequaln(seedLineage,terminal.payload.state.runMetadata.caseMetadata.latePhysicalBoxBranch);
  a.numericalSpecificationExact=numerical_specification_equal(seedConfig,terminal.payload.state.config);
 else
  a.lineageExact=strcmp(branch.lineage.parentCheckpoint,r.sourceFile)&&branch.lineage.boxFraction==fractions(k)&& ...
   branch.lineage.absolutePhysicalEpoch==r.absolutePhysicalEpoch&&~branch.lineage.parentHistoryInherited&& ...
   ~strcmp(branch.newCaseId,r.parentCaseId);
 end
 a.passed=a.passed&&a.freshHistoryStartsAtZero&&a.finitePhysicalFields&&a.physicalEndpointStop&&a.withinAdditionalStepCap&& ...
  a.historyPrefixExact&&a.runtimeReferencesExact&&a.lineageExact&&a.numericalSpecificationExact&& ...
  a.physicalEndpointError<=r.limits.physicalEndpointAbsolute&&(a.productionResolutionPassed||~r.large);
 branch.audit=a;observation=observe(result,r.canonicalCovarianceFactor);clear result terminal seedHistory seedReferences seedLineage seedConfig
 branches{k}=branch;observations{k}=observation;
 save(fullfile(directory,'observation.mat'),'observation','-v7.3');save(fullfile(directory,'branch_report.mat'),'branch','-v7.3');
 write_json(fullfile(directory,'branch_report.json'),branch);clear observation
 fprintf('THREE_BOX_CASE %d/3 %s passed=%d steps=%d endpoint=%.3e prefix=%d core=%.8g/%.8g\n', ...
  k,labels{k},a.passed,branch.segmentAcceptedSteps,a.physicalEndpointError,a.historyPrefixExact,a.coreCells);
 if ~a.passed
  write_failure(outputDirectory,'branch_gate',labels{k},a);error('ipm:ThreeBoxBranchGate','Branch failed; dependent cases were not started.');
 end
 if k==1
  covariance=paired_covariance(referenceObservation,observations{1},r.limits);
  write_json(fullfile(outputDirectory,'actual_endpoint_covariance.json'),covariance);
  if ~covariance.passed
   write_failure(outputDirectory,'actual_endpoint_covariance',labels{k},covariance);error('ipm:ThreeBoxCovariance','Actual endpoint covariance failed; crops not advanced.');
  end
  clear referenceObservation
 elseif k==2
  comparisons.fullVsSmall=mesh_compare_fresh_box_data(observations{1},observations{2},r.physicalWindow,r.limits);
  write_json(fullfile(outputDirectory,'full_vs_small.json'),comparisons.fullVsSmall);
  if ~comparisons.fullVsSmall.passed
   write_failure(outputDirectory,'full_vs_small',labels{k},comparisons.fullVsSmall);error('ipm:ThreeBoxComparison','Long small/full comparison failed; middle not started.');
  end
 else
  comparisons.fullVsMiddle=mesh_compare_fresh_box_data(observations{1},observations{3},r.physicalWindow,r.limits);
  comparisons.middleVsSmall=mesh_compare_fresh_box_data(observations{3},observations{2},r.physicalWindow,r.limits);
 end
end
comparisons.passed=comparisons.fullVsSmall.passed&&comparisons.fullVsMiddle.passed&&comparisons.middleVsSmall.passed;
ids=cellfun(@(b)b.newCaseId,branches,'UniformOutput',false);assert(numel(unique(ids))==3);
report=struct('registration',registration,'branches',{branches},'actualEndpointCovariance',covariance, ...
 'naturalEpsilonLadder',comparisons,'completed',true,'passed',covariance.passed&&comparisons.passed, ...
 'nativeCovarianceCheckedAtTarget',true,'commonEpsilonControlExtended',false,'dynamicBoxConvergenceEstablished',false, ...
 'productionDomainModified',false,'interpretation','Three natural-epsilon boxes at the actual longer native reference endpoint; prior common-epsilon control is not extended. No infinite-box or earlier-history claim.');
save(fullfile(outputDirectory,'three_box_report.mat'),'report','-v7.3');write_json(fullfile(outputDirectory,'three_box_report.json'),report);
fprintf('THREE_BOX_COMPLETE covariance=%d ladder=%d productionPromotion=0\n',covariance.passed,comparisons.passed);
end

function a=paired_audit(result,cp)
s=cp.payload.state;
a=struct('entireHistoryTrusted',all(ipm.output.trustedMask(result)), ...
 'nativeEntireHistoryTrusted',all(ipm.output.trustedMask(cp.payload.log.history,s.config)), ...
 'nativeFullHistoryExact',isequaln(result.history,cp.payload.log.history), ...
 'nativeTerminalScalesExact',exp(s.scale.logC_l)==result.scale.Cx&&exp(s.scale.logC_l)==result.scale.Cy&& ...
 exp(s.scale.logC_omega)==result.scale.Comega&&s.scale.X_shift==result.scale.Xshift, ...
 'nativeTerminalPaired',isequal(s.rho,result.state.rho)&&isequal(s.x(:),result.grid.x(:))&&isequal(s.y(:),result.grid.y(:))&& ...
 s.step==result.state.steps&&s.scale.physicalTime==result.state.physicalTime&&s.scale.canonicalTime==result.state.canonicalTime&& ...
 strcmp(char(s.runMetadata.caseId),char(result.metadata.caseId)));
a.passed=a.entireHistoryTrusted&&a.nativeEntireHistoryTrusted&&a.nativeTerminalPaired&&a.nativeFullHistoryExact&&a.nativeTerminalScalesExact;
end
function passed=numerical_specification_equal(a,b)
% Deserialized anonymous handles have distinct MATLAB identities. Compare
% their complete FUNCTIONS records, including captured numeric arrays.
passed=true;
for domain={'grid','physics','transport','elliptic','scaling','remesh'}
 x=a.(domain{1});y=b.(domain{1});passed=passed&&isequal(fieldnames(x),fieldnames(y));
 for name=fieldnames(x)'
  u=x.(name{1});v=y.(name{1});
  if isa(u,'function_handle')&&isa(v,'function_handle'),u=functions(u);v=functions(v);end
  passed=passed&&isequaln(u,v);
 end
end
end
function passed=references_equal(a,b)
names=fieldnames(a);names=names(startsWith(names,'reference')|startsWith(names,'adaptiveTarget')|startsWith(names,'safety')| ...
 ismember(names,{'pinX','strainTarget','peakTrackingHalfWidth'}));passed=true;
for k=1:numel(names),passed=passed&&isfield(b,names{k})&&isequaln(a.(names{k}),b.(names{k}));end
end
function o=observe(result,rateScale)
p=result.physical;h=result.history;
o=struct('x',p.x,'y',p.y,'rho',p.rho,'rhoX',p.rhoX,'rhoY',p.rhoY, ...
 'cLParentUnits',result.scale.rates.cx/rateScale,'cOmegaParentUnits',result.scale.rates.comega/rateScale, ...
 'physicalGradientMaximum',h.common.physicalRhoXInf(end),'physicalQuadraticPeak',h.common.physicalQuadraticPeak(end), ...
 'physicalCoreWidthX',h.mesh.trackedWallCoreWidth(end)/result.scale.Cx,'physicalCoreWidthY',h.mesh.trackedVerticalCoreWidth(end)/result.scale.Cy);
end
function c=paired_covariance(a,b,l)
assert(isequal(size(a.rho),size(b.rho)));
c=struct('coordinateRelativeInf',max(relative_inf(a.x,b.x),relative_inf(a.y,b.y)), ...
 'rhoRelativeInf',relative_inf(a.rho,b.rho),'rhoXRelativeInf',relative_inf(a.rhoX,b.rhoX),'rhoYRelativeInf',relative_inf(a.rhoY,b.rhoY), ...
 'cLRelativeDifference',relative_inf(a.cLParentUnits,b.cLParentUnits),'cOmegaRelativeDifference',relative_inf(a.cOmegaParentUnits,b.cOmegaParentUnits));
c.passed=c.coordinateRelativeInf<=l.covarianceCoordinateRelativeInf&& ...
 max([c.rhoRelativeInf,c.rhoXRelativeInf,c.rhoYRelativeInf])<=l.covarianceFieldRelativeInf&& ...
 max(c.cLRelativeDifference,c.cOmegaRelativeDifference)<=l.covarianceRateRelative;
end
function v=relative_inf(a,b)
v=max(abs(a-b),[],'all')/max(max(abs(a),[],'all'),realmin);
end
function write_failure(directory,stage,label,audit)
value=struct('stage',stage,'label',label,'audit',audit,'dependentsStopped',true,'productionDomainModified',false);
write_json(fullfile(directory,'failure.json'),value);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,PrettyPrint=true),'char');
end
