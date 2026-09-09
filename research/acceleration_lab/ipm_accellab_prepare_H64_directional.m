function report=ipm_accellab_prepare_H64_directional(registrationFile)
% Pure, explicitly non-native operator inputs from one immutable actual CP.
reg=jsondecode(fileread(registrationFile));assert(maxNumCompThreads==10&&~isfolder(reg.outputDirectory));
mkdir(reg.outputDirectory);profile clear;profile on;timer=tic;
try
 cp=ipm.output.readCheckpoint(reg.checkpointFile);s=cp.payload.state;m=s.runMetadata.autonomousMesh;
 assert(s.step==reg.sourceStep&&strcmp(s.runMetadata.caseId,reg.caseId)&&m.version==2&&m.currentLevelId==4&& ...
  isequal([numel(s.x),numel(s.y)],[641,321])&&all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
 p=s.config.remesh.autonomousMesh;draft=p;draft.version=3;draft.nodeFamily=struct('maximumTotalNodes',310000);
 draft=ipm.config.autonomousMeshPolicy(draft);
 assert(isequaln(rmfield(draft,{'version','nodeFamily'}),rmfield(p,{'version','nodeFamily'})));
 family=ipm.remesh.referenceAxisFamily(m.initialization.selectedBaseX,m.initialization.selectedBaseY,m.referenceFamily.anchor,draft);
 assert(isequaln(family.members(1:4),m.referenceFamily.members)&& ...
  isequaln(family.rootX,m.referenceFamily.rootX)&&isequaln(family.rootY,m.referenceFamily.rootY));
 assert(isequal(vertcat(family.members.cellFactors),[1,1;2,1;1,2;2,2;3,1;1,3;3,2;2,3]));
 Dx=ipm.mesh.fdMatrix(s.x,1,7);Omega=s.rho*Dx';
 view=struct('rho',s.rho,'x',s.x,'y',s.y,'Dx',Dx,'source',Omega,'trusted',true);
 feature=ipm.diagnostics.meshFeatureIntervals(view);
 snapshot=struct('datasetIndex',1,'sourceLabel',reg.checkpointFile,'trusted',true,'rho',s.rho,'x',s.x,'y',s.y, ...
  'scale',struct('Cx',exp(s.scale.logC_l),'Cy',exp(s.scale.logC_l),'Comega',exp(s.scale.logC_omega)), ...
  'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime);
 dataset=struct('schemaVersion',2,'kind','ipm_frozen_grid_dataset','tag','actual_H64_directional_operator_input', ...
  'snapshotCount',1,'snapshots',snapshot,'xLimits',s.x([1,end]),'yLimits',s.y([1,end])');
 limits=rmfield(p.qualityLimits,{'minWeightToControlWidth','maxWeightToControlWidth'});
 report=struct('kind','actual_H64_directional_operator_preparation_v1','registration',reg,'sourceStep',s.step, ...
  'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime,'sourceNodeCount',[numel(s.x),numel(s.y)], ...
  'sourceCore',feature.actualCoreCells,'sourceFront',feature.leftFrontCells,'sourceConfigUnchanged',true, ...
  'nativeSignatureValidated',true,'allHistoryTrusted',true,'oldFourMembersBitwise',true, ...
  'familyGeneratedFromActualInitialSelectedBase',true,'sourcePolicyVersion',2,'experimentalPolicyVersion',3, ...
  'nativeTransferPerformed',false,'checkpointWritten',false,'noLU',true,'noPDE',true,'directions',{{}});
 for id=[7,8]
  member=family.members(id);assert(member.resourceAdmitted&&member.qualityPassed);
  dpath=fullfile(reg.outputDirectory,sprintf('member_%d_%dx%d',id,member.nodeCount));mkdir(dpath);
  [candidates,axisReport]=ipm.remesh.plannedAxisPairs(view,struct('x',member.baseX,'y',member.baseY),family.anchor,draft);
  save(fullfile(dpath,'axis_proposals.mat'),'candidates','axisReport','feature','-v7.3');
  record=struct('memberId',id,'cellFactors',member.cellFactors,'nodeCount',member.nodeCount, ...
   'plannerStatus',axisReport.status,'candidateCount',numel(candidates),'selectedIndex',1, ...
   'readyForOperatorResourceProbe',false,'LUResourceQualified',false,'nativeTransferQualified',false,'PDEQualified',false);
  if isempty(candidates)
   record.status='rejected_registered_axis_family';write_json(fullfile(dpath,'report.json'),record);
   report.directions{end+1}=record;continue
  end
  c=candidates(1);assert(~c.unchanged&&numel(c.x)==member.nodeCount(1)&&numel(c.y)==member.nodeCount(2));
  score=ipm_gridlab_score_frozen_pair(dataset,c.x,c.y,struct('trustedOnly',true,'meshLimits',limits, ...
   'minimumXCoreCells',p.targetCoreCells(1),'minimumYCoreCells',p.targetCoreCells(2),'minimumXFrontCells',p.minimumFrontCells));
  assert(score.admissible);z=score.aggregate;
  rx=ipm.remesh.interpolate(s.x,s.rho,c.x,struct('sampleDimension',2,'conservation','constant','stencilWidth',6));
  rhoFrozen=ipm.remesh.interpolate(s.y,rx,c.y,struct('sampleDimension',1,'conservation','constant','stencilWidth',6));
  newDx=ipm.mesh.fdMatrix(c.x,1,7);nf=ipm.diagnostics.meshFeatureIntervals(struct('rho',rhoFrozen,'x',c.x,'y',c.y, ...
   'Dx',newDx,'source',rhoFrozen*newDx','trusted',true));
  gates=struct('wholeMesh',score.admissible,'sourceCore',all([z.minimumXCoreCells,z.minimumYCoreCells]>=p.targetCoreCells-1e-6), ...
   'sourceFront',z.minimumXLeftFrontCells>=p.minimumFrontCells-1e-6,'interpolatedCore',all(nf.actualCoreCells>=p.transactionMinimumCoreCells), ...
   'interpolatedFront',nf.leftFrontCells>=p.minimumFrontCells,'peak',z.worst.rhoXMaximumRelativeChange<=p.maximumSinglePeakJump, ...
   'cumulativePeak',m.cumulativeAbsolutePeakJump+z.worst.rhoXMaximumRelativeChange<=p.maximumCumulativeAbsolutePeakJump, ...
   'mass',z.worst.conservationRelativeDefect<=p.maximumMassRelativeDefect, ...
   'range',z.worst.relativeRangeViolation<=p.maximumRelativeRangeViolation);
  ready=all(structfun(@(v)islogical(v)&&isscalar(v)&&v,gates));
  parameters=struct('xIndex',c.xIndex,'yIndex',c.yIndex,'xTrial',axisReport.xTrials(c.xIndex),'yTrial',axisReport.yTrials(c.yIndex));
  bundle=struct('kind','actual_H64_directional_operator_only_input_v1','sourceFile',reg.checkpointFile,'sourceStep',s.step, ...
   'sourceScale',s.scale,'sourceNormalizedTime',s.normalizedTime,'sourceRemeshCount',s.remeshCount,'sourceConfig',s.config, ...
   'sourceReferences',s.rescaling,'originalReferenceFamily',m.referenceFamily,'futureDraftPolicy',draft, ...
   'futureReferenceFamily',family,'targetMemberId',id,'candidate',c,'candidateParameters',parameters,'rhoFrozen',rhoFrozen, ...
   'score',score,'interpolatedFeature',nf,'inputGates',gates,'readyForOperatorResourceProbe',ready, ...
   'nativeTransferQualified',false,'LUResourceQualified',false,'PDEQualified',false, ...
   'interpolation','Six-point constant-conservative X then Y; deliberately NOT native Y bubble transfer.');
  save(fullfile(dpath,'resource_bundle.mat'),'bundle','-v7.3');
  record.inputGates=gates;record.readyForOperatorResourceProbe=ready;record.interpolatedCore=nf.actualCoreCells;
  record.interpolatedFront=nf.leftFrontCells;record.pairAggregate=z;record.candidateParameters=parameters;
  record.resourceBundleFile=fullfile(dpath,'resource_bundle.mat');
  if ready,record.status='input_ready_resource_unmeasured';else,record.status='rejected_frozen_input_gate';end
  write_json(fullfile(dpath,'report.json'),record);report.directions{end+1}=record;
 end
 assert(isequaln(s,cp.payload.state));
 profile off;profileInfo=profile('info');files={profileInfo.FunctionTable.FileName};
 for suffix={'+mesh/build.m','+evolve/flow.m','+field/velocity.m','+field/poisson.m','+remesh/transfer.m', ...
   '+evolve/advance.m','+output/restoreCheckpoint.m','+ipm/solve.m'}
  assert(~any(endsWith(files,suffix{1})),['Forbidden call: ',suffix{1}]);
 end
 report.callGraphNoLUVerified=true;report.wallSeconds=toc(timer);report.completed=true;
 save(fullfile(reg.outputDirectory,'report.mat'),'report','profileInfo','family','-v7.3');
 write_json(fullfile(reg.outputDirectory,'report.json'),report);
 fprintf('H64_DIRECTIONAL_PREP step=%d directions=%d ready=%d noLU=1 noPDE=1\n',s.step,numel(report.directions), ...
  sum(cellfun(@(d)d.readyForOperatorResourceProbe,report.directions)));
catch e
 profile off;failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
 write_json(fullfile(reg.outputDirectory,'failure.json'),failure);rethrow(e)
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
