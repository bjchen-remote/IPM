function report=ipm_perflab_profile_frozen_state(state,outputDirectory)
%IPM_PERFLAB_PROFILE_FROZEN_STATE Two original flows on an existing LU cache.
% No restore, advance, matrix construction, candidate kernel or checkpoint.
assert(maxNumCompThreads==10&&~isfolder(outputDirectory));
assert(isfield(state,'rhsCache')&&strcmp(state.ops.dynamicScaleGeometry,'isotropic')&& ...
 strcmp(state.ops.rescaling.lengthGauge,'transport_anchor')&& ...
 strcmp(state.ops.rescaling.cOmegaGauge,'wall_omega_quadratic_peak'));
status=profile('status');assert(strcmp(status.ProfilerStatus,'off'),'ipm:ProfilerBusy','Refusing to interrupt an active profiler.');
mkdir(outputDirectory);before=state;oldProfile=profile('info');
if ~isempty(oldProfile.FunctionTable),save(fullfile(outputDirectory,'prior_profiler_data.mat'),'oldProfile','-v7.3');end
clear oldProfile
timer=tic;[rhoRate,flow]=ipm.evolve.flow(state.rho,state.ops,state.scale);unprofiledSeconds=toc(timer);
first=cache_audit(state,rhoRate,flow);clear rhoRate flow
assert(first.passed,'ipm:FrozenProfileCache','Unprofiled flow did not match the existing accepted-state cache.');
profile clear;
cleanup=onCleanup(@()profile('off'));
profile('on','-timer','performance','-detail','builtin','-nohistory');
timer=tic;[rhoRate,flow]=ipm.evolve.flow(state.rho,state.ops,state.scale);profiledSeconds=toc(timer);
profile off;profileInfo=profile('info');
second=cache_audit(state,rhoRate,flow);clear rhoRate flow
stateUnchanged=isequaln(before,state);clear before
assert(second.passed&&stateUnchanged,'ipm:FrozenProfileMutation','Profiled flow/cache/state equality failed.');
ft=profileInfo.FunctionTable;
keys={'flow','rhs','velocity','poisson','greenBoundary','isotropicGauge','assembleRhs','trackFeatures'};
tails={'/+evolve/flow.m','/+evolve/rhs.m','/+field/velocity.m','/+field/poisson.m', ...
 '/+field/greenBoundary.m','/+evolve/isotropicGauge.m','/+evolve/assembleRhs.m','/+diagnostics/trackFeatures.m'};
indices=struct();totals=struct();calls=struct();
for k=1:numel(keys)
 index=primary_index(ft,tails{k},keys{k});indices.(keys{k})=index;
 totals.(keys{k})=ft(index).TotalTime;calls.(keys{k})=ft(index).NumCalls;
end
index=primary_index(ft,'/decomposition.m','mldivide');indices.decompositionBackslash=index;
totals.decompositionBackslash=ft(index).TotalTime;calls.decompositionBackslash=ft(index).NumCalls;
assert(all(structfun(@(n)n==1,calls)),'ipm:FrozenProfileCallGraph','Expected one quadratic/anchor flow and one cached solve/assembly.');
names={'cachedDecompositionBackslash','greenBoundary','remainingPoisson','remainingVelocity', ...
 'transportAssembly','remainingIsotropicGauge','featureTracking','remainingRhs','flowWrapper'};
seconds=[totals.decompositionBackslash,totals.greenBoundary, ...
 totals.poisson-totals.decompositionBackslash-totals.greenBoundary, ...
 totals.velocity-totals.poisson,totals.assembleRhs,totals.isotropicGauge-totals.assembleRhs, ...
 totals.trackFeatures,totals.rhs-totals.velocity-totals.isotropicGauge-totals.trackFeatures, ...
 totals.flow-totals.rhs];
partition=struct('names',{names},'seconds',seconds,'fractionOfProfiledFlow',seconds/totals.flow, ...
 'sumSeconds',sum(seconds),'profiledFlowSeconds',totals.flow, ...
 'interpretation','Disjoint call-graph differences, under the verified single-flow/single-solve/single-assembly path. Includes profiler overhead. Cached decomposition backslash is inclusive of its internal native solve.');
functionRows=cell(1,numel(ft));lineRows=struct('functionName',{},'file',{},'line',{},'executions',{},'seconds',{});
for k=1:numel(ft)
 childSeconds=0;if ~isempty(ft(k).Children),childSeconds=sum([ft(k).Children.TotalTime]);end
 functionRows{k}=struct('index',k,'name',ft(k).FunctionName,'file',ft(k).FileName,'type',ft(k).Type, ...
  'calls',ft(k).NumCalls,'totalSeconds',ft(k).TotalTime,'childSeconds',childSeconds, ...
  'selfSecondsByChildDifference',ft(k).TotalTime-childSeconds,'partialData',ft(k).PartialData);
 lines=ft(k).ExecutedLines;
 for j=1:size(lines,1)
  lineRows(end+1)=struct('functionName',ft(k).FunctionName,'file',ft(k).FileName, ...
   'line',lines(j,1),'executions',lines(j,2),'seconds',lines(j,3)); %#ok<AGROW>
 end
end
[~,order]=sort(cellfun(@(a)a.selfSecondsByChildDifference,functionRows),'descend');functionRows=functionRows(order);
if ~isempty(lineRows),[~,order]=sort([lineRows.seconds],'descend');lineRows=lineRows(order);end
report=struct('kind','existing_cache_exclusive_original_flow_profile','nodes',[state.ops.nx,state.ops.ny], ...
 'acceptedStep',state.step,'canonicalTime',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
 'sourceFlow',which('ipm.evolve.flow'),'sourceRhs',which('ipm.evolve.rhs'), ...
 'sourceAssemble',which('ipm.evolve.assembleRhs'),'sourcePoisson',which('ipm.field.poisson'), ...
 'threads',10,'originalFlowCalls',2,'unprofiledSeconds',unprofiledSeconds,'profiledSeconds',profiledSeconds, ...
 'profileWallOverheadRatio',profiledSeconds/unprofiledSeconds,'profileTimer','performance','profileDetail','builtin', ...
 'firstCacheAudit',first,'secondCacheAudit',second,'stateUnchangedExact',stateUnchanged, ...
 'verifiedProfileCallCounts',calls,'primaryFunctionIndices',indices,'partition',partition, ...
 'functionsBySelfSeconds',{functionRows},'linesBySeconds',lineRows,'profiledSourceChanged',any([ft.PartialData]), ...
 'restores',0,'operatorBuilds',0,'poissonFactorizations',0,'pdeSteps',0,'checkpointWrites',0, ...
 'passed',first.passed&&second.passed&&stateUnchanged, ...
 'interpretation','Two original flow evaluations using caller-owned factor/cache. Full rhoRate, derived scaleRate, flow and cache match exactly; no state update. Profiler function/line times are diagnostic, not a speed benchmark.');
save(fullfile(outputDirectory,'profile_raw.mat'),'profileInfo','-v7.3');
save(fullfile(outputDirectory,'profile_report.mat'),'report','-v7.3');
fid=fopen(fullfile(outputDirectory,'profile_report.json'),'w');assert(fid>=0);fileCleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
fprintf('FROZEN_PROFILE_PASS nodes=%dx%d step=%d unprofiled=%.6g profiled=%.6g cachedSolveFraction=%.6g\n', ...
 report.nodes,state.step,unprofiledSeconds,profiledSeconds,partition.fractionOfProfiledFlow(1));
clear cleanup
end
function a=cache_audit(s,rhoRate,flow)
cache=ipm.evolve.makeRhsCache(s.rho,rhoRate,flow,s.ops,s.scale);
a=struct('rhoRateExact',isequaln(rhoRate,s.rhsCache.rhoRate),'scaleRateExact',isequaln(cache.scaleRate,s.rhsCache.scaleRate), ...
 'entireFlowExact',isequaln(flow,s.flow),'entireCacheExact',isequaln(cache,s.rhsCache));
a.passed=a.rhoRateExact&&a.scaleRateExact&&a.entireFlowExact&&a.entireCacheExact;
end
function index=primary_index(ft,tail,name)
keep=false(size(ft));
for k=1:numel(ft)
 short=regexp(ft(k).FunctionName,'[^.>/]+$','match','once');
 keep(k)=endsWith(strrep(ft(k).FileName,'\','/'),tail)&&strcmp(short,name);
end
index=find(keep);assert(isscalar(index),'ipm:ProfilerFunctionIndex','Missing or ambiguous primary profile function: %s',name);
end
