function report=ipm_perflab_weno_actual_abba(checkpointFile,outputDirectory,options)
%IPM_PERFLAB_WENO_ACTUAL_ABBA One native restore, frozen-state paired timing.
% Generated perflab_b16/perflab_b32 packages contain only the nine numerical
% call-chain functions needed to redirect WENO; neither has solve/restore.
if nargin<3,options=struct();end
d=struct('blocks',4,'seed',20260908,'allowLarge',false);
assert(isstruct(options)&&isscalar(options)&&isempty(setdiff(fieldnames(options),fieldnames(d))));
for n=fieldnames(d)',if ~isfield(options,n{1}),options.(n{1})=d.(n{1});end,end
assert(maxNumCompThreads==10&&~isfolder(outputDirectory));
validateattributes(options.blocks,{'numeric'},{'scalar','integer','>=',2,'<=',6});
mkdir(outputDirectory);stream=RandStream('mt19937ar','Seed',options.seed);
names={'rho_x','rho_y','unit_x','unit_y','entire_rhs','entire_advance'};
blockSizes=[16,32];orders=cell(2,6);
for j=1:2
 for k=1:6
  order=repmat([1,2,2,1],options.blocks,1);
  flip=rand(stream,options.blocks,1)>.5;order(flip,:)=3-order(flip,:);orders{j,k}=order;
 end
end
registration=struct('checkpointFile',checkpointFile,'options',options,'threads',10,'blockSizes',blockSizes, ...
 'kernels',{names},'orders',{orders},'nativeRestoreCount',1,'stepSamplesReturnToSameFrozenInput',true, ...
 'stableWholeGainRule','Bitwise all kernels/RHS/advance and each paired ABBA block speedup >=1.05 for entire RHS or advance', ...
 'sameRestoredAcceptedRhsCache',true,'checkpointWrites',0,'solverPathModified',false);
save(fullfile(outputDirectory,'registration.mat'),'registration','-v7.3');write_json(fullfile(outputDirectory,'registration.json'),registration);
timer=tic;
[state,log,~]=ipm.output.restoreCheckpoint(checkpointFile,struct('finalTime',Inf,'physicalFinalTime',Inf, ...
 'saveResults',false,'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false,'checkpoint',struct('enabled',false)));
restoreSeconds=toc(timer);
assert(options.allowLarge||numel(state.rho)<=20000);
assert(all(ipm.output.trustedMask(log.history,state.config))&&strcmp(state.ops.dynamicScaleGeometry,'isotropic')&& ...
 strcmp(state.config.time.timeIntegrator,'ssprk54')&&strcmp(state.ops.spatialDiscretization,'high_order')&& ...
 strcmp(state.ops.transportScheme,'weno5_fd')&&~state.config.remesh.adaptiveRemesh);clear log
q=state.rho;unit=ones(size(q),'like',q);u=state.flow.transportU1;v=state.flow.transportU2;o=state.ops;
ko=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate','extrapolationDegree',5,'epsilon',o.wenoEpsilon);
inputs={{q,u,o.metricX,o.computationalSpacingX,ko}, ...
 {q',v',o.metricY',o.computationalSpacingY,ko}, ...
 {unit,u,o.metricX,o.computationalSpacingX,ko}, ...
 {unit',v',o.metricY',o.computationalSpacingY,ko}};
records=cell(2,6);candidates=cell(1,2);
for j=1:2
 prefix=sprintf('perflab_b%d',blockSizes(j));candidateExact=true;
 for k=1:6
  if k<=4
   args=inputs{k};a=@()ipm.field.weno5FluxDerivative(args{:});
   method=str2func([prefix '.field.weno5FluxDerivative']);b=@()method(args{:});
  elseif k==5
   a=@()call_rhs(@ipm.evolve.rhs,state);method=str2func([prefix '.evolve.rhs']);b=@()call_rhs(method,state);
  else
   a=@()call_advance(@ipm.evolve.advance,state);method=str2func([prefix '.evolve.advance']);b=@()call_advance(method,state);
  end
  record=paired_measure(a,b,orders{j,k});record.blockSize=blockSizes(j);record.operation=names{k};records{j,k}=record;
  candidateExact=candidateExact&&record.allOutputsBitwise;
  write_json(fullfile(outputDirectory,sprintf('block%d_%s.json',blockSizes(j),names{k})),record);
  fprintf('WENO_ACTUAL block=%d op=%s bitwise=%d speedup=%.6g blockRange=%.6g/%.6g\n', ...
   blockSizes(j),names{k},record.allOutputsBitwise,record.medianSpeedup,min(record.pairedBlockSpeedup),max(record.pairedBlockSpeedup));
  if ~record.allOutputsBitwise,break;end
 end
 gain=false;
 if candidateExact
  gain=records{j,5}.stableFivePercentGain||records{j,6}.stableFivePercentGain;
 end
 candidates{j}=struct('blockSize',blockSizes(j),'allKernelRhsStepBitwise',candidateExact,'eligibleForSeparateMatchedPhysicalTest',candidateExact&&gain);
end
report=struct('registration',registration,'shape',[o.ny,o.nx],'acceptedStep',state.step, ...
 'canonicalTime',state.scale.canonicalTime,'localPhysicalTime',state.scale.physicalTime, ...
 'restoreSeconds',restoreSeconds,'records',{records},'candidates',{candidates},'completed',true, ...
 'anyStableWholeGain',any(cellfun(@(a)a.eligibleForSeparateMatchedPhysicalTest,candidates)), ...
 'checkpointWritten',false,'productionChanged',false,'matlabVersion',version, ...
 'interpretation','Shared-machine frozen-input ABBA timing, with a single native LU cache. No committed trajectory or checkpoint; stable >=5 percent whole gain only permits a separate matched-physical-time test.');
save(fullfile(outputDirectory,'benchmark_report.mat'),'report','-v7.3');write_json(fullfile(outputDirectory,'benchmark_report.json'),report);
fprintf('WENO_ACTUAL_COMPLETE bitwise16=%d bitwise32=%d stableWholeGain=%d restores=1\n', ...
 candidates{1}.allKernelRhsStepBitwise,candidates{2}.allKernelRhsStepBitwise,report.anyStableWholeGain);
end
function r=paired_measure(a,b,order)
reference=a();candidate=b();warmExact=isequaln(reference,candidate);clear candidate
seconds=zeros(size(order));exact=true(size(order));
if warmExact
 for i=1:size(order,1)
  for j=1:4
   timer=tic;if order(i,j)==1,value=a();else,value=b();end
   seconds(i,j)=toc(timer);exact(i,j)=isequaln(value,reference);clear value
  end
 end
else
 exact(:)=false;seconds(:)=NaN;
end
baseline=seconds(order==1);candidate=seconds(order==2);ratio=NaN(size(order,1),1);
for i=1:size(order,1),ratio(i)=median(seconds(i,order(i,:)==1))/median(seconds(i,order(i,:)==2));end
r=struct('order',order,'seconds',seconds,'warmOutputsBitwise',warmExact,'sampleBitwise',exact, ...
 'allOutputsBitwise',warmExact&&all(exact,'all'),'baselineMedianSeconds',median(baseline), ...
 'candidateMedianSeconds',median(candidate),'medianSpeedup',median(baseline)/median(candidate), ...
 'pairedBlockSpeedup',ratio,'stableFivePercentGain',warmExact&&all(exact,'all')&&all(ratio>=1.05));
end
function value=call_rhs(method,s)
[rhs,flow]=method(s.rho,s.ops,s.scale);value=struct('rhs',rhs,'flow',flow);
end
function value=call_advance(method,s)
[advanced,stop]=method(s);assert(isempty(stop)&&advanced.step==s.step+1);
% Both algorithms share the same immutable OPS/decomposition handle; compare
% every remaining state field, including rho, flow, scale and RHS cache.
value=struct('state',rmfield(advanced,'ops'),'stopReason',stop);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,PrettyPrint=true),'char');
end
