function report=mesh_audit_controller_complexity(nativeFiles,outDir)
%MESH_AUDIT_CONTROLLER_COMPLEXITY Read-only native telemetry and pure timings.
% Benchmark states and saved payloads below are synthetic, not checkpoints.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
report=struct('kind','controller_window_complexity_nolu_v1','noLU',true,'noPDE',true, ...
    'sampleCounts',[32,128,512,1024,4096,16384,65536], ...
    'actualNative',struct([]),'microbenchmarks',struct([]));
for name={'ipm.remesh.controllerTelemetry','ipm.evolve.planAutonomousMesh',mfilename}
    issues=checkcode(which(name{1}),'-id');assert(isempty(issues),jsonencode(issues));
end
for k=1:numel(nativeFiles)
    checkpoint=ipm.output.readCheckpoint(nativeFiles{k});s=checkpoint.payload.state;h=checkpoint.payload.log.history;
    times=h.common.canonicalTau(:);epochs=h.mesh.remeshCount(:);currentEpoch=s.remeshCount;
    selected=find(epochs==currentEpoch & times>=s.scale.canonicalTime-.35);
    rawDt=NaN;if isfield(h.common,'acceptedCanonicalDt'),rawDt=h.common.acceptedCanonicalDt(end);end
    steps=[];if isfield(h.common,'acceptedStep'),steps=h.common.acceptedStep(:);end
    countLower=NaN;countUpper=NaN;firstStep=NaN;lastSteps=[];
    if ~isempty(steps) && ~isempty(selected) && all(isfinite(steps(selected)))
        firstStep=steps(selected(1));countLower=s.step-firstStep+1;lastSteps=steps(max(1,end-11):end)';
        before=selected(1)-1;if before>=1 && isfinite(steps(before)),countUpper=s.step-steps(before);else,countUpper=s.step+1;end
    end
    stored=0;storedSpan=NaN;storedFirst=NaN;
    if isfield(s.runMetadata,'autonomousMesh')
        w=s.runMetadata.autonomousMesh.window;stored=numel(w.time);storedSpan=w.time(end)-w.time(1);storedFirst=w.step(1);
    end
    row=struct('file',nativeFiles{k},'strictNativeValidationPassed',true,'step',s.step, ...
        'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime, ...
        'nativeWindowUnit','native_canonical','currentEpoch',currentEpoch,'historyRecords',numel(times), ...
        'sampledCurrentWindowRecords',numel(selected),'firstSampledWindowStep',firstStep, ...
        'everyAcceptedStepWindowLowerBound',countLower,'everyAcceptedStepWindowUpperBound',countUpper, ...
        'lastRecordedSteps',lastSteps,'lastAcceptedCanonicalDt',rawDt, ...
        'steadyDtEstimateFullPoint35Window',floor(.35/rawDt)+1,'minDt',s.config.time.minDt, ...
        'formalMinDtPoint35Bound',floor(.35/s.config.time.minDt)+2,'maxSteps',s.config.time.maxSteps, ...
        'storedControllerWindowRecords',stored,'storedControllerWindowSpan',storedSpan,'storedControllerFirstStep',storedFirst);
    report.actualNative=append_row(report.actualNative,row);clear checkpoint s h;
end
policy=ipm.config.autonomousMeshPolicy(struct());
for m=report.sampleCounts
    t=linspace(0,.35,m)';dt=.35/(m-1);core=[40*exp(-.1*t),44*exp(-.2*t)];
    w=struct('time',t,'step',(0:m-1)','remeshCount',zeros(m,1),'coreCells',core,'safety',.3*ones(m,1));
    memory=struct('version',1,'policy',policy,'transactions',struct([]), ...
        'cumulativeAbsolutePeakJump',0,'window',w,'initialization',struct(),'lastDecision',struct(),'lastFailure',struct());
    state=struct('config',struct('remesh',struct('autonomousMesh',policy,'maxRemeshes',300)), ...
        'step',m,'scale',struct('canonicalTime',.35+dt,'physicalTime',.35+dt), ...
        'ops',struct('remeshCount',0),'runMetadata',struct('autonomousMesh',memory), ...
        'flow',struct('coreGridPoints',40*exp(-.1*(.35+dt)), ...
        'verticalCoreGridPoints',44*exp(-.2*(.35+dt)),'safetyFactor',.3));
    repetitions=max(3,min(40,floor(3e5/m)));
    telemetrySeconds=median_time(@()ipm.remesh.controllerTelemetry(memory,state,policy),repetitions);
    plannerSeconds=median_time(@()ipm.evolve.planAutonomousMesh(state,false),repetitions);
    [planned,plan]=ipm.evolve.planAutonomousMesh(state,false);assert(~plan.requested);
    payload=mock_payload(memory);signature=ipm.output.checkpointSignature(payload);
    signatureSeconds=median_time(@()ipm.output.checkpointSignature(payload),repetitions);
    other=signature;names=fieldnames(w);
    for j=1:numel(names)
        value=other.runMetadata.autonomousMesh.window.(names{j});
        other.runMetadata.autonomousMesh.window.(names{j})=value+zeros(size(value));
    end
    compareSeconds=median_time(@()isequaln(signature,other),repetitions);
    windowInfo=whos('w');payloadInfo=whos('payload');signatureInfo=whos('signature');
    artifact=struct('kind','synthetic_payload_storage_only','payload',payload,'signature',signature);
    file=fullfile(outDir,sprintf('synthetic_storage_%d.mat',m));timer=tic;save(file,'artifact','-v7.3');saveSeconds=toc(timer);fileInfo=dir(file);
    row=struct('records',m,'repetitions',repetitions,'windowScalarDoubles',6*m, ...
        'numericWindowBytes',48*m,'windowWhosBytes',windowInfo.bytes, ...
        'payloadWhosBytes',payloadInfo.bytes,'signatureWhosBytes',signatureInfo.bytes, ...
        'telemetryMedianSeconds',telemetrySeconds,'plannerIncludingTelemetryMedianSeconds',plannerSeconds, ...
        'signatureConstructionMedianSeconds',signatureSeconds,'signatureIsequalnMedianSeconds',compareSeconds, ...
        'syntheticSaveSeconds',saveSeconds,'syntheticCompressedFileBytes',fileInfo.bytes, ...
        'resultingRecords',numel(planned.runMetadata.autonomousMesh.window.time),'decayEstimate',plan.decayEstimate);
    report.microbenchmarks=append_row(report.microbenchmarks,row);
    fprintf('CONTROLLER_COST m=%d telemetry=%.6g plan=%.6g save=%.6g bytes=%d\n',m,telemetrySeconds,plannerSeconds,saveSeconds,fileInfo.bytes);
end
% The OLS slope is a positive weighted average of consecutive secants.
% This numeric check is supporting evidence, not a substitute for the proof.
rng(20260909);maximumOvershoot=0;
for k=1:200
    t=cumsum(.001+rand(37,1));v=randn(37,2);
    for axis=1:2
        fit=polyfit(t-t(end),v(:,axis),1);secants=diff(v(:,axis))./diff(t);
        maximumOvershoot=max([maximumOvershoot,fit(1)-max(secants),min(secants)-fit(1)]);
    end
end
report.olsSecantConvexHull=struct('testedFits',400,'maximumOvershoot',maximumOvershoot,'passed',maximumOvershoot<=1e-12);
report.allPassed=all([report.actualNative.strictNativeValidationPassed])&&report.olsSecantConvexHull.passed;
write_json(fullfile(outDir,'report.json'),report);save(fullfile(outDir,'report.mat'),'report');assert(report.allPassed);
fprintf('CONTROLLER_COMPLEXITY_COMPLETE noLU=1 cases=%d benchmarks=%d\n',numel(nativeFiles),numel(report.sampleCounts));
end
function seconds=median_time(f,n)
f();values=zeros(n,1);for k=1:n,timer=tic;f();values(k)=toc(timer);end;seconds=median(values);
end
function payload=mock_payload(memory)
s=struct('config',struct(),'runMetadata',struct('autonomousMesh',memory),'rho',zeros(9), ...
    'x',-4:4,'y',(0:8)','baseX',-4:4,'baseY',(0:8)','rescaling',struct(), ...
    'remeshCount',0,'scale',struct(),'normalizedTime',0,'step',0,'mass0',0,'rhoRange0',[0,0]);
log=struct('history',struct(),'snapshotRho',{{}},'snapshotNormalizedTime',[]);
payload=struct('state',s,'log',log,'cursor',struct());
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
