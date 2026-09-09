function report=screen_profile_mode(files,outputDirectory)
%SCREEN_PROFILE_MODE Strict-native inner-shape forecasts with held-out times.
% The fixed experiment uses the first four of six frames for training.
assert(iscell(files) && numel(files)==6 && ~isfolder(outputDirectory));
assert(maxNumCompThreads==10);
mkdir(outputDirectory);
registration=struct('files',{files},'trainingCount',4,'observer','linear C1 tensor Hermite', ...
    'window',[-1,1,0,2],'nodeCounts',[49,33;97,65],'rateBounds',[.02,4], ...
    'normalization','Original recorded quadratic P and strict widths; no fitted amplitude adjustment.', ...
    'pdeSteps',0,'nativeCheckpointWritten',false);
save(fullfile(outputDirectory,'registration.mat'),'registration');
data=cell(1,2);axes=cell(1,2);weights=cell(1,2);
for level=1:2
    xi=linspace(-1,1,registration.nodeCounts(level,1));
    eta=linspace(0,2,registration.nodeCounts(level,2));
    [XI,ETA]=meshgrid(xi,eta);axes{level}=struct('XI',XI,'ETA',ETA);
    weights{level}=trap_weights(eta)*trap_weights(xi)';
    data{level}=zeros(numel(files),numel(XI));
end
records=cell(1,numel(files));
for k=1:numel(files)
    cp=ipm.output.readCheckpoint(files{k});s=cp.payload.state;h=cp.payload.log.history;
    assert(all(ipm.output.trustedMask(h,s.config)) && ...
        strcmp(s.config.transport.spatialDiscretization,'high_order') && ...
        strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic'));
    Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);omega=s.rho*Dx';
    a=h.gauge.omegaGaugeQuadraticPeakX(end);P=h.gauge.omegaGaugeQuadraticPeakValue(end);
    wx=h.mesh.trackedWallCoreWidth(end);wy=h.mesh.trackedVerticalCoreWidth(end);
    assert(all(isfinite([a,P,wx,wy])) && all([P,wx,wy]>0));
    for level=1:2
        points=axes{level};
        v=ipm_accellab_tensor_hermite(omega,s.x,s.y,Dx,Dy, ...
            a+wx*points.XI,wy*points.ETA);
        data{level}(k,:)=reshape(v.value/P,1,[]);
    end
    records{k}=struct('file',files{k},'tau',s.scale.canonicalTime, ...
        'physicalTime',s.scale.physicalTime,'remeshCount',s.remeshCount, ...
        'a',a,'P',P,'widths',[wx,wy],'caseId',char(s.runMetadata.caseId));
    clear cp s omega h
end
records=[records{:}];assert(isscalar(unique({records.caseId})));
fits=cell(1,2);summaries=cell(1,2);
for level=1:2
    fits{level}=fit_profile_mode([records.tau],data{level},weights{level}(:),4);
    summaries{level}=rmfield(fits{level},'candidateLimit');
end
report=struct('registration',registration,'records',records,'fits',{summaries}, ...
    'rateRelativeObservationSensitivity',abs(fits{2}.rate-fits{1}.rate)/fits{2}.rate, ...
    'interpretation',['Forecasts are tested against subsequent actual native frames. ' ...
      'The inferred limit concerns shape after translation, width and amplitude normalization. ' ...
      'It does not establish original-gauge convergence, a singular time, or a limit solution.']);
save(fullfile(outputDirectory,'data_and_fits.mat'),'report','data','axes','weights','fits','-v7.3');
fid=fopen(fullfile(outputDirectory,'report.json'),'w');assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('PROFILE_MODE %s\n',jsonencode(report));
end

function w=trap_weights(x)
d=diff(x(:));w=[d(1)/2;(d(1:end-1)+d(2:end))/2;d(end)/2];
end
