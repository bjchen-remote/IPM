function report=screen_peak_template_history(checkpointFile,outputFile)
%SCREEN_PEAK_TEMPLATE_HISTORY Read-only native sampled peak-drift audit.
% A changed recorded center detects at least one selector event in an
% interval. Unchanged endpoints do not exclude intermediate switches.
assert(~isfile(outputFile),'Refusing to overwrite a diagnostic.');
assert(maxNumCompThreads==10,'Strict historical signature requires 10 threads.');
cp=ipm.output.readCheckpoint(checkpointFile);
h=cp.payload.log.history;
assert(all(ipm.output.trustedMask(h,cp.payload.state.config)));
g=h.gauge;c=h.common;m=h.mesh;
P=g.omegaGaugeQuadraticPeakValue(:);
center=g.omegaGaugeQuadraticStencilCenterX(:);
tau=c.canonicalTau(:);step=c.acceptedStep(:);
remesh=m.remeshCount(:);
valid=isfinite(P) & P>0 & isfinite(center) & isfinite(tau);
pair=valid(1:end-1) & valid(2:end);
dP=diff(P);relative=dP./P(1:end-1);
meshChanged=diff(remesh)~=0;
centerChanged=diff(center)~=0;
groups={'remeshRecorded','sameMeshCenterChanged','sameMeshCenterUnchanged'};
masks={pair & meshChanged,pair & ~meshChanged & centerChanged, ...
    pair & ~meshChanged & ~centerChanged};
summary=struct;
for k=1:numel(groups)
    values=relative(masks{k});
    summary.(groups{k})=struct('intervalCount',numel(values), ...
        'signedRelativeJumpSum',sum(values), ...
        'absoluteRelativeJumpSum',sum(abs(values)), ...
        'maxAbsoluteRelativeJump',max([0;abs(values)]));
end
rate=g.omegaGaugeQuadraticForcing(:)+c.c_omega(:).*P;
denominator=abs(g.omegaGaugeQuadraticForcing(:))+abs(c.c_omega(:).*P);
rateValid=valid & isfinite(rate) & isfinite(denominator) & denominator>0;
indices=find(pair);[~,order]=sort(abs(relative(indices)),'descend');
indices=indices(order(1:min(20,numel(order))));
events=cell(1,numel(indices));
for k=1:numel(indices)
    j=indices(k);
    events{k}=struct('rowInterval',[j,j+1],'stepInterval',step(j:j+1)', ...
        'tauInterval',tau(j:j+1)','P',P(j:j+1)', ...
        'recordedCenterX',center(j:j+1)','remeshCount',remesh(j:j+1)', ...
        'relativePeakChange',relative(j),'remeshRecorded',meshChanged(j), ...
        'centerChangeRecorded',centerChanged(j));
end
report=struct('checkpointFile',checkpointFile,'nativeSignatureValidated',true, ...
    'fullHistoryTrusted',true,'sampleCount',numel(P),'groups',summary, ...
    'largestIntervals',[events{:}], ...
    'maximumInstantaneousRateResidual',max(abs(rate(rateValid))), ...
    'maximumNormalizedInstantaneousRateResidual',max(abs(rate(rateValid))./denominator(rateValid)), ...
    'totalRelativePeakDrift',(P(end)-P(find(valid,1)))/P(find(valid,1)), ...
    'interpretation',['Sampled drift accounting, not event-resolved attribution. ' ...
      'An instantaneous zero rate does not imply conservation across selector changes. ' ...
      'Unchanged sampled centers cannot rule out unrecorded intermediate switches.']);
fid=fopen(outputFile,'w');assert(fid>=0);cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('PEAK_TEMPLATE_HISTORY %s\n',jsonencode(rmfield(report,'largestIntervals')));
end
