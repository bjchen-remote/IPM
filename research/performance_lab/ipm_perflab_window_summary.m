function s=ipm_perflab_window_summary(m)
% Research diagnostics. Neither fit is the original exact-window raw fit.
b=m.buckets;n=m.bucketCount;
last=m.raw(m.rawCount,:);cutoff=last(1)-m.windowWidth;
rows=[m.raw(1:m.rawCount,:);b.last(1:n,:)];
[~,keep]=unique(rows(:,2),'sorted');rows=rows(keep,:);
[~,order]=sort(rows(:,2));rows=rows(order,:);
[bound,whichBucket]=max(b.maximumDecay(1:n,:),[],1);
from=zeros(2,6);to=zeros(2,6);
has=bound>0;
fromNames={'fromX','fromY'};toNames={'toX','toY'};
for axis=1:2
    if has(axis)
        from(axis,:)=b.(fromNames{axis})(whichBucket(axis),:);
        to(axis,:)=b.(toNames{axis})(whichBucket(axis),:);
    else
        from(axis,:)=last;to(axis,:)=last;
    end
end
fitRepresentative=[NaN NaN];
if size(rows,1)>=2
    for axis=1:2
        p=polyfit(rows(:,1)-last(1),log(rows(:,3+axis)),1);
        fitRepresentative(axis)=-p(1);
    end
end
% Chan/Welford merge in fixed chronological bucket order. Full buckets are
% retained: their oldest samples may precede the exact moving-window cutoff.
count=0;meanTime=0;meanLog=[0 0];sumTime2=0;sumTimeLog=[0 0];
for k=1:n
    combined=count+b.count(k);
    dt=b.meanTime(k)-meanTime;dy=b.meanLog(k,:)-meanLog;
    weight=count*b.count(k)/combined;
    sumTime2=sumTime2+b.sumTime2(k)+dt*dt*weight;
    sumTimeLog=sumTimeLog+b.sumTimeLog(k,:)+dt*dy*weight;
    meanTime=meanTime+dt*(b.count(k)/combined);
    meanLog=meanLog+dy*(b.count(k)/combined);
    count=combined;
end
fitBuckets=[NaN NaN];
if sumTime2>0
    fitBuckets=-sumTimeLog/sumTime2;
end
% A research-only decision projection makes serialized decision replay
% checkable. It does not call, replace, or qualify the native controller.
trendAvailable=count>=4 && last(1)>b.firstTime(1);
decay=[0 0];
if trendAvailable,decay=bound;end
predicted=last(4:5).*exp(-.2*decay);
reason='none';
if last(6)>.70
    reason='safety_stop';
elseif any(last(4:5)<26)
    reason='current_core_trigger';
elseif trendAvailable && any(predicted<22 & last(4:5)<32)
    reason='forecast_core_trigger';
end
s=struct('lastObservation',last,'samples',rows,'sampleCount',size(rows,1), ...
    'rawCount',m.rawCount,'bucketCount',n,'expandedObservationCount',count, ...
    'exactWindowCutoff',cutoff,'effectiveBucketStart',b.firstTime(1), ...
    'bucketTimeOverhang',max(0,cutoff-b.firstTime(1)), ...
    'sampleCoverage',last(1)-rows(1,1),'bucketCoverage',last(1)-b.firstTime(1), ...
    'hasFullWindowCoverage',b.firstTime(1)<=cutoff, ...
    'maximumDecay',bound,'hasDecayEvidence',has, ...
    'evidenceFrom',from,'evidenceTo',to, ...
    'evidenceOverhang',max(0,cutoff-from(:,1))', ...
    'evidenceRightEndpointExpired',(to(:,1)<cutoff)', ...
    'evidenceCrossesCutoff',(from(:,1)<cutoff & to(:,1)>=cutoff)', ...
    'fitRepresentative',fitRepresentative,'fitWholeBuckets',fitBuckets, ...
    'trendAvailable',trendAvailable,'predictedCoreCells',predicted, ...
    'researchDecision',reason);
end
