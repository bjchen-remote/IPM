function memory = ipm_perflab_window_update(memory,row)
% Research only: bounded real observations and conservative secant evidence.
% row = [canonicalTime acceptedStep remeshEpoch coreX coreY safety].
% Ordinary observations must be consecutive accepted steps. A remesh reset
% must follow its old-epoch source observation at the SAME time and step.
validateattributes(row,{'double'},{'real','finite','size',[1 6]});
assert(row(1)>=0 && all(row(2:3)>=0) && all(row(2:3)<flintmax) && ...
    all(row(2:3)==fix(row(2:3))) && all(row(4:5)>0) && row(6)>=0, ...
    'ipm:ResearchWindowObservation','Invalid observation.');
if isempty(memory)
    assert(isequal(row(1:3),[0 0 0]),'ipm:ResearchWindowInitial', ...
        'A new research stream starts at time, step and epoch zero.');
    memory=initialize(row);
else
    assert(isstruct(memory) && isscalar(memory) && ...
        isequal(memory.kind,'research_bounded_window_v1') && ...
        isequal(memory.version,1),'ipm:ResearchWindowVersion','Wrong state.');
    previous=memory.raw(memory.rawCount,:);
    if row(3)~=previous(3)
        assert(row(3)==previous(3)+1 && isequal(row(1:2),previous(1:2)), ...
            'ipm:ResearchWindowEpoch','Reset requires the next epoch at its actual source time and step.');
        memory=initialize(row);
    elseif row(2)==previous(2)
        assert(isequal(row,previous),'ipm:ResearchWindowDuplicate', ...
            'A repeated accepted step must be an exact duplicate.');
        return
    else
        assert(row(2)==previous(2)+1 && row(1)>previous(1), ...
            'ipm:ResearchWindowOrder','Accepted observations must be consecutive and time increasing.');
    end
end
previous=[];
if memory.rawCount>0
    previous=memory.raw(memory.rawCount,:);
end
bucketIndex=floor((row(1)-memory.origin)/memory.bucketWidth);
assert(isfinite(bucketIndex) && bucketIndex>=0 && bucketIndex<2^48, ...
    'ipm:ResearchWindowResolution','Bucket index exceeds the registered arithmetic range.');
cutoff=row(1)-memory.windowWidth;
b=memory.buckets;
n=memory.bucketCount;
% A whole bucket expires only after its last real observation has expired.
expired=find(b.last(1:n,1)<cutoff,1,'last');
if ~isempty(expired)
    names=fieldnames(b);
    for k=1:numel(names)
        a=b.(names{k});
        a(1:n-expired,:)=a(expired+1:n,:);
        a(n-expired+1:n,:)=0;
        b.(names{k})=a;
    end
    n=n-expired;
end
if n==0 || b.index(n)~=bucketIndex
    assert(n==0 || bucketIndex>b.index(n),'ipm:ResearchWindowOrder','Bucket order changed.');
    assert(n<memory.bucketCapacity,'ipm:ResearchWindowCapacity','Bucket budget exhausted.');
    n=n+1;
    b.index(n)=bucketIndex;
    b.firstTime(n)=row(1);
    b.firstStep(n)=row(2);
end
b.last(n,:)=row;
% Stable streaming sufficient statistics, with time relative to epoch origin.
offset=row(1)-memory.origin;
logs=log(row(4:5));
oldCount=b.count(n);
newCount=oldCount+1;
assert(newCount<flintmax,'ipm:ResearchWindowCount','Bucket count lost integer resolution.');
deltaTime=offset-b.meanTime(n);
deltaLog=logs-b.meanLog(n,:);
b.meanTime(n)=b.meanTime(n)+deltaTime/newCount;
b.meanLog(n,:)=b.meanLog(n,:)+deltaLog/newCount;
b.sumTime2(n)=b.sumTime2(n)+deltaTime*(offset-b.meanTime(n));
b.sumTimeLog(n,:)=b.sumTimeLog(n,:)+deltaTime*(logs-b.meanLog(n,:));
b.count(n)=newCount;
if ~isempty(previous)
    decay=-(logs-log(previous(4:5)))/(row(1)-previous(1));
    assert(all(isfinite(decay)),'ipm:ResearchWindowRate','Nonfinite secant is never silently discarded.');
    fromNames={'fromX','fromY'}; toNames={'toX','toY'};
    for axis=1:2
        if decay(axis)>b.maximumDecay(n,axis)
            b.maximumDecay(n,axis)=decay(axis);
            b.hasEvidence(n,axis)=true;
            b.(fromNames{axis})(n,:)=previous;
            b.(toNames{axis})(n,:)=row;
        end
    end
end
assert(all(isfinite([b.meanTime(n),b.meanLog(n,:),b.sumTime2(n),b.sumTimeLog(n,:)])), ...
    'ipm:ResearchWindowStatistics','Nonfinite sufficient statistics.');
memory.buckets=b;
memory.bucketCount=n;
raw=memory.raw(1:memory.rawCount,:);
raw=[raw;row];
raw=raw(max(1,size(raw,1)-memory.rawCapacity+1):end,:);
raw=raw(raw(:,1)>=cutoff,:);
memory.raw(:)=0;
memory.raw(1:size(raw,1),:)=raw;
memory.rawCount=size(raw,1);
end

function m=initialize(row)
m=struct('kind','research_bounded_window_v1','version',1, ...
    'windowWidth',.35,'rawCapacity',32,'bucketCapacity',480, ...
    'bucketWidth',.35/478,'origin',row(1),'epoch',row(3), ...
    'raw',zeros(32,6),'rawCount',0,'bucketCount',0);
n=480;
b=struct('index',zeros(n,1),'count',zeros(n,1), ...
    'firstTime',zeros(n,1),'firstStep',zeros(n,1),'last',zeros(n,6), ...
    'meanTime',zeros(n,1),'meanLog',zeros(n,2),'sumTime2',zeros(n,1), ...
    'sumTimeLog',zeros(n,2),'maximumDecay',zeros(n,2), ...
    'hasEvidence',false(n,2),'fromX',zeros(n,6),'toX',zeros(n,6), ...
    'fromY',zeros(n,6),'toY',zeros(n,6));
m.buckets=b;
end
