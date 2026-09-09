function report=ipm_perflab_test_bounded_window(outDir)
% Synthetic, pure MATLAB research regression. No native checkpoint or PDE.
assert(~isfolder(outDir),'ipm:ResearchWindowOutput','Use a new output directory.');
mkdir(outDir);mkdir(fullfile(outDir,'checkpoints'));
cases=registered_cases();
registration=struct('kind','research_bounded_window_tests_v1', ...
    'windowWidth',.35,'rawCapacity',32,'bucketCapacity',480, ...
    'bucketWidth',.35/478,'maximumRepresentativeRows',512, ...
    'maximumEvidenceEndpointRows',1920, ...
    'epochRule','next_epoch_at_same_step_and_time_after_source_observation', ...
    'secantOracle','both_true_adjacent_endpoints_in_current_exact_window', ...
    'expiration','drop_whole_bucket_only_if_last_time_strictly_before_cutoff', ...
    'fitRole','diagnostic_only_never_reduces_or_replaces_secant_bound', ...
    'caseNames',{{cases.name}},'caseObservationCounts',arrayfun(@(c)size(c.rows,1),cases), ...
    'nativeCheckpoint',false,'PDE',false,'LU',false);
save(fullfile(outDir,'registration_and_inputs.mat'),'registration','cases','-v7');
write_json(fullfile(outDir,'registration.json'),registration);
records=struct([]);fitRows=struct([]);terminalWindows=cell(1,numel(cases));
for c=1:numel(cases)
    rows=cases(c).rows;N=size(rows,1);
    oracle=secant_oracle(rows,.35);
    checkpoints=unique([max(1,round(linspace(3,N-1,6))), ...
        find(diff(rows(:,3))~=0)',find(diff(rows(:,3))~=0)'+1]);
    queries=unique([1 min(4,N) round(linspace(min(4,N),N,65)) checkpoints ...
        cases(c).extraQueries]);
    memory=[];replay=[];rawMax=0;bucketMax=0;representativeMax=0;
    maximumExcess=[0 0];checkpointBytes=0;maximumBytes=0;
    updateSeconds=0;replaySeconds=0;saveSeconds=0;loadSeconds=0;
    serializedSplits=0;replayedUpdates=0;summaryComparisons=0;
    tCase=tic;
    for k=1:N
        tick=tic;memory=ipm_perflab_window_update(memory,rows(k,:));
        updateSeconds=updateSeconds+toc(tick);
        if ~isempty(replay)
            tick=tic;replay=ipm_perflab_window_update(replay,rows(k,:));
            assert(isequaln(memory,replay),'ipm:ResearchWindowReplay','Restart state differs.');
            replaySeconds=replaySeconds+toc(tick);replayedUpdates=replayedUpdates+1;
        end
        assert(isequal(memory.raw(memory.rawCount,:),rows(k,:)), ...
            'ipm:ResearchWindowLast','Last real observation changed.');
        assert(memory.rawCount<=32 && memory.bucketCount<=480 && ...
            memory.rawCount+memory.bucketCount<=512,'ipm:ResearchWindowBudget','Unbounded record pool.');
        bound=max(memory.buckets.maximumDecay(1:memory.bucketCount,:),[],1);
        assert(all(bound>=oracle(k,:)), ...
            'ipm:ResearchWindowEnvelope','An in-window true adjacent secant was lost.');
        maximumExcess=max(maximumExcess,bound-oracle(k,:));
        rawMax=max(rawMax,memory.rawCount);bucketMax=max(bucketMax,memory.bucketCount);
        % The raw steps are consecutive: overlap with bucket representatives
        % is exactly the bucket-last steps at or after the first raw step.
        exactRepresentativeCount=memory.rawCount+memory.bucketCount- ...
            nnz(memory.buckets.last(1:memory.bucketCount,2)>=memory.raw(1,2));
        representativeMax=max(representativeMax,exactRepresentativeCount);
        assert(exactRepresentativeCount<=512,'ipm:ResearchWindowBudget','Representative budget exceeded.');
        if k>1 && rows(k,3)~=rows(k-1,3)
            assert(memory.rawCount==1 && memory.bucketCount==1 && ...
                all(bound==0),'ipm:ResearchWindowEpoch','Remesh retained cross-epoch evidence.');
        end
        if ismember(k,queries)
            ipm_perflab_window_validate(memory);
            s=ipm_perflab_window_summary(memory);
            assert(s.sampleCount==exactRepresentativeCount, ...
                'ipm:ResearchWindowCount','Representative count shortcut disagrees with actual union.');
            info=whos('memory');maximumBytes=max(maximumBytes,info.bytes);
            assert(s.sampleCount<=512 && all(ismember(s.samples,rows(1:k,:),'rows')), ...
                'ipm:ResearchWindowRealSamples','A representative is not a real observation.');
            for axis=1:2
                if s.hasDecayEvidence(axis)
                    assert(ismember(s.evidenceFrom(axis,:),rows(1:k,:),'rows') && ...
                        ismember(s.evidenceTo(axis,:),rows(1:k,:),'rows'), ...
                        'ipm:ResearchWindowRealEvidence','A secant endpoint is not a real observation.');
                end
            end
            if ~isempty(replay)
                assert(isequaln(s,ipm_perflab_window_summary(replay)), ...
                    'ipm:ResearchWindowReplayDecision','Serialized decision or fit differs.');
                summaryComparisons=summaryComparisons+1;
            end
            selected=rows(1:k,3)==rows(k,3) & rows(1:k,1)>=s.exactWindowCutoff;
            expanded=rows(1:k,3)==rows(k,3) & rows(1:k,1)>=s.effectiveBucketStart;
            exactFit=raw_fit(rows(selected,:));expandedFit=raw_fit(rows(expanded,:));
            assert(sum(expanded)==s.expandedObservationCount, ...
                'ipm:ResearchWindowDoubleCount','Bucket OLS statistics double-count observations.');
            f=struct('caseName',cases(c).name,'observationIndex',k,'time',rows(k,1), ...
                'step',rows(k,2),'epoch',rows(k,3),'rawWindowCount',sum(selected), ...
                'expandedCount',sum(expanded),'expiredRowsRetained',sum(expanded & ~selected), ...
                'representativeCount',s.sampleCount,'bucketCount',s.bucketCount, ...
                'exactWindowRawFit',exactFit,'representativeFit',s.fitRepresentative, ...
                'wholeBucketFit',s.fitWholeBuckets,'expandedRawFit',expandedFit, ...
                'representativeBias',s.fitRepresentative-exactFit, ...
                'wholeBucketBias',s.fitWholeBuckets-exactFit, ...
                'sufficientStatisticRoundoff',s.fitWholeBuckets-expandedFit, ...
                'oracleMaximumSecant',oracle(k,:),'bound',s.maximumDecay, ...
                'boundExcess',s.maximumDecay-oracle(k,:), ...
                'bucketTimeOverhang',s.bucketTimeOverhang,'evidenceOverhang',s.evidenceOverhang, ...
                'evidenceRightEndpointExpired',s.evidenceRightEndpointExpired, ...
                'evidenceCrossesCutoff',s.evidenceCrossesCutoff, ...
                'sampleCoverage',s.sampleCoverage,'bucketCoverage',s.bucketCoverage, ...
                'researchDecision',s.researchDecision);
            if isempty(fitRows),fitRows=f;else,fitRows(end+1)=f;end %#ok<AGROW>
        end
        if ismember(k,checkpoints)
            researchCheckpoint=struct('kind','research_only_not_native_ipm_checkpoint', ...
                'registration',registration,'window',memory, ...
                'summary',ipm_perflab_window_summary(memory));
            filename=fullfile(outDir,'checkpoints',sprintf('%s_at_%06d.mat',cases(c).name,k));
            tick=tic;save(filename,'researchCheckpoint','-v7');saveSeconds=saveSeconds+toc(tick);
            tick=tic;loaded=load(filename,'researchCheckpoint');loadSeconds=loadSeconds+toc(tick);
            ipm_perflab_window_validate(loaded.researchCheckpoint.window);
            assert(isequaln(loaded.researchCheckpoint,researchCheckpoint), ...
                'ipm:ResearchWindowSerialization','MAT roundtrip changed data.');
            replay=loaded.researchCheckpoint.window;
            assert(isequaln(ipm_perflab_window_update(replay,rows(k,:)),replay), ...
                'ipm:ResearchWindowIdempotence','Duplicate input is not idempotent.');
            entry=dir(filename);checkpointBytes=max(checkpointBytes,entry.bytes);
            serializedSplits=serializedSplits+1;
        end
    end
    terminalWindows{c}=memory;
    r=struct('caseName',cases(c).name,'observations',N,'envelopeChecks',N, ...
        'maximumRawRows',rawMax,'maximumBucketRecords',bucketMax, ...
        'maximumRepresentativeRowsAllUpdates',representativeMax, ...
        'maximumStateBytesAtQueries',maximumBytes,'maximumCheckpointBytes',checkpointBytes, ...
        'maximumBoundExcess',maximumExcess,'serializedSplits',serializedSplits, ...
        'bitwiseReplayedUpdates',replayedUpdates,'bitwiseSummaryComparisons',summaryComparisons, ...
        'updateSeconds',updateSeconds,'updateMicrosecondsPerObservation',1e6*updateSeconds/N, ...
        'replayAndEqualitySeconds',replaySeconds,'saveSeconds',saveSeconds,'loadSeconds',loadSeconds, ...
        'totalTestSeconds',toc(tCase),'passed',true);
    if isempty(records),records=r;else,records(end+1)=r;end %#ok<AGROW>
    fprintf('%s: %d observations, raw %d / buckets %d, %.1f us/update, all envelopes/replay passed\n', ...
        r.caseName,N,rawMax,bucketMax,r.updateMicrosecondsPerObservation);
end
negative=negative_cases();
% Explicit boundary witnesses must be observed, not merely listed as cases.
overhang=max(vertcat(fitRows.evidenceOverhang),[],'all');
assert(overhang>.35/478,'ipm:ResearchWindowWitness','Missing long-crossing-secant witness.');
assert(any([fitRows.expiredRowsRetained]>0),'ipm:ResearchWindowWitness','Missing partial expired bucket witness.');
assert(any(vertcat(fitRows.evidenceRightEndpointExpired),'all'), ...
    'ipm:ResearchWindowWitness','Missing retained expired maximum witness.');
summaryTable=fit_summary(cases,fitRows);
costState=terminalWindows{1};summaryTimes=zeros(30,1);validateTimes=zeros(5,1);
for k=1:numel(summaryTimes)
    tick=tic;ipm_perflab_window_summary(costState);summaryTimes(k)=toc(tick);
end
for k=1:numel(validateTimes)
    tick=tic;ipm_perflab_window_validate(costState);validateTimes(k)=toc(tick);
end
boundedCost=struct('activeBuckets',costState.bucketCount,'summaryRepeats',numel(summaryTimes), ...
    'summaryMedianSeconds',median(summaryTimes),'validationRepeats',numel(validateTimes), ...
    'validationMedianSeconds',median(validateTimes),'summaryRawSeconds',summaryTimes, ...
    'validationRawSeconds',validateTimes);
report=struct('kind','research_bounded_window_report_v1','passed',true, ...
    'registration',registration,'cases',records,'negativeTests',negative, ...
    'fitSummary',summaryTable,'fitQueryCount',numel(fitRows), ...
    'maximumEvidenceOverhang',overhang,'boundedCost',boundedCost,'nativeOrPDEQualification',false, ...
    'limitations',{{'Fits differ from original exact-window raw polyfit.', ...
    'Whole-bucket fitting includes expired observations without a directional bias guarantee.', ...
    'A retained secant can begin arbitrarily earlier than the cutoff after a large time gap.', ...
    'Fixed controller state does not bound the full solver history or transaction ledger.', ...
    'Finite stored aggregate statistics cannot authenticate discarded input without external provenance.'}});
save(fullfile(outDir,'report.mat'),'report','fitRows','terminalWindows','-v7');
write_json(fullfile(outDir,'report.json'),report);
write_json(fullfile(outDir,'fit_rows.json'),fitRows);
end

function cases=registered_cases()
cases=struct('name',{},'rows',{},'extraQueries',{});
n=25001;t=(0:n-1)'*2e-5;
cases(end+1)=make_case('constant_25001',t,zeros(n,1),[40+0*t 44+0*t],[]);
n=100001;dt=[ones(1000,1)*5e-4;logspace(-5,-8,n-1001)'];t=[0;cumsum(dt)];
cases(end+1)=make_case('exponential_nonuniform_100001',t,zeros(n,1), ...
    [40*exp(-.6*t) 44*exp(-.35*t)],[]);
n=10001;t=(0:n-1)'*1e-4;jump=zeros(n,1);jump(5001:5004)=-.025;
cases(end+1)=make_case('diagnostic_template_jump',t,zeros(n,1), ...
    [40*exp(-.4*t+jump) 44*exp(-.25*t+.6*jump)],4999:5010);
cases(end).extraQueries=[cases(end).extraQueries 8499:8515];
n=10001;dt=1e-5*(1+.8*sin((1:n-1)'*.019));dt(2000)=.6;t=[0;cumsum(dt)];
core=[40*exp(-.2*t) 44*exp(-.1*t)];core(2001:end,:)=core(2001:end,:)*exp(-.1);
cases(end+1)=make_case('large_gap_then_tiny_dt',t,zeros(n,1),core,1999:2010);
n=4001;t=(0:n-1)'*.0001;core=[40-5*t 44-3*t];
c=make_case('same_step_remesh_epoch_reset',t,zeros(n,1),core,[2000 2001 2002]);
before=c.rows(1:2001,:);after=c.rows(2001:end,:);after(:,3)=1;after(:,4:5)=after(:,4:5)+[8 9];
c.rows=[before;after];cases(end+1)=c;
h=.35/478;
t=unique(sort([0;h;2*h;.35-eps(.35);.35;.35+eps(.35); ...
    .35+h-eps(.35+h);.35+h;.35+h+eps(.35+h);.7;.7+eps(.7)]));
cases(end+1)=make_case('exact_and_ulp_boundaries',t,zeros(size(t)), ...
    [40+0*t 44+0*t],1:numel(t));
% Spike ends just before the cutoff while its bucket still has a retained row.
t=[0;.1;.10001;.10002;.1005;.45002;.45003;.4505;.4508];
core=[40+0*t 44+0*t];core(3:4,:)=core(3:4,:)*.999;
cases(end+1)=make_case('expired_spike_in_partial_bucket',t,zeros(size(t)),core,1:numel(t));
end

function c=make_case(name,t,epoch,core,extra)
n=numel(t);c=struct('name',name,'rows',[t (0:n-1)' epoch core .3+zeros(n,1)], ...
    'extraQueries',extra);
end

function bound=secant_oracle(rows,W)
% Independent monotone queues retain exact-window true adjacent secants.
N=size(rows,1);bound=zeros(N,2);queue=zeros(N,2);values=zeros(N,2);head=[1 1];tail=[0 0];
for k=2:N
    if rows(k,3)~=rows(k-1,3)
        head=[1 1];tail=[0 0];continue
    end
    d=-(log(rows(k,4:5))-log(rows(k-1,4:5)))/(rows(k,1)-rows(k-1,1));
    for a=1:2
        while tail(a)>=head(a) && values(tail(a),a)<=d(a),tail(a)=tail(a)-1;end
        tail(a)=tail(a)+1;queue(tail(a),a)=k;values(tail(a),a)=d(a);
        while tail(a)>=head(a) && rows(queue(head(a),a)-1,1)<rows(k,1)-W
            head(a)=head(a)+1;
        end
        if tail(a)>=head(a),bound(k,a)=max(0,values(head(a),a));end
    end
end
end

function p=raw_fit(rows)
p=[NaN NaN];if size(rows,1)<2,return,end
for a=1:2
    q=polyfit(rows(:,1)-rows(end,1),log(rows(:,3+a)),1);p(a)=-q(1);
end
end

function summary=fit_summary(cases,rows)
summary=struct([]);
for c=1:numel(cases)
    r=rows(strcmp({rows.caseName},cases(c).name));
    s=struct('caseName',cases(c).name,'queries',numel(r), ...
        'maximumAbsoluteRepresentativeBias',max(abs(vertcat(r.representativeBias)),[],1,'omitnan'), ...
        'maximumAbsoluteWholeBucketBias',max(abs(vertcat(r.wholeBucketBias)),[],1,'omitnan'), ...
        'maximumSufficientStatisticRoundoff',max(abs(vertcat(r.sufficientStatisticRoundoff)),[],1,'omitnan'), ...
        'maximumExpiredRowsRetained',max([r.expiredRowsRetained]), ...
        'maximumBucketTimeOverhang',max([r.bucketTimeOverhang]), ...
        'maximumEvidenceOverhang',max(vertcat(r.evidenceOverhang),[],1));
    if isempty(summary),summary=s;else,summary(end+1)=s;end %#ok<AGROW>
end
end

function results=negative_cases()
row=[0 0 0 40 44 .3];m=ipm_perflab_window_update([],row);
tests=struct('name',{},'callback',{});
tests(end+1)=neg('NaN_core',@()ipm_perflab_window_update(m,[.1 1 0 NaN 44 .3]));
tests(end+1)=neg('zero_core',@()ipm_perflab_window_update(m,[.1 1 0 0 44 .3]));
tests(end+1)=neg('nonfinite_rate',@()ipm_perflab_window_update(m,[realmin 1 0 1e-300 44 .3]));
tests(end+1)=neg('skipped_step',@()ipm_perflab_window_update(m,[.1 2 0 40 44 .3]));
tests(end+1)=neg('same_time_next_step',@()ipm_perflab_window_update(m,[0 1 0 40 44 .3]));
tests(end+1)=neg('changed_duplicate',@()ipm_perflab_window_update(m,[0 0 0 39 44 .3]));
tests(end+1)=neg('skipped_epoch',@()ipm_perflab_window_update(m,[0 0 2 40 44 .3]));
tests(end+1)=neg('epoch_at_wrong_clock',@()ipm_perflab_window_update(m,[.1 1 1 40 44 .3]));
reset=ipm_perflab_window_update(m,[0 0 1 48 53 .3]);
tests(end+1)=neg('backward_epoch',@()ipm_perflab_window_update(reset,row));
tests(end+1)=neg('bucket_precision_limit',@()ipm_perflab_window_update(m,[2^49*(.35/478) 1 0 40 44 .3]));
bad=m;bad.extra=zeros(1000,1);tests(end+1)=neg('hidden_unbounded_payload',@()ipm_perflab_window_validate(bad));
bad2=m;bad2.buckets.meanLog(1,1)=NaN;tests(end+1)=neg('nonfinite_statistic',@()ipm_perflab_window_validate(bad2));
bad3=m;bad3.buckets.count(1)=2;tests(end+1)=neg('aggregate_count_corruption',@()ipm_perflab_window_validate(bad3));
bad4=ipm_perflab_window_update(m,[.1 1 0 39 43 .3]);bad4.buckets.maximumDecay(2,1)=9;
tests(end+1)=neg('secant_evidence_corruption',@()ipm_perflab_window_validate(bad4));
results=struct([]);
for k=1:numel(tests)
    caught=false;identifier='';
    try
        tests(k).callback();
    catch err
        caught=true;identifier=err.identifier;
    end
    assert(caught,'ipm:ResearchWindowNegative','Negative case was accepted: %s',tests(k).name);
    r=struct('name',tests(k).name,'rejected',caught,'identifier',identifier);
    if isempty(results),results=r;else,results(end+1)=r;end %#ok<AGROW>
end
end

function t=neg(name,callback)
t=struct('name',name,'callback',callback);
end

function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
