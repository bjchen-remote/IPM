function report = ipm_perflab_signature_probe(checkpointFile)
%IPM_PERFLAB_SIGNATURE_PROBE Identify exact fields affected by thread count.
%   Read the existing checkpoint only at 10 threads. Recompute its existing
%   signature from unchanged in-memory payload at 10/1/2/4/10 threads. No
%   checkpoint is changed or written, and the 10-thread setting is restored.

maxNumCompThreads(10);
threadCleanup = onCleanup(@()maxNumCompThreads(10));
checkpoint = ipm.output.readCheckpoint(checkpointFile);
reference = checkpoint.signature;
records = struct('threads',{},'exact',{},'differences',{});
for threads = [10,1,2,4,10]
    maxNumCompThreads(threads);
    actual = ipm.output.checkpointSignature(checkpoint.payload);
    differences = compare_values(actual,reference,'signature');
    record = struct('threads',maxNumCompThreads, ...
        'exact',isequaln(actual,reference),'differences',differences);
    records(end+1) = record; %#ok<AGROW>
    fprintf('signature threads=%d exact=%d differingFields=%d\n', ...
        record.threads,record.exact,numel(differences));
    for index = 1:numel(differences)
        difference = differences(index);
        fprintf('  %s reference=%.17g actual=%.17g delta=%.17g rel=%.6g\n', ...
            difference.path,difference.reference,difference.actual, ...
            difference.absoluteDifference,difference.relativeDifference);
    end
end
maxNumCompThreads(10);
report = struct('schemaVersion',1,'kind','ipm_checkpoint_signature_thread_probe', ...
    'checkpointFile',char(checkpointFile),'matlabVersion',version, ...
    'computer',computer,'records',records,'restoredThreads',maxNumCompThreads, ...
    'checkpointWritten',false,'payloadChanged',false);
clear threadCleanup;
end

function differences = compare_values(actual,reference,path)
differences = struct('path',{},'reference',{},'actual',{}, ...
    'absoluteDifference',{},'relativeDifference',{});
if isequaln(actual,reference)
    return;
end
if isstruct(reference)
    names = fieldnames(reference);
    for element = 1:numel(reference)
        for index = 1:numel(names)
            name = names{index};
            nested = sprintf('%s(%d).%s',path,element,name);
            child = compare_values(actual(element).(name), ...
                reference(element).(name),nested);
            differences = [differences,child]; %#ok<AGROW>
        end
    end
elseif iscell(reference)
    for index = 1:numel(reference)
        child = compare_values(actual{index},reference{index}, ...
            sprintf('%s{%d}',path,index));
        differences = [differences,child]; %#ok<AGROW>
    end
else
    assert(isnumeric(reference) && isnumeric(actual) && ...
        isscalar(reference) && isscalar(actual), ...
        'Only scalar numeric reductions may differ between these signatures.');
    delta = actual-reference;
    differences = struct('path',path,'reference',reference,'actual',actual, ...
        'absoluteDifference',delta,'relativeDifference', ...
        abs(delta)/max(abs(reference),realmin));
end
end
