function report=verify_shadow_native_equivalence(mainCheckpoint,shadowCheckpoint,sourceStep,outputDirectory)
%VERIFY_SHADOW_NATIVE_EQUIVALENCE Match signed original and isolated native
% histories on identical accepted steps, including a real remesh ledger.
assert(~isfolder(outputDirectory));
main=ipm.output.readCheckpoint(mainCheckpoint);
shadow=ipm.output.readCheckpoint(shadowCheckpoint);
a=main.payload.log.history;b=shadow.payload.log.history;
[matched,where]=ismember(b.common.acceptedStep,a.common.acceptedStep);
% The short branch's final accepted step is clipped to its artificial
% finalTime, so compare only preceding identical free-running steps.
use=find(matched & b.common.acceptedStep>=sourceStep & ...
    b.common.acceptedStep<b.common.acceptedStep(end));
assert(numel(use)>=2 && all(where(use)>0));
groups={'common','mesh','gauge'};checked=0;mismatches={};
for k=1:numel(groups)
    group=groups{k};left=a.(group);right=b.(group);
    names=intersect(fieldnames(left),fieldnames(right));
    for j=1:numel(names)
        name=names{j};lv=left.(name);rv=right.(name);
        if strcmp(group,'common') && strcmp(name,'canonicalEndpointDtLimit')
            % This reports each run's configured terminal horizon.
            continue;
        end
        if ~(isnumeric(lv) && isnumeric(rv) && ...
                size(lv,1)==numel(a.common.acceptedStep) && ...
                size(rv,1)==numel(b.common.acceptedStep) && ...
                size(lv,2)==size(rv,2))
            continue;
        end
        checked=checked+1;
        if ~isequaln(lv(where(use),:),rv(use,:))
            mismatches{end+1}=[group,'.',name]; %#ok<AGROW>
        end
    end
end
na=shadow.payload.state.remeshCount;
ledgerEqual=isequaln( ...
    main.payload.state.runMetadata.autonomousMesh.transactions(1:na), ...
    shadow.payload.state.runMetadata.autonomousMesh.transactions);
assert(isempty(mismatches) && ledgerEqual, ...
    'ipm:ShadowNativeEquivalence','Shadow branch changed a native recorded value.');
report=struct('kind','shadow_native_no_feedback_equivalence_v1', ...
    'sourceStep',sourceStep,'matchedSteps',b.common.acceptedStep(use), ...
    'artificialHorizonStepExcluded',b.common.acceptedStep(end), ...
    'configuredHorizonDiagnosticExcluded','common.canonicalEndpointDtLimit', ...
    'checkedNumericHistoryFields',checked, ...
    'mismatchedFields',{mismatches}, ...
    'remeshLedgerEqual',ledgerEqual, ...
    'terminalRemeshCount',na,'nativeNumericalOutputsExactlyEqual',true);
mkdir(outputDirectory);
fid=fopen(fullfile(outputDirectory,'report.json'),'w');assert(fid>=0);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('SHADOW_EQUIVALENCE steps=%d fields=%d remesh=%d\n', ...
    numel(use),checked,na);
end
