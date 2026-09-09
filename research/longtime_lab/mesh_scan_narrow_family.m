function search=mesh_scan_narrow_family(currentInput,rootResultFile,outputDirectory)
%MESH_SCAN_NARROW_FAMILY Frozen six-candidate scan; no LU or PDE integration.
% Complete all targets/sigmas, retaining failed designs and exact-anchor
% second-pass scores. This does not change an active driver's registered set.
assert(~isfolder(outputDirectory),'Use a new directory for every scan.');
mkdir(outputDirectory);
originalThreads=maxNumCompThreads(10);
threadCleanup=onCleanup(@()maxNumCompThreads(originalThreads));
trials=[21,.35;24,.35;21,.25;24,.25;21,.18;24,.18];
search=struct('kind','frozen_narrow_family_scan','sourceResult',currentInput, ...
    'rootResult',rootResultFile,'trials',trials,'completedTrials',{{}}, ...
    'pdeAdvanced',false,'poissonOperatorBuilt',false, ...
    'productionDriverModified',false,'completed',false);
write_search(outputDirectory,search);
for k=1:size(trials,1)
    target=trials(k,1); sigma=trials(k,2);
    controls=struct('targetXCoreCells',target,'targetYCoreCells',target, ...
        'minimumXFrontCells',20,'xEqualizationSigma',sigma);
    directory=fullfile(outputDirectory,sprintf('trial_%02d',k));
    item=struct('controls',controls,'directory',directory, ...
        'status','not_started','candidateReady',false, ...
        'candidateFile','','errorIdentifier','','errorMessage','');
    timer=tic;
    try
        report=mesh_design_stage(currentInput,rootResultFile,directory,controls);
        item.status=report.status;
        if report.transactionReady
            candidateFile=fullfile(directory,'mesh_anchor1_exact_candidate.mat');
            candidate=mesh_balance_y_candidate(report.designFile, ...
                report.referenceFile,candidateFile,1);
            item.candidateFile=candidateFile;
            item.candidateReady=candidate.transactionReady;
            item.exactAnchors=any(candidate.candidateX==1) && any(candidate.candidateX==-1);
            assert(item.exactAnchors,'The independent candidate must retain exact x=+/-1.');
            item.aggregate=candidate.pairScore.aggregate;
            item.quality=candidate.pairScore.quality;
            item.rejectionReasons=candidate.pairScore.rejectionReasons;
            item.xEqualizationInfo=candidate.xEqualizationInfo;
            item.yEqualizationInfo=candidate.yEqualizationInfo;
            clear candidate
        end
    catch exception
        item.status='exception';
        item.errorIdentifier=exception.identifier;
        item.errorMessage=exception.message;
        item.errorStack=exception.stack;
    end
    item.wallSeconds=toc(timer);
    search.completedTrials{end+1}=item;
    write_search(outputDirectory,search);
    fprintf('NARROW_MESH_TRIAL target=%g sigma=%g status=%s ready=%d seconds=%.3f\n', ...
        target,sigma,item.status,item.candidateReady,item.wallSeconds);
    clear report
end
search.completed=true;
search.admittedTrialIndices=find(cellfun(@(r)r.candidateReady,search.completedTrials));
write_search(outputDirectory,search);
fprintf('NARROW_MESH_FINISHED admitted=%s\n',mat2str(search.admittedTrialIndices));
end

function write_search(directory,search)
save(fullfile(directory,'mesh_search_report.mat'),'search','-v7.3');
fid=fopen(fullfile(directory,'mesh_search_report.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(search,PrettyPrint=true));
end
