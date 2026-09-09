function report=mesh_scan_fresh_anchor_family(branchReportFile,initialDataFile,outDir)
%MESH_SCAN_FRESH_ANCHOR_FAMILY Registered no-LU 24-case constant-patch scan.
% Only completed, audited fresh branch outputs may be consumed.
assert(~isfolder(outDir));branch=jsondecode(fileread(branchReportFile));
assert(branch.audit.passed && isfile(branch.checkpointFile) && isfile(branch.resultFile));
mkdir(outDir);oldThreads=maxNumCompThreads(10);cleanup=onCleanup(@()maxNumCompThreads(oldThreads));
targets=[24,32];fineCounts=[96,128,160,192,224,256];fractions=[.90,.95];
rows=cell(1,numel(targets)*numel(fineCounts)*numel(fractions));index=0;
design=[];
for target=targets
 for fine=fineCounts
  for fraction=fractions
   index=index+1;controls=struct('targetXCoreCells',target,'targetYCoreCells',target, ...
       'positiveFineCells',fine,'anchorFineCellFraction',fraction,'minimumXFrontCells',20);
   trial=fullfile(outDir,sprintf('target%d_fine%d_fraction%02d',target,fine,round(100*fraction)));
   candidateFile=fullfile(trial,'candidate.mat');
   try
    if isempty(design)
        design=mesh_fresh_grid_design(branch.checkpointFile,branch.resultFile,initialDataFile,trial,controls);
    else
        mkdir(trial);mesh_anchor_grid_candidate(design.dataset,design.reference,design.anchor,candidateFile,controls);
    end
    if ~isfile(candidateFile)
        loaded=load(fullfile(trial,'failure.mat'),'failure');error('ipm:FreshScanConstruction','%s',loaded.failure.message);
    end
    loaded=load(candidateFile,'candidate');c=loaded.candidate;
    row=struct('controls',controls,'candidateFile',candidateFile,'passed',c.frozenCandidatePassed, ...
        'meshAdmissible',c.pairScore.admissible,'rejectionReasons',{c.pairScore.rejectionReasons}, ...
        'localQuadraturePassed',c.localQuadraturePassed,'resolutionPassed',c.resolutionPassed, ...
        'transferPassed',c.transferPassed,'aggregate',c.pairScore.aggregate, ...
        'quality',c.pairScore.quality,'constructionFailure',struct());
   catch exception
    failure=struct('identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(trial,'scan_failure.mat'),'failure','controls','-v7.3');
    row=struct('controls',controls,'candidateFile',candidateFile,'passed',false, ...
        'meshAdmissible',false,'rejectionReasons',{{}},'localQuadraturePassed',false, ...
        'resolutionPassed',false,'transferPassed',false,'aggregate',struct(), ...
        'quality',struct(),'constructionFailure',failure);
   end
   rows{index}=row;write_json(fullfile(trial,'scan_trial.json'),row);
   fprintf('FRESH_ANCHOR_SCAN %d/24 target=%d fine=%d fraction=%.2f passed=%d\n',index,target,fine,fraction,row.passed);
   clear loaded c row
  end
 end
end
report=struct('kind','completed_native_fresh_constant_patch_family_scan', ...
    'branchReportFile',branchReportFile,'initialDataFile',initialDataFile,'sourceCheckpoint',branch.checkpointFile, ...
    'targets',targets,'fineCounts',fineCounts,'anchorFractions',fractions,'rows',{rows}, ...
    'completed',true,'anyPassed',any(cellfun(@(r)r.passed,rows)), ...
    'nativeTransactionPerformed',false,'pdeAdvanced',false,'ellipticOperatorsBuilt',false);
save(fullfile(outDir,'scan_report.mat'),'report','-v7.3');write_json(fullfile(outDir,'scan_report.json'),report);
fprintf('FRESH_ANCHOR_SCAN_COMPLETE anyPassed=%d noLU=1\n',report.anyPassed);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
