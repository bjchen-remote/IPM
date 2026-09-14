function report=audit_outer_cutoff_frozen_old(oldRunDirectory,outputFile)
%AUDIT_OUTER_CUTOFF_FROZEN_OLD Late old-C fields normalized by outer density.
%   This is a frozen-data shape screen only. It does not calculate the new
%   gauge's RHS or create a new-C orbit.
steps=[12593,14751,17950,19135];
assert(isfolder(oldRunDirectory) && ~isfile(outputFile));
manifest=splitlines(string(fileread(fullfile(oldRunDirectory, ...
    'checkpoint_manifest.jsonl'))));
paths=cell(1,numel(steps));
for k=1:numel(manifest)
    if strlength(manifest(k))==0,continue;end
    entry=jsondecode(manifest(k));
    index=find(steps==entry.step);
    if ~isempty(index),paths{index}=entry.checkpointFile;end
end
assert(all(~cellfun(@isempty,paths)), ...
    'ipm:OuterCutoffFrozenFiles','Missing signed late old-C checkpoints.');
records=struct([]);pairs=struct([]);refinement=struct([]);
referenceWindowRms=NaN;
previous=[];
for k=1:numel(paths)
    cp=ipm.output.readCheckpoint(paths{k});
    saved=cp.payload.state;
    assert(cp.trustedRecord && strcmp(saved.config.scaling.cOmegaGauge, ...
        'wall_omega_quadratic_peak'));
    x=saved.x(:)';y=saved.y(:);
    Dx=ipm.mesh.fdMatrix(x,1,7);
    Rx=saved.rho*Dx';
    positive=find(x>=0);
    dx=diff(x(positive));
    wx=zeros(size(x));
    wx(positive)=[dx(1),dx(1:end-1)+dx(2:end),dx(end)]/2;
    z=(x-2)/0.5;
    bump=zeros(size(x));inside=abs(z)<1;
    bump(inside)=cos(pi*z(inside)/2).^2;
    weight=wx.*bump;weight=weight/sum(weight);
    M=sum(weight.*saved.rho(1,:).^2);
    if k==1,referenceWindowRms=sqrt(M);end
    amplitudeFactor=referenceWindowRms/sqrt(M);
    native=struct('x',x,'y',y,'Rx',amplitudeFactor*Rx, ...
        'tau',saved.scale.canonicalTime);
    row=struct('checkpointFile',paths{k},'step',saved.step, ...
        'tau',native.tau,'nodeCount',[numel(x),numel(y)], ...
        'oldWindowRms',sqrt(M), ...
        'counterfactualAmplitudeFactor',amplitudeFactor);
    if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
    if k>1
        pair=outer_cutoff_pair_norms(previous,native);
        if isempty(pairs),pairs=pair;else,pairs(end+1)=pair;end %#ok<AGROW>
        fine=outer_cutoff_pair_norms(previous,native,[],[],[801,401]);
        coarseNorms=[pair.rows.incrementNorm];
        fineNorms=[fine.rows.incrementNorm];
        row=struct('fromTau',pair.fromTau,'toTau',pair.toTau, ...
            'coarseNodes',pair.observationNodes, ...
            'fineNodes',fine.observationNodes, ...
            'relativeDifference',abs(coarseNorms-fineNorms)./ ...
                max(abs(fineNorms),realmin), ...
            'maximumRelativeDifference',max(abs(coarseNorms-fineNorms)./ ...
                max(abs(fineNorms),realmin)));
        if isempty(refinement),refinement=row; ...
        else,refinement(end+1)=row;end %#ok<AGROW>
    end
    previous=native;
end
report=struct('kind','old_peak_c_frozen_outer_normalized_cutoff_screen_v1', ...
    'referenceWindowRms',referenceWindowRms, ...
    'windowCenter',[2,0],'windowRadius',0.5, ...
    'records',records,'pairs',pairs, ...
    'observationGridRefinement',refinement, ...
    'newGaugeTrajectory',false,'pdeAdvanced',false, ...
    'interpretation', ...
    ['Late old-gauge fields are renormalized by their own outer density ' ...
     'window. This screens candidate cutoff radii and p values only; ' ...
     'the active new-C trajectory must be measured separately.']);
fid=fopen(outputFile,'w');assert(fid>=0);
closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
