function report = ipm_accellab_assess_peak_time(fineDirectory,coarseDirectory)
%IPM_ACCELLAB_ASSESS_PEAK_TIME Saved-field reference check, no LU or RHS.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(coarseDirectory,['fine_reference_',token]); mkdir(destination);
entries = cell(2,2); modes = {'native_quadratic','research_hermite_C1'}; grids = [49,65];
for g = 1:2
    for m = 1:2
        referenceFile = fullfile(fineDirectory,sprintf('%s_n%d_level3.mat',modes{m},grids(g)));
        reference = load(referenceFile,'run'); fine = reference.run;
        check = load(fullfile(fineDirectory,sprintf('%s_n%d_level2.mat',modes{m},grids(g))),'run');
        sensitivity = field_error(check.run,fine);
        errors = cell(3,1);
        for k = 1:3
            sourceFile = fullfile(coarseDirectory,sprintf('%s_n%d_level%d.mat',modes{m},grids(g),k));
            source = load(sourceFile,'run'); errors{k} = field_error(source.run,fine);
            errors{k}.sourceFile = sourceFile; errors{k}.dt = source.run.dt;
        end
        rho = cell2mat(cellfun(@(e)e.rhoAbsolute,errors,'UniformOutput',false));
        omega = cell2mat(cellfun(@(e)e.omegaAbsolute,errors,'UniformOutput',false));
        entries{g,m} = struct('nx',grids(g),'mode',modes{m},'referenceFile',referenceFile, ...
            'referenceDt',fine.dt,'fineReferenceSensitivity',sensitivity,'coarseErrors',{errors}, ...
            'rhoReferenceOrders',log2(rho(1:2,:)./rho(2:3,:)), ...
            'omegaReferenceOrders',log2(omega(1:2,:)./omega(2:3,:)), ...
            'coarseRhoInfinityOverReferenceSensitivity',rho(:,2)/max(sensitivity.rhoAbsolute(2),realmin));
    end
end
report = struct('status','saved_field_assessment_complete','fineDirectory',fineDirectory, ...
    'coarseDirectory',coarseDirectory,'outputDirectory',destination,'entries',{entries}, ...
    'poissonBuilds',0,'rhsEvaluations',0,'pdeSteps',0, ...
    'interpretation','The finest registered solution is a numerical reference, whose adjacent-level sensitivity is retained. Small errors or peak cancellation do not by themselves certify fourth-order convergence.');
save(fullfile(destination,'report.mat'),'report'); fid = fopen(fullfile(destination,'report.json'),'w');
cleanup = onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); clear cleanup;
fprintf('HERMITE_PEAK_TIME_REFERENCE %s\n',destination);
for k = 1:numel(entries)
    e = entries{k}; fprintf('PEAK_REFERENCE n%d %s rhoOrders [%g %g; %g %g]\n', ...
        e.nx,e.mode,e.rhoReferenceOrders(1,:),e.rhoReferenceOrders(2,:));
end
end

function e = field_error(a,b)
assert(isequal(a.x,b.x) && isequal(a.y,b.y) && abs(a.z(5)-b.z(5)) < 1e-12);
q = a.rho-b.rho; w = b.integrationWeights; f = q*b.Dx';
[~,index] = max(abs(q(:))); [j,i] = ind2sub(size(q),index);
e = struct('rhoAbsolute',norms(q,w),'omegaAbsolute',norms(f,w), ...
    'rhoMaximumErrorLocation',[b.x(i),b.y(j)],'rhoMaximumErrorSigned',q(j,i), ...
    'allFiveScaleErrors',a.z-b.z);
end

function v = norms(a,w)
v = [sqrt(sum(a.^2.*w,'all')/sum(w,'all')),max(abs(a),[],'all')];
end
