function report = ipm_accellab_continuous_sampling_from_cache(fieldsFile,reportFile,outputRoot)
%IPM_ACCELLAB_CONTINUOUS_SAMPLING_FROM_CACHE Separate mask and sampling errors.
% Uses the already saved small observer arrays; no global field is loaded.
data = load(fieldsFile,'details'); prior = load(reportFile,'report');
assert(isfield(prior.report,'residual')); reference = prior.report.residual.gauss;
rows = cell(1,numel(data.details.sampleFields));
for k = 1:numel(rows)
    values = data.details.sampleFields{k}; x = values.xi; y = values.eta;
    fullWeights = weights(y)'*weights(x);
    cx = abs(x)<=1; cy = y<=1.5;
    coreX = zeros(size(x)); coreY = zeros(size(y));
    coreX(cx) = weights(x(cx)); coreY(cy) = weights(y(cy));
    coreWeights = coreY'*coreX; holdoutWeights = fullWeights-coreWeights;
    [XI,ETA] = meshgrid(x,y); mask = abs(XI)<=1 & ETA<=1.5;
    assert(abs(sum(coreWeights,'all')-3)<1e-12 && abs(sum(holdoutWeights,'all')-9)<1e-12 && ...
        min(holdoutWeights,[],'all')>=0);
    row = struct('nodeCount',[numel(x),numel(y)], ...
        'originalMaskedCoreEffectiveArea',sum(fullWeights(mask)), ...
        'clippedCoreArea',sum(coreWeights,'all'),'clippedHoldoutArea',sum(holdoutWeights,'all'));
    for name = {'G','U','amplitudeNormalizedEulerian','translationOnly'}
        key = name{1}; f = values.(key);
        row.(key) = struct('coreL2',norm2(f,coreWeights),'holdoutL2',norm2(f,holdoutWeights), ...
            'fullLocalL2',norm2(f,fullWeights));
    end
    row.originalMaskedG = prior.report.residual.sampledMetrics{k}.G;
    for name = {'coreL2','holdoutL2','fullLocalL2'}
        key = name{1}; row.GRelativeDifferenceFromGauss.(key) = row.G.(key)/reference.G.(key)-1;
        row.GOverAmplitudeNormalizedRhs.(key) = row.G.(key)/row.amplitudeNormalizedEulerian.(key);
    end
    rows{k} = row;
end
report = struct('kind','same_rectangle_boundary_clipped_sampling_audit','sourceFields',fieldsFile, ...
    'sourceReport',reportFile,'gaussReference',reference,'sampling',{rows}, ...
    'originalMeasurementsPreserved',true,'globalFieldLoads',0,'pdeSteps',0,'rhsEvaluations',0,'poissonOperatorBuilds',0, ...
    'interpretation','Original full-grid weights followed by a core mask had different effective areas. Boundary-clipped tensor trapezoids use the same rectangles as Gauss; their remaining differences are observation quadrature sensitivity, not a PDE error bound.');
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['continuous_sampling_',token]); mkdir(destination);
save(fullfile(destination,'report.mat'),'report');
fid=fopen(fullfile(destination,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); fprintf('CONTINUOUS_SAMPLING %s\n',destination);
end
function w = weights(x)
x=x(:)';h=diff(x);w=[h(1)/2,(h(1:end-1)+h(2:end))/2,h(end)/2];
end
function v = norm2(f,w)
v=sqrt(sum(w.*f.^2,'all')/sum(w,'all'));
end
