function report = layout()
%LAYOUT Ensure the new packages resolve locally and the solver is independent.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
manifest = jsondecode(fileread(fullfile(root,'tests','source_map.json')));
oldNames = fieldnames(manifest.functions);
mappedFunctions = 0;
retiredFunctions = 0;
for index = 1:numel(oldNames)
    name = manifest.functions.(oldNames{index});
    if isempty(name)
        retiredFunctions = retiredFunctions+1;
        continue;
    end
    assert((ischar(name) && isrow(name)) || ...
        (isstring(name) && isscalar(name) && ~ismissing(name)), ...
        'ipm:SourceMapFunction', ...
        'Function mapping %s must be package text or JSON null.',oldNames{index});
    name = char(name);
    parts = strsplit(name,'.');
    folders = cellfun(@(part)['+',part],parts(1:end-1), ...
        'UniformOutput',false);
    parent = root;
    if strcmp(parts{1},'ipmtests')
        parent = fullfile(root,'tests');
    end
    expected = fullfile(parent,folders{:},[parts{end},'.m']);
    assert(strcmp(which(name),expected),'ipm:PackageResolution', ...
        '%s does not resolve to the independent copy: %s',name,expected);
    mappedFunctions = mappedFunctions+1;
end

dependencies = matlab.codetools.requiredFilesAndProducts( ...
    fullfile(root,'+ipm','solve.m'));
for index = 1:numel(dependencies)
    file = dependencies{index};
    localRuntime = startsWith(file,[fullfile(root,'+ipm'),filesep]);
    matlabRuntime = startsWith(file,[matlabroot,filesep]);
    assert(localRuntime || matlabRuntime,'ipm:ExternalSolverDependency', ...
        'The solver depends on a file outside its runtime packages: %s',file);
end
report = struct('passed',true,'mappedFunctions',mappedFunctions, ...
    'retiredFunctions',retiredFunctions, ...
    'runtimeDependencies',numel(dependencies));
fprintf(['Package layout passed: %d mapped functions, %d retired, ' ...
    '%d dependencies.\n'],report.mappedFunctions,report.retiredFunctions, ...
    report.runtimeDependencies);
end
