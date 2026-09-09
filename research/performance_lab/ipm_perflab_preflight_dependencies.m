function report=ipm_perflab_preflight_dependencies(sourceRoot,entries)
%IPM_PERFLAB_PREFLIGHT_DEPENDENCIES Conservative registered project closure.
% MATLAB's analyzer is supplemented by explicit project-call scanning so a
% missing transitive .m file cannot disappear from the required-files list.
prefix='(?:ipm(?:\.[A-Za-z]\w*){1,}|(?:ipm_gridlab_|ipm_accellab_|mesh_|fresh_profile_|run_fresh_profile_|continuous_inner_)\w+|fourth_order_smooth_peak_grid)';
queue=entries(:)';seen={};files={};issues=struct('functionName',{},'reason',{},'resolvedFile',{});
edges=struct('caller',{},'callee',{});dynamicFiles={};
while ~isempty(queue)
    name=queue{1};queue(1)=[];if ismember(name,seen),continue;end;seen{end+1}=name; %#ok<AGROW>
    file=which(name);
    if isempty(file)||~isfile(file)||~startsWith(canonical(file),[sourceRoot,filesep])
        issues(end+1)=struct('functionName',name,'reason','missing_from_registered_source_root','resolvedFile',file); %#ok<AGROW>
        continue
    end
    file=canonical(file);files{end+1}=file; %#ok<AGROW>
    code=fileread(file);
    % MATLAB strings/comments may conservatively add dependencies. This is
    % intentional: a false negative is safer than silently omitting a file.
    calls=regexp(code,['\b(',prefix,')\s*(?=\()'],'tokens');
    handles=regexp(code,['@\s*(',prefix,')'],'tokens');calls=[calls,handles]; %#ok<AGROW>
    names=unique(cellfun(@(x)x{1},calls,'UniformOutput',false),'stable');
    for k=1:numel(names)
        edges(end+1)=struct('caller',name,'callee',names{k}); %#ok<AGROW>
    end
    queue=[queue,names]; %#ok<AGROW>
    if ~isempty(regexp(code,'\b(?:eval|evalin|str2func|feval)\s*\(','once'))
        dynamicFiles{end+1}=file; %#ok<AGROW>
    end
end
files=unique(files,'stable');analyzerFiles={};analyzerProducts=[];analyzerFailure='';outside={};
if isempty(issues)
    try
        [analyzerFiles,analyzerProducts]=matlab.codetools.requiredFilesAndProducts(files);
        analyzerFiles=cellstr(analyzerFiles);
        for k=1:numel(analyzerFiles)
            f=canonical(analyzerFiles{k});
            if ~startsWith(f,[sourceRoot,filesep])&&~startsWith(f,[matlabroot,filesep])
                outside{end+1}=f; %#ok<AGROW>
            end
        end
    catch exception
        analyzerFailure=exception.message;
    end
end
allFiles=unique([files(:);analyzerFiles(:)],'stable');
manifest=struct('file',{},'sha256',{});
for k=1:numel(allFiles)
    f=canonical(allFiles{k});
    if startsWith(f,[sourceRoot,filesep])
        manifest(end+1)=struct('file',f,'sha256',file_sha256(f)); %#ok<AGROW>
    end
end
report=struct('entries',{entries},'projectFunctions',{seen},'projectFiles',{files}, ...
    'edges',edges,'unavailable',issues,'analyzerFiles',{analyzerFiles}, ...
    'analyzerProducts',analyzerProducts,'analyzerFailure',analyzerFailure, ...
    'dependenciesOutsideSource',{outside},'dynamicCallSites',{unique(dynamicFiles)},'manifest',manifest, ...
    'passed',isempty(issues)&&isempty(analyzerFailure)&&isempty(outside), ...
    'scope','Complete conservative closure of registered project entry points plus MATLAB requiredFilesAndProducts. Runtime-generated function names and undeclared callbacks require explicit registration; no arbitrary dynamic-code completeness claim.');
end
function hash=file_sha256(file)
fid=fopen(file,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
function value=canonical(value)
value=char(java.io.File(value).getCanonicalPath());
end
