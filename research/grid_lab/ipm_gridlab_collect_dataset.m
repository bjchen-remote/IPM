function dataset = ipm_gridlab_collect_dataset(inputs,options)
%IPM_GRIDLAB_COLLECT_DATASET Freeze fields from version-2 solver results.
%   DATASET = IPM_GRIDLAB_COLLECT_DATASET(INPUTS) accepts one result, one
%   MAT-file path, or a cell array of those inputs.  Stored snapshots are
%   copied with their paired grids, physical reconstruction, scaling, and
%   provenance; a result without stored snapshots contributes only its
%   terminal accepted field.
%
%   OPTIONS.targetCanonicalTimes selects the nearest available samples in
%   each result. OPTIONS.trustedOnly discards samples outside the solver's
%   continuous trusted prefix. OPTIONS.outputFile saves DATASET with -v7.3.

if nargin < 2 || isempty(options)
    options = struct();
end
options = resolve_options(options);
[items,labels] = normalize_inputs(inputs);

template = struct( ...
    'datasetIndex',0,'sourceIndex',0,'sourceLabel','', ...
    'snapshotIndex',0,'historyIndex',0,'terminalFallback',false, ...
    'trusted',false,'rho',[],'physicalRho',[], ...
    'x',[],'y',[],'physicalX',[],'physicalY',[], ...
    'scale',struct(), ...
    'normalizedTime',NaN,'canonicalTime',NaN,'physicalTime',NaN, ...
    'signature',struct());
snapshots = repmat(template,0,1);
sourceSummary = repmat(struct('label','','storedSnapshotCount',0, ...
    'selectedSnapshotCount',0,'stopReason','','nodeCount',[0,0], ...
    'metadata',struct(),'config',struct(),'gridSignature',struct()), ...
    numel(items),1);

for sourceIndex = 1:numel(items)
    result = ipm.output.validate(items{sourceIndex});
    allSnapshots = ipm.output.snapshotAt(result,'all');
    terminalFallback = isempty(result.snapshots.rho);
    selected = select_indices(allSnapshots,options.targetCanonicalTimes);
    trustedPrefix = ipm.output.continuousTrustedPrefix( ...
        ipm.output.trustedMask(result));
    historyTau = result.history.common.canonicalTau(:);
    kept = 0;
    for localIndex = selected(:)'
        snapshot = allSnapshots(localIndex);
        [~,historyIndex] = min(abs(historyTau-snapshot.canonicalTime));
        isTrusted = trustedPrefix(historyIndex);
        if options.trustedOnly && ~isTrusted
            continue
        end
        rho = double(snapshot.rho);
        physicalRho = double(snapshot.physicalRho);
        x = double(snapshot.x(:)');
        y = double(snapshot.y(:));
        physicalX = double(snapshot.physicalX(:)');
        physicalY = double(snapshot.physicalY(:));
        record = template;
        record.datasetIndex = numel(snapshots)+1;
        record.sourceIndex = sourceIndex;
        record.sourceLabel = labels{sourceIndex};
        record.snapshotIndex = snapshot.index;
        record.historyIndex = historyIndex;
        record.terminalFallback = terminalFallback;
        record.trusted = isTrusted;
        record.rho = rho;
        record.physicalRho = physicalRho;
        record.x = x;
        record.y = y;
        record.physicalX = physicalX;
        record.physicalY = physicalY;
        record.scale = scale_at(result,historyIndex);
        record.normalizedTime = snapshot.normalizedTime;
        record.canonicalTime = snapshot.canonicalTime;
        record.physicalTime = snapshot.physicalTime;
        record.signature = field_signature(rho,x,y);
        snapshots(end+1,1) = record; %#ok<AGROW>
        kept = kept+1;
    end
    sourceSummary(sourceIndex) = struct( ...
        'label',labels{sourceIndex}, ...
        'storedSnapshotCount',numel(result.snapshots.rho), ...
        'selectedSnapshotCount',kept, ...
        'stopReason',char(result.state.stopReason), ...
        'nodeCount',[numel(result.grid.x),numel(result.grid.y)], ...
        'metadata',result.metadata,'config',result.config, ...
        'gridSignature',grid_signature(result.grid.x,result.grid.y));
end
if isempty(snapshots)
    error('ipm:gridlab:EmptyDataset', ...
        'No snapshots survived the requested selection.');
end

leftEndpoints = arrayfun(@(value) value.x(1),snapshots);
rightEndpoints = arrayfun(@(value) value.x(end),snapshots);
referenceLimits = [snapshots(1).x(1),snapshots(1).x(end)];
domainTolerance = 64*eps(max(1,max(abs(referenceLimits))));
if any(abs(leftEndpoints-referenceLimits(1)) > domainTolerance) || ...
        any(abs(rightEndpoints-referenceLimits(2)) > domainTolerance)
    error('ipm:gridlab:DatasetDomain', ...
        'All frozen fields must share the same x-domain endpoints.');
end
yLowerEndpoints = arrayfun(@(value) value.y(1),snapshots);
yUpperEndpoints = arrayfun(@(value) value.y(end),snapshots);
yReferenceLimits = [snapshots(1).y(1),snapshots(1).y(end)];
yDomainTolerance = 64*eps(max(1,max(abs(yReferenceLimits))));
if any(abs(yLowerEndpoints-yReferenceLimits(1)) > yDomainTolerance) || ...
        any(abs(yUpperEndpoints-yReferenceLimits(2)) > yDomainTolerance)
    error('ipm:gridlab:DatasetDomain', ...
        'All frozen fields must share the same y-domain endpoints.');
end

dataset = struct();
dataset.schemaVersion = 2;
dataset.kind = 'ipm_frozen_grid_dataset';
dataset.createdUtc = char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
dataset.tag = options.tag;
dataset.snapshotCount = numel(snapshots);
dataset.xLimits = referenceLimits;
dataset.yLimits = yReferenceLimits;
dataset.snapshots = snapshots;
dataset.sources = sourceSummary;
dataset.selection = rmfield(options,'outputFile');

if ~isempty(options.outputFile)
    outputFile = char(options.outputFile);
    outputFolder = fileparts(outputFile);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        error('ipm:gridlab:DatasetFolder', ...
            'The output folder does not exist: %s',outputFolder);
    end
    save(outputFile,'dataset','-v7.3');
end
end

function options = resolve_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:DatasetOptions', ...
        'options must be a scalar structure.');
end
defaults = struct('targetCanonicalTimes',[], ...
    'trustedOnly',false,'outputFile','','tag','');
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:DatasetOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
end
if ~isempty(options.targetCanonicalTimes)
    validateattributes(options.targetCanonicalTimes,{'numeric'}, ...
        {'vector','real','finite','nonnegative'},mfilename, ...
        'options.targetCanonicalTimes');
elseif ~isnumeric(options.targetCanonicalTimes)
    error('ipm:gridlab:DatasetOptions', ...
        'options.targetCanonicalTimes must be numeric.');
end
if ~(islogical(options.trustedOnly) && isscalar(options.trustedOnly))
    error('ipm:gridlab:DatasetOptions', ...
        'options.trustedOnly must be a logical scalar.');
end
for field = {'outputFile','tag'}
    value = options.(field{1});
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('ipm:gridlab:DatasetOptions', ...
            'options.%s must be text.',field{1});
    end
    options.(field{1}) = char(value);
end
end

function [items,labels] = normalize_inputs(inputs)
if ~iscell(inputs)
    if isstring(inputs) && ~isscalar(inputs)
        items = cellstr(inputs(:));
    else
        items = {inputs};
    end
else
    items = inputs(:);
end
if isempty(items)
    error('ipm:gridlab:DatasetInput','At least one result is required.');
end
labels = cell(size(items));
for index = 1:numel(items)
    item = items{index};
    if ischar(item) || (isstring(item) && isscalar(item))
        labels{index} = char(item);
    elseif isstruct(item) && isscalar(item)
        labels{index} = sprintf('memory_result_%d',index);
    else
        error('ipm:gridlab:DatasetInput', ...
            'Each input must be a scalar result or MAT-file path.');
    end
end
end

function selected = select_indices(snapshots,targetTimes)
if isempty(targetTimes)
    selected = 1:numel(snapshots);
    return
end
available = [snapshots.canonicalTime];
selected = zeros(size(targetTimes));
for index = 1:numel(targetTimes)
    [~,selected(index)] = min(abs(available-targetTimes(index)));
end
selected = unique(selected,'stable');
end

function signature = field_signature(rho,x,y)
signature = struct( ...
    'size',size(rho), ...
    'minimum',min(rho,[],'all'), ...
    'maximum',max(rho,[],'all'), ...
    'sum',sum(rho,'all'), ...
    'frobeniusNorm',norm(rho,'fro'), ...
    'xEndpoints',[x(1),x(end)], ...
    'yEndpoints',[y(1),y(end)], ...
    'xWeightedSum',weighted_sum(x), ...
    'yWeightedSum',weighted_sum(y), ...
    'xSpacingNorm',norm(diff(x)), ...
    'ySpacingNorm',norm(diff(y)));
end

function scale = scale_at(result,historyIndex)
common = result.history.common;
scale = struct( ...
    'Cx',history_value(common,{'C_x','C_l'},historyIndex, ...
        result.scale.Cx), ...
    'Cy',history_value(common,{'C_y'},historyIndex,result.scale.Cy), ...
    'Comega',history_value(common,{'C_omega'},historyIndex, ...
        result.scale.Comega), ...
    'Xshift',history_value(common,{'X_shift'},historyIndex, ...
        result.scale.Xshift));
if any(~isfinite([scale.Cx,scale.Cy,scale.Comega,scale.Xshift])) || ...
        any([scale.Cx,scale.Cy,scale.Comega] <= 0)
    error('ipm:gridlab:DatasetScale', ...
        'Frozen scaling factors must be finite and positive.');
end
end

function value = history_value(history,names,index,fallback)
value = [];
for nameIndex = 1:numel(names)
    name = names{nameIndex};
    if isfield(history,name) && numel(history.(name)) >= index
        value = history.(name)(index);
        break
    end
end
if isempty(value)
    value = fallback;
end
value = double(value);
end

function signature = grid_signature(x,y)
x = double(x(:)');
y = double(y(:));
signature = struct('nodeCount',[numel(x),numel(y)], ...
    'xEndpoints',x([1,end]),'yEndpoints',y([1,end]), ...
    'xWeightedSum',weighted_sum(x), ...
    'yWeightedSum',weighted_sum(y), ...
    'xSpacingNorm',norm(diff(x)), ...
    'ySpacingNorm',norm(diff(y)));
end

function value = weighted_sum(axis)
axis = double(axis(:));
weights = (1:numel(axis))';
value = sum(weights.*axis);
end
