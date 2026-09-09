function study = ipm_gridlab_evaluate_candidates( ...
        dataset,models,nodeCount,options)
%IPM_GRIDLAB_EVALUATE_CANDIDATES Build and rank grids on frozen fields.
%   Failed or inadmissible candidates are retained with infinite objective,
%   so a parameter sweep has a complete audit trail.

if nargin < 4 || isempty(options)
    options = struct();
end
options = resolve_options(options);
models = normalize_models(models);
validateattributes(nodeCount,{'numeric'}, ...
    {'scalar','integer','>=',9},mfilename,'nodeCount');
if mod(nodeCount,2) ~= 1
    error('ipm:gridlab:CandidateNodeCount', ...
        'nodeCount must be odd for the symmetric double-odd grid.');
end
if ~isfield(dataset,'xLimits') || numel(dataset.xLimits) ~= 2
    error('ipm:gridlab:CandidateDataset', ...
        'Dataset xLimits are missing or invalid.');
end

profiles = ipm_gridlab_feature_profiles(dataset,options.featureOptions);
scoreOptions = options.scoreOptions;
scoreOptions.profiles = profiles;
axisOptions = options.axisOptions;
if ~isfield(axisOptions,'referenceAxis')
    nodeCounts = arrayfun(@(value) numel(value.x),dataset.snapshots);
    distance = abs(nodeCounts-nodeCount);
    referenceIndex = find(distance == min(distance),1,'last');
    axisOptions.referenceAxis = dataset.snapshots(referenceIndex).x;
end
template = struct('index',0,'name','','model',struct(), ...
    'status','','admissible',false,'objective',inf, ...
    'axis',[],'axisInfo',struct(),'report',struct(), ...
    'errorIdentifier','','errorMessage','');
candidates = repmat(template,numel(models),1);

for index = 1:numel(models)
    model = models{index};
    name = candidate_name(model,index);
    record = template;
    record.index = index;
    record.name = name;
    record.model = model;
    try
        [axis,axisInfo] = ipm_gridlab_axis( ...
            dataset.xLimits,nodeCount,model,axisOptions);
        report = ipm_gridlab_score_frozen(dataset,axis,scoreOptions);
        record.status = 'scored';
        record.admissible = axisInfo.admissible && report.admissible;
        record.objective = report.objective;
        if ~record.admissible
            record.objective = inf;
        end
        record.axis = axis;
        record.axisInfo = axisInfo;
        record.report = report;
    catch exception
        record.status = 'failed';
        record.errorIdentifier = exception.identifier;
        record.errorMessage = exception.message;
    end
    candidates(index) = record;
end

admissible = [candidates.admissible]';
objectives = [candidates.objective]';
indices = (1:numel(candidates))';
[~,order] = sortrows([~admissible,objectives,indices],[1,2,3]);
study = struct('schemaVersion',1,'kind','ipm_frozen_grid_study', ...
    'nodeCount',nodeCount,'datasetTag',dataset.tag, ...
    'datasetSnapshotCount',dataset.snapshotCount, ...
    'axisReferenceMode',reference_mode(axisOptions), ...
    'profiles',profiles,'candidates',candidates,'order',order, ...
    'admissibleCount',nnz(admissible),'options',options);
end

function mode = reference_mode(axisOptions)
if isfield(axisOptions,'referenceAxis') && ...
        ~isempty(axisOptions.referenceAxis)
    mode = 'dataset_grid_multiplier';
else
    mode = 'uniform_physical_coordinate';
end
end

function options = resolve_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:CandidateOptions', ...
        'options must be a scalar structure.');
end
defaults = struct('featureOptions',struct(), ...
    'axisOptions',struct(),'scoreOptions',struct());
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:CandidateOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    elseif ~isstruct(options.(name)) || ~isscalar(options.(name))
        error('ipm:gridlab:CandidateOptions', ...
            'options.%s must be a scalar structure.',name);
    end
end
end

function models = normalize_models(models)
if isstruct(models)
    models = arrayfun(@(value) value,models(:),'UniformOutput',false);
elseif iscell(models)
    models = models(:);
else
    error('ipm:gridlab:CandidateModels', ...
        'models must be a structure array or cell array of structures.');
end
if isempty(models) || any(~cellfun(@(value) ...
        isstruct(value) && isscalar(value),models))
    error('ipm:gridlab:CandidateModels', ...
        'Every candidate model must be a scalar structure.');
end
end

function name = candidate_name(model,index)
if isfield(model,'name') && ...
        (ischar(model.name) || (isstring(model.name) && isscalar(model.name)))
    name = char(model.name);
else
    name = sprintf('candidate_%03d',index);
end
end
