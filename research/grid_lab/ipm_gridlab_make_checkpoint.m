function checkpoint = ipm_gridlab_make_checkpoint(acceptedState,options)
%IPM_GRIDLAB_MAKE_CHECKPOINT Pack one accepted solver state for persistence.
%   Derived matrices, decompositions, and flow are intentionally omitted.
%   The exact grid and rescaling references are sufficient to rebuild them.
%   This function defines a contract; the production solver does not yet
%   call it automatically.

if nargin < 2 || isempty(options)
    options = struct();
end
options = resolve_options(options);
validate_state(acceptedState);
ops = acceptedState.ops;

payload = struct( ...
    'config',acceptedState.config, ...
    'runMetadata',acceptedState.runMetadata, ...
    'rho',double(acceptedState.rho), ...
    'x',double(ops.x(:)'), ...
    'y',double(ops.y(:)), ...
    'baseX',double(ops.baseX(:)'), ...
    'baseY',double(ops.baseY(:)), ...
    'rescaling',ops.rescaling, ...
    'remeshCount',ops.remeshCount, ...
    'scale',acceptedState.scale, ...
    'normalizedTime',acceptedState.normalizedTime, ...
    'step',acceptedState.step, ...
    'mass0',acceptedState.mass0, ...
    'rhoRange0',acceptedState.rhoRange0);

checkpoint = struct();
checkpoint.schemaVersion = 1;
checkpoint.kind = 'ipm_accepted_step_checkpoint';
checkpoint.createdUtc = char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
checkpoint.tag = options.tag;
checkpoint.acceptedStep = true;
checkpoint.payload = payload;
checkpoint.signature = state_signature(payload);
checkpoint.resumePolicy = struct( ...
    'rebuildOperators',true,'recomputeFlow',true, ...
    'historyMode','new_segment', ...
    'requiresIpmSolveInjectionHook',true);

if ~isempty(options.outputFile)
    outputFile = char(options.outputFile);
    outputFolder = fileparts(outputFile);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        error('ipm:gridlab:CheckpointFolder', ...
            'The output folder does not exist: %s',outputFolder);
    end
    save(outputFile,'checkpoint','-v7.3');
end
end

function options = resolve_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:CheckpointOptions', ...
        'options must be a scalar structure.');
end
defaults = struct('outputFile','','tag','');
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:CheckpointOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
    value = options.(name);
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('ipm:gridlab:CheckpointOptions', ...
            'options.%s must be text.',name);
    end
    options.(name) = char(value);
end
end

function validate_state(state)
required = {'config','runMetadata','ops','rho','flow','scale', ...
    'normalizedTime','step','mass0','rhoRange0'};
if ~isstruct(state) || ~isscalar(state)
    error('ipm:gridlab:CheckpointState', ...
        'acceptedState must be a scalar solver-state structure.');
end
missing = required(~isfield(state,required));
if ~isempty(missing)
    error('ipm:gridlab:CheckpointState', ...
        'acceptedState is missing: %s.',strjoin(missing,', '));
end
requiredOps = {'x','y','baseX','baseY','rescaling','remeshCount'};
missingOps = requiredOps(~isfield(state.ops,requiredOps));
if ~isempty(missingOps)
    error('ipm:gridlab:CheckpointState', ...
        'acceptedState.ops is missing: %s.',strjoin(missingOps,', '));
end
if ~ipm.evolve.isFinite(state)
    error('ipm:gridlab:CheckpointFinite', ...
        'Only a complete finite accepted state may be checkpointed.');
end
if ~isequal(size(state.rho),[numel(state.ops.y),numel(state.ops.x)])
    error('ipm:gridlab:CheckpointSize', ...
        'rho must match the accepted x/y grid.');
end
end

function signature = state_signature(payload)
rho = payload.rho;
signature = struct('rhoSize',size(rho), ...
    'rhoMinimum',min(rho,[],'all'),'rhoMaximum',max(rho,[],'all'), ...
    'rhoSum',sum(rho,'all'),'rhoFrobeniusNorm',norm(rho,'fro'), ...
    'xEndpoints',payload.x([1,end]), ...
    'yEndpoints',payload.y([1,end]), ...
    'normalizedTime',payload.normalizedTime, ...
    'canonicalTime',payload.scale.canonicalTime, ...
    'physicalTime',payload.scale.physicalTime,'step',payload.step);
end
