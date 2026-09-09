function state = ipm_gridlab_restore_checkpoint(inputCheckpoint)
%IPM_GRIDLAB_RESTORE_CHECKPOINT Rebuild one accepted state from its payload.
%   This reconstructs operators and flow but does not advance the PDE.  A
%   future production injection hook can pass the returned state to the
%   single maintained ipm.solve loop.

checkpoint = read_checkpoint(inputCheckpoint);
validate_checkpoint(checkpoint);
payload = checkpoint.payload;
verify_signature(payload,checkpoint.signature);

config = payload.config;
grid = config.grid;
grid.nx = numel(payload.x);
grid.ny = numel(payload.y);
grid.customX = payload.x;
grid.customY = payload.y;
ops = ipm.mesh.build(config,grid);
ops.baseX = payload.baseX;
ops.baseY = payload.baseY;
ops.rescaling = payload.rescaling;
ops.remeshCount = payload.remeshCount;
[~,flow] = ipm.evolve.flow(payload.rho,ops,payload.scale);

state = struct('config',config,'runMetadata',payload.runMetadata, ...
    'ops',ops,'rho',payload.rho,'flow',flow,'scale',payload.scale, ...
    'normalizedTime',payload.normalizedTime,'step',payload.step, ...
    'mass0',payload.mass0,'rhoRange0',payload.rhoRange0);
if ~ipm.evolve.isFinite(state)
    error('ipm:gridlab:RestoredStateFinite', ...
        'The reconstructed checkpoint state is not finite.');
end
end

function checkpoint = read_checkpoint(inputCheckpoint)
if ischar(inputCheckpoint) || ...
        (isstring(inputCheckpoint) && isscalar(inputCheckpoint))
    fileName = char(inputCheckpoint);
    if ~isfile(fileName)
        error('ipm:gridlab:CheckpointMissing', ...
            'Checkpoint file does not exist: %s',fileName);
    end
    loaded = load(fileName);
    if ~isfield(loaded,'checkpoint')
        error('ipm:gridlab:CheckpointFile', ...
            'Checkpoint MAT-file must contain a variable named checkpoint.');
    end
    checkpoint = loaded.checkpoint;
else
    checkpoint = inputCheckpoint;
end
end

function validate_checkpoint(checkpoint)
required = {'schemaVersion','kind','acceptedStep','payload','signature'};
if ~isstruct(checkpoint) || ~isscalar(checkpoint)
    error('ipm:gridlab:CheckpointContract', ...
        'Checkpoint must be a scalar structure.');
end
missing = required(~isfield(checkpoint,required));
if ~isempty(missing) || checkpoint.schemaVersion ~= 1 || ...
        ~strcmp(checkpoint.kind,'ipm_accepted_step_checkpoint') || ...
        ~(islogical(checkpoint.acceptedStep) && checkpoint.acceptedStep)
    error('ipm:gridlab:CheckpointContract', ...
        'Checkpoint does not satisfy the accepted-step version-1 contract.');
end
requiredPayload = {'config','runMetadata','rho','x','y','baseX','baseY', ...
    'rescaling','remeshCount','scale','normalizedTime','step', ...
    'mass0','rhoRange0'};
missingPayload = requiredPayload(~isfield(checkpoint.payload,requiredPayload));
if ~isempty(missingPayload)
    error('ipm:gridlab:CheckpointContract', ...
        'Checkpoint payload is missing: %s.',strjoin(missingPayload,', '));
end
payload = checkpoint.payload;
validateattributes(payload.x,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'checkpoint.payload.x');
validateattributes(payload.y,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'checkpoint.payload.y');
if ~isnumeric(payload.rho) || ~isreal(payload.rho) || ...
        any(~isfinite(payload.rho),'all') || ...
        ~isequal(size(payload.rho),[numel(payload.y),numel(payload.x)])
    error('ipm:gridlab:CheckpointSize', ...
        'Checkpoint rho must be finite and match its x/y axes.');
end
end

function verify_signature(payload,signature)
actual = struct('rhoSize',size(payload.rho), ...
    'rhoMinimum',min(payload.rho,[],'all'), ...
    'rhoMaximum',max(payload.rho,[],'all'), ...
    'rhoSum',sum(payload.rho,'all'), ...
    'rhoFrobeniusNorm',norm(payload.rho,'fro'), ...
    'xEndpoints',payload.x([1,end]), ...
    'yEndpoints',payload.y([1,end]), ...
    'normalizedTime',payload.normalizedTime, ...
    'canonicalTime',payload.scale.canonicalTime, ...
    'physicalTime',payload.scale.physicalTime,'step',payload.step);
names = fieldnames(actual);
missing = names(~isfield(signature,names));
if ~isempty(missing)
    error('ipm:gridlab:CheckpointSignature', ...
        'Checkpoint signature is incomplete.');
end
for index = 1:numel(names)
    name = names{index};
    left = actual.(name);
    right = signature.(name);
    if ~isequal(size(left),size(right)) || ...
            any(abs(double(left)-double(right)) > ...
            100*eps(max(1,max(abs(double(left)),[],'all'))),'all')
        error('ipm:gridlab:CheckpointSignature', ...
            'Checkpoint signature mismatch in %s.',name);
    end
end
end
