function [checkpoint,audit,fileName,memoryInfo] = ipm_perflab_regrid_checkpoint( ...
        inputCheckpoint,candidateX,candidateY,options)
%IPM_PERFLAB_REGRID_CHECKPOINT Isolated old-factor-release candidate.
%   Research copy of ipm_gridlab_regrid_checkpoint, 2026-09-08.
%   Only releases the obsolete Poisson factor before copying oldState and
%   before the maintained transfer builds the new factor. All numerical
%   transfer, flow, audit and checkpoint code below remains unchanged.
%   The production checkpoint is restored without advancing time, its field
%   is transferred with the maintained remesh path, and all derived
%   operators/flow diagnostics are rebuilt.  Scaling clocks, normalization
%   references, initial invariants, history, and output cursors are retained.
%   The old terminal history row is replaced at the identical accepted time
%   so the new grid's diagnostics are the segment restart datum.

if nargin < 4 || isempty(options)
    options = struct();
end
options = resolve_options(options);
sourceCheckpoint = ipm.output.readCheckpoint(inputCheckpoint);
[state,log,cursor] = ipm.output.restoreCheckpoint( ...
    sourceCheckpoint,struct());
% restoreCheckpoint correctly marks an actual continuation.  This operation
% is instead a zero-time grid transaction, so start its provenance from the
% immutable source metadata and record the regrid explicitly below.
state.runMetadata = sourceCheckpoint.payload.state.runMetadata;

[candidateX,candidateY] = validate_axes( ...
    candidateX,candidateY,state.ops,options.meshLimits);
% No downstream old-grid transfer/audit consumer reads ops.poisson. Remove
% its only live owner BEFORE oldState is copied; removing it afterward from
% state alone would leave the LU alive through oldState.ops.poisson.
% Source checkpoints contain no ops or factorization object.
memoryInfo = struct('schemaVersion',1, ...
    'kind','ipm_old_poisson_factor_release', ...
    'oldPoissonClass',class(state.ops.poisson), ...
    'releasedFields',{{'poisson'}}, ...
    'releasePoint','after axis validation, before oldState copy and transfer');
state.ops = rmfield(state.ops,'poisson');
oldState = state;
memoryInfo.oldStateHasPoisson = isfield(oldState.ops,'poisson');
memoryInfo.stateHasPoissonAtTransfer = isfield(state.ops,'poisson');
newConfig = regridded_config( ...
    state.config,candidateX,candidateY);
proposal = struct('x',candidateX,'y',candidateY);
[rhoNew,opsNew,transferMetrics] = ipm.remesh.transfer( ...
    state.rho,state.ops,newConfig,proposal);
memoryInfo.transferMetrics = transferMetrics;
memoryInfo.newPoissonClass = class(opsNew.poisson);
opsNew.remeshCount = state.ops.remeshCount+1;
[~,flowNew] = ipm.evolve.flow(rhoNew,opsNew,state.scale);

state.config = newConfig;
state.ops = opsNew;
state.rho = rhoNew;
state.flow = flowNew;
audit = audit_transaction(oldState,state,transferMetrics,options);
assert_audit(audit);
state.runMetadata = append_regrid_metadata( ...
    state.runMetadata,audit,options.tag);

log = remove_terminal_record(log,state.config.output.storeSnapshots);
[log,~] = ipm.output.record(log,state);
cursor.lastCheckpointStep = state.step;
cursor.lastCheckpointCanonicalTime = state.scale.canonicalTime;
checkpoint = ipm.output.makeCheckpoint(state,log,cursor);
fileName = "";
if ~isempty(options.outputFile)
    [checkpoint,fileName] = ipm.output.writeCheckpoint( ...
        checkpoint,options.outputFile);
end
end

function options = resolve_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('ipm:gridlab:RegridOptions', ...
        'options must be a scalar structure.');
end
defaultLimits = struct('maxAdjacentCellRatio',1.08, ...
    'maxLogSpacingCurvature',0.03,'minStencilRcond',1e-9, ...
    'minQuadratureWeightRatio',1e-4);
defaults = struct('outputFile','','tag','', ...
    'minimumXCoreCells',16,'minimumYCoreCells',16, ...
    'maximumRhoXRelativeChange',5e-3, ...
    'maximumMassRelativeDefect',5e-12, ...
    'maximumRelativeRangeViolation',2e-4, ...
    'meshLimits',defaultLimits);
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(options),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:RegridOptions', ...
        'Unknown option(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
end
for field = {'outputFile','tag'}
    value = options.(field{1});
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('ipm:gridlab:RegridOptions', ...
            'options.%s must be text.',field{1});
    end
    options.(field{1}) = char(value);
end
for field = {'minimumXCoreCells','minimumYCoreCells'}
    validateattributes(options.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename, ...
        ['options.' field{1}]);
end
for field = {'maximumRhoXRelativeChange', ...
        'maximumMassRelativeDefect','maximumRelativeRangeViolation'}
    validateattributes(options.(field{1}),{'numeric'}, ...
        {'scalar','real','finite','nonnegative'},mfilename, ...
        ['options.' field{1}]);
end
options.meshLimits = merge_limits(options.meshLimits,defaultLimits);
end

function limits = merge_limits(limits,defaults)
if ~isstruct(limits) || ~isscalar(limits)
    error('ipm:gridlab:RegridOptions', ...
        'options.meshLimits must be a scalar structure.');
end
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(limits),names,'stable');
if ~isempty(unexpected)
    error('ipm:gridlab:RegridOptions', ...
        'Unknown meshLimits field(s): %s.',strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(limits,name)
        limits.(name) = defaults.(name);
    end
    validateattributes(limits.(name),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename, ...
        ['options.meshLimits.' name]);
end
end

function [x,y] = validate_axes(x,y,ops,limits)
x = double(x(:)');
y = double(y(:));
validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'candidateX');
validateattributes(y,{'numeric'}, ...
    {'vector','real','finite','increasing','nonnegative'}, ...
    mfilename,'candidateY');
if numel(x) ~= ops.nx || numel(y) ~= ops.ny
    error('ipm:gridlab:RegridNodeCount', ...
        'This qN continuation transaction preserves both node counts.');
end
tolerance = 1e3*eps(max([1,max(abs(x)),max(abs(y))]));
if mod(numel(x),2) ~= 1 || max(abs(x+fliplr(x))) > tolerance || ...
        abs(y(1)) > tolerance
    error('ipm:gridlab:RegridAxis', ...
        'candidateX must be odd symmetric and candidateY must start at 0.');
end
if ~same_endpoints(x,ops.x) || ~same_endpoints(y,ops.y)
    error('ipm:gridlab:RegridDomain', ...
        'A continuation regrid may not change the computational box.');
end
admit_axis(x,limits,'x');
admit_axis(y,limits,'y');
end

function admit_axis(axis,limits,label)
weights = ipm.mesh.quadrature(axis);
quality = ipm.mesh.quality(axis,7,weights);
reasons = {};
if quality.maximumAdjacentCellRatio > limits.maxAdjacentCellRatio
    reasons{end+1} = 'adjacent ratio';
end
if quality.maximumLogSpacingCurvature > limits.maxLogSpacingCurvature
    reasons{end+1} = 'log-spacing curvature';
end
if quality.minimumStencilRcond < limits.minStencilRcond
    reasons{end+1} = 'stencil conditioning';
end
if quality.minimumQuadratureWeightRatio < ...
        limits.minQuadratureWeightRatio || any(weights <= 0)
    reasons{end+1} = 'quadrature margin';
end
if ~isempty(reasons)
    error('ipm:gridlab:RegridMesh', ...
        'Candidate %s axis fails: %s.',label,strjoin(reasons,', '));
end
end

function config = regridded_config(config,x,y)
schema = ipm.config.schema();
flat = struct();
for domainIndex = 1:numel(schema.domainNames)
    domain = config.(schema.domainNames{domainIndex});
    names = fieldnames(domain);
    for nameIndex = 1:numel(names)
        flat.(names{nameIndex}) = domain.(names{nameIndex});
    end
end
flat.nx = numel(x);
flat.ny = numel(y);
flat.customX = x;
flat.customY = y;
config = ipm.config.resolve(flat);
end

function audit = audit_transaction(oldState,newState,metrics,options)
oldWeights = oldState.ops.integrationWeights;
newWeights = newState.ops.integrationWeights;
oldMass = sum(oldState.rho.*oldWeights,'all');
newMass = sum(newState.rho.*newWeights,'all');
massDefect = abs(newMass-oldMass)/max(abs(oldMass),eps);

oldRhoX = oldState.rho*oldState.ops.Dx';
newRhoX = newState.rho*newState.ops.Dx';
oldMaximum = max(abs(oldRhoX),[],'all');
newMaximum = max(abs(newRhoX),[],'all');
rhoXChange = abs(newMaximum-oldMaximum)/max(oldMaximum,eps);
scaleFactor = exp(oldState.scale.logC_l-oldState.scale.logC_omega);

oldMinimum = min(oldState.rho,[],'all');
oldMaximumRho = max(oldState.rho,[],'all');
rangeScale = max([oldMaximumRho-oldMinimum,abs(oldMinimum), ...
    abs(oldMaximumRho),eps]);
rangeViolation = max([oldMinimum-min(newState.rho,[],'all'), ...
    max(newState.rho,[],'all')-oldMaximumRho,0])/rangeScale;
xWeights = ipm.mesh.quadrature(newState.ops.x);
yWeights = ipm.mesh.quadrature(newState.ops.y);
xQuality = ipm.mesh.quality(newState.ops.x,7,xWeights);
yQuality = ipm.mesh.quality(newState.ops.y,7,yWeights);

passed = metrics.xCore >= options.minimumXCoreCells && ...
    metrics.yCore >= options.minimumYCoreCells && ...
    rhoXChange <= options.maximumRhoXRelativeChange && ...
    massDefect <= options.maximumMassRelativeDefect && ...
    rangeViolation <= options.maximumRelativeRangeViolation;
audit = struct('schemaVersion',1, ...
    'kind','ipm_gridlab_checkpoint_regrid_audit', ...
    'passed',passed,'step',newState.step, ...
    'normalizedTime',newState.normalizedTime, ...
    'canonicalTime',newState.scale.canonicalTime, ...
    'physicalTime',newState.scale.physicalTime, ...
    'oldNodeCount',[oldState.ops.nx,oldState.ops.ny], ...
    'newNodeCount',[newState.ops.nx,newState.ops.ny], ...
    'oldRemeshCount',oldState.ops.remeshCount, ...
    'newRemeshCount',newState.ops.remeshCount, ...
    'xCoreCells',metrics.xCore,'yCoreCells',metrics.yCore, ...
    'rescaledRhoXMaximumOld',oldMaximum, ...
    'rescaledRhoXMaximumNew',newMaximum, ...
    'physicalRhoXMaximumOld',scaleFactor*oldMaximum, ...
    'physicalRhoXMaximumNew',scaleFactor*newMaximum, ...
    'rhoXMaximumRelativeChange',rhoXChange, ...
    'massRelativeDefect',massDefect, ...
    'relativeRangeViolation',rangeViolation, ...
    'xSymmetryRelativeDefect',norm( ...
        newState.rho-fliplr(newState.rho),'fro') / ...
        max(norm(newState.rho,'fro'),eps), ...
    'xQuality',xQuality,'yQuality',yQuality, ...
    'limits',rmfield(options,{'outputFile','tag','meshLimits'}), ...
    'meshLimits',options.meshLimits);
end

function assert_audit(audit)
if audit.xCoreCells < audit.limits.minimumXCoreCells || ...
        audit.yCoreCells < audit.limits.minimumYCoreCells
    error('ipm:gridlab:RegridResolution', ...
        'Regridded x/y core %.6g/%.6g misses the required %.6g/%.6g.', ...
        audit.xCoreCells,audit.yCoreCells, ...
        audit.limits.minimumXCoreCells,audit.limits.minimumYCoreCells);
end
if audit.rhoXMaximumRelativeChange > ...
        audit.limits.maximumRhoXRelativeChange
    error('ipm:gridlab:RegridRhoX', ...
        'Regridding changed max|rho_x| by %.6g.', ...
        audit.rhoXMaximumRelativeChange);
end
if audit.massRelativeDefect > audit.limits.maximumMassRelativeDefect
    error('ipm:gridlab:RegridMass', ...
        'Regridding has relative mass defect %.6g.', ...
        audit.massRelativeDefect);
end
if audit.relativeRangeViolation > ...
        audit.limits.maximumRelativeRangeViolation
    error('ipm:gridlab:RegridRange', ...
        'Regridding has relative range violation %.6g.', ...
        audit.relativeRangeViolation);
end
if ~audit.passed
    error('ipm:gridlab:RegridAudit', ...
        'The checkpoint regrid audit did not pass.');
end
end

function metadata = append_regrid_metadata(metadata,audit,tag)
summary = rmfield(audit,{'xQuality','yQuality','limits','meshLimits'});
summary.tag = tag;
if ~isfield(metadata,'gridLabRegrids') || isempty(metadata.gridLabRegrids)
    metadata.gridLabRegrids = summary;
else
    existing = metadata.gridLabRegrids;
    if ~isequal(sort(fieldnames(existing)),sort(fieldnames(summary)))
        error('ipm:gridlab:RegridMetadata', ...
            'Existing grid-lab regrid provenance has an incompatible schema.');
    end
    metadata.gridLabRegrids(end+1,1) = summary;
end
metadata.latestGridLabRegridTag = string(tag);
end

function log = remove_terminal_record(log,storeSnapshots)
groups = fieldnames(log.history);
for groupIndex = 1:numel(groups)
    group = log.history.(groups{groupIndex});
    names = fieldnames(group);
    for nameIndex = 1:numel(names)
        values = group.(names{nameIndex});
        if ~isempty(values)
            group.(names{nameIndex}) = values(1:end-1);
        end
    end
    log.history.(groups{groupIndex}) = group;
end
if storeSnapshots
    log.snapshotRho(end) = [];
    log.snapshotNormalizedTime(end) = [];
    log.snapshotX(end) = [];
    log.snapshotY(end) = [];
end
end

function answer = same_endpoints(first,second)
first = double(first(:));
second = double(second(:));
scale = max([1,abs(first(1)),abs(first(end)), ...
    abs(second(1)),abs(second(end))]);
answer = max(abs([first(1)-second(1), ...
    first(end)-second(end)])) <= 128*eps(scale);
end
