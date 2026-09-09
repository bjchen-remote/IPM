function report = autonomousMeshConfig(oldReferenceFile)
%IPMTESTS.BASELINE.AUTONOMOUSMESHCONFIG Configuration-only, without operators.
if nargin < 1,oldReferenceFile = '';end
schema = ipm.config.schema();descriptor = schema.options.autonomousMesh;
assert(~descriptor.hasDefault && strcmp(descriptor.domain,'remesh') && strcmp(descriptor.kind,'autonomous_mesh'));
base = ipm.config.resolve(struct());assert(~isfield(base.remesh,'autonomousMesh'));
p = ipm.config.autonomousMeshPolicy(struct());
assert(p.version == 1 && p.enabled && isequal(p.targetCoreCells,[32,32]) && ...
    isequal(p.transactionMinimumCoreCells,[31,31]) && isequal(p.regridCoreTrigger,[26,26]) && ...
    isequal(p.predictedCoreBuffer,[22,22]) && isequal(p.endpointCoreFloor,[20,20]) && p.historyCoreFloor == 14);
targets = [32,32;42,42;56,56;32,56;256,256];
for k = 1:size(targets,1)
    target = targets(k,:);q = ipm.config.autonomousMeshPolicy(struct('targetCoreCells',target'));
    assert(isequal(q.targetCoreCells,target) && isequal(q.transactionMinimumCoreCells,target-1) && ...
        isequal(q.regridCoreTrigger,(26/32)*target) && isequal(q.predictedCoreBuffer,(22/32)*target) && ...
        isequal(q.endpointCoreFloor,(20/32)*target) && isequaln(q,ipm.config.autonomousMeshPolicy(q)));
end
q = orderfields(p,numel(fieldnames(p)):-1:1);q.timeUnit = "NATIVE_CANONICAL";q.enabled = 1;
q.qualityLimits = orderfields(q.qualityLimits);q.search = orderfields(q.search);
assert(isequaln(p,ipm.config.autonomousMeshPolicy(q)));
options = struct('rescalingMode','dynamic','dynamicScaleGeometry','isotropic', ...
    'symmetryMode','double_odd_omega','lengthGauge','transport_anchor', ...
    'cOmegaGauge','wall_omega_quadratic_peak','spatialDiscretization','high_order', ...
    'transportScheme','weno5_fd','timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'adaptiveRemesh',true,'initialAnalyticRemesh',true,'autonomousMesh',struct());
c = ipm.config.resolve(options);
assert(isequaln(c.remesh.autonomousMesh,p) && isequaln(ipm.config.resolve(flatten(c)),c) && ...
    isequaln(ipm.config.autonomousMeshPolicy(p,flatten(c)),p));
without = c;without.remesh = rmfield(without.remesh,'autonomousMesh');
assert(isequaln(without,ipm.config.resolve(rmfield(options,'autonomousMesh'))));
disabled = ipm.config.resolve(struct('autonomousMesh',struct('enabled',0)));
assert(~disabled.remesh.autonomousMesh.enabled);disabled.remesh = rmfield(disabled.remesh,'autonomousMesh');
assert(isequaln(disabled,base));
alternatives = {struct('symmetryMode','half_plane'),struct('rescalingMode','dynamic','dynamicScaleGeometry','anisotropic')};
for k = 1:numel(alternatives)
    other = alternatives{k};other.autonomousMesh = struct('enabled',false);
    otherResolved = ipm.config.resolve(other);assert(~otherResolved.remesh.autonomousMesh.enabled);
end
off = ipm.config.autonomousMeshPolicy(struct('enabled',false),struct());assert(~off.enabled);
id = 'ipm:BadAutonomousMeshPolicy';negativeCount = 0;
invalid = {[],true,repmat(struct(),1,2),struct('unknown',1),struct('enabled',2),struct('version',2), ...
    struct('targetCoreCells',32),struct('targetCoreCells',[31,32]),struct('targetCoreCells',[32,257]), ...
    struct('targetCoreCells',[32.5,42]),struct('targetCoreCells',[32,NaN]),struct('targetCoreCells',[32,Inf]), ...
    struct('timeUnit','parent_equivalent_canonical'),struct('qualityLimits',[]), ...
    struct('qualityLimits',struct('unknown',1)),struct('search',struct('unknown',1))};
for k = 1:numel(invalid)
    expect(@()ipm.config.autonomousMeshPolicy(invalid{k}),id);negativeCount = negativeCount+1;
end
derived = {'transactionMinimumCoreCells','regridCoreTrigger','predictedCoreBuffer','endpointCoreFloor'};
for k = 1:numel(derived)
    q = struct('targetCoreCells',[42,56]);q.(derived{k}) = p.(derived{k});
    expect(@()ipm.config.autonomousMeshPolicy(q),id);negativeCount = negativeCount+1;
end
for group = {'qualityLimits','search'}
    fields = fieldnames(p.(group{1}));
    for k = 1:numel(fields)
        q = p;name = fields{k};factor = 2;
        if strcmp(group{1},'qualityLimits') && startsWith(name,'min'),factor = .5;end
        q.(group{1}).(name) = factor*p.(group{1}).(name);
        expect(@()ipm.config.autonomousMeshPolicy(q),id);negativeCount = negativeCount+1;
    end
end
fixed = {'minimumFrontCells','historyCoreFloor','maximumSafety','maximumSinglePeakJump', ...
    'maximumCumulativeAbsolutePeakJump','maximumMassRelativeDefect','maximumRelativeRangeViolation', ...
    'trendWindow','maximumReviewInterval'};
for k = 1:numel(fixed)
    q = p;q.(fixed{k}) = 2*p.(fixed{k});expect(@()ipm.config.autonomousMeshPolicy(q),id);negativeCount = negativeCount+1;
end
unsupported = {'rescalingMode','physical';'dynamicScaleGeometry','anisotropic';'symmetryMode','half_plane'; ...
    'lengthGauge','local_strain';'cOmegaGauge','wall_omega_peak';'spatialDiscretization','legacy_second_order'; ...
    'transportScheme','muscl_minmod';'timeIntegrator','ssprk3';'remeshTransferScheme','pchip'; ...
    'adaptiveRemesh',false;'initialAnalyticRemesh',false;'adaptiveLevels',[.1,.5,.8]};
for k = 1:size(unsupported,1)
    other = options;other.(unsupported{k,1}) = unsupported{k,2};
    expect(@()ipm.config.resolve(other),'ipm:AutonomousMeshCombination');negativeCount = negativeCount+1;
end
expect(@()ipm.config.autonomousMeshPolicy(p,struct()),'ipm:AutonomousMeshCombination');negativeCount = negativeCount+1;
old = struct('supplied',false,'configurationsCompared',0,'checkpointRead',false,'operatorRestored',false);
if ~isempty(oldReferenceFile)
    assert(maxNumCompThreads == 10);
    r = load(oldReferenceFile,'inputs','configs','checkpointFile','oldSignature');
    for k = 1:numel(r.inputs)
        c = ipm.config.resolve(r.inputs{k});assert(isequaln(c,r.configs{k}) && ~isfield(c.remesh,'autonomousMesh'));
    end
    checkpoint = ipm.output.readCheckpoint(r.checkpointFile);c = checkpoint.payload.state.config;
    assert(isequaln(checkpoint.signature,r.oldSignature) && ~isfield(c.remesh,'autonomousMesh') && ...
        isequaln(ipm.config.resolve(flatten(c)),c));
    old = struct('supplied',true,'configurationsCompared',numel(r.inputs),'checkpointRead',true, ...
        'checkpointFile',r.checkpointFile,'step',checkpoint.payload.state.step,'signatureExactlyUnchanged',true, ...
        'storedConfigExactlyReresolved',true,'newFieldAbsent',true,'operatorRestored',false);
end
report = struct('allPassed',true,'optionalNoDefault',true,'targetsTested',targets,'normalizedRereadExact',true, ...
    'unrelatedConfigExact',true,'disabledAlternativeCombinations',3,'negativeCasesPassed',negativeCount, ...
    'oldEvidence',old,'configOnly',true,'productionAutonomousRuntimeCertified',false);
fprintf('AUTONOMOUS_MESH_CONFIG_PASS targets=%d negatives=%d oldConfigs=%d oldCheckpoint=%d noLU=1\n', ...
    size(targets,1),negativeCount,old.configurationsCompared,old.checkpointRead);
end
function flat = flatten(config)
flat = struct();schema = ipm.config.schema();names = schema.domainNames;
for k = 1:numel(names)
    fields = fieldnames(config.(names{k}));
    for j = 1:numel(fields),flat.(fields{j}) = config.(names{k}).(fields{j});end
end
end
function expect(action,id)
try
    action();
catch exception
    assert(strcmp(exception.identifier,id),'ipm:AutonomousMeshUnexpectedError', ...
        'Expected %s, received %s: %s',id,exception.identifier,exception.message);return;
end
error('ipm:AutonomousMeshMissingError','Expected rejection %s.',id);
end
