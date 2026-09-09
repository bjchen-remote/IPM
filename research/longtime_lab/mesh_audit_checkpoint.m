function audit = mesh_audit_checkpoint(checkpointFile,outputFile)
%MESH_AUDIT_CHECKPOINT Read-only accepted-state/mesh audit; advances no PDE.
% Run with the ipm_structured root and this directory on MATLAB's path.
% The optional JSON report is a new research artifact, never a restart.
if nargin < 2
    outputFile = '';
end
checkpoint = ipm.output.readCheckpoint(checkpointFile);
state = checkpoint.payload.state;
history = checkpoint.payload.log.history;
trusted = ipm.output.trustedMask(history,state.config);
prefix = ipm.output.continuousTrustedPrefix(trusted);
audit = struct('checkpointFile',char(checkpointFile), ...
    'schemaVersion',checkpoint.schemaVersion, ...
    'configSchemaVersion',state.config.schemaVersion, ...
    'steps',state.step,'canonicalTime',state.scale.canonicalTime, ...
    'physicalTime',state.scale.physicalTime, ...
    'nx',numel(state.x),'ny',numel(state.y), ...
    'lengthGauge',state.config.scaling.lengthGauge, ...
    'cOmegaGauge',state.config.scaling.cOmegaGauge, ...
    'frozenTime',state.config.time, ...
    'continuousTrustedPrefix',all(prefix), ...
    'trustedRecords',nnz(prefix),'historyRecords',numel(prefix));
names = {'physicalRhoXInf','rhoXInf','dt','dtPhysical','c_l', ...
    'c_omega','canonicalTau','physicalTime','timeSpeed'};
audit.common = terminal_fields(history.common,names);
audit.mesh = terminal_fields(history.mesh,fieldnames(history.mesh));
audit.gauge = terminal_fields(history.gauge,fieldnames(history.gauge));
audit.x = ipm.mesh.quality(state.x,7,ipm.mesh.quadrature(state.x));
audit.y = ipm.mesh.quality(state.y,7,ipm.mesh.quadrature(state.y));
tau = history.common.canonicalTau(:);
tail = tau >= tau(end)-0.35;
audit.tailWindow = [tau(find(tail,1)),tau(end)];
audit.tailCore = struct();
for item = {'coreGridPoints','verticalCoreGridPoints','safetyFactor'}
    name = item{1};
    values = history.mesh.(name)(:);
    fit = polyfit(tau(tail)-tau(end),values(tail),1);
    audit.tailCore.(name) = struct('minimum',min(values(tail)), ...
        'maximum',max(values(tail)),'terminal',values(end), ...
        'linearSlopePerTau',fit(1));
end
fprintf('MESH_AUDIT_JSON %s\n',jsonencode(audit));
if ~isempty(outputFile)
    if exist(outputFile,'file')
        error('ipm:LongtimeAuditExistingOutput', ...
            'Refusing to replace an existing audit report.');
    end
    fid = fopen(outputFile,'w');
    if fid < 0
        error('ipm:LongtimeAuditOutput','Cannot open audit report.');
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(audit,PrettyPrint=true));
end
end

function selected = terminal_fields(group,names)
selected = struct();
for index = 1:numel(names)
    name = names{index};
    if isfield(group,name) && isnumeric(group.(name)) && ...
            ~isempty(group.(name))
        selected.(name) = group.(name)(end);
    end
end
end
