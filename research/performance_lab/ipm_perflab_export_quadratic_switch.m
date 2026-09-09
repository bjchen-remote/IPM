function evidence = ipm_perflab_export_quadratic_switch(referenceFile,outputDirectory)
%IPM_PERFLAB_EXPORT_QUADRATIC_SWITCH Replay native snapshots around a switch.
% No replacement peak functional or alternate PDE path. Frames are research
% evidence, not checkpoints. Require exact original terminal/history replay.
assert(~exist(outputDirectory,'dir'),'Pass a new evidence directory.');
mkdir(outputDirectory);
loaded = load(referenceFile,'reference');
reference = loaded.reference;
assert(reference.config.grid.nx<=65 && reference.config.grid.ny<=65, ...
    'Small-grid evidence only.');
schema = ipm.config.schema();
opts = struct();
for k=1:numel(schema.domainNames)
    domain = reference.config.(schema.domainNames{k});
    names = fieldnames(domain);
    for j=1:numel(names)
        opts.(names{j}) = domain.(names{j});
    end
end
opts.storeSnapshots = true;
opts.saveResults = false;
result = ipm.solve(opts);
assert(isequaln(result.state.rho,reference.state.rho) && ...
    isequaln(result.history,reference.history), ...
    'Snapshot replay changed the original numerical endpoint or history.');
assert(all(ipm.output.trustedMask(result)),'Replay must be entirely trusted.');
center = result.history.gauge.omegaGaugeQuadraticStencilCenterX;
switchRows = find(diff(center)~=0)+1;
assert(~isempty(switchRows),'No template switch in this reference.');
rows = unique([switchRows-1;switchRows;min(switchRows+1,numel(center))]);
state = ipm.evolve.initialize(opts);
frames = cell(numel(rows),1);
for k=1:numel(rows)
    row = rows(k);
    rho = result.snapshots.rho{row};
    [rhs,flow] = ipm.evolve.flow(rho,state.ops,state.scale);
    [baseRhs,~] = ipm.evolve.assembleRhs(rho,flow.u1,flow.u2, ...
        flow.c_l,flow.c_l,0,flow.c_r,state.ops,state.ops.transportBoundaryMode);
    [physicalRhs,~] = ipm.evolve.assembleRhs(rho,flow.u1,flow.u2, ...
        0,0,0,0,state.ops,state.ops.transportBoundaryMode);
    assert(isequaln(rhs,baseRhs+flow.c_omega*rho), ...
        'Stored base RHS must reproduce the maintained full RHS exactly.');
    P = flow.omegaGaugeQuadraticPeakValue;
    forcing = flow.omegaGaugeQuadraticForcing;
    assert(isequaln(P,result.history.gauge.omegaGaugeQuadraticPeakValue(row)), ...
        'Frame flow does not exactly replay the stored peak functional.');
    frames{k} = struct('row',row,'rho',rho,'rhoX',rho*state.ops.Dx.', ...
        'rhs',rhs,'baseRhs',baseRhs,'physicalRhs',physicalRhs,'flow',flow, ...
        'x',result.snapshots.x{row},'y',result.snapshots.y{row}, ...
        'canonicalTime',result.history.common.canonicalTau(row), ...
        'physicalTime',result.history.common.physicalTime(row), ...
        'center',center(row),'P',P,'forcing',forcing, ...
        'cOmega',flow.c_omega,'PPrime',forcing+flow.c_omega*P);
end
P = result.history.gauge.omegaGaugeQuadraticPeakValue;
forcing = result.history.gauge.omegaGaugeQuadraticForcing;
time = result.history.common.canonicalTau;
evidence = struct('schemaVersion',1,'kind','native_quadratic_peak_switch_frames', ...
    'sourceReferenceFile',referenceFile,'terminalRhoBitwise',true, ...
    'historyBitwise',true,'allTrusted',true,'frameRows',rows,'switchRows',switchRows, ...
    'switchTimes',time(switchRows),'centerBefore',center(switchRows-1), ...
    'centerAfter',center(switchRows),'PBefore',P(switchRows-1),'PAfter',P(switchRows), ...
    'PJump',P(switchRows)-P(switchRows-1), ...
    'jumpDividedByRecordedTimeInterval', ...
        (P(switchRows)-P(switchRows-1))./(time(switchRows)-time(switchRows-1)), ...
    'maximumRawPPrime',max(abs(forcing+result.history.common.c_omega.*P)), ...
    'interpretation','The algebraic within-stencil derivative excludes the finite functional jump.');
Dx = state.ops.Dx;
Dy = state.ops.Dy;
config = result.config;
save(fullfile(outputDirectory,'switch_frames.mat'),'frames','evidence','Dx','Dy','config','-v7.3');
save(fullfile(outputDirectory,'native_replay_with_snapshots.mat'),'result','-v7.3');
fid=fopen(fullfile(outputDirectory,'switch_evidence.json'),'w');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(evidence,PrettyPrint=true),'char');
end
