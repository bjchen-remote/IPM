function report=mesh_probe_autonomous_runtime(outDir)
%MESH_PROBE_AUTONOMOUS_RUNTIME Registered bounded solver-integration probes.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
opts=struct('nx',321,'ny',161,'xlim',[-8,8],'ymax',4, ...
    'gridMode','uniform','gridStretchAutomatic',false,'gridStretch',[0,0], ...
    'initialCondition','degenerate_primitive','degeneratePower',8, ...
    'rescalingMode','dynamic','dynamicScaleGeometry','isotropic', ...
    'symmetryMode','double_odd_omega','lengthGauge','transport_anchor', ...
    'transportAnchorX',1,'cOmegaGauge','wall_omega_quadratic_peak', ...
    'spatialDiscretization','high_order','transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54','remeshTransferScheme','high_order', ...
    'wallTransportMode','conservative_flux','transportBoundaryMode','open', ...
    'adaptiveRemesh',true,'initialAnalyticRemesh',true,'autonomousMesh',struct(), ...
    'finalTime',4e-5,'maxDt',1e-5,'minDt',1e-12,'maxSteps',4, ...
    'outputEvery',1e-5,'storeSnapshots',true,'saveResults',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false,'verbose',false);
report=struct('kind','bounded_autonomous_solver_integration_probe_v1', ...
    'longtimeQualification',false,'tests',struct([]));
save(fullfile(outDir,'registration.mat'),'opts','report');
for k=1:2
    o=opts;
    if k==1,label='degenerate_primitive';
    else
        label='smooth_broad_fixture';
        o.initialCondition=@(X,Y)(1-exp(-(X/1.5).^4)).*exp(-Y.^4);
    end
    row=struct('label',label,'passed',false,'errorIdentifier','','errorMessage','', ...
        'errorStack',struct([]),'seconds',NaN,'coreCells',[NaN,NaN],'steps',NaN,'stopReason','','initialization',struct());
    timer=tic;
    try
        result=ipm.solve(o);
        row.seconds=toc(timer);row.steps=result.state.steps;row.stopReason=result.state.stopReason;
        m=result.metadata.autonomousMesh;
        row.initialization=m.initialization;row.coreCells=m.window.coreCells(end,:);
        assert(m.initialization.originalPhysicalTime==0 && m.initialization.originalCanonicalTime==0 && ...
            m.initialization.originalStep==0 && ~m.initialization.newPhysicalEpochCreated);
        assert(result.history.common.physicalTime(1)==0 && result.history.common.canonicalTau(1)==0);
        assert(isequal(result.snapshots.rho{1},o_sample(o,result.snapshots.x{1},result.snapshots.y{1})));
        row.passed=result.state.steps==4 && strcmp(result.state.stopReason,'final_time');
        save(fullfile(outDir,[label '.mat']),'result','o','-v7.3');
    catch e
        row.seconds=toc(timer);row.errorIdentifier=e.identifier;row.errorMessage=e.message;row.errorStack=e.stack;
        save(fullfile(outDir,[label '_failure.mat']),'e','o');
    end
    if isempty(report.tests),report.tests=row;else,report.tests(end+1)=row;end
    write_json(fullfile(outDir,'report.json'),report);
    fprintf('AUTONOMOUS_RUNTIME_PROBE %s passed=%d steps=%g reason=%s error=%s seconds=%.3f\n', ...
        label,row.passed,row.steps,row.stopReason,row.errorIdentifier,row.seconds);
    clear result m
end
report.allPassed=all([report.tests.passed]);write_json(fullfile(outDir,'report.json'),report);
save(fullfile(outDir,'report.mat'),'report');
end
function rho=o_sample(opts,x,y)
[X,Y]=meshgrid(x,y);
physics=struct('initialCondition',opts.initialCondition,'degeneratePower',opts.degeneratePower);
ops=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y),'symmetryMode',opts.symmetryMode);
rho=ipm.field.initialDensity(ops,physics);
end
function write_json(file,data)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data,PrettyPrint=true));
end
