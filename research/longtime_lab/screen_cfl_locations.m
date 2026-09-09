function report = screen_cfl_locations(resultFiles,outputFile)
%SCREEN_CFL_LOCATIONS Locate the existing CFL constraint without a PDE/LU.
% Reconstruct the exact fourth-order control widths and rescaled velocities
% from paired v2 endpoint data. Local CFL values are diagnostics, not a
% proposal to relax the global explicit stability condition.
assert(iscell(resultFiles) && ~isempty(resultFiles));
assert(~isfile(outputFile),'Refusing to overwrite an existing screen.');
records=cell(1,numel(resultFiles));
for k=1:numel(resultFiles)
    r=ipm.output.validate(resultFiles{k});
    assert(all(ipm.output.trustedMask(r)));
    assert(strcmp(r.config.transport.spatialDiscretization,'high_order'));
    x=r.grid.x(:)'; y=r.grid.y(:);
    hx=zeros(size(x)); hy=zeros(size(y));
    hx([1,end])=[x(2)-x(1),x(end)-x(end-1)]/2;
    hy([1,end])=[y(2)-y(1);y(end)-y(end-1)]/2;
    hx(2:end-1)=(x(3:end)-x(1:end-2))/2;
    hy(2:end-1)=(y(3:end)-y(1:end-2))/2;
    rates=r.scale.rates;
    v1=r.state.velocity.x+rates.cx*x+rates.cr;
    v2=r.state.velocity.y+rates.cy*y;
    rx=abs(v1)./hx; ry=abs(v2)./hy; total=rx+ry;
    % Check assembly against the maintained selector; its result does not
    % alter an accepted step or the solver's history.
    ops=struct('spatialDiscretization','high_order','hx',hx,'hy',hy);
    flow=struct('transportU1',v1,'transportU2',v2);
    scale=struct('logC_omega',log(r.scale.Comega),'logC_l',log(r.scale.Cx), ...
        'canonicalTime',r.state.canonicalTime,'physicalTime',r.state.physicalTime);
    [~,~,selection]=ipm.evolve.selectTimestep(flow,scale,ops,r.config.time);
    assert(isequal(selection.rate,max(total,[],'all')));
    g=r.history.gauge; m=r.history.mesh;
    a=g.omegaGaugeQuadraticPeakX(end);
    wx=m.trackedWallCoreWidth(end); wy=m.trackedVerticalCoreWidth(end);
    [X,Y]=meshgrid(x,y);
    near=(abs(abs(X)-a)<=3*wx) & Y<=3*wy;
    lab=abs(X)<=2 & Y<=2;
    sourceMaximum=max(abs(r.state.omega),[],'all');
    [sorted,indices]=maxk(total(:),30);
    [iy,ix]=ind2sub(size(total),indices);
    cells=struct('x',x(ix)','y',y(iy),'rate',sorted,'xPart',rx(indices), ...
        'yPart',ry(indices),'v1',v1(indices),'v2',v2(indices), ...
        'sourceRelativeToPeak',abs(r.state.omega(indices))/sourceMaximum, ...
        'horizontalCoreDistance',abs(abs(x(ix)')-a)/wx, ...
        'verticalCoreDistance',y(iy)/wy);
    records{k}=struct('resultFile',resultFiles{k},'tau',r.state.canonicalTime, ...
        'physicalTime',r.state.physicalTime,'globalRate',selection.rate, ...
        'globalCflDt',selection.dtCfl,'maxDt',r.config.time.maxDt, ...
        'strictCoreWidths',[wx,wy],'peakX',a, ...
        'nearPeakRate',max(total(near)),'nearPeakCflDt',r.config.time.cfl/max(total(near)), ...
        'wallRate',max(total(1,:)),'wallCflDt',r.config.time.cfl/max(total(1,:)), ...
        'labRate',max(total(lab)),'labCflDt',r.config.time.cfl/max(total(lab)), ...
        'xOnlyRate',max(rx,[],'all'),'yOnlyRate',max(ry,[],'all'), ...
        'topCells',cells,'maintainedSelectorRateExact',true, ...
        'recordedAcceptedStepRate',r.history.common.transportRate(end));
end
report=struct('kind','endpoint_cfl_location_screen','pdeAdvanced',false, ...
    'operatorBuilt',false,'records',[records{:}], ...
    'interpretation',['All local rates are descriptive. An explicit step must still satisfy ' ...
    'the global constraint. Recorded accepted-step rate belongs to its step selection, ' ...
    'whereas this screen recomputes the saved endpoint rate.']);
save(outputFile,'report');
[folder,name]=fileparts(outputFile);
fid=fopen(fullfile(folder,[name,'.json']),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
