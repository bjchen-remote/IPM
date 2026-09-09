function report=ipm_perflab_low_order_preconditioned_poisson(out)
%IPM_PERFLAB_LOW_ORDER_PRECONDITIONED_POISSON Same high-order A, low-order M.
% Research only. No solver integration, PDE step, checkpoint or large matrix.
assert(maxNumCompThreads==10&&~isfile(fullfile(out,'report.json')));
clock=tic; rows={}; details={}; grids={};
for n=[49,65,129]
    for g=1:3
        H=8; stretch=0;
        if g>=2,stretch=(n-1)*log(1.07)/2;end
        if g==3,H=1e6;end
        [x,y]=axes_for(n,H,stretch);
        qx=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));
        qy=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
        assert(quality_pass(qx)&&quality_pass(qy),'perflab:GridQuality','Registered ordinary axis rejected.');
        config=ipm.config.resolve(struct('nx',n,'ny',n,'xlim',[-H,H],'ymax',H, ...
            'gridMode','uniform','spatialDiscretization','high_order', ...
            'transportScheme','weno5_fd','timeIntegrator','ssprk54', ...
            'remeshTransferScheme','high_order','rescalingMode','physical', ...
            'farBoundaryMode','green','transportBoundaryMode','open', ...
            'wallTransportMode','conservative_flux','adaptiveRemesh',false, ...
            'initialAnalyticRemesh',false,'saveResults',false,'makePlots',false, ...
            'livePlot',false,'writeVideo',false,'verbose',false));
        override=config.grid;override.customX=x;override.customY=y;
        timer=tic;ops=ipm.mesh.build(config,override);nativeBuild=toc(timer);
        [pre,preInfo]=make_preconditioner(x,y);
        A=ops.A;[gridInfo,luProxy]=matrix_info(A,preInfo);
        gridInfo.nodes=[n,n];gridInfo.H=H;gridInfo.stretch=stretch;
        gridInfo.productionAxisQualityPassed=true;gridInfo.xQuality=qx;gridInfo.yQuality=qy;
        gridInfo.nativeWholeBuildSeconds=nativeBuild;gridInfo.rcondNative=rcond(ops.poisson);
        grids{end+1}=gridInfo; %#ok<AGROW>
        [U,V]=meshgrid(linspace(0,1,n),linspace(0,1,n));
        known=sin(pi*U).*sin(pi*V)+.1*sin(7*pi*U).*sin(5*pi*V);
        known=known(2:end-1,2:end-1);b=A*known(:);
        ref=ops.poisson\b;boundary=zeros(n);
        [rr,dd]=evaluate(A,b,ref,known(:),boundary,ops,pre,'known_discrete',gridInfo);
        rows=[rows,rr];details=[details,dd]; %#ok<AGROW>
        exact=mms(ops.X,ops.Y,H);
        bd=struct('left',exact.psi(:,1),'right',exact.psi(:,end), ...
            'bottom',exact.psi(1,:),'top',exact.psi(end,:));
        [psi,b]=ipm.field.poisson(exact.source,ops,1,bd);
        [rr,dd]=evaluate(A,b(:),reshape(psi(2:end-1,2:end-1),[],1),[],psi,ops,pre,'nonzero_boundary_mms',gridInfo);
        for j=1:numel(rr)
            rr{j}.continuumPsiError=relative(dd{j}.psi,exact.psi);
            [u1,u2]=velocity(dd{j}.psi,ops);
            rr{j}.continuumVelocityError=max(relative(u1,exact.u1),relative(u2,exact.u2));
        end
        rows=[rows,rr];details=[details,dd]; %#ok<AGROW>
        rho=-exp(-(ops.X/H).^2-(ops.Y/H).^2);source=rho*ops.Dx';
        [psi,b]=ipm.field.poisson(source,ops);
        [rr,dd]=evaluate(A,b(:),reshape(psi(2:end-1,2:end-1),[],1),[],psi,ops,pre,'unchanged_image_green',gridInfo);
        rows=[rows,rr];details=[details,dd]; %#ok<AGROW>
        savedOps=rmfield(ops,'poisson');
        save(fullfile(out,sprintf('grid_%d_%d.mat',n,g)),'gridInfo','luProxy','savedOps','pre','A','-v7.3');
        clear savedOps;
        clear ops pre A luProxy;
        assert(toc(clock)<300,'perflab:TimeBudget','Small-screen time budget exceeded.');
    end
end
% Deliberately inadmissible stress axes, never mislabeled as valid PDE grids.
for stretch=[12,16]
    n=65;H=1e6;[x,y]=axes_for(n,H,stretch);
    xx=-ipm.mesh.fdMatrix(x,2,7);yy=-ipm.mesh.fdMatrix(y,2,7);
    A=kron(speye(n-2),yy(2:end-1,2:end-1))+kron(xx(2:end-1,2:end-1),speye(n-2));
    timer=tic;native=decomposition(A,'lu');nativeBuild=toc(timer);
    [pre,preInfo]=make_preconditioner(x,y);[gridInfo,luProxy]=matrix_info(A,preInfo);
    gridInfo.nodes=[n,n];gridInfo.H=H;gridInfo.stretch=stretch;
    gridInfo.xQuality=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));
    gridInfo.yQuality=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
    gridInfo.productionAxisQualityPassed=quality_pass(gridInfo.xQuality)&&quality_pass(gridInfo.yQuality);
    assert(~gridInfo.productionAxisQualityPassed);
    gridInfo.nativeWholeBuildSeconds=nativeBuild;gridInfo.rcondNative=rcond(native);
    grids{end+1}=gridInfo; %#ok<AGROW>
    [U,V]=meshgrid(linspace(0,1,n),linspace(0,1,n));
    known=sin(pi*U).*sin(pi*V)+.1*sin(7*pi*U).*sin(5*pi*V);
    known=known(2:end-1,2:end-1);b=A*known(:);ref=native\b;
    ops=struct('nx',n,'ny',n,'Dx',ipm.mesh.fdMatrix(x,1,7), ...
        'Dy',ipm.mesh.fdMatrix(y,1,7),'poisson',native);
    [rr,dd]=evaluate(A,b,ref,known(:),zeros(n),ops,pre,'inadmissible_raw_axis_known_discrete',gridInfo);
    rows=[rows,rr];details=[details,dd]; %#ok<AGROW>
    save(fullfile(out,sprintf('raw_%d.mat',stretch)),'gridInfo','luProxy','x','y','A','b','known','-v7.3');
    clear ops pre native A luProxy;
end
report=struct('kind','same_high_order_A_low_order_preconditioned_gmres_v1', ...
    'grids',{grids},'rows',{rows},'wallSeconds',toc(clock), ...
    'largestNodes',[129,129],'PDEExecuted',false,'checkpointReadOrWritten',false, ...
    'bitwiseReplacementClaim',false,'productionEligible',false, ...
    'nativeBoundaryAndOriginalAUnchanged',true,'maximumCorrectionSolves',2, ...
    'gmresRestart',30,'gmresTolerance',1e-12,'gmresMaximumOuterIterations',10);
save(fullfile(out,'report.mat'),'report','details','-v7.3');
write_json(fullfile(out,'report.json'),report);
end
function [x,y]=axes_for(n,H,s)
x=linspace(-1,1,n);y=linspace(0,1,n)';
if s>0,x=sinh(s*x)/sinh(s);y=sinh(s*y)/sinh(s);end
x=H*x;y=H*y;
end
function yes=quality_pass(q)
yes=q.quadratureWeightsStrictlyPositive&&q.maximumAdjacentCellRatio<=1.08&& ...
 q.maximumLogSpacingCurvature<=.01&&q.minimumStencilRcond>=1e-9&& ...
 q.minimumQuadratureWeightRatio>=1e-8&&q.minimumQuadratureWeightToControlWidthRatio>=.35&& ...
 q.maximumQuadratureWeightToControlWidthRatio<=1.65;
end
function [p,info]=make_preconditioner(x,y)
timer=tic;Mx=low_order(x);My=low_order(y);nx=numel(x)-2;ny=numel(y)-2;
M=kron(speye(nx),My)+kron(Mx,speye(ny));
hx=(x(3:end)-x(1:end-2))/2;hy=(y(3:end)-y(1:end-2))/2;w=kron(hx',hy);
Q=spdiags(w,0,numel(w),numel(w))*M;Q=(Q+Q')/2;
[L,flag,perm]=chol(Q,'lower','vector');assert(flag==0,'perflab:PreconditionerSPD','Low-order Cholesky failed.');
p=struct('L',L,'permutation',perm,'weights',w);
item=whos('p');info=struct('setupSeconds',toc(timer),'lowOrderNnz',nnz(M), ...
 'weightedMatrixNnz',nnz(Q),'factorNnz',nnz(L),'cacheBytes',item.bytes);
end
function M=low_order(z)
z=z(:);n=numel(z)-2;hL=diff(z(1:end-1));hR=diff(z(2:end));
a=-2./(hL.*(hL+hR));d=2./(hL.*hR);c=-2./(hR.*(hL+hR));
M=sparse([(2:n)';(1:n)';(1:n-1)'],[(1:n-1)';(1:n)';(2:n)'],[a(2:end);d;c(1:end-1)],n,n);
end
function z=apply_preconditioner(p,r)
v=p.weights.*r;z=zeros(size(r));z(p.permutation)=p.L'\(p.L\v(p.permutation));
end
function [x,info]=candidate(A,b,p,corrections)
timer=tic;calls=0;[x,flag,relres,iter,resvec]=gmres(A,b,30,1e-12,10,@pre);
records=struct('flag',flag,'reportedPreconditionedRelativeResidual',relres, ...
 'iterations',iter,'residualHistory',resvec);
for j=1:corrections
    residual=b-A*x;
    if ~all(isfinite(x))||~any(residual),break;end
    [delta,flag,relres,iter,resvec]=gmres(A,residual,30,1e-12,10,@pre);
    x=x+delta;
    records(end+1)=struct('flag',flag,'reportedPreconditionedRelativeResidual',relres, ...
        'iterations',iter,'residualHistory',resvec); %#ok<AGROW>
end
info=struct('seconds',toc(timer),'preconditionerApplications',calls,'solves',records);
    function v=pre(r)
        calls=calls+1;v=apply_preconditioner(p,r);
    end
end
function [rows,data]=evaluate(A,b,ref,known,boundary,ops,pre,name,g)
rows=cell(1,3);data=cell(1,3);solutions=cell(1,3);infos=cell(1,3);
timer=tic;solutions{1}=ops.poisson\b;assert(isequaln(solutions{1},ref));
infos{1}=struct('seconds',toc(timer),'preconditionerApplications',0,'solves',struct([]));
[solutions{2},infos{2}]=candidate(A,b,pre,0);[solutions{3},infos{3}]=candidate(A,b,pre,2);
refPsi=boundary;refPsi(2:end-1,2:end-1)=reshape(ref,ops.ny-2,ops.nx-2);
[refU1,refU2]=velocity(refPsi,ops);methods={'native_lu','low_order_gmres','low_order_gmres_plus_two_defects'};
for j=1:3
    v=solutions{j};psi=boundary;psi(2:end-1,2:end-1)=reshape(v,ops.ny-2,ops.nx-2);
    [u1,u2]=velocity(psi,ops);r=A*v-b;rowMax=full(max(abs(A),[],2));
    rr=struct('problem',name,'method',methods{j},'nodes',g.nodes,'H',g.H,'stretch',g.stretch, ...
        'productionAxisQualityPassed',g.productionAxisQualityPassed,'finite',all(isfinite(v)), ...
        'originalGlobalRelativeResidual',norm(r,inf)/max(norm(b,inf),realmin), ...
        'componentwiseBackwardError',max(abs(r)./max(abs(A)*abs(v)+abs(b),realmin)), ...
        'rowNormalizedResidual',norm(r./rowMax,inf)/max(norm(v,inf),realmin), ...
        'relativePsiDifferenceToNative',relative(v,ref), ...
        'relativeVelocityDifferenceToNative',max(relative(u1,refU1),relative(u2,refU2)), ...
        'knownDiscreteError',NaN,'nativeKnownDiscreteError',NaN,'nativeOracleReliable',true, ...
        'boundaryExact',isequaln(psi([1,end],:),boundary([1,end],:))&&isequaln(psi(:,[1,end]),boundary(:,[1,end])), ...
        'timingAndConvergence',infos{j},'bitwiseNative',isequaln(v,ref));
    if ~isempty(known)
        rr.knownDiscreteError=relative(v,known);rr.nativeKnownDiscreteError=relative(ref,known);
        rr.nativeOracleReliable=rr.nativeKnownDiscreteError<=1e-8;
    end
    rr.screenPassed=rr.finite&&rr.boundaryExact&&rr.nativeOracleReliable&& ...
        rr.componentwiseBackwardError<=1e-12&&rr.rowNormalizedResidual<=1e-10&& ...
        rr.relativePsiDifferenceToNative<=1e-9&&rr.relativeVelocityDifferenceToNative<=1e-8&& ...
        (isempty(known)||rr.knownDiscreteError<=1e-8);
    rows{j}=rr;data{j}=struct('solution',v,'psi',psi,'native',ref,'b',b,'known',known);
    fprintf('LOW_ORDER_PC %s n%d H%.0g s%.4g %s pass%d %.4gsec back%.3g dx%.3g du%.3g\n', ...
        name,g.nodes(1),g.H,g.stretch,methods{j},rr.screenPassed,infos{j}.seconds, ...
        rr.componentwiseBackwardError,rr.relativePsiDifferenceToNative,rr.relativeVelocityDifferenceToNative);
end
end
function [g,proxy]=matrix_info(A,p)
% Public LU is only a fill/storage proxy; never used to solve or as native
% decomposition memory. Its work is not included in either candidate setup.
timer=tic;[L,U,pr,pc,D]=lu(A,'vector');proxy=struct('factorNnz',[nnz(L),nnz(U)], ...
 'L',L,'U',U,'rowPermutation',pr,'columnPermutation',pc,'D',D);
i=whos('proxy');rowMax=full(max(abs(A),[],2));
g=struct('unknowns',size(A,1),'originalNnz',nnz(A),'rowScaleOrders',log10(max(rowMax)/min(rowMax)), ...
 'preconditioner',p,'publicLUFactorNnz',[nnz(L),nnz(U)],'publicLUCacheProxyBytes',i.bytes, ...
 'publicLUProxySeconds',toc(timer),'publicLUProxyIsNativeCacheMeasurement',false, ...
 'gmresThirtyVectorLowerBoundBytes',8*size(A,1)*31);
end
function e=mms(X,Y,H)
u=(X+H)/(2*H);v=Y/H;s=sin(pi*u).*sin(pi*v);c=.03*cos(2*pi*u).*(1+v);
e=struct('psi',s+c,'source',1.25*(pi/H)^2*s+(pi/H)^2*c, ...
 'u1',-pi/H*sin(pi*u).*cos(pi*v)-.03/H*cos(2*pi*u), ...
 'u2',pi/(2*H)*cos(pi*u).*sin(pi*v)-.03*pi/H*sin(2*pi*u).*(1+v));
end
function [u1,u2]=velocity(psi,ops)
u1=-(ops.Dy*psi);u2=psi*ops.Dx';
end
function r=relative(a,b)
r=max(abs(a-b),[],'all')/max(max(abs(b),[],'all'),realmin);
end
function write_json(file,value)
f=fopen(file,'w');assert(f>=0);c=onCleanup(@()fclose(f));fwrite(f,jsonencode(value,PrettyPrint=true),'char');
end
