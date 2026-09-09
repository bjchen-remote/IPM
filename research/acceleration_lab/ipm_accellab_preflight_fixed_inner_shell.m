function report=ipm_accellab_preflight_fixed_inner_shell(initialAuditFile,out)
%IPM_ACCELLAB_PREFLIGHT_FIXED_INNER_SHELL Pure axes/analytic initial fields.
% Uses audited A initialization arrays, not a migrated state or fresh epoch.
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
d=load(initialAuditFile,'data');a=d.data{1};config=a.config;
x=a.initialization.selectedBaseX;y=a.initialization.selectedBaseY;
assert(isequal(x([1,end]),[-8,8]) && isequal(y([1,end]),[0;4]));
xp=continue_axis(x(x>=0),16);yp=continue_axis(y',8)';
xx=[-fliplr(xp(2:end)),xp];yy=yp;
[foundX,ix]=ismember(x,xx);[foundY,iy]=ismember(y,yy);
assert(all(foundX)&&all(foundY));
registration=struct('kind','A_axes_fixed_inner_H8_H16_shell_geometry_v1', ...
    'sourceAuditedInitialization',initialAuditFile,'innerBox',[8,4],'outerBox',[16,8], ...
    'allOriginalAxesPreserved',true,'analyticDatumIndependentlySampled',true, ...
    'noAdaptiveInitialReselection',true,'noPDE',true,'noLU',true, ...
    'construction','n=ceil(extension/last_spacing); log h_i=log h0+s0*i+c*(i/n)^2; unique scalar c matches endpoint; no quality-threshold change.', ...
    'productionAdmissionClaim',false,'greenIntegralEvaluated',false);
write_json(fullfile(out,'registration.json'),registration);
profile clear;profile on;
Dx=ipm.mesh.fdMatrix(x,1,7);Dy=ipm.mesh.fdMatrix(y,1,7);
DX=ipm.mesh.fdMatrix(xx,1,7);DY=ipm.mesh.fdMatrix(yy,1,7);
wx=ipm.mesh.quadrature(x);wy=ipm.mesh.quadrature(y);
WX=ipm.mesh.quadrature(xx);WY=ipm.mesh.quadrature(yy);
qx=ipm.mesh.quality(xx,7,WX);qy=ipm.mesh.quality(yy,7,WY);
[X,Y]=meshgrid(x,y);[XX,YY]=meshgrid(xx,yy);
rho=ipm.field.initialDensity(struct('X',X,'Y',Y,'symmetryMode',config.physics.symmetryMode),config.physics);
RHO=ipm.field.initialDensity(struct('X',XX,'Y',YY,'symmetryMode',config.physics.symmetryMode),config.physics);
omega=rho*Dx';OMEGA=RHO*DX';
cx=find(x>=-2 & x<=2);cy=find(y>=0 & y<=2);
derivativeRowsExact=isequaln(Dx(cx,:),DX(ix(cx),ix)) && isequaln(Dy(cy,:),DY(iy(cy),iy));
ex=setdiff(1:numel(xx),ix);ey=setdiff(1:numel(yy),iy);
outsideStencilZero=nnz(DX(ix(cx),ex))==0 && nnz(DY(iy(cy),ey))==0;
quadratureInnerExact=isequaln(wx(cx),WX(ix(cx))) && isequaln(wy(cy),WY(iy(cy)));
analyticCommonExact=isequaln(rho,RHO(iy,ix));
pairedInner=OMEGA(iy(cy),ix(cx));oldInner=omega(cy,cx);
report=registration;report.nodes=[numel(xx),numel(yy)];report.totalNodes=numel(xx)*numel(yy);
report.originalResourceCap=config.remesh.autonomousMesh.nodeFamily.maximumTotalNodes;
report.originalResourceCapPassed=report.totalNodes<=report.originalResourceCap;
report.xQuality=qx;report.yQuality=qy;
report.qualityReasons=[quality_reasons(qx,config.remesh.autonomousMesh.qualityLimits,'x'), ...
    quality_reasons(qy,config.remesh.autonomousMesh.qualityLimits,'y')];
report.qualityPassed=isempty(report.qualityReasons);
report.commonInnerBox=[-2,2,0,2];report.derivativeRowsExact=derivativeRowsExact;
report.outsideStencilZero=outsideStencilZero;report.quadratureInnerExact=quadratureInnerExact;
report.analyticCommonExact=analyticCommonExact;report.pairedInnerOmegaExact=isequaln(pairedInner,oldInner);
report.pairedInnerOmegaRelativeInf=max(abs(pairedInner-oldInner),[],'all')/max(abs(oldInner),[],'all');
assert(derivativeRowsExact&&outsideStencilZero&&quadratureInnerExact&&analyticCommonExact);
report.diagnosticGeometryReady=report.qualityPassed;
report.productionReady=false;
profile off;p=profile('info');names={p.FunctionTable.FunctionName};
for key={'ipm.mesh.build','ipm.field.poisson','ipm.field.velocity','ipm.evolve.flow','ipm.solve','ipm.evolve.advance'}
    assert(~any(strcmp(names,key{1})));
end
save(fullfile(out,'geometry.mat'),'report','x','y','xx','yy','ix','iy','config','p','-v7.3');
write_json(fullfile(out,'report.json'),report);disp(jsonencode(report));
end

function extended=continue_axis(axis,target)
axis=axis(:)';h=diff(axis);h0=h(end);s0=log(h(end)/h(end-1));
n=ceil((target-axis(end))/h0);j=1:n;f=(j/n).^2;gap=target-axis(end);
value=@(c)sum(exp(log(h0)+s0*j+c*f));
lower=-100;upper=100;assert(value(lower)<gap && value(upper)>gap);
for k=1:80
    middle=(lower+upper)/2;
    if value(middle)>gap,upper=middle;else,lower=middle;end
end
steps=exp(log(h0)+s0*j+((lower+upper)/2)*f);
added=axis(end)+cumsum(steps);added(end)=target;
extended=[axis,added];assert(all(diff(extended)>0));
end

function reasons=quality_reasons(q,p,axis)
bad=[q.maximumAdjacentCellRatio>p.maxAdjacentCellRatio,q.maximumLogSpacingCurvature>p.maxLogSpacingCurvature, ...
    q.minimumStencilRcond<p.minStencilRcond,q.minimumQuadratureWeightRatio<p.minQuadratureWeightRatio, ...
    q.minimumQuadratureWeightToControlWidthRatio<p.minWeightToControlWidth, ...
    q.maximumQuadratureWeightToControlWidthRatio>p.maxWeightToControlWidth,~q.quadratureWeightsStrictlyPositive];
names={'adjacent','curvature','rcond','weight_ratio','weight_control_min','weight_control_max','positive_weights'};
reasons=cellfun(@(s)[axis,':',s],names(bad),'UniformOutput',false);
end

function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
