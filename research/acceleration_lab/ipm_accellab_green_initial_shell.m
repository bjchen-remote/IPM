function report=ipm_accellab_green_initial_shell(geometryFile,initialReportFile,out)
%IPM_ACCELLAB_GREEN_INITIAL_SHELL Analytic t0 Green integrals, never a flow.
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
g=load(geometryFile,'xx','yy');native=jsondecode(fileread(initialReportFile));
registration=struct('kind','analytic_k8_t0_anchor_Green_shell_v1', ...
    'geometryFile',geometryFile,'nativeInitialReport',initialReportFile, ...
    'anchor',1,'boxes',[8,4;16,8],'source','xi^7/(1+xi^8+eta^8)', ...
    'kernel','4*xi*eta/(pi*((1-xi)^2+eta^2)*((1+xi)^2+eta^2))', ...
    'cartesianOrders',[2,4,8],'polarOrders',[16,32,64,128], ...
    'relativeShellAgreementTolerance',1e-9,'absoluteIntegralAgreementTolerance',1e-10, ...
    'physicalTime',0,'noFDVelocityTruthClaim',true,'noPDE',true,'noLU',true, ...
    'notApplicableToUnknownEvolvedOuterField',true);
write_json(fullfile(out,'registration.json'),registration);
timer=tic;
xs=g.xx(g.xx>=0);ys=g.yy;
cart=zeros(1,3);
for k=1:3
    n=registration.cartesianOrders(k);
    cart(k)=rectangle(xs(xs>=8),ys,n)+rectangle(xs(xs<=8),ys(ys>=4),n);
end
orders=registration.polarOrders;shell=zeros(size(orders));tail=zeros(numel(orders),2);
inside=tail;leading=zeros(size(orders));momentError=zeros(size(orders));
for k=1:numel(orders)
    n=orders(k);[z,w]=rule(n);
    momentError(k)=max(abs([sum(w)-2,sum(w.*z),sum(w.*z.^2)-2/3,sum(w.*z.^6)-2/7]));
    shell(k)=origin_polar(8,n,'shell');
    tail(k,:)=[origin_polar(8,n,'tail'),origin_polar(16,n,'tail')];
    inside(k,:)=[anchor_polar(8,n),anchor_polar(16,n)];
    [theta,wt]=axis_rule([0,atan(.5),pi/2],n);
    cs=cos(theta);sn=sin(theta);
    leading(k)=sum(wt.*(4*cs.^8.*sn./(pi*(cs.^8+sn.^8))).*max(cs,2*sn));
end
report=registration;
report.cartesianShell=cart;report.polarShell=shell;
report.analyticInsideCL=inside;report.analyticTailCL=tail;
report.infiniteSourceCLFromInsidePlusTail=inside+tail;
report.leadingTailCoefficient=leading;report.leadingShellApproximation=leading/16;
report.gaussMomentError=momentError;report.wallSeconds=toc(timer);
report.lastCartesianChange=abs(cart(end)-cart(end-1));
report.lastPolarShellChange=abs(shell(end)-shell(end-1));
report.coordinateShellDifference=abs(cart(end)-shell(end));
report.lastInsideChange=abs(inside(end,:)-inside(end-1,:));
report.lastTailChange=abs(tail(end,:)-tail(end-1,:));
report.shellVsInsideDifference=abs(shell(end)-(inside(end,2)-inside(end,1)));
report.shellVsTailDifference=abs(shell(end)-(tail(end,1)-tail(end,2)));
relativeBudget=registration.relativeShellAgreementTolerance*abs(shell(end));
report.quadratureAgreementPassed=max([report.lastCartesianChange,report.lastPolarShellChange, ...
    report.coordinateShellDifference,report.shellVsTailDifference])<=relativeBudget && ...
    max([report.lastInsideChange,report.lastTailChange,report.shellVsInsideDifference])<=1e-10 && ...
    max(momentError)<=5e-13;
cn=[native.cases(1).initialObservation.physicalRateCL,native.cases(2).initialObservation.physicalRateCL];
report.nativeInitialCL=cn;
report.nativeBMinusA=cn(2)-cn(1);
report.analyticShellFractionOfNativeDifference=shell(end)/report.nativeBMinusA;
report.nativeMinusAnalyticInside=cn-inside(end,:);
report.differenceNotExplainedByAnalyticShell=report.nativeBMinusA-shell(end);
report.decompositionMeaning='Algebraic comparison with exact-source continuum Green integrals, not a pure FD/BC error attribution. Native gaps include derivative, quadrature, cutoff/compression, boundary and elliptic discretization, and different selected grids.';
report.integrationApproved=false;report.atTimeGreaterThanZeroValidated=false;
save(fullfile(out,'report.mat'),'report');write_json(fullfile(out,'report.json'),report);
disp(jsonencode(report));
end

function v=rectangle(xs,ys,n)
[qx,wx]=axis_rule(xs,n);[z,w]=rule(n);v=0;
% Stream one y cell, bounding temporary memory independently of grid size.
for j=1:numel(ys)-1
    qy=(ys(j)+ys(j+1))/2+(ys(j+1)-ys(j))*z/2;
    wy=(ys(j+1)-ys(j))*w/2;
    X=qx+zeros(numel(qy),1);Y=qy+zeros(1,numel(qx));
    v=v+sum(integrand(X,Y).*(wy*wx),'all');
end
end

function v=origin_polar(H,n,kind)
[theta,wt]=axis_rule([0,atan(.5),pi/2],n);[z,w]=rule(n);
cs=cos(theta);sn=sin(theta);r0=H./max(cs,2*sn);
if strcmp(kind,'shell')
    r=r0.*(1.5+.5*z);wr=.5*w*r0;
    v=sum(integrand(r.*cs,r.*sn).*r.*wr.*wt,'all');
else
    invr=(.5+.5*z)./r0;winvr=.5*w./r0;r=1./invr;
    v=sum(integrand(r.*cs,r.*sn).*r./invr.^2.*winvr.*wt,'all');
end
end

function v=anchor_polar(H,n)
height=H/2;a=1;
edges=unique([linspace(0,pi,9),atan2(height,H-a),atan2(height,-a)]);
[phi,wp]=axis_rule(edges,n);[z,w]=rule(n);
cs=cos(phi);sn=sin(phi);rb=Inf(size(cs));
rb(cs>0)=(H-a)./cs(cs>0);rb(cs<0)=-a./cs(cs<0);
rb=min(rb,height./sn);
r=(.5+.5*z).*rb;wr=.5*w*rb;X=a+r.*cs;Y=r.*sn;
% The r Jacobian cancels the integrable kernel singularity at (a,0).
value=4*X.*sn./(pi*((X+a).^2+Y.^2)).*omega(X,Y);
v=sum(value.*wr.*wp,'all');
end

function v=integrand(x,y)
v=4*x.*y./(pi*((1-x).^2+y.^2).*((1+x).^2+y.^2)).*omega(x,y);
end
function v=omega(x,y)
v=x.^7./(1+x.^8+y.^8);
end
function [q,w]=axis_rule(edges,n)
edges=edges(:);[z,weights]=rule(n);
q=(edges(1:end-1)+edges(2:end))/2+diff(edges)*z'/2;
weights=diff(edges)*weights'/2;q=q(:)';w=weights(:)';
end
function [x,w]=rule(n)
j=(1:n-1)';b=j./sqrt(4*j.^2-1);[Q,T]=eig(diag(b,1)+diag(b,-1),'vector');
[x,index]=sort(T);w=2*Q(1,index)'.^2;
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
