function [values,info]=ipm_accellab_green_full_tail_trace(query,H,options)
%IPM_ACCELLAB_GREEN_FULL_TAIL_TRACE Full exterior Green potential, no PDE.
% Positive source region is [0,infinity)^2 minus [0,H]x[0,H/2].
% Source is W/C times the analytic physical source evaluated at S/C:
% x^7/(x^8+y^8), optionally plus eps*d_x[x^2/(x^2+y^2)^(3/2)].
% These exterior manufactured fields are NOT an unforced IPM solution.
arguments
    query (:,2) double
    H (1,1) double {mustBePositive,mustBeFinite}
    options.method (1,:) char {mustBeMember(options.method,{'radial_analytic','target_polar'})}='radial_analytic'
    options.order (1,1) double {mustBeInteger,mustBePositive}=128
    options.Cx (1,1) double {mustBePositive,mustBeFinite}=1
    options.Comega (1,1) double {mustBePositive,mustBeFinite}=1
    options.perturbation (1,1) double {mustBeFinite}=0
    options.sourceOnly (1,:) char {mustBeMember(options.sourceOnly,{'leading','perturbation','sum','exact_initial','initial_remainder'})}='leading'
end
assert(all(isfinite(query),'all')&&isreal(query)&&all(query(:,2)>=0));
assert(all(abs(query(:,1))<=H)&&all(query(:,2)<=H/2));
if strcmp(options.method,'radial_analytic')
    assert(strcmp(options.sourceOnly,'leading')&&options.perturbation==0, ...
        'ipm:TailTraceScope','Analytic radial primitive is registered only for the leading source.');
end
[u,w]=gauss01(options.order);values=zeros(size(query,1),1);timer=tic;
for k=1:size(query,1)
    q=query(k,:)/H;
    if q(1)==0 || q(2)==0,continue;end
    if strcmp(options.method,'radial_analytic')
        values(k)=options.Comega*H*radial_analytic(q,u,w);
    else
        values(k)=sign(q(1))*target_polar(abs(q),H,u,w,options);
    end
end
info=struct('method',options.method,'order',options.order,'queryCount',size(query,1), ...
    'secondsExcludingRuleConstruction',toc(timer),'H',H,'Cx',options.Cx,'Comega',options.Comega, ...
    'sourceOnly',options.sourceOnly,'perturbation',options.perturbation, ...
    'fullGreenKernel',true,'innerPointMultipoleApproximationUsed',false,'LUCount',0,'PDECount',0);
end

function value=radial_analytic(q,u,w)
% With source f(theta)/r the radial measure is f(theta) dr. The four
% log-quadratic antiderivatives are combined before the limit at infinity.
qt=atan2(q(2),abs(q(1)));edges=unique([0,atan(.5),qt,pi/2]);value=0;
for p=1:numel(edges)-1
    lo=edges(p);hi=edges(p+1);d=hi-lo;
    if lo==qt
        theta=lo+d*u.^2;wt=2*d*u.*w;
    elseif hi==qt
        theta=hi-d*u.^2;wt=2*d*u.*w;
    else
        theta=lo+d*u;wt=d*w;
    end
    c=cos(theta);s=sin(theta);r0=1./max(c,2*s);
    bd=q(1)*c+q(2)*s;by=q(1)*c-q(2)*s;
    % Cross-products avoid cancellation in sqrt(|q|^2-b^2) near alignment.
    kd=abs(q(1)*s-q(2)*c);ky=abs(q(1)*s+q(2)*c);
    fy=primitive(r0,by,ky)+primitive(r0,-by,ky);
    fd=primitive(r0,bd,kd)+primitive(r0,-bd,kd);
    radial=(2*pi*(ky-kd)-fy+fd)/(4*pi);
    f=c.^7./(c.^8+s.^8);
    value=value+sum(wt.*f.*radial);
end
end
function v=primitive(r,b,k)
t=r-b;v=t.*log(t.^2+k.^2)-2*t+2*k.*atan2(t,k);
end

function value=target_polar(q,H,u,w,options)
% Independent coordinates S=q+s*(cos phi,sin phi). For a target inside/on
% the rectangle, each ray's exterior interval is [exit_rectangle,exit_quadrant].
% s=s0+L*u^2 or s=s0+(u/(1-u))^2 regularizes s*log(s) at a boundary target.
corners=[0,0;1,0;0,.5;1,.5];vectors=corners-q;
active=any(vectors~=0,2);angles=mod(atan2(vectors(active,2),vectors(active,1)),2*pi);
edges=unique([0;pi/2;pi;3*pi/2;2*pi;angles]);value=0;
for j=1:numel(edges)-1
    phi=edges(j)+(edges(j+1)-edges(j))*u';
    wp=(edges(j+1)-edges(j))*w';c=cos(phi);s=sin(phi);
    rx=Inf(size(c));ry=rx;qx=rx;qy=rx;
    pos=c>0;neg=c<0;rx(pos)=(1-q(1))./c(pos);rx(neg)=-q(1)./c(neg);qx(neg)=rx(neg);
    pos=s>0;neg=s<0;ry(pos)=(.5-q(2))./s(pos);ry(neg)=-q(2)./s(neg);qy(neg)=ry(neg);
    first=max(0,min(rx,ry));last=min(qx,qy);
    active=last>first;
    if ~any(active),continue;end
    c=c(active);s=s(active);wp=wp(active);first=first(active);last=last(active);
    radius=zeros(numel(u),numel(first));dr=radius;finite=isfinite(last);
    radius(:,finite)=first(finite)+(last(finite)-first(finite)).*u.^2;
    dr(:,finite)=2*(last(finite)-first(finite)).*u;
    radius(:,~finite)=first(~finite)+(u./(1-u)).^2;
    dr(:,~finite)=2*u./(1-u).^3+zeros(1,nnz(~finite));
    x=q(1)+radius.*c;y=q(2)+radius.*s;
    assert(all(x>=-32*eps & y>=-32*eps,'all'));
    % Source evaluation uses physical coordinates, independently of the
    % analytic leading-amplitude cancellation used by radial_analytic.
    xp=H*x/options.Cx;yp=H*y/options.Cx;
    source=zeros(size(x));
    if strcmp(options.sourceOnly,'exact_initial')
        source=xp.^7./(1+xp.^8+yp.^8);
    elseif strcmp(options.sourceOnly,'initial_remainder')
        powerSum=xp.^8+yp.^8;
        source=-xp.^7./(powerSum.*(1+powerSum));
    elseif ~strcmp(options.sourceOnly,'perturbation')
        source=xp.^7./(xp.^8+yp.^8);
    end
    if any(strcmp(options.sourceOnly,{'perturbation','sum'}))
        r2=xp.^2+yp.^2;
        extra=2*xp./r2.^1.5-3*xp.^3./r2.^2.5;
        if strcmp(options.sourceOnly,'perturbation'),source=extra;
        else,source=source+options.perturbation*extra;end
    end
    source=options.Comega/options.Cx*source;
    % Exact image-kernel identity, avoiding subtraction of distant logs.
    % The direct squared distance is exactly radius^2 in these coordinates.
    opposite=(q(1)+x).^2+(q(2)+y).^2;
    green=log1p(16*q(1)*q(2)*x.*y./(radius.^2.*opposite))/(4*pi);
    value=value+H^2*sum(green.*source.*radius.*dr.*(w*wp),'all');
end
end
function [u,w]=gauss01(n)
j=(1:n-1)';b=j./sqrt(4*j.^2-1);[Q,T]=eig(diag(b,1)+diag(b,-1),'vector');
[z,index]=sort(T);u=(z+1)/2;w=Q(1,index)'.^2;
end
