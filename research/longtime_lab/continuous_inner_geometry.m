function report=continuous_inner_geometry(omega,x,y,Dx,Dy,window)
%CONTINUOUS_INNER_GEOMETRY Paired C1 peak and connected 0.9 level geometry.
% Frozen geometry only: zero forcing passed to the peak helper is a dummy,
% never a measurement of the actual PDE peak derivative.
x=x(:)';y=y(:);
p=ipm_accellab_hermite_peak(omega(1,:),zeros(size(x)),x,Dx,window);
assert(p.valid && p.smoothUniquePeak,'A unique nondegenerate interior C1 peak is required.');
P=p.value;a=p.x;threshold=.9*P;
wallSlopes=omega(1,:)*Dx';
[locations,slopes]=level_roots(omega(1,:),wallSlopes,x,threshold,window);
left=find(locations<a,1,'last');right=find(locations>a,1,'first');
assert(~isempty(left) && ~isempty(right) && slopes(left)>0 && slopes(right)<0, ...
    'Connected wall level roots must be transverse and bracket the peak.');
trace=ipm_accellab_tensor_hermite(omega,x,y,Dx,Dy,a+zeros(size(y)),y);
assert(abs(trace.value(1)-P)<=1e-12*max(1,P));
[vertical,verticalSlopes]=level_roots(trace.value,trace.derivativeY,y,threshold,[y(1),y(end)],true);
assert(~isempty(vertical) && vertical(1)>0 && verticalSlopes(1)<0, ...
    'The first vertical level root must be a transverse downcrossing.');
report=struct('kind','continuous_C1_inner_geometry_only','peakX',a,'peak',P, ...
    'peakCell',p.selectedCell,'peakFraction',p.selectedFraction, ...
    'wallCurvature',p.selectedSecondDerivative,'level',.9, ...
    'wallCrossings',[locations(left),locations(right)], ...
    'wallCrossingSlopes',[slopes(left),slopes(right)], ...
    'wallWidth',locations(right)-locations(left), ...
    'verticalWidth',vertical(1),'verticalCrossingSlope',verticalSlopes(1), ...
    'window',window,'pdeDerivativeEvaluated',false,'nativeGaugeChanged',false);
end

function [locations,slopes]=level_roots(values,jets,x,target,window,firstOnly)
if nargin<6,firstOnly=false;end
values=values(:)';jets=jets(:)';x=x(:)';locations=[];slopes=[];
for k=1:numel(x)-1
    if x(k+1)<window(1) || x(k)>window(2),continue;end
    h=x(k+1)-x(k);v0=values(k);v1=values(k+1);d0=jets(k);d1=jets(k+1);
    c=[2*v0-2*v1+h*(d0+d1),-3*v0+3*v1-h*(2*d0+d1),h*d0,v0-target];
    scale=max(abs(c));assert(scale>0,'A whole cell lies at the threshold.');
    points=roots(c/scale);
    points=sort(real(points(abs(imag(points))<=1e-10 & real(points)>=-1e-12 & real(points)<=1+1e-12)));
    for pointValue=points(:)'
        t=min(1,max(0,pointValue));z=x(k)+h*t;
        if z<window(1) || z>window(2),continue;end
        if ~isempty(locations) && abs(z-locations(end))<=64*eps(max(1,abs(z))),continue;end
        locations(end+1)=z;slopes(end+1)=(3*c(1)*t^2+2*c(2)*t+c(3))/h; %#ok<AGROW>
        if firstOnly,return;end
    end
end
end
