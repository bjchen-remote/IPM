function report=ipm_accellab_smooth_inner_rates(omega,forcing,x,y,Dx,Dy,window,theta)
%IPM_ACCELLAB_SMOOTH_INNER_RATES C1 connected positive-part integral widths.
% Primary theta=.5; .3/.7 are registered holdouts. No fitted derivative.
if nargin<8,theta=.5;end
assert(ismember(theta,[.3,.5,.7]));m=4;x=x(:)';y=y(:);
peak=ipm_accellab_hermite_peak(omega(1,:),forcing(1,:),x,Dx,window);
report=struct('kind','connected_C1_smooth_amplitude_integral_width','valid',false, ...
    'theta',theta,'power',m,'peak',peak,'invalidReasons',{{}},'PPrimeNeverAssignedZero',true);
if ~peak.valid || ~peak.smoothUniquePeak
    report.invalidReasons={'invalid_or_nonunique_peak'};return
end
a=peak.x;P=peak.value;PPrime=peak.PPrime;
wall=omega(1,:);wallD=wall*Dx';force=forcing(1,:);forceD=force*Dx';
v=ipm_accellab_hermite_line(wall,wallD,x,a,peak.selectedCell);
f=ipm_accellab_hermite_line(force,forceD,x,a,peak.selectedCell);
distance=min(abs(x-a));
if distance<=peak.activePositionTolerance
    report.invalidReasons={'peak_phase_at_native_interface_requires_one_sided_analysis'};return
end
aPrime=-f.derivative/v.secondDerivative;
wallRoots=ipm_accellab_hermite_level_roots(wall,wallD,x,theta*P,window);
left=find(wallRoots.locations<a,1,'last');right=find(wallRoots.locations>a,1,'first');
if isempty(left) || isempty(right)
    report.invalidReasons={'wall_connected_component_not_bounded'};return
end
indices=[left,right];bounds=wallRoots.locations(indices);
if ~regular_component(wallRoots,x,bounds,indices,[1,-1])
    report.invalidReasons={'wall_component_topology_or_threshold_degeneracy'};return
end
trace=ipm_accellab_tensor_hermite(omega,x,y,Dx,Dy,a+zeros(size(y)),y);
traceF=ipm_accellab_tensor_hermite(forcing,x,y,Dx,Dy,a+zeros(size(y)),y);
traceRate=traceF.value+aPrime*trace.derivativeX;
traceRateY=traceF.derivativeY+aPrime*trace.derivativeXY;
verticalRoots=ipm_accellab_hermite_level_roots(trace.value,trace.derivativeY,y,theta*P,[0,y(end)]);
j=find(verticalRoots.locations>0,1,'first');
if isempty(j) || ~regular_component(verticalRoots,y,[0,verticalRoots.locations(j)],j,-1)
    report.invalidReasons={'vertical_component_topology_or_threshold_degeneracy'};return
end
horizontal=integral_width(wall,wallD,force,forceD,x,bounds,P,PPrime,theta,m);
vertical=integral_width(trace.value,trace.derivativeY,traceRate,traceRateY,y,[0,verticalRoots.locations(j)],P,PPrime,theta,m);
report.valid=horizontal.value>0 && vertical.value>0 && all(isfinite([horizontal.value,horizontal.prime,vertical.value,vertical.prime,aPrime]));
report.peakPrime=PPrime;report.translationRate=aPrime;
report.rawWallWidth=horizontal.value;report.rawWallWidthPrime=horizontal.prime;
report.rawVerticalWidth=vertical.value;report.rawVerticalWidthPrime=vertical.prime;
report.logScaleXRate=horizontal.prime/horizontal.value;report.logScaleYRate=vertical.prime/vertical.value;
report.horizontal=horizontal;report.vertical=vertical;
report.wallThresholdRoots=wallRoots;report.verticalThresholdRoots=verticalRoots;
report.signature=struct('peakCell',peak.selectedCell,'wallRootCells',wallRoots.cells(indices),'verticalRootCell',verticalRoots.cells(j),'theta',theta);
report.verticalTraceTranslationIncluded=true;
report.interpretation='Read-only new width observation. Threshold roots delimit the connected component but no root slope or root time derivative enters the integral derivative. Boundary integrand is zero. Actual PPrime and moving-peak trace forcing are retained. Topology changes and degenerate peak phase are guarded, not claimed differentiable.';
end
function passed=regular_component(roots,x,bounds,indices,signs)
x=x(:)';flat=roots.flatCells;
passed=~any(x(flat)<=bounds(2) & x(flat+1)>=bounds(1)) && ...
    all(roots.relativeResiduals(indices)<1e-11) && all(roots.slopes(indices).*signs>0) && ...
    all(abs(roots.normalizedSlopes(indices))>1e-10);
end
function out=integral_width(value,jets,forcing,forcingJets,x,bounds,P,PPrime,theta,m)
x=x(:)';edges=unique([bounds,x(x>bounds(1)&x<bounds(2))]);
z=[-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498, ...
    .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363];
w=[.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620, ...
    .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763];
q=(edges(1:end-1)'+edges(2:end)')/2+diff(edges)'*z/2;weights=diff(edges)'*w/2;
v=ipm_accellab_hermite_line(value,jets,x,q);f=ipm_accellab_hermite_line(forcing,forcingJets,x,q);
s=(v.value/P-theta)/(1-theta);assert(min(s,[],'all')>-1e-9,'ipm:SmoothComponent','Selected interval contains a subthreshold region.');
s=max(0,s);normalizedRate=f.value/P-v.value*PPrime/P^2;
integrand=s.^m;rateIntegrand=m*s.^(m-1).*normalizedRate/(1-theta);
out=struct('value',sum(weights.*integrand,'all'),'prime',sum(weights.*rateIntegrand,'all'), ...
    'bounds',bounds,'cellCount',numel(edges)-1,'gaussNodesPerCell',8,'minimumPositivePartArgument',min(s,[],'all'), ...
    'endpointWeight',[0,0],'rootSlopeDivisionUsed',false,'maximumIntegrandPolynomialDegree',3*m, ...
    'quadratureMeaning','Gauss8 integrates the degree-12 value and degree-12 true directional derivative within each native cubic cell, up to roundoff.');
end
