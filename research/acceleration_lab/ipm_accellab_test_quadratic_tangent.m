function report=ipm_accellab_test_quadratic_tangent(outFile)
% Pure three-value MMS of the peak Hessian identity, NOT an IPM RHS test.
assert(~isfile(outFile));x=[-.5,0,.4];a=.07;rho=2-2*(x-a).^2;
p=ipm.evolve.quadraticPeakFunctional(rho,x,2);
v=.3+.7*(x-a)+.2*(x-a).^2;w=-.2+.9*(x-a)-.1*(x-a).^2;
expectedHessian=-.7*.9/(-4);epsilons=[1e-3,5e-4,2.5e-4];measured=zeros(size(epsilons));
for j=1:numel(epsilons)
 e=epsilons(j);plus=ipm.evolve.quadraticPeakFunctional(rho+e*w,x,2);
 minus=ipm.evolve.quadraticPeakFunctional(rho-e*w,x,2);
 measured(j)=(sum(plus.weights.*v)-sum(minus.weights.*v))/(2*e);
end
direction=w-rho*sum(p.weights.*w)/p.value;
F=tangent_field(rho,x,v);directionPoly=polyfit(x,direction,2);forcePoly=polyfit(x,F,2);
hessianTerm=-polyval(polyder(directionPoly),p.x)*polyval(polyder(forcePoly),p.x)/(2*p.curvature);
dpJ=zeros(size(epsilons));closure=zeros(size(epsilons));
for j=1:numel(epsilons)
 e=epsilons(j);Jv=(tangent_field(rho+e*direction,x,v)-tangent_field(rho-e*direction,x,v))/(2*e);
 dpJ(j)=sum(p.weights.*Jv);closure(j)=dpJ(j)+hessianTerm;
end
report=struct('kind','synthetic_quadratic_peak_tangent_Hessian_MMS_not_IPM', ...
 'epsilon',epsilons,'expectedHessian',expectedHessian,'measuredHessian',measured, ...
 'hessianError',abs(measured-expectedHessian),'DPDirection',sum(p.weights.*direction), ...
 'DPF',sum(p.weights.*F),'D2PDirectionF',hessianTerm,'DPJDirection',dpJ, ...
 'differentiatedConstraintResidual',closure,'actualPDERHSCalls',0,'operatorBuilds',0, ...
 'interpretation','Synthetic projected polynomial vector field verifies the differential identity only. Nonzero DPJ explicitly disproves a generic invariant tangent-space assumption.');
report.passed=max(report.hessianError)<1e-7&&max(abs(closure))<1e-7&& ...
 abs(report.DPDirection)<1e-14&&abs(report.DPF)<1e-14&&min(abs(dpJ))>1e-2;
save(outFile,'report');assert(report.passed);disp(report);
end
function F=tangent_field(rho,x,v)
p=ipm.evolve.quadraticPeakFunctional(rho,x,2);F=v-rho*sum(p.weights.*v)/p.value;
end
