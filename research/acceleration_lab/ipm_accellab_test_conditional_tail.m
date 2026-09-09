function report=ipm_accellab_test_conditional_tail(out)
%IPM_ACCELLAB_TEST_CONDITIONAL_TAIL Pure quadrature MMS, not an IPM solution.
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
reg=struct('kind','conditional_far_tail_manufactured_integral_v1', ...
    'H',[8,16,32,64,128],'epsilon',[0,.1,.2],'orders',[64,128,256], ...
    'physicalAnchor',1,'aspectHeightOverWidth',.5, ...
    'manufacturedDensity','rho0+epsilon*x^2/(1+x^2+y^2)^(3/2)', ...
    'manufacturedSolutionOfUnforcedIPM',false,'LUCount',0,'PDECount',0, ...
    'absoluteQuadratureTolerance',1e-11,'covarianceRelativeTolerance',1e-11, ...
    'CxValues',[.5,1,2],'ComegaValues',[.5,1,2], ...
    'covariancePhysicalH',16,'covariancePhysicalEpsilon',.2);
write_json(fullfile(out,'registration.json'),reg);
v=zeros(numel(reg.H),numel(reg.epsilon),numel(reg.orders));
for j=1:numel(reg.H)
    for k=1:numel(reg.epsilon)
        for n=1:numel(reg.orders)
            v(j,k,n)=tail(reg.H(j),1,reg.epsilon(k),reg.orders(n),1,1);
        end
    end
end
report=reg;report.tail=v;report.lastQuadratureChange=abs(v(:,:,end)-v(:,:,end-1));
report.quadraturePassed=max(report.lastQuadratureChange,[],'all')<=reg.absoluteQuadratureTolerance;
% The manufactured perturbation has g(theta)/r^2+O(r^-4),
% g=2*cos(theta)-3*cos(theta)^3, independently of the k8 base density.
[th,wt]=theta_rule(reg.orders(end));c=cos(th);s=sin(th);m=max(c,2*s);
report.analyticPerturbationCoefficientB=sum(wt.*(2/pi*c.*s.*(2*c-3*c.^3).*m.^2));
report.perturbationTail=v(:,2:end,end)-v(:,1,end);
report.HSquaredPerturbationPerEpsilon=reg.H(:).^2.*report.perturbationTail./reg.epsilon(2:end);
report.coefficientError=abs(report.HSquaredPerturbationPerEpsilon-report.analyticPerturbationCoefficientB);
report.errorRatios=report.coefficientError(1:end-1,:)./report.coefficientError(2:end,:);
report.linearityError=max(abs(report.perturbationTail(:,2)-2*report.perturbationTail(:,1)));
report.manufacturedCoefficientApproachPassed=all(diff(report.coefficientError)<0,'all') && ...
    all(report.errorRatios(end,:)>=3.8 & report.errorRatios(end,:)<=4.2);
covariance=struct([]);physical=tail(16,1,.2,256,1,1);
for C=reg.CxValues
    for W=reg.ComegaValues
        canonical=tail(16*C,C,.2,256,C,W);
        recoveredPhysical=C/W*canonical;
        row=struct('Cx',C,'Comega',W,'computationalBoxH',16*C, ...
            'computationalAnchor',C,'canonicalTailCL',canonical, ...
            'recoveredPhysicalCL',recoveredPhysical,'physicalCL',physical, ...
            'relativeError',abs(recoveredPhysical-physical)/abs(physical));
        if isempty(covariance),covariance=row;else,covariance(end+1)=row;end
    end
end
report.covariance=covariance;report.covariancePassed=max([covariance.relativeError])<=reg.covarianceRelativeTolerance;
report.allPureIntegralChecksPassed=report.quadraturePassed&&report.covariancePassed&& ...
    report.manufacturedCoefficientApproachPassed&&report.linearityError<=reg.absoluteQuadratureTolerance;
report.claim='Verifies algebra and quadrature for a chosen smooth manufactured density and coordinate changes only. It neither propagates the far-field bounds for IPM nor validates exterior data at t>0.';
save(fullfile(out,'report.mat'),'report');write_json(fullfile(out,'report.json'),report);
disp(jsonencode(struct('B',report.analyticPerturbationCoefficientB, ...
    'HSquaredPerturbationPerEpsilon',report.HSquaredPerturbationPerEpsilon, ...
    'errorRatios',report.errorRatios,'maximumQuadratureChange',max(report.lastQuadratureChange,[],'all'), ...
    'maximumCovarianceRelativeError',max([covariance.relativeError]), ...
    'allPureIntegralChecksPassed',report.allPureIntegralChecksPassed)));
end
function v=tail(H,a,epsilon,n,C,W)
[theta,wt]=theta_rule(n);[z,w]=rule(n);
c=cos(theta);s=sin(theta);r0=H./max(c,2*s);
zrad=(.5+.5*z)./r0;wz=.5*w./r0;r=1./zrad;
x=r.*c;y=r.*s;xp=x/C;yp=y/C;
omega0=xp.^7./(1+xp.^8+yp.^8);
den=1+xp.^2+yp.^2;
omegaPerturbation=2*xp./den.^1.5-3*xp.^3./den.^2.5;
source=W/C*(omega0+epsilon*omegaPerturbation);
kernel=4*x.*y./(pi*((a-x).^2+y.^2).*((a+x).^2+y.^2));
v=sum(kernel.*source.*r./zrad.^2.*wz.*wt,'all');
end
function [theta,wt]=theta_rule(n)
[z,w]=rule(n);e=[0;atan(.5);pi/2];
theta=(e(1:end-1)+e(2:end))/2+diff(e)*z'/2;wt=diff(e)*w'/2;
theta=theta(:)';wt=wt(:)';
end
function [z,w]=rule(n)
j=(1:n-1)';b=j./sqrt(4*j.^2-1);[Q,T]=eig(diag(b,1)+diag(b,-1),'vector');
[z,index]=sort(T);w=2*Q(1,index)'.^2;
end
function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
