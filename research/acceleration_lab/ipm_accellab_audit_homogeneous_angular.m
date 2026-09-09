function report=ipm_accellab_audit_homogeneous_angular(registrationFile)
% Independent fixed-sector quadrature audit; no checkpoint or PDE data.
reg=jsondecode(fileread(registrationFile));assert(~isfolder(reg.outputDirectory));mkdir(reg.outputDirectory);
assert(reg.gaussOrder==32&&reg.thetaCount==129&&reg.maximumYSlope==.25&& ...
 reg.absoluteTolerance==1e-13&&reg.relativeTolerance==1e-13&&reg.maximumDifference==1e-11);
theta=linspace(0,asin(reg.maximumYSlope),reg.thetaCount)';
% Exactly the registered model's 32-point Gaussian construction and formulas.
k=(1:31)';b=k./sqrt(4*k.^2-1);[q,l]=eig(diag(b,1)+diag(b,-1));
[nodes,i]=sort(diag(l));weights=2*q(1,i)'.^2;
ang=.5*(nodes+1)*theta';co=cos(ang);si=sin(ang);ff=co.^7./(co.^8+si.^8);
ic=(.5*theta').*(weights'*(co.*ff));is=(.5*theta').*(weights'*(si.*ff));
c=cos(theta);z=sin(theta);p=z.*(pi/4-ic')+c.*is';dp=c.*(pi/4-ic')-z.*is';
u1=-z.*p-c.*dp;u2=c.*p-z.*dp;
icReference=zeros(size(theta));isReference=icReference;
fc=@(a)cos(a).^8./(cos(a).^8+sin(a).^8);
fs=@(a)sin(a).*cos(a).^7./(cos(a).^8+sin(a).^8);
for j=2:numel(theta)
 icReference(j)=integral(fc,0,theta(j),'AbsTol',reg.absoluteTolerance,'RelTol',reg.relativeTolerance);
 isReference(j)=integral(fs,0,theta(j),'AbsTol',reg.absoluteTolerance,'RelTol',reg.relativeTolerance);
end
pr=z.*(pi/4-icReference)+c.*isReference;dpr=c.*(pi/4-icReference)-z.*isReference;
u1r=-z.*pr-c.*dpr;u2r=c.*pr-z.*dpr;
gauss=[ic',is',p,dp,u1,u2];adaptive=[icReference,isReference,pr,dpr,u1r,u2r];
error=abs(gauss-adaptive);maximumAbsoluteDifference=max(error,[],1);
report=struct('kind','fixed_near_wall_sector_homogeneous_angular_quadrature_audit', ...
 'thetaRange',[theta(1),theta(end)],'thetaCount',numel(theta),'gaussOrder',32, ...
 'referenceMethod','MATLAB integral adaptive Gauss-Kronrod, AbsTol=RelTol=1e-13, separate Ic and Is', ...
 'fieldNames',{{'Ic','Is','p','pPrime','U1','U2'}},'maximumAbsoluteDifference',maximumAbsoluteDifference, ...
 'maximumAcrossFields',max(maximumAbsoluteDifference),'registeredGate',reg.maximumDifference, ...
 'passed',all(isfinite(gauss),'all')&&all(isfinite(adaptive),'all')&&all(maximumAbsoluteDifference<=reg.maximumDifference), ...
 'minimumSampledP',min(p),'wallValues',gauss(1,:), ...
 'nativeCPRead',false,'futureDataRead',false,'LUCalls',0,'PDECalls',0,'modelParametersChanged',false, ...
 'scope','Only theta in [0,asin(.25)] on the fixed129-point grid. This is observed angular quadrature agreement, not a full-angle or full-exterior error bound, nor a proof of evolved IPM closure.');
save(fullfile(reg.outputDirectory,'angular_audit.mat'),'report','theta','gauss','adaptive','error','reg');
fid=fopen(fullfile(reg.outputDirectory,'report.json'),'w');guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));assert(report.passed);disp(report);
end
