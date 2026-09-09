function report = ipm_accellab_test_recycle_algebra(outputDirectory)
%IPM_ACCELLAB_TEST_RECYCLE_ALGEBRA Manufactured array tests; no LU or PDE.
% Manufactured X and B=A*X are not native solver samples or an IPM trajectory.
assert(~isfolder(outputDirectory),'ipm:ResearchOutputExists','Use a new directory.');
mkdir(outputDirectory);
registration = struct('kind','manufactured_recycle_algebra_v1', ...
    'nodeCount',24,'sampleCount',9,'memory',4,'LUCount',0,'PDECount',0, ...
    'isNativeData',false,'integrationApproved',false);
write_json(fullfile(outputDirectory,'registration.json'),registration);
n=24; x=(1:n)'/(n+1); t=(0:8)/8;
T=spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);
A=spdiags(10.^linspace(-12,0,n)',0,n,n)*T;
V=sin(pi*x*(1:4));
X=ones(n,1)+V*[t;t.^2;t.^3;t.^4];
B=A*X;
velocityMap=spdiags([-ones(n,1),ones(n,1)],[-1,1],n,n)/2;
positive=ipm_accellab_recycle_shadow(A,B,X,velocityMap);
assert(all([positive.entries.screenPassed]));
Xbad=X; Xbad(:,6)=Xbad(:,6)+1e-3*sin(7*pi*x);
negative=ipm_accellab_recycle_shadow(A,A*Xbad,Xbad,velocityMap);
assert(~negative.entries(1).screenPassed);

% Both raw and componentwise residuals can miss near-nullspace error.
epsilon=1e-14;
Aill=[1,1;1,1+epsilon]; exact=[1;2]; wrong=[2;1]; bill=Aill*exact;
r=Aill*wrong-bill;
trap=struct('relativeGlobalResidual',norm(r,inf)/norm(bill,inf), ...
    'componentwiseBackwardError',max(abs(r)./(abs(Aill)*abs(wrong)+abs(bill))), ...
    'velocityRelativeInf',norm(wrong-exact,inf)/norm(exact,inf));
assert(trap.relativeGlobalResidual<5e-12 && ...
    trap.componentwiseBackwardError<5e-14 && trap.velocityRelativeInf>.1);
report=struct('registration',registration,'positive',positive,'negative',negative, ...
    'residualOnlyFalsePositive',trap,'allManufacturedTestsPassed',true, ...
    'IPMValidationPassed',false,'integrationApproved',false);
save(fullfile(outputDirectory,'report.mat'),'report');
summary=report;
for k=1:numel(summary.positive.entries), summary.positive.entries(k).candidate=[]; end
for k=1:numel(summary.negative.entries), summary.negative.entries(k).candidate=[]; end
write_json(fullfile(outputDirectory,'report.json'),summary);
end

function write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
