function report=ipm_perflab_green_reuse_probe(resultFile,outputDirectory)
%IPM_PERFLAB_GREEN_REUSE_PROBE No LU: reconstruct only Green geometry.
assert(maxNumCompThreads==10&&~isfolder(outputDirectory));mkdir(outputDirectory);
cases=cell(1,48);count=0;
for sizeIndex=1:2
 nx=32*sizeIndex+1;ny=16*sizeIndex+1;
 for stretch=[0,2]
  sx=linspace(-1,1,nx);sy=linspace(0,1,ny)';
  if stretch>0,sx=sinh(stretch*sx)/sinh(stretch);sy=sinh(stretch*sy)/sinh(stretch);end
  for H=[4,1000]
   for family=1:3
    for mode={'double_odd_omega','half_plane'}
     count=count+1;o=geometry(H*sx,H*sy,mode{1},1e-12,32);
     source=sin(o.X).*exp(-(o.X/3).^2-(o.Y/2).^2);
     if family==2,source=source+cos(o.Y).*o.X/100;o.greenMaxSources=10000;end
     if family==3,source(:)=0;o.farBoundaryMode='dirichlet_zero';end
     a=ipm.field.greenBoundary(source,o);b=ipm_perflab_green_candidate(source,o);
     cases{count}=struct('nodes',[nx,ny],'stretch',stretch,'H',H,'family',family,'symmetry',mode{1}, ...
      'exact',isequaln(a,b),'maximumAbsoluteDifference',boundary_difference(a,b));
    end
   end
  end
 end
end
cases=cases(1:count);assert(all(cellfun(@(a)a.exact,cases)),'ipm:GreenTinyExact','Tiny Green reuse failed bitwise; real data not evaluated.');
result=ipm.output.validate(resultFile);c=result.config;
assert(strcmp(c.transport.spatialDiscretization,'high_order')&&strcmp(result.scale.geometry,'isotropic'));
o=geometry(result.grid.x,result.grid.y,c.physics.symmetryMode,c.elliptic.greenSourceTolerance,c.elliptic.greenMaxSources);
o.farBoundaryMode=c.elliptic.farBoundaryMode;source=result.state.omega;
a=ipm.field.greenBoundary(source,o);b=ipm_perflab_green_candidate(source,o);
expectedPsiBoundary=zeros(size(source));expectedPsiBoundary(:,1)=a.left;expectedPsiBoundary(:,end)=a.right;
expectedPsiBoundary(1,:)=a.bottom;expectedPsiBoundary(end,:)=a.top;
mask=false(size(source));mask(:,[1,end])=true;mask([1,end],:)=true;
boundaryMatchesSavedPsi=isequal(expectedPsiBoundary(mask),result.state.psi(mask));
assert(boundaryMatchesSavedPsi,'ipm:GreenSavedPair','No-LU Green reconstruction must match the saved native boundary bitwise.');
realExact=isequaln(a,b);assert(realExact,'ipm:GreenActualExact','Frozen real Green candidate failed bitwise.');clear result
stream=RandStream('mt19937ar','Seed',20260909);order=repmat([1,2,2,1],6,1);
flip=rand(stream,6,1)>.5;order(flip,:)=3-order(flip,:);seconds=zeros(size(order));exact=true(size(order));
for i=1:6
 for j=1:4
  timer=tic;if order(i,j)==1,value=ipm.field.greenBoundary(source,o);else,value=ipm_perflab_green_candidate(source,o);end
  seconds(i,j)=toc(timer);exact(i,j)=isequaln(value,a);
 end
end
ratios=zeros(6,1);for i=1:6,ratios(i)=median(seconds(i,order(i,:)==1))/median(seconds(i,order(i,:)==2));end
report=struct('kind','Green_within_call_exact_reuse_no_LU','tinyCases',{cases},'tinyCount',count, ...
 'resultFile',resultFile,'nodes',[o.nx,o.ny],'boundaryMatchesSavedPsiBitwise',boundaryMatchesSavedPsi, ...
 'actualCandidateBitwise',realExact,'allTimingOutputsBitwise',all(exact,'all'), ...
 'order',order,'seconds',seconds,'baselineMedianSeconds',median(seconds(order==1)), ...
 'candidateMedianSeconds',median(seconds(order==2)),'pairedBlockSpeedup',ratios, ...
 'medianSpeedup',median(seconds(order==1))/median(seconds(order==2)), ...
 'threads',10,'operatorBuilds',0,'poissonSolves',0,'pdeSteps',0, ...
 'interpretation','No-LU geometry reconstructed using original quadrature and paired to the saved psi boundary. Timing is the complete Green function, not an exclusive whole-RHS profile.');
save(fullfile(outputDirectory,'green_reuse_report.mat'),'report','-v7.3');
fid=fopen(fullfile(outputDirectory,'green_reuse_report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
fprintf('GREEN_REUSE_NO_LU_PASS tiny=%d actualBitwise=%d seconds=%.6g/%.6g speedup=%.6g range=%.6g/%.6g\n', ...
 count,report.allTimingOutputsBitwise,report.baselineMedianSeconds,report.candidateMedianSeconds,report.medianSpeedup,min(ratios),max(ratios));
end
function o=geometry(x,y,symmetry,tolerance,maxSources)
x=x(:)';y=y(:);o=struct('x',x,'y',y,'nx',numel(x),'ny',numel(y), ...
 'X',repmat(x,numel(y),1),'Y',repmat(y,1,numel(x)), ...
 'integrationWeights',ipm.mesh.quadrature(y)*ipm.mesh.quadrature(x), ...
 'symmetryMode',symmetry,'farBoundaryMode','green','greenSourceTolerance',tolerance,'greenMaxSources',maxSources);
end
function v=boundary_difference(a,b)
v=0;for n={'left','right','top','bottom'},v=max(v,max(abs(a.(n{1})-b.(n{1})),[],'all'));end
end
