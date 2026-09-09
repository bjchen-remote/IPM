function report=ipm_perflab_public_lu_tiny(outputFile)
%IPM_PERFLAB_PUBLIC_LU_TINY Public LU triangular cache on the same tiny A.
% No internal MATLAB API or change to refinement/warning policy of native LU.
assert(maxNumCompThreads==10&&~isfile(outputFile));cases=cell(1,6);k=0;
for nx=[33,65]
 ny=(nx+1)/2;
 for stretch=[0,4,12]
  k=k+1;x=linspace(-1,1,nx);y=linspace(0,1,ny)';
  if stretch>0,x=sinh(stretch*x)/sinh(stretch);y=sinh(stretch*y)/sinh(stretch);end
  xx=-ipm.mesh.fdMatrix(x,2,7);yy=-ipm.mesh.fdMatrix(y,2,7);
  A=kron(speye(nx-2),yy(2:end-1,2:end-1))+kron(xx(2:end-1,2:end-1),speye(ny-2));
  native=decomposition(A,'lu');
  [L,U,p,q,D]=lu(A,'vector');
  % The public permutation/scaling contract is P*(D\A)*Q=L*U.
  % L/U condition numbers are not substituted for the native A condition.
  cache=struct('L',decomposition(L,'triangular','lower','CheckCondition',false), ...
   'U',decomposition(U,'triangular','upper','CheckCondition',false),'p',p,'q',q,'D',D);
  [X,Y]=meshgrid(linspace(0,1,nx-2),linspace(0,1,ny-2));
  oracle=sin(pi*(X+.03)).*cos(2*pi*(Y+.07));b=A*oracle(:);
  a=native\b;c=triangular(cache,b);
  refError=norm(a-oracle(:),inf)/norm(oracle(:),inf);candidateError=norm(c-oracle(:),inf)/norm(oracle(:),inf);
  order=repmat([1,2,2,1],4,1);order(2:2:end,:)=3-order(2:2:end,:);times=zeros(size(order));
  for i=1:4
   for j=1:4
    timer=tic;if order(i,j)==1,value=native\b;else,value=triangular(cache,b);end
    times(i,j)=toc(timer);
    if order(i,j)==1,assert(isequal(value,a));else,assert(isequal(value,c));end
   end
  end
  cases{k}=struct('nodes',[nx,ny],'stretch',stretch,'rcondA',rcond(native), ...
   'factorNnz',[nnz(L),nnz(U)],'bitwise',isequal(a,c),'maximumAbsoluteDifference',norm(a-c,inf), ...
   'relativeForwardDifference',norm(a-c,inf)/max(norm(a,inf),realmin), ...
   'nativeMmsForwardError',refError,'candidateMmsForwardError',candidateError, ...
   'nativeComponentwiseBackwardError',componentwise(A,a,b), ...
   'candidateComponentwiseBackwardError',componentwise(A,c,b), ...
   'order',order,'seconds',times,'nativeMedianSeconds',median(times(order==1)), ...
   'candidateMedianSeconds',median(times(order==2)), ...
   'speedup',median(times(order==1))/median(times(order==2)));
  fprintf('PUBLIC_LU_TINY n%d stretch%d exact=%d speedup=%.5g diff=%.3g backward=%.3g/%.3g\n', ...
   nx,stretch,cases{k}.bitwise,cases{k}.speedup,cases{k}.relativeForwardDifference, ...
   cases{k}.nativeComponentwiseBackwardError,cases{k}.candidateComponentwiseBackwardError);
 end
end
report=struct('kind','public_sparse_LU_same_matrix_triangular_cache_tiny','cases',{cases}, ...
 'allBitwise',all(cellfun(@(a)a.bitwise,cases)),'productionEligible',false,'largeLuBuilds',0,'pdeSteps',0, ...
 'reference','https://www.mathworks.com/help/matlab/ref/lu.html', ...
 'interpretation','Identical A and b in real arithmetic. Public triangular solves need not reproduce the native UMFPACK solve/refinement rounding. Any nonbitwise result rejects this as an exact replacement; tiny speed is not a large-grid claim.');
fid=fopen(outputFile,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
end
function x=triangular(c,b)
scaled=c.D\b;z=c.U\(c.L\scaled(c.p));x=zeros(size(b));x(c.q)=z;
end
function e=componentwise(A,x,b)
e=max(abs(A*x-b)./max(abs(A)*abs(x)+abs(b),realmin));
end
