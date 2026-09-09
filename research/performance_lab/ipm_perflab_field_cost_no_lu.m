function report=ipm_perflab_field_cost_no_lu(resultFile,outputFile)
%IPM_PERFLAB_FIELD_COST_NO_LU Exact saved-field derivatives and A multiply.
assert(maxNumCompThreads==10&&~isfile(outputFile));r=ipm.output.validate(resultFile);
x=r.grid.x(:)';y=r.grid.y(:);Dx=ipm.mesh.fdMatrix(x,1,7);Dy=ipm.mesh.fdMatrix(y,1,7);
nx=numel(x);ny=numel(y);psi=r.state.psi;rho=r.state.rho;
source=rho*Dx';u1=-(Dy*psi);u2=psi*Dx';u2(1,:)=0;
assert(isequal(source,r.state.omega)&&isequal(u1,r.state.velocity.x)&&isequal(u2,r.state.velocity.y));
xx=-ipm.mesh.fdMatrix(x,2,7);yy=-ipm.mesh.fdMatrix(y,2,7);
Tx=xx(2:end-1,2:end-1);Ty=yy(2:end-1,2:end-1);
A=kron(speye(nx-2),Ty)+kron(Tx,speye(ny-2));
boundaryX=[psi(2:end-1,1),psi(2:end-1,end)];boundaryY=[psi(1,2:end-1);psi(end,2:end-1)];
rhs=source(2:end-1,2:end-1)-boundaryX*xx(2:end-1,[1,end])'-yy(2:end-1,[1,end])*boundaryY;
interior=psi(2:end-1,2:end-1);residual=norm(A*interior(:)-rhs(:),inf)/max(norm(rhs(:),inf),eps);
assert(residual==r.elliptic.solveInfo.relativeResidual);
operations={@()rho*Dx',@()-(Dy*psi),@()psi*Dx',@()A*interior(:)};
names={'source_rho_Dx','velocity_Dy_psi','velocity_psi_Dx','poisson_A_psi'};
seconds=zeros(8,4);for k=1:4,op=operations{k};op();end
for i=1:8
 order=1:4;if mod(i,2)==0,order=4:-1:1;end
 for k=order,op=operations{k};timer=tic;op();seconds(i,k)=toc(timer);end
end
report=struct('kind','saved_field_exact_derivative_and_residual_cost_no_LU','resultFile',resultFile, ...
 'nodes',[nx,ny],'names',{names},'seconds',seconds,'medianSeconds',median(seconds,1), ...
 'sourceAndVelocityBitwise',true,'poissonResidualBitwise',true,'poissonResidual',residual, ...
 'sparseOperatorNnz',nnz(A),'twoDimensionalOperatorAssemblies',1,'luBuilds',0,'poissonSolves',0,'pdeSteps',0, ...
 'interpretation','Separate no-LU component timings, exactly paired to saved source/velocity/residual. Not an exclusive full-RHS timing profile.');
fid=fopen(outputFile,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,PrettyPrint=true),'char');
fprintf('FIELD_COST_NO_LU_PASS medians_ms=%.6g/%.6g/%.6g/%.6g\n',1000*report.medianSeconds);
end
