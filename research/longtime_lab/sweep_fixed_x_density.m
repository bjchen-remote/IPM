function report=sweep_fixed_x_density(distributionFile,outputFile)
%SWEEP_FIXED_X_DENSITY Same-N rounded-density spacing/shape experiment.
assert(isfile(distributionFile)&&~isfile(outputFile));
d=jsondecode(fileread(distributionFile));f=d.feature;
nx=d.nodeCount(1);H=d.x.nodes(end);anchor=1;
spacing0=min(diff(f.coreInterval)/32,diff(f.frontInterval)/20);
factors=[.95,1,1.05,1.10];
rows=struct([]);best=struct();bestRatio=Inf;
for factor=factors
 for fine=[32,64,96]
  for roundCells=[8,16,32]
   for fraction=[.4,.5,.65]
    try
     [axis,info]=ipm.remesh.corePatchAxis(anchor,f.coreCenter,H,nx, ...
        spacing0/factor,fine,fraction,roundCells,.25);
     q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
     pos=axis(axis>=0);
     core=interval_count(pos,f.coreInterval);
     front=interval_count(pos,f.frontInterval);
     valid=q.maximumAdjacentCellRatio<=1.08 && ...
         q.maximumLogSpacingCurvature<=.01 && q.minimumStencilRcond>=1e-9 && ...
         q.minimumQuadratureWeightRatio>=1e-8 && ...
         q.minimumQuadratureWeightToControlWidthRatio>=.35 && ...
         q.maximumQuadratureWeightToControlWidthRatio<=1.65 && ...
         core>=32 && front>=20 && any(axis==1) && any(axis==-1);
     if valid && q.maximumAdjacentCellRatio<bestRatio
         bestRatio=q.maximumAdjacentCellRatio;
         best=struct('factor',factor,'fine',fine,'roundCells',roundCells, ...
             'fraction',fraction,'axis',axis,'quality',q, ...
             'coreCells',core,'frontCells',front,'construction',info);
     end
     rows=append(rows,struct('factor',factor,'fine',fine, ...
         'roundCells',roundCells,'fraction',fraction,'valid',valid, ...
         'ratio',q.maximumAdjacentCellRatio,'curvature',q.maximumLogSpacingCurvature, ...
         'coreCells',core,'frontCells',front));
    catch e
     if ~strcmp(e.identifier,'MATLAB:assertion:failed')
         rethrow(e);
     end
    end
   end
  end
 end
end
report=struct('kind','same_n_rounded_density_factor_sweep_v1', ...
    'source',distributionFile,'baselineRatio',d.x.maximumAdjacentRatio, ...
    'nodeCount',nx,'bestRatio',bestRatio,'best',best, ...
    'evaluated',numel(rows),'rows',rows,'geometryOnly',true);
fid=fopen(outputFile,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end
function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
