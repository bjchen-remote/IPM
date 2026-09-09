function report=compare_space_observers(firstCheckpoint,secondCheckpoint,window,outputFile)
%COMPARE_SPACE_OBSERVERS Bidirectional fixed-lab observation sensitivity.
% No LU/PDE. Linear and C1 Hermite reports coexist; neither changes a gate.
assert(~isfile(outputFile),'Refusing to overwrite an observation report.');
assert(maxNumCompThreads==10,'Strict historical signatures require 10 threads.');
validateattributes(window,{'double'},{'finite','real','numel',3});
assert(window(1)<window(2) && window(3)>0);
a=read_view(firstCheckpoint);b=read_view(secondCheckpoint);
assert(abs(a.physicalTime-b.physicalTime)<=1e-10,'Physical endpoints must match.');
report=struct('firstCheckpoint',firstCheckpoint,'secondCheckpoint',secondCheckpoint, ...
    'window',window,'physicalTimeDifference',b.physicalTime-a.physicalTime, ...
    'firstNativeValidated',true,'secondNativeValidated',true, ...
    'secondToFirstNodes',direction(a,b,window), ...
    'firstToSecondNodes',direction(b,a,window), ...
    'interpretation',['Two fixed linear observation operators on the same fields. ' ...
    'Differences measure spatial plus transfer/evolution sensitivity and observation error. ' ...
    'No decision gate is changed and no observer is treated as ground truth.']);
fid=fopen(outputFile,'w');assert(fid>=0);cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('SPACE_OBSERVERS %s\n',jsonencode(report));
end

function v=read_view(file)
cp=ipm.output.readCheckpoint(file);s=cp.payload.state;
assert(all(ipm.output.trustedMask(cp.payload.log.history,s.config)));
assert(strcmp(s.config.transport.spatialDiscretization,'high_order') && ...
    strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic'));
Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);
Cx=exp(s.scale.logC_l);Comega=exp(s.scale.logC_omega);
v=struct('x',s.x(:)','y',s.y(:),'Dx',Dx,'Dy',Dy,'Cx',Cx, ...
    'shift',s.scale.X_shift,'physicalTime',s.scale.physicalTime, ...
    'px',(s.x(:)'-s.scale.X_shift)/Cx,'py',s.y(:)/Cx, ...
    'fields',struct('rho',s.rho/Comega,'omega',(s.rho*Dx')*Cx/Comega, ...
    'rhoY',(Dy*s.rho)*Cx/Comega));
end

function output=direction(reference,candidate,window)
xi=find(reference.px>=window(1) & reference.px<=window(2));
yi=find(reference.py>=0 & reference.py<=window(3));
assert(numel(xi)>=2 && numel(yi)>=2,'Insufficient reference nodes in window.');
[X,Y]=meshgrid(reference.px(xi),reference.py(yi));
QX=X*candidate.Cx+candidate.shift;QY=Y*candidate.Cx;
[RX,RY]=meshgrid(reference.x(xi),reference.y(yi));
w=widths(reference.py)*widths(reference.px)';w=w(yi,xi);
names={'rho','omega','rhoY'};
output=struct('referenceNodeCount',[numel(reference.x),numel(reference.y)], ...
    'windowNodeCount',[numel(xi),numel(yi)],'fields',struct);
for k=1:numel(names)
    name=names{k};native=reference.fields.(name);values=native(yi,xi);
    linear=interp2(candidate.x,candidate.y,candidate.fields.(name),QX,QY,'linear',NaN);
    hermite=ipm_accellab_tensor_hermite(candidate.fields.(name),candidate.x, ...
        candidate.y,candidate.Dx,candidate.Dy,QX,QY);
    identity=ipm_accellab_tensor_hermite(reference.fields.(name),reference.x, ...
        reference.y,reference.Dx,reference.Dy,RX,RY);
    identityError=max(abs(identity.value-values),[],'all')/max(abs(values),[],'all');
    assert(identityError<=1e-12,'Hermite native-node identity check failed.');
    output.fields.(name)=struct('linear',norms(linear-values,values,w), ...
        'hermite',norms(hermite.value-values,values,w), ...
        'betweenObservers',norms(hermite.value-linear,values,w), ...
        'referenceHermiteNodalIdentityRelativeInf',identityError);
end
end

function value=norms(error,reference,w)
assert(all(isfinite(error),'all') && all(w>0,'all'));
denL2=sum(reference.^2.*w,'all');denInf=max(abs(reference),[],'all');
assert(denL2>0 && denInf>0);
value=struct('relativeL2',sqrt(sum(error.^2.*w,'all')/denL2), ...
    'relativeInf',max(abs(error),[],'all')/denInf);
end

function w=widths(axis)
d=diff(axis(:));w=[d(1)/2;(d(1:end-1)+d(2:end))/2;d(end)/2];
end
