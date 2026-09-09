function report=ipm_accellab_compare_ABC_inner(registrationFile)
%IPM_ACCELLAB_COMPARE_ABC_INNER Fixed actual endpoint shape observation only.
reg=jsondecode(fileread(registrationFile));assert(maxNumCompThreads==10&&~isfolder(reg.outputDirectory));
mkdir(reg.outputDirectory);profile clear;profile on;
old=load(reg.comparatorMat,'report');prior=old.report;q=load(prior.factorialRegistration,'registration','plans');
assert(isequal(reg.wallInterval(:)',[-2,2])&&isequal(reg.verticalInterval(:)',[0,3]));
report=struct('kind','ABC_same_physical_time_C1_inner_shape_observer_v1', ...
    'registrationFile',registrationFile,'cases',struct([]),'pairs',struct([]), ...
    'samePhysicalEndpoint',prior.physicalFinalTime,'oldFixedPhysicalGatesReplaced',false, ...
    'oldAllScientificGatesPassed',prior.allScientificGatesPassed,'noLU',true,'noPDE',true, ...
    'pdeRateEvaluated',false,'nativeGaugeChanged',false,'trajectoryAdvanced',false);
views=cell(1,3);curves=cell(3,2);
for k=1:3
    oldrow=prior.cases(k);assert(strcmp(oldrow.status,'native_qualified'));
    assert(strcmp(oldrow.label,reg.labels{k})&&strcmp(oldrow.checkpointFile,reg.checkpointFiles{k}));
    plan=q.plans(k);
    if k==1&&plan.reuseExistingBaseline
        b=load(fullfile(q.registration.baselineDirectory,'registration.mat'),'o');expected=ipm.config.resolve(b.o);
    else
        expected=plan.config;
    end
    audit=mesh_audit_box_space_case(oldrow.resultFile,oldrow.checkpointFile,prior.physicalFinalTime,expected);
    assert(audit.passed&&isequaln(audit,oldrow.nativeAudit));
    cp=ipm.output.readCheckpoint(oldrow.checkpointFile);s=cp.payload.state;
    Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);omega=s.rho*Dx';
    Cx=exp(s.scale.logC_l);W=exp(s.scale.logC_omega);
    geometry=continuous_inner_geometry(omega,s.x,s.y,Dx,Dy,prior.physicalPeakWindow*Cx+s.scale.X_shift);
    assert(isequaln(geometry,oldrow.continuousGeometry.detail));
    v=struct('x',s.x,'y',s.y,'Dx',Dx,'Dy',Dy,'omega',omega,'g',geometry);
    physical=struct('peakX',(geometry.peakX-s.scale.X_shift)/Cx, ...
        'peak',geometry.peak*Cx/W,'wallWidth',geometry.wallWidth/Cx,'verticalWidth',geometry.verticalWidth/Cx);
    row=struct('label',oldrow.label,'step',s.step,'canonicalTime',s.scale.canonicalTime, ...
        'physicalTime',s.scale.physicalTime,'nodeCount',[numel(s.x),numel(s.y)], ...
        'nativeAudit',audit,'geometryExactlyMatchesPriorComparator',true, ...
        'nativeGeometry',geometry,'physicalGeometry',physical,'observerDifferences',struct([]));
    for direction=1:2
        interval=domain(reg,direction);n=reg.curvePointCounts(direction);
        query=linspace(interval(1),interval(2),n)';hermite=sample(v,query,direction,'hermite');linear=sample(v,query,direction,'linear');
        fine=linspace(interval(1),interval(2),2*n-1)';hf=sample(v,fine,direction,'hermite');lf=sample(v,fine,direction,'linear');
        [gauss,weights]=gauss_queries(v,[],interval,direction);hg=sample(v,gauss,direction,'hermite');lg=sample(v,gauss,direction,'linear');
        obs=norms(lg,hg,weights,linear,hermite,lf,hf);obs.direction=direction;
        row.observerDifferences=append(row.observerDifferences,obs);
        curves{k,direction}=struct('coordinate',query,'hermite',hermite,'linear',linear,'fineCoordinate',fine, ...
            'fineHermite',hf,'fineLinear',lf,'gaussQuery',gauss,'gaussWeights',weights, ...
            'gaussHermite',hg,'gaussLinear',lg);
        writematrix([query,hermite,linear],fullfile(reg.outputDirectory,sprintf('%s_%d_curve.csv',oldrow.label,direction)));
    end
    views{k}=v;report.cases=append(report.cases,row);clear cp s omega Dx Dy
end
assert(max([report.cases.physicalTime])-min([report.cases.physicalTime])<=1e-12);
for ids=[1,2;1,3;2,3]'
    a=ids(1);b=ids(2);row=struct('first',reg.labels{a},'second',reg.labels{b},'directions',struct([]));
    for direction=1:2
        interval=domain(reg,direction);[query,weights]=gauss_queries(views{a},views{b},interval,direction);
        entry=struct('direction',direction,'interval',interval,'linear',struct(),'hermite',struct(), ...
            'integration','Gauss4 on union of both native pullback cell boundaries; exact for represented squared polynomials in exact arithmetic');
        for method={'linear','hermite'}
            name=method{1};av=sample(views{a},query,direction,name);bv=sample(views{b},query,direction,name);
            ca=curves{a,direction};cb=curves{b,direction};fineName=['fine',upper(name(1)),name(2:end)];
            entry.(name)=norms(bv,av,weights,cb.(name),ca.(name),cb.(fineName),ca.(fineName));
        end
        row.directions=append(row.directions,entry);
    end
    report.pairs=append(report.pairs,row);
end
report.samplingRefinementMaximumChange=0;
for row=report.pairs
    for d=row.directions
        report.samplingRefinementMaximumChange=max([report.samplingRefinementMaximumChange, ...
            d.linear.sampledRelativeInfRefinementChange,d.hermite.sampledRelativeInfRefinementChange]);
    end
end
report.finiteTimeShapeObservationCompleted=true;
report.interpretation='Each case uses its own continuous peak and connected 0.9 widths. Linear and Hermite use that same geometry and peak; their difference isolates interpolation of the field, not redefinition of geometry. No claim of stationary limit, closed evolved tail, or disappearance of original failed fixed-physical-window errors.';
profile off;profileInfo=profile('info');files={profileInfo.FunctionTable.FileName};
for suffix={'+mesh/build.m','+evolve/flow.m','+field/velocity.m','+field/poisson.m', ...
        '+evolve/initializeScaling.m','+remesh/transfer.m','+evolve/advance.m','+output/restoreCheckpoint.m'}
    assert(~any(endsWith(files,suffix{1})),['Forbidden call: ',suffix{1}]);
end
report.callGraphNoLUVerified=true;save(fullfile(reg.outputDirectory,'profiles.mat'),'report','curves','profileInfo','-v7.3');
fid=fopen(fullfile(reg.outputDirectory,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
disp(jsonencode(struct('completed',true,'pairs',report.pairs,'sampledInfRefinementMaximumChange',report.samplingRefinementMaximumChange)));
end
function y=sample(v,z,direction,method)
g=v.g;
if direction==1,x=g.peakX+g.wallWidth*z;yq=zeros(size(z));
else,x=g.peakX+zeros(size(z));yq=g.verticalWidth*z;end
assert(all(x>=v.x(1)&x<=v.x(end)&yq>=v.y(1)&yq<=v.y(end)));
if strcmp(method,'hermite')
    f=ipm_accellab_tensor_hermite(v.omega,v.x,v.y,v.Dx,v.Dy,x,yq);y=f.value/g.peak;
else
    y=interp2(v.x,v.y,v.omega,x,yq,'linear')/g.peak;
end
assert(all(isfinite(y)));
end
function [z,w]=gauss_queries(a,b,interval,direction)
knots=[interval(:);own_knots(a,interval,direction)];
if ~isempty(b),knots=[knots;own_knots(b,interval,direction)];end
knots=unique(knots);mid=(knots(1:end-1)+knots(2:end))/2;half=diff(knots)/2;
t=[-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7), ...
    sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
q=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
z=reshape(mid+half*t,[],1);w=reshape(half*q,[],1);assert(all(w>0));
end
function z=own_knots(v,interval,direction)
if direction==1,z=(v.x(:)-v.g.peakX)/v.g.wallWidth;else,z=v.y(:)/v.g.verticalWidth;end
z=z(z>interval(1)&z<interval(2));
end
function o=norms(b,a,w,bc,ac,bf,af)
e=b-a;coarse=max(abs(bc-ac))/max(abs(ac));fine=max(abs(bf-af))/max(abs(af));
o=struct('absoluteL2',sqrt(sum(w.*e.^2)),'relativeL2ToFirst',sqrt(sum(w.*e.^2)/sum(w.*a.^2)), ...
    'relativeL2ToSecond',sqrt(sum(w.*e.^2)/sum(w.*b.^2)), ...
    'sampledAbsoluteInf',max(abs(bf-af)),'sampledRelativeInfToFirst',fine, ...
    'sampledRelativeInfToSecond',max(abs(bf-af))/max(abs(bf)), ...
    'coarseSampledRelativeInfToFirst',coarse,'sampledRelativeInfRefinementChange',abs(fine-coarse), ...
    'infIsUniformSampleObservationNotPolynomialExtremumProof',true);
end
function interval=domain(reg,direction)
if direction==1,interval=reg.wallInterval(:)';else,interval=reg.verticalInterval(:)';end
end
function a=append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
