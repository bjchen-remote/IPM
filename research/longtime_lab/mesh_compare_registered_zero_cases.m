function report=mesh_compare_registered_zero_cases(factorialFile,timeFile,outDir)
%MESH_COMPARE_REGISTERED_ZERO_CASES Actual A/B/C/D/E only; no LU or PDE.
assert(maxNumCompThreads==10&&~isfolder(outDir));mkdir(outDir);
issues=checkcode(which(mfilename),'-id');assert(isempty(issues),jsonencode(issues));
q=load(factorialFile,'registration','plans');t=load(timeFile,'registration','plans');reg=q.registration;
assert(reg.physicalFinalTime==t.registration.physicalFinalTime&& ...
    isequal(reg.physicalWindows,t.registration.physicalWindows)&&isequaln(reg.limits,t.registration.limits));
assert(isequal(reg.physicalWindows,[0,.8,.25;0,2,1])&&isequal(reg.physicalPeakWindow,[.1,.6])&& ...
    isequal(reg.commonQueryShape,[321,161]));
expectedLabels={'A_H8_T32','B_H16_T32','C_H8_T42','D_H16_T42'};
assert(numel(q.plans)==4&&isequal({q.plans.label},expectedLabels)&&isscalar(t.plans)&& ...
    strcmp(t.plans.label,'E_H8_T32_HALF_TIME'),'ipm:ZeroCaseLabels','Exact registered A/B/C/D/E order is required.');
assert(isequal(t.registration.physicalPeakWindow,reg.physicalPeakWindow), ...
    'ipm:ZeroCaseWindow','The time-control peak window must match the factorial registration.');
parent=load(t.registration.baselineFactorialRegistration,'registration');
assert(isequaln(parent.registration,reg),'ipm:ZeroCaseParent','Time control must reference the identical factorial observation registration.');
if isfield(t.registration,'commonQueryShape')
    assert(isequal(t.registration.commonQueryShape,reg.commonQueryShape),'ipm:ZeroCaseWindow','Common query shapes must match.');
end
plans=[q.plans,t.plans];
registration=struct('kind','registered_original_zero_ABCDE_observer_v1','factorialRegistration',factorialFile, ...
    'timeRegistration',timeFile,'physicalFinalTime',reg.physicalFinalTime,'physicalWindows',reg.physicalWindows, ...
    'physicalPeakWindow',reg.physicalPeakWindow,'commonQueryShape',reg.commonQueryShape,'limits',reg.limits, ...
    'physicalRateFormula','[history.common.c_l(end), history.common.c_omega(end)] * exp(state.scale.logC_l-logC_omega)', ...
    'sourceStepPhysicalClockSpeedUsed',false,'missingCasesRemainPending',true, ...
    'nativeResultExpectedConfigRequired',true,'legacyAndAuxiliaryRemainSeparate',true, ...
    'continuousObserverCaveat','Absolute/covariance tests passed; per-phase refinement ratio >=8 did not all pass.', ...
    'interpretation','Box, resolution, time step and automatic mesh-event sensitivity; no pure-box or convergence-order claim.', ...
    'noLU',true,'noFlow',true,'noPDE',true);
write_json(fullfile(outDir,'registration.json'),registration);
profile clear;profile on;
report=registration;report.cases=struct([]);views=cell(1,5);queries=cell(1,5);
for k=1:5
    p=plans(k);row=struct('label',p.label,'status','pending','resultFile','','checkpointFile','', ...
        'nativeAudit',struct(),'observables',struct(),'coverage',false(2,1), ...
        'continuousGeometry',struct(),'ownQueryObserverDifferences',struct([]),'failure',struct());
    if k==1 && p.reuseExistingBaseline
        resultFile=p.baselineResultFile;checkpointFile=p.baselineCheckpointFile;
        a=load(fullfile(reg.baselineDirectory,'registration.mat'),'o');expected=ipm.config.resolve(a.o);
    else
        resultFile=p.options.resultFile;checkpointFile='';expected=p.config;
    end
    row.resultFile=resultFile;
    if ~isfile(resultFile)
        report.cases=append_row(report.cases,row);continue
    end
    try
        r=ipm.output.validate(resultFile);
        if isempty(checkpointFile)&&isfield(r.metadata,'latestCheckpointFile'),checkpointFile=char(r.metadata.latestCheckpointFile);end
        row.checkpointFile=checkpointFile;
        if isempty(checkpointFile)||~isfile(checkpointFile)
            row.failure=struct('identifier','missing_paired_native_checkpoint','message','Result exists; immutable paired native checkpoint is pending.');
            report.cases=append_row(report.cases,row);continue
        end
        audit=mesh_audit_box_space_case(resultFile,checkpointFile,reg.physicalFinalTime,expected);
        row.nativeAudit=audit;
        assert(audit.passed,'ipm:ZeroCaseNativeRejected','Actual endpoint/native/expected-config gates failed.');
        cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;
        [v,observation]=read_view(s,cp.payload.log.history,reg.physicalPeakWindow);
        v.label=p.label;v.caseId=s.runMetadata.caseId;v.resultFile=resultFile;v.checkpointFile=checkpointFile;
        row.observables=observation;row.continuousGeometry=v.geometry;row.status='native_qualified';
        query=cell(1,2);
        for w=1:2
            box=reg.physicalWindows(w,:);row.coverage(w)=covered(v,box);
            if row.coverage(w)
                query{w}=query_fields(v,box,reg.commonQueryShape);
                row.ownQueryObserverDifferences=append_row(row.ownQueryObserverDifferences,query{w}.observerDifference);
            else
                row.ownQueryObserverDifferences=append_row(row.ownQueryObserverDifferences, ...
                    struct('window',box,'rho',struct(),'rhoX',struct(),'rhoY',struct()));
            end
        end
        save(fullfile(outDir,[p.label,'_observations.mat']),'observation','query','-v7.3');
        views{k}=v;queries{k}=query;
    catch e
        views{k}=[];queries{k}=[];
        row.status='rejected';row.failure=struct('identifier',e.identifier,'message',e.message);
    end
    report.cases=append_row(report.cases,row);
end
% A self-comparison is explicit; A is never relabelled as B/C/D/E.
report.selfTest=struct('status','pending','caseLabel',plans(1).label,'comparisons',struct([]));
if ~isempty(views{1})
    selfRows=struct([]);maximum=0;
    for w=1:2
        c=pair_comparison(views{1},views{1},reg.physicalWindows(w,:),reg.limits);
        selfRows=append_row(selfRows,c);assert(strcmp(c.status,'compared'));
        for direction={'secondToFirst','firstToSecond'}
            fields=c.hermite.(direction{1}).fields;
            for name={'rho','rhoX','rhoY'}
                z=fields.(name{1});maximum=max([maximum,z.linear.relativeInf,z.hermite.relativeInf]);
            end
        end
        assert(c.originalLinearLegacyPassed&&c.auxiliaryHermiteContinuousPassed);
    end
    assert(maximum<=1e-12,'ipm:ObserverIdentity','Actual A self-comparison exceeds nodal identity tolerance.');
    report.selfTest=struct('status','passed','caseLabel',plans(1).label,'comparisons',selfRows, ...
        'maximumRelativeNodalIdentityError',maximum,'identityTolerance',1e-12, ...
        'commonQueryLinearHermiteDifferenceIsSeparate',true,'otherCasesFabricated',false);
end
pairs=[1,2;3,4;1,3;2,4;1,5];edgeNames={'box_T32_B_minus_A','box_T42_D_minus_C', ...
    'target_H8_C_minus_A','target_H16_D_minus_B','time_E_minus_A'};
report.edges=struct([]);
for k=1:size(pairs,1)
    ids=pairs(k,:);edge=struct('name',edgeNames{k},'first',plans(ids(1)).label,'second',plans(ids(2)).label, ...
        'status','pending','windows',struct([]));
    if all(~cellfun(@isempty,views(ids)))
        edge.status='compared';
        for w=1:2
            edge.windows=append_row(edge.windows,pair_comparison(views{ids(1)},views{ids(2)},reg.physicalWindows(w,:),reg.limits));
        end
    elseif any(strcmp({report.cases(ids).status},'rejected'))
        edge.status='blocked_by_rejected_native_case';
    end
    report.edges=append_row(report.edges,edge);
end
report.interaction=struct('status','pending','formula','D-C-B+A on the same fixed physical query nodes', ...
    'windows',struct([]),'scalar',struct());
if all(~cellfun(@isempty,views(1:4)))
    report.interaction.status='compared';interactionFields=cell(1,2);
    for w=1:2
        if ~all(arrayfun(@(v)v.coverage(w),report.cases(1:4)))
            report.interaction.status='coverage_failed';continue
        end
        grid=queries{1}{w};entry=struct('window',reg.physicalWindows(w,:),'linear',struct(),'hermite',struct());
        fields=struct();
        for observer={'linear','hermite'}
            method=observer{1};fields.(method)=struct();
            for name={'rho','rhoX','rhoY'}
                key=name{1};a=queries{1}{w}.(method).(key);b=queries{2}{w}.(method).(key);
                c=queries{3}{w}.(method).(key);d=queries{4}{w}.(method).(key);
                signed=d-c-b+a;fields.(method).(key)=signed;
                entry.(method).(key)=norms(signed,a,grid.weights);
            end
        end
        interactionFields{w}=fields;report.interaction.windows=append_row(report.interaction.windows,entry);
    end
    for name={'physicalGradInf','physicalRhoXInf','physicalQuadraticPeak','physicalRateCL','physicalRateCOmega', ...
            'canonicalRateCL','canonicalRateCOmega','legacyPhysicalWidthX','legacyPhysicalWidthY', ...
            'continuousPhysicalWidthX','continuousPhysicalWidthY'}
        key=name{1};v=arrayfun(@(c)c.observables.(key),report.cases(1:4));
        z=v(4)-v(3)-v(2)+v(1);report.interaction.scalar.(key)=struct('signed',z,'relativeToA',z/max(abs(v(1)),realmin));
    end
    save(fullfile(outDir,'signed_interaction_fields.mat'),'interactionFields','-v7.3');
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for forbidden={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,forbidden{1})));
end
report.observationCompleted=true;report.allNativeEndpointsAvailable=all(strcmp({report.cases.status},'native_qualified'));
report.allScientificGatesPassed=false;
if report.allNativeEndpointsAvailable
    gates=true;
    for edge=report.edges
        for item=edge.windows
            gates=gates&&strcmp(item.status,'compared')&&item.originalLinearLegacyPassed&&item.auxiliaryHermiteContinuousPassed;
        end
    end
    report.allScientificGatesPassed=gates;
end
save(fullfile(outDir,'report.mat'),'report','profileInfo','-v7.3');write_json(fullfile(outDir,'report.json'),report);
fprintf('ZERO_CASE_OBSERVER_COMPLETE self=%s actual=%d pending=%d allScientific=%d noLU=1\n', ...
    report.selfTest.status,nnz(strcmp({report.cases.status},'native_qualified')),nnz(strcmp({report.cases.status},'pending')),report.allScientificGatesPassed);
end

function [v,o]=read_view(s,h,peakWindow)
Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);
Cx=exp(s.scale.logC_l);Comega=exp(s.scale.logC_omega);amplitude=exp(s.scale.logC_l-s.scale.logC_omega);
omega=s.rho*Dx';x=(s.x-s.scale.X_shift)/Cx;y=s.y/Cx;
o=struct('physicalTime',s.scale.physicalTime,'canonicalTime',s.scale.canonicalTime,'step',s.step, ...
    'physicalGradInf',h.common.physicalGradInf(end),'physicalRhoXInf',h.common.physicalRhoXInf(end), ...
    'physicalQuadraticPeak',h.common.physicalQuadraticPeak(end), ...
    'physicalRateCL',h.common.c_l(end)*amplitude,'physicalRateCOmega',h.common.c_omega(end)*amplitude, ...
    'canonicalRateCL',h.common.canonicalCL(end),'canonicalRateCOmega',h.common.canonicalCOmega(end), ...
    'legacyPhysicalWidthX',h.mesh.trackedWallCoreWidth(end)/Cx, ...
    'legacyPhysicalWidthY',h.mesh.trackedVerticalCoreWidth(end)/Cx, ...
    'continuousPhysicalWidthX',NaN,'continuousPhysicalWidthY',NaN, ...
    'nodeCount',[numel(x),numel(y)],'physicalDomain',[x(1),x(end),y(end)],'Cx',Cx,'Comega',Comega, ...
    'cfl',s.config.time.cfl,'maxDt',s.config.time.maxDt,'remeshCount',s.remeshCount, ...
    'configuredWenoEpsilon',s.config.transport.wenoEpsilon, ...
    'epsilonInPhysicalRhoAmplitudeUnitsOnly',s.config.transport.wenoEpsilon/Comega^2, ...
    'cumulativeAbsolutePeakJump',s.runMetadata.autonomousMesh.cumulativeAbsolutePeakJump);
v=struct('x',x,'y',y,'nativeX',s.x,'nativeY',s.y,'Dx',Dx,'Dy',Dy,'Cx',Cx,'shift',s.scale.X_shift, ...
    'rho',s.rho/Comega,'rhoX',omega*amplitude,'rhoY',(Dy*s.rho)*amplitude, ...
    'physicalGradientMaximum',o.physicalRhoXInf,'physicalQuadraticPeak',o.physicalQuadraticPeak, ...
    'physicalCoreWidthX',o.legacyPhysicalWidthX,'physicalCoreWidthY',o.legacyPhysicalWidthY, ...
    'cLParentUnits',o.physicalRateCL,'cOmegaParentUnits',o.physicalRateCOmega,'observation',o, ...
    'geometry',struct('available',false,'detail',struct(),'failure',struct()));
try
    g=continuous_inner_geometry(omega,s.x,s.y,Dx,Dy,peakWindow*Cx+s.scale.X_shift);
    v.geometry.available=true;v.geometry.detail=g;
    o.continuousPhysicalWidthX=g.wallWidth/Cx;o.continuousPhysicalWidthY=g.verticalWidth/Cx;
catch e
    v.geometry.failure=struct('identifier',e.identifier,'message',e.message);
end
v.observation=o;
end
function yes=covered(v,w)
yes=v.x(1)<=w(1)&&v.x(end)>=w(2)&&v.y(1)<=0&&v.y(end)>=w(3);
end
function c=pair_comparison(a,b,w,limits)
c=struct('window',w,'status','coverage_failed','original',struct(),'hermite',struct(), ...
    'scalar',struct(),'originalLinearLegacyPassed',false,'auxiliaryHermiteContinuousPassed',false);
if ~covered(a,w)||~covered(b,w),return;end
assert(abs(a.observation.physicalTime-b.observation.physicalTime)<=2e-12);
c.status='compared';ab=mesh_compare_fresh_box_data(a,b,w,limits);ba=mesh_compare_fresh_box_data(b,a,w,limits);
c.original=struct('secondRelativeToFirst',ab,'firstRelativeToSecond',ba,'bothDirectionsPassed',ab.passed&&ba.passed);
c.originalLinearLegacyPassed=c.original.bothDirectionsPassed;
f=direction(a,b,w);r=direction(b,a,w);c.hermite=struct('secondToFirst',f,'firstToSecond',r);
physicalRates=two_rate_gate([a.cLParentUnits,a.cOmegaParentUnits],[b.cLParentUnits,b.cOmegaParentUnits],limits);
canonicalRates=two_rate_gate([a.observation.canonicalRateCL,a.observation.canonicalRateCOmega], ...
    [b.observation.canonicalRateCL,b.observation.canonicalRateCOmega],limits);
peak=both_relative([a.physicalGradientMaximum,a.physicalQuadraticPeak],[b.physicalGradientMaximum,b.physicalQuadraticPeak]);
legacy=both_relative([a.physicalCoreWidthX,a.physicalCoreWidthY],[b.physicalCoreWidthX,b.physicalCoreWidthY]);
continuous=both_relative([a.observation.continuousPhysicalWidthX,a.observation.continuousPhysicalWidthY], ...
    [b.observation.continuousPhysicalWidthX,b.observation.continuousPhysicalWidthY]);
c.scalar=struct('physicalRhoXInfAndQuadraticPeak',peak,'physicalGradInfDiagnosticOnly', ...
    both_relative(a.observation.physicalGradInf,b.observation.physicalGradInf), ...
    'legacyWidths',legacy,'continuousWidths',continuous,'physicalRateGates',physicalRates, ...
    'canonicalRateGatesReportedSeparately',canonicalRates,'continuousGeometryAvailable',a.geometry.available&&b.geometry.available);
fieldPass=true;
for name={'rho','rhoX','rhoY'}
    key=name{1};bound=limits.cropGradientRelativeL2Inf;if strcmp(key,'rho'),bound=limits.cropRhoRelativeL2Inf;end
    fieldPass=fieldPass&&max([f.fields.(key).hermite.relativeL2,f.fields.(key).hermite.relativeInf, ...
        r.fields.(key).hermite.relativeL2,r.fields.(key).hermite.relativeInf])<=bound;
end
c.auxiliaryHermiteContinuousPassed=fieldPass&&a.geometry.available&&b.geometry.available&& ...
    all(peak.maximum<=limits.cropPeakRelative)&&all(continuous.maximum<=limits.cropPhysicalWidthRelative)&&physicalRates.bothPassed;
end
function d=direction(a,b,w)
xi=find(a.x>=w(1)&a.x<=w(2));yi=find(a.y>=0&a.y<=w(3));assert(numel(xi)>=5&&numel(yi)>=5);
[X,Y]=meshgrid(a.x(xi),a.y(yi));weights=widths(a.y(yi))*widths(a.x(xi))';
[NX,NY]=meshgrid(a.nativeX(xi),a.nativeY(yi));QX=X*b.Cx+b.shift;QY=Y*b.Cx;
if isequal(a.nativeX,b.nativeX)&&isequal(a.nativeY,b.nativeY)&&a.Cx==b.Cx&&a.shift==b.shift
    QX=NX;QY=NY;
end
d=struct('referenceNodeCount',[numel(a.x),numel(a.y)],'windowNodeCount',[numel(xi),numel(yi)],'fields',struct());
for name={'rho','rhoX','rhoY'}
    key=name{1};values=a.(key)(yi,xi);linear=interp2(b.x,b.y,b.(key),X,Y,'linear',NaN);
    h=ipm_accellab_tensor_hermite(b.(key),b.nativeX,b.nativeY,b.Dx,b.Dy,QX,QY);
    identity=ipm_accellab_tensor_hermite(a.(key),a.nativeX,a.nativeY,a.Dx,a.Dy,NX,NY);
    identityError=max(abs(identity.value-values),[],'all')/max(abs(values),[],'all');assert(identityError<=1e-12);
    d.fields.(key)=struct('linear',norms(linear-values,values,weights),'hermite',norms(h.value-values,values,weights), ...
        'betweenObservers',norms(h.value-linear,values,weights),'nativeHermiteIdentityRelativeInf',identityError);
end
end
function q=query_fields(v,w,shape)
x=linspace(w(1),w(2),shape(1));y=linspace(0,w(3),shape(2))';[X,Y]=meshgrid(x,y);
q=struct('x',x,'y',y,'weights',widths(y)*widths(x)','linear',struct(),'hermite',struct(), ...
    'observerDifference',struct('window',w,'rho',struct(),'rhoX',struct(),'rhoY',struct()));
for name={'rho','rhoX','rhoY'}
    key=name{1};l=interp2(v.x,v.y,v.(key),X,Y,'linear',NaN);
    h=ipm_accellab_tensor_hermite(v.(key),v.nativeX,v.nativeY,v.Dx,v.Dy,X*v.Cx+v.shift,Y*v.Cx);
    q.linear.(key)=l;q.hermite.(key)=h.value;
    q.observerDifference.(key)=norms(h.value-l,l,q.weights);
end
end
function r=two_rate_gate(a,b,limits)
difference=abs(b-a);forward=max(limits.cropRateAbsolute,limits.cropRateRelative*abs(a));
reverse=max(limits.cropRateAbsolute,limits.cropRateRelative*abs(b));
r=struct('absoluteDifference',difference,'secondRelativeToFirstLimit',forward,'firstRelativeToSecondLimit',reverse, ...
    'secondRelativeToFirstPassed',all(difference<=forward),'firstRelativeToSecondPassed',all(difference<=reverse), ...
    'bothPassed',all(difference<=forward)&all(difference<=reverse));
end
function r=both_relative(a,b)
r=struct('secondRelativeToFirst',abs(b-a)./max(abs(a),realmin),'firstRelativeToSecond',abs(b-a)./max(abs(b),realmin));
r.maximum=max(r.secondRelativeToFirst,r.firstRelativeToSecond);
end
function n=norms(error,reference,w)
assert(all(isfinite(error),'all')&&all(isfinite(reference),'all')&&all(w>0,'all'));
den=sum(reference.^2.*w,'all');di=max(abs(reference),[],'all');assert(den>0&&di>0);
n=struct('relativeL2',sqrt(sum(error.^2.*w,'all')/den),'relativeInf',max(abs(error),[],'all')/di, ...
    'absoluteL2',sqrt(sum(error.^2.*w,'all')),'absoluteInf',max(abs(error),[],'all'));
end
function w=widths(x)
d=diff(x(:));w=[d(1)/2;(d(1:end-1)+d(2:end))/2;d(end)/2];
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
