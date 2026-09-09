function report=ipm_accellab_audit_physical_far_coverage(registrationFile)
%IPM_ACCELLAB_AUDIT_PHYSICAL_FAR_COVERAGE Saved paired fields only; no flow.
reg=jsondecode(fileread(registrationFile));out=reg.outputDirectory;
assert(maxNumCompThreads==10 && ~isfolder(out));mkdir(out);
profile clear;profile on;
original=load(reg.initialABMat,'data','report');data=cell(1,3);
report=struct('kind','saved_native_physical_far_coverage_v1','registrationFile',registrationFile, ...
    'cases',{{}},'leadingA',reg.leadingA,'farPrimaryBands',reg.primaryPhysicalBands, ...
    'auxiliaryFractionalBands',reg.auxiliaryFractionalBands,'conditionalModelOnly',true, ...
    'addedToNativeRate',false,'globalFarFieldBoundProven',false,'LUCount',0,'PDECount',0,'restoreCount',0);
for k=1:3
    item=reg.cases(k);cp=ipm.output.readCheckpoint(item.checkpointFile);
    loaded=load(item.resultFile,'result');r=loaded.result;s=cp.payload.state;
    assert(isequaln(r.history,cp.payload.log.history)&&isequaln(r.state.rho,s.rho));
    assert(isequaln(r.grid.x,s.x)&&isequaln(r.grid.y,s.y)&&isequaln(r.config,s.config));
    [trusted,~]=ipm.output.trustedMask(r.history,s.config);assert(all(trusted));
    if item.historyIndex==1
        h=r.history.common;m=r.history.mesh;C=h.C_l(1);W=h.C_omega(1);shift=h.X_shift(1);
        assert(h.acceptedStep(1)==0 && h.t(1)==0 && h.canonicalTau(1)==0 && h.physicalTime(1)==0);
        assert(isequaln([C,h.C_y(1),W,shift],[1,1,1,0]));
        assert(r.snapshots.normalizedTime(1)==0 && r.snapshots.physicalTime(1)==0 && r.snapshots.canonicalTime(1)==0);
        x=r.snapshots.x{1};y=r.snapshots.y{1};rho=r.snapshots.rho{1};
        init=original.data{k}.initialization;
        assert(isequaln(x,init.selectedBaseX)&&isequaln(y,init.selectedBaseY));
        [X,Y]=meshgrid(x,y);temp=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y),'symmetryMode',r.config.physics.symmetryMode);
        assert(isequaln(rho,ipm.field.initialDensity(temp,r.config.physics)));
        Dx=ipm.mesh.fdMatrix(x,1,7);Omega=rho*Dx';
        assert(max(abs(Omega),[],'all')==h.physicalRhoXInf(1));
        assert(isequaln([h.c_l(1),h.c_omega(1)], ...
            [original.report.cases{k}.initialObservation.canonicalCL,original.report.cases{k}.initialObservation.canonicalCOmega]));
        v=[];physicalTime=0;canonicalTime=0;step=0;nativeCL=h.c_l(1);nativeCW=h.c_omega(1);
        sourceMeaning='Saved actual zero-time snapshot, exactly paired with original selected axes and analytic initial resampling. Omega recomputed by the maintained paired seven-point derivative. No t0 velocity array was saved.';
        nodes=[numel(x),numel(y)];core=[m.coreGridPoints(1),m.verticalCoreGridPoints(1)];
    else
        h=r.history.common;assert(item.historyIndex==-1 && s.scale.canonicalTime==8);
        C=exp(s.scale.logC_l);W=exp(s.scale.logC_omega);shift=s.scale.X_shift;
        assert(r.scale.Cx==C && r.scale.Cy==C && r.scale.Comega==W && shift==0);
        assert(isequaln(r.grid.physicalX,s.x/C)&&isequaln(r.grid.physicalY,s.y/C));
        Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);
        assert(isequaln(r.state.omega,r.state.rho*Dx'));
        factor=exp(s.scale.logC_l-s.scale.logC_omega);
        assert(isequaln(r.physical.rhoX,factor*r.state.omega));
        inverseW=exp(-s.scale.logC_omega);
        assert(isequaln(r.physical.u1,inverseW*r.state.velocity.x)&&isequaln(r.physical.u2,inverseW*r.state.velocity.y));
        x=r.grid.physicalX;y=r.grid.physicalY;rho=r.physical.rho;Omega=r.physical.rhoX;
        u1=r.state.velocity.x;u2=r.state.velocity.y;
        v=struct('u1',r.physical.u1,'u2',r.physical.u2, ...
            'u1x',factor*(u1*Dx'),'u1y',factor*(Dy*u1), ...
            'u2x',factor*(u2*Dx'),'u2y',factor*(Dy*u2), ...
            'u2xy',exp(2*s.scale.logC_l-s.scale.logC_omega)*((Dy*u2)*Dx'));
        assert(all(v.u2(1,:)==0));
        physicalTime=s.scale.physicalTime;canonicalTime=s.scale.canonicalTime;step=s.step;
        nativeCL=h.c_l(end);nativeCW=h.c_omega(end);nodes=[numel(x),numel(y)];
        core=[r.history.mesh.coreGridPoints(end),r.history.mesh.verticalCoreGridPoints(end)];
        sourceMeaning='Saved endpoint physical Omega/u1/u2, exact reconstruction pairing to the native computational fields; velocity derivatives use the maintained native axes and physical derivative scale factors.';
    end
    Hc=r.config.grid.xlim(2);Hp=Hc/C;assert(x(1)==-Hp && x(end)==Hp && y(end)==Hp/2);
    [X,Y]=meshgrid(x,y);R=hypot(X,Y);OmegaLeading=X.^7./(X.^8+Y.^8);
    support=false(size(X));support(1:end-3,4:end-3)=true;
    row=struct('label',item.label,'resultFile',item.resultFile,'checkpointFile',item.checkpointFile, ...
        'strictNativeValidated',true,'completeResultNativePairing',true,'allHistoryTrusted',true, ...
        'observedStep',step,'canonicalTime',canonicalTime,'physicalTime',physicalTime, ...
        'Cx',C,'Comega',W,'Xshift',shift,'Hcanonical',Hc,'Hphysical',Hp,'YmaxPhysical',Hp/2, ...
        'maximumPhysicalRadiusInBox',hypot(Hp,Hp/2), ...
        'nativeCanonicalCL',nativeCL,'nativeCanonicalCOmega',nativeCW, ...
        'nativePhysicalRateCL',C/W*nativeCL,'nativePhysicalRateCOmega',C/W*nativeCW, ...
        'conditionalLeadingTailScalePhysical',reg.leadingA/Hp, ...
        'conditionalLeadingTailScaleCanonical',W*reg.leadingA/Hc, ...
        'physicalTailScaleDividedByNativePhysicalCL',(reg.leadingA/Hp)/(C/W*nativeCL), ...
        'leadingScaleValidAsCorrection',false,'velocityAtObservedTimeAvailable',~isempty(v), ...
        'sourceMeaning',sourceMeaning,'nodes',nodes,'nativeCore',core,'bands',{{}});
    windows=[reg.primaryPhysicalBands;reg.auxiliaryFractionalBands*Hp];
    masks=cell(1,4);
    for j=1:4
        band=windows(j,:);q=band_stats(band,reg,x,y,X,Y,R,Omega,OmegaLeading,support,v);
        if j<=2,q.kind='fixed_physical_primary';else,q.kind='fraction_of_current_box_auxiliary';end
        row.bands{j}=q;masks{j}=R>=band(1)&R<=band(2)&X>0&Y<=reg.maximumYSlope*R&support;
    end
    row.primaryFarCoveragePassed=all(cellfun(@(q)q.eligibleFarWindow,row.bands(1:2)));
    row.globalFarFieldClaim=false;report.cases{k}=row;
    data{k}=struct('x',x,'y',y,'rho',rho,'Omega',Omega,'OmegaLeading',OmegaLeading,'velocity',v,'masks',{masks});
    clear cp loaded r s Dx Dy rho Omega v X Y R
end
profile off;profileInfo=profile('info');names={profileInfo.FunctionTable.FunctionName};
for key={'ipm.mesh.build','ipm.evolve.flow','ipm.field.velocity','ipm.field.poisson', ...
        'ipm.output.restoreCheckpoint','ipm.evolve.advance','ipm.solve','ipm.remesh.transfer'}
    assert(~any(strcmp(names,key{1})),['Forbidden call: ',key{1}]);
end
report.noLUCallGraphChecked=true;
save(fullfile(out,'report.mat'),'report','data','profileInfo','-v7.3');write_json(fullfile(out,'report.json'),report);
for k=1:3
    q=report.cases{k};disp(jsonencode(struct('label',q.label,'physicalTime',q.physicalTime,'Cx',q.Cx, ...
        'Hphysical',q.Hphysical,'tailScale',q.conditionalLeadingTailScalePhysical, ...
        'nativePhysicalCL',q.nativePhysicalRateCL,'primaryFarCoveragePassed',q.primaryFarCoveragePassed, ...
        'velocityAvailable',q.velocityAtObservedTimeAvailable,'bands',{q.bands})));
end
end

function s=band_stats(band,reg,x,y,X,Y,R,Omega,leading,support,v)
mask=R>=band(1)&R<=band(2)&X>0&Y<=reg.maximumYSlope*R&support;
[row,col]=find(mask);s=struct('physicalRadiusBand',band,'pointCount',nnz(mask), ...
    'distinctX',numel(unique(col)),'distinctY',numel(unique(row)), ...
    'minimumDeclaredRadiusAtLeastFour',band(1)>=reg.minimumFarPhysicalRadius, ...
    'velocityAvailable',~isempty(v),'fullSpaceBoundClaim',false);
[qr,qt]=meshgrid(linspace(band(1),band(2),reg.coverageRadialSamples), ...
    linspace(0,asin(reg.maximumYSlope),reg.coverageAngularSamples));
qx=qr.*cos(qt);qy=qr.*sin(qt);
covered=qx>=x(4)&qx<=x(end-3)&qy>=0&qy<=y(end-3);
s.nominalQueryCoverageFraction=nnz(covered)/numel(covered);
s.eligibleFarWindow=s.minimumDeclaredRadiusAtLeastFour && all(covered,'all') && ...
    s.distinctX>=reg.minimumDistinctX && s.distinctY>=reg.minimumDistinctY;
s.finiteNativeSourceData=all(isfinite(Omega(mask)))&&all(isfinite(leading(mask)));
s.sampledRadiusRange=[NaN,NaN];s.sourceRelativeToLeadingInf=NaN;s.rSourceDifferenceInf=NaN;
s.rSquaredSourceDifferenceInf=NaN;s.maximumNativeSource=NaN;
s.velocitySpeedInf=NaN;s.rVelocityGradientInf=NaN;s.rAbsU2OverYInf=NaN;s.rSquaredAbsU2XOverYInf=NaN;
s.wallRatiosUseMaintainedOneSidedDerivativeLimit=false;
if ~any(mask,'all'),s.coverageStatus='none';return;end
if all(covered,'all'),s.coverageStatus='full_registered_query_support';else,s.coverageStatus='partial';end
s.sampledRadiusRange=[min(R(mask)),max(R(mask))];delta=Omega-leading;
s.sourceRelativeToLeadingInf=max(abs(delta(mask))./abs(leading(mask)));
s.rSourceDifferenceInf=max(R(mask).*abs(delta(mask)));
s.rSquaredSourceDifferenceInf=max(R(mask).^2.*abs(delta(mask)));
s.maximumNativeSource=max(abs(Omega(mask)));
if isempty(v),return;end
gradient=max(cat(3,abs(v.u1x),abs(v.u1y),abs(v.u2x),abs(v.u2y)),[],3);
s.velocitySpeedInf=max(hypot(v.u1(mask),v.u2(mask)));
s.rVelocityGradientInf=max(R(mask).*gradient(mask));
first=R.*abs(v.u2)./Y;second=R.^2.*abs(v.u2x)./Y;
first(1,:)=R(1,:).*abs(v.u2y(1,:));second(1,:)=R(1,:).^2.*abs(v.u2xy(1,:));
s.rAbsU2OverYInf=max(first(mask));s.rSquaredAbsU2XOverYInf=max(second(mask));
s.wallRatiosUseMaintainedOneSidedDerivativeLimit=true;
s.finiteNativeVelocityData=all(isfinite([v.u1(mask);v.u2(mask);gradient(mask);first(mask);second(mask)]));
end
function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
