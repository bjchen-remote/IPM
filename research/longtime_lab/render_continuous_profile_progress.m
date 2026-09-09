function report=render_continuous_profile_progress(checkpointFiles,outputDirectory)
%RENDER_CONTINUOUS_PROFILE_PROGRESS Actual amplitudes and paired native axes.
% Fresh clocks are mapped only through their validated direct native parent.
assert(iscell(checkpointFiles) && numel(checkpointFiles)>=3 && ~isfolder(outputDirectory));
assert(maxNumCompThreads==10);mkdir(outputDirectory);
n=numel(checkpointFiles);curves=cell(1,n);records=cell(1,n);
xi=linspace(-2,2,401);eta=linspace(0,3,401);
for k=1:n
    cp=ipm.output.readCheckpoint(checkpointFiles{k});s=cp.payload.state;h=cp.payload.log.history;
    assert(all(ipm.output.trustedMask(h,s.config)) && strcmp(s.config.scaling.dynamicScaleGeometry,'isotropic'));
    factor=1;cx0=1;epoch=0;tau0=0;fresh=false;
    if isfield(s.runMetadata.caseMetadata,'latePhysicalBoxBranch')
        l=s.runMetadata.caseMetadata.latePhysicalBoxBranch;parent=ipm.output.readCheckpoint(l.parentCheckpoint);p=parent.payload.state;
        assert(~isfield(p.runMetadata.caseMetadata,'latePhysicalBoxBranch') && ...
            strcmp(p.runMetadata.caseId,l.parentCaseId) && p.step==l.parentStep && ...
            p.scale.physicalTime==l.absolutePhysicalEpoch && p.scale.canonicalTime==l.parentCanonicalTime && ...
            exp(p.scale.logC_l)==l.parentCx && exp(p.scale.logC_omega)==l.parentComega && ...
            l.canonicalCovarianceFactor==l.parentCx/l.parentComega && ...
            ~l.parentHistoryInherited && ~l.parentCheckpointModified);
        factor=l.canonicalCovarianceFactor;cx0=l.parentCx;epoch=l.absolutePhysicalEpoch;tau0=l.parentCanonicalTime;fresh=true;
    end
    Dx=ipm.mesh.fdMatrix(s.x,1,7);Dy=ipm.mesh.fdMatrix(s.y,1,7);W=s.rho*Dx';
    g=continuous_inner_geometry(W,s.x,s.y,Dx,Dy,s.config.scaling.transportAnchorX*[.8,1.15]);
    center=ipm_accellab_tensor_hermite(s.rho,s.x,s.y,Dx,Dy,g.peakX,0);
    wall=ipm_accellab_tensor_hermite(W,s.x,s.y,Dx,Dy,g.peakX+g.wallWidth*xi,zeros(size(xi)));
    vertical=ipm_accellab_tensor_hermite(W,s.x,s.y,Dx,Dy,g.peakX+zeros(size(eta)),g.verticalWidth*eta);
    Cx=exp(s.scale.logC_l);Cw=exp(s.scale.logC_omega);
    curves{k}=struct('physicalCenteredX',(s.x-g.peakX)/Cx, ...
        'physicalDensity',(s.rho(1,:)-center.value)/Cw,'physicalGradient',W(1,:)*Cx/Cw, ...
        'parentCoordinateX',cx0*s.x,'normalizedWall',wall.value/g.peak, ...
        'normalizedVertical',vertical.value/g.peak,'parentNormalizedGradient',W(1,:)/g.peak);
    records{k}=struct('checkpointFile',checkpointFiles{k},'caseId',char(s.runMetadata.caseId), ...
        'freshCase',fresh,'step',s.step,'tau',tau0+factor*s.scale.canonicalTime, ...
        'absolutePhysicalTime',epoch+s.scale.physicalTime,'physicalPeak',g.peak*Cx/Cw, ...
        'physicalGridGradientMaximum',h.common.physicalRhoXInf(end), ...
        'physicalWidths',[g.wallWidth,g.verticalWidth]/Cx,'physicalPeakX',g.peakX/Cx, ...
        'cLParentUnits',h.common.c_l(end)/factor,'cOmegaParentUnits',h.common.c_omega(end)/factor, ...
        'coreCells',[h.mesh.coreGridPoints(end),h.mesh.verticalCoreGridPoints(end)], ...
        'safety',h.mesh.safetyFactor(end),'geometry',g,'fullHistoryTrusted',true);
end
records=[records{:}];tau=[records.tau];assert(all(diff(tau)>0));
report=struct('kind','actual_native_continuous_profile_progress','records',records, ...
    'nativeSignatureValidated',true,'geometryUniversalFourthOrderClaim',false,'pdeAdvanced',false, ...
    'interpretation','Actual native finite-time profiles, paired with their own axes. Center and additive density level are aligned; physical density and gradient amplitudes are retained. Inner normalization is a shape diagnostic, not an evolved extrapolated state or a singularity proof.');
colors=turbo(n);labels=compose('tau %.3f',tau);
fig=figure('Visible','off','Color','w','Position',[60,60,1500,920]);clean=onCleanup(@()close(fig));
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;hold on
for k=1:n,plot(curves{k}.physicalCenteredX,curves{k}.physicalDensity,'Color',colors(k,:),'LineWidth',1.5);end
xlim([-.008,.008]);xlabel('physical x - current peak x');ylabel('rho - rho at peak');
title('Physical wall density; amplitude retained');grid on;legend(labels,'Location','best');
nexttile;hold on
for k=1:n,plot(curves{k}.physicalCenteredX,curves{k}.physicalGradient,'Color',colors(k,:),'LineWidth',1.5);end
xlim([-.004,.004]);xlabel('physical x - current peak x');ylabel('physical rho_x');
title('Resolved gradient sharpening');grid on
nexttile;hold on
for k=1:n,plot(xi,curves{k}.normalizedWall,'Color',colors(k,:),'LineWidth',1.5);end
xlabel('(X - continuous peak X) / wall width');ylabel('R_X / continuous peak');
title('Wall shape at measured inner scale');grid on
nexttile;hold on
for k=1:n,plot(eta,curves{k}.normalizedVertical,'Color',colors(k,:),'LineWidth',1.5);end
xlabel('Y / continuous vertical width');ylabel('R_X at continuous peak / peak');
title('Vertical shape at the same peak');grid on
nexttile;widths=vertcat(records.physicalWidths);semilogy(tau,widths(:,1),'-o',tau,widths(:,2),'-s','LineWidth',1.5);
xlabel('parent-equivalent tau');ylabel('physical 0.9 core width');title('Measured core contraction');grid on;legend('wall','vertical','Location','best');
nexttile;plot(tau,[records.cLParentUnits],'-o',tau,[records.cOmegaParentUnits],'-s','LineWidth',1.5);
xlabel('parent-equivalent tau');ylabel('original gauge rate in parent units');title('Rates retain their measured drift');grid on;legend('c_l','c_omega','Location','best');
sgtitle(sprintf('IPM actual profiles: tau %.4f | physical t %.9f | gradient %.4f | core %.1f x %.1f', ...
    tau(end),records(end).absolutePhysicalTime,records(end).physicalGridGradientMaximum,records(end).coreCells));
exportgraphics(fig,fullfile(outputDirectory,'profile_evidence.png'),'Resolution',160);
save(fullfile(outputDirectory,'profiles.mat'),'report','curves','xi','eta','-v7.3');
fid=fopen(fullfile(outputDirectory,'profile_report.json'),'w');assert(fid>=0);fcloseClean=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
