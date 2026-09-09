function report=mesh_plot_native_progress(checkpointFile,outDir)
%MESH_PLOT_NATIVE_PROGRESS Plot the actual signed from-zero controller history.
assert(maxNumCompThreads==10 && ~isfolder(outDir));mkdir(outDir);
cp=ipm.output.readCheckpoint(checkpointFile);s=cp.payload.state;h=cp.payload.log.history;
assert(s.config.remesh.autonomousMesh.version==2 && all(ipm.output.trustedMask(h,s.config)));
c=h.common;m=h.mesh;memory=s.runMetadata.autonomousMesh;
assert(c.acceptedStep(1)==0 && c.canonicalTau(1)==0 && c.physicalTime(1)==0 && ...
    memory.initialization.originalPhysicalTime==0 && ~memory.initialization.newPhysicalEpochCreated);
tau=c.canonicalTau;transactions=memory.transactions;
growth=transactions(arrayfun(@(v)any(v.targetNodeCount~=v.sourceNodeCount),transactions));
growthTau=arrayfun(@(v)v.controllerDecision.sourceCanonicalTime,growth);
report=struct('kind','actual_from_zero_signed_history_progress_v1', ...
    'checkpointFile',checkpointFile,'caseId',s.runMetadata.caseId,'step',s.step, ...
    'canonicalTime',s.scale.canonicalTime,'physicalTime',s.scale.physicalTime, ...
    'currentNodeCount',[numel(s.x),numel(s.y)],'actualRemeshCount',numel(transactions), ...
    'actualGrowthCount',numel(growth),'growthCanonicalTimes',growthTau, ...
    'historyRecords',numel(tau),'entireNativeHistoryTrusted',true, ...
    'physicalRhoXMaximum',c.physicalRhoXInf(end), ...
    'physicalGradientInfinity',c.physicalGradInf(end), ...
    'coreCells',[m.coreGridPoints(end),m.verticalCoreGridPoints(end)], ...
    'noLU',true,'noPDE',true,'asymptoticOrErrorQualification',false);
fig=figure('Visible','off','Color','w','Position',[50,50,1250,900]);
layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
axesList=gobjects(1,6);
axesList(1)=nexttile;stairs(tau,c.nodeCountX,'LineWidth',1.6);hold on;
stairs(tau,c.nodeCountY,'LineWidth',1.6);ylabel('Nodes');title('Automatic directional growth');
legend({'N_x','N_y'},'Location','northwest');
axesList(2)=nexttile;plot(tau,[m.coreGridPoints,m.verticalCoreGridPoints],'LineWidth',1.3);hold on;
yline(s.config.remesh.autonomousMesh.regridCoreTrigger(1),'--','Current-core trigger', ...
    'LabelHorizontalAlignment','left','Color',[.25,.25,.25]);
ylim([min([m.coreGridPoints;m.verticalCoreGridPoints;s.config.remesh.autonomousMesh.regridCoreTrigger(1)])-1, ...
    max([m.coreGridPoints;m.verticalCoreGridPoints])+1]);
ylabel('Cells across 90% core');title('Actual paired-grid resolution');legend({'x core','y core'},'Location','northwest');
axesList(3)=nexttile;semilogy(tau,c.physicalRhoXInf,'LineWidth',1.6);ylabel('Physical max |rho_x|');title('Gradient growth');
axesList(4)=nexttile;semilogy(tau,[m.trackedWallCoreWidth./c.C_x,m.trackedVerticalCoreWidth./c.C_y],'LineWidth',1.4);
ylabel('Physical core width');title('Physical concentration');legend({'wall','vertical'},'Location','northeast');
axesList(5)=nexttile;plot(tau,[c.canonicalCL,c.canonicalCOmega],'LineWidth',1.3);
ylabel('Canonical scale rate');title('Unmodified instantaneous gauge');legend({'c_l','c_{omega}'},'Location','best');
axesList(6)=nexttile;plot(tau,c.physicalTime,'LineWidth',1.6);ylabel('Physical time');title('Physical clock, with no epoch reset');
for k=1:numel(axesList)
    ax=axesList(k);set(ax,'Color','w','XColor','k','YColor','k','FontSize',11,'LineWidth',.8);
    ax.Title.Color='k';ax.XLabel.Color='k';ax.YLabel.Color='k';grid(ax,'on');xlabel(ax,'Canonical time tau');
    xlim(ax,[0,tau(end)]);
    for time=growthTau,xline(ax,time,':','Color',[.7,.2,.15],'HandleVisibility','off');end
end
set(findall(fig,'Type','legend'),'Color','w','TextColor','k','EdgeColor',[.7,.7,.7]);
title(layout,sprintf('Uninterrupted from physical t = 0 | tau = %.4f | %d x %d nodes', ...
    s.scale.canonicalTime,numel(s.x),numel(s.y)),'Color','k','FontSize',16);
subtitle(layout,'Signed native history; dotted lines: actual node growth. Finite-box run, not a singularity certificate.', ...
    'Color',[.25,.25,.25],'FontSize',11);
exportgraphics(fig,fullfile(outDir,'native_progress.png'),'Resolution',160,'BackgroundColor','white');
savefig(fig,fullfile(outDir,'native_progress.fig'));close(fig);
save(fullfile(outDir,'native_progress.mat'),'h','report','growth','-v7.3');
fid=fopen(fullfile(outDir,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('NATIVE_PROGRESS_PLOT_PASS step=%d tau=%.12g\n',s.step,s.scale.canonicalTime);
end
