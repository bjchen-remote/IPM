function ipm_accellab_render_one_sided_profile_index(indicesFile,outputFile)
% Presentation-only rendering of already-complete immutable observations.
assert(~isfile(outputFile));d=load(indicesFile,'report');report=d.report;
s=report.samples;colors=lines(numel(s));fig=figure('Visible','off','Color','w','Position',[50,50,1300,850]);
clean=onCleanup(@()close(fig));tiledlayout(2,2,'TileSpacing','compact');labels=compose('tau %.3f',[s.parentEquivalentTau]);
nexttile;hold on
for k=1:numel(s),semilogx(s(k).physicalRadius,s(k).effectiveSpatialIndex,'Color',colors(k,:),'LineWidth',1.4);end
set(gca,'XScale','log');xlabel('fixed physical radius r');ylabel('a_{eff}(r)');title('Admitted positive-side finite-time indices');grid on;legend(labels,'Location','best');
nexttile;hold on
for k=1:numel(s),v=s(k).physicalDensityIncrement;v(~s(k).valid)=NaN;loglog(s(k).physicalRadius,v,'Color',colors(k,:),'LineWidth',1.4);end
set(gca,'XScale','log','YScale','log');xlabel('fixed physical radius r');ylabel('rho(x_{peak}+r)-rho(x_{peak})');title('Same-polynomial physical density increments');grid on
nexttile;hold on
for k=1:numel(s),semilogx(s(k).physicalRadius,double(s(k).valid)+.12*(k-1),'Color',colors(k,:),'LineWidth',1.3);end
set(gca,'XScale','log');xlabel('fixed physical radius r');ylabel('valid (traces offset by 0.12)');title('All registered-radius inclusion masks');grid on
nexttile;t=report.temporalScaleSlopes;mid=.5*([t.startTau]+[t.endTau]);
plot(mid,[t.amplitudeWidthLogSlopeContinuousPeak],'-o',mid,[t.amplitudeWidthLogSlopeNativeGridMaximum],'--s','LineWidth',1.3);
xlabel('adjacent interval midpoint in parent-equivalent tau');ylabel('Delta log(G w) / Delta log(w)');
title('Separate temporal scale slopes; remesh jumps retained');grid on;legend('continuous peak','native grid maximum','Location','best');
sgtitle('Finite-time profile observations; no limiting regularity or singularity claim','Color',[.1,.1,.1]);
set(findall(fig,'Type','axes'),'Color','w','XColor',[.15,.15,.15],'YColor',[.15,.15,.15]);
set(findall(fig,'Type','text'),'Color',[.1,.1,.1]);
set(findall(fig,'Type','legend'),'Color','w','TextColor',[.1,.1,.1],'EdgeColor',[.6,.6,.6]);
exportgraphics(fig,outputFile,'Resolution',160);
end
