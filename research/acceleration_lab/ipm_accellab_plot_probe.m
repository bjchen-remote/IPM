function ipm_accellab_plot_probe(experimentFile,outputFile)
%IPM_ACCELLAB_PLOT_PROBE Plot saved evidence without an operator/PDE load.
assert(~isfile(outputFile),'ipm:AccelPlotOverwrite','Plot output already exists.');
loaded = load(experimentFile,'experiment');
e = loaded.experiment;
s = e.sourceShapeDiagnostics;
figureHandle = figure('Visible','off','Color','w','Position',[50,50,1200,800]);
cleanup = onCleanup(@()close(figureHandle));
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout);
plot(ax,s.canonicalTime,s.wallCoreWidth/s.wallCoreWidth(1),'-o', ...
    s.canonicalTime,s.verticalCoreWidth/s.verticalCoreWidth(1),'-s','LineWidth',1.4);
xlabel(ax,'Canonical time'); ylabel(ax,'Width / first sampled width');
title(ax,'Both inner scales keep contracting');
legend(ax,{'Wall 90% core','Vertical 90% core'},'Location','southwest');
grid(ax,'on');
colors = parula(numel(s.canonicalTime));
labels = compose('%.3f',s.canonicalTime);
ax = nexttile(layout); hold(ax,'on');
for k = 1:numel(s.canonicalTime)
    plot(ax,s.wallInnerCoordinates,s.wallInnerProfiles(:,k),'Color',colors(k,:),'LineWidth',1.3);
end
xlabel(ax,'(X - peak X) / wall core width'); ylabel(ax,'R_X / P');
title(ax,'Wall shape in the shrinking inner coordinate');
legend(ax,labels,'Location','south','NumColumns',3); grid(ax,'on');
ax = nexttile(layout); hold(ax,'on');
for k = 1:numel(s.canonicalTime)
    plot(ax,s.verticalInnerCoordinates,s.verticalInnerProfiles(:,k),'Color',colors(k,:),'LineWidth',1.3);
end
xlabel(ax,'Y / vertical core width'); ylabel(ax,'R_X / P');
title(ax,'Vertical shape in the shrinking inner coordinate');
legend(ax,labels,'Location','southwest'); grid(ax,'on');
ax = nexttile(layout); hold(ax,'on');
markers = {'o','s','^'};
handles = gobjects(numel(e.trials),1);
trialLabels = cell(numel(e.trials),1);
for k = 1:numel(e.trials)
    trial = e.trials{k};
    xValues = []; yValues = [];
    for j = 1:numel(trial.attempts)
        attempt = trial.attempts(j);
        if isfield(attempt.audit,'coreGradientRhsRelativeInf')
            xValues(end+1) = attempt.residualNorm/trial.baselineNorm; %#ok<AGROW>
            yValues(end+1) = attempt.audit.coreGradientRhsRelativeInf/ ...
                trial.baselineAudit.coreGradientRhsRelativeInf; %#ok<AGROW>
        end
    end
    handles(k) = scatter(ax,xValues,yValues,55,markers{k},'filled');
    trialLabels{k} = sprintf('Memory %d',trial.options.memory);
end
xline(ax,1,':','Color',[0.4,0.4,0.4]); yline(ax,1,':','Color',[0.4,0.4,0.4]);
xlim(ax,[0.78,1.02]); ylim(ax,[0.92,1.92]);
xlabel(ax,'Combined L2 RHS norm / baseline');
ylabel(ax,'Core gradient RHS infinity norm / baseline');
title(ax,'Every evaluated secant candidate is rejected');
text(ax,0.79,1.85,'Mean residual improves; local maximum worsens','FontSize',9);
legend(ax,handles,trialLabels,'Location','southeast'); grid(ax,'on');
title(layout,'q512 exact-gauge probe: contracting inner shape and rejected extrapolation');
exportgraphics(layout,outputFile,'Resolution',180);
clear cleanup;
end
