function [figureHandle,outputFile] = ...
    ipm_gaugelab_plot_tournament(products,outputFile)
%IPM_GAUGELAB_PLOT_TOURNAMENT Compare trusted q256 gauge histories.

if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(products.outputDirectory, ...
        'q256_gauge_tournament_rates.png');
end
if isfile(outputFile)
    error('ipm:GaugeLabPlotCollision', ...
        'Refusing to overwrite an existing tournament plot: %s',outputFile);
end
plottable = find(arrayfun(@(run)~isempty(run.resultFile),products.runs));
curveTemplate = struct('t',[],'cl',[],'cw',[],'rhoX',[], ...
    'label','','eligible',false,'negativeControl',false);
curves = repmat(curveTemplate,numel(plottable),1);
curveCount = 0;
for j = 1:numel(plottable)
    run = products.runs(plottable(j));
    result = ipm.output.validate(run.resultFile);
    trusted = ipm.output.continuousTrustedPrefix( ...
        ipm.output.trustedMask(result));
    if ~any(trusted)
        continue;
    end
    h = result.history.common;
    if ~isempty(run.errorIdentifier)
        status = 'error';
    elseif strcmp(run.role,'negative_control')
        status = 'negative control';
    elseif isfield(run,'productionEligible') && run.productionEligible
        status = 'production eligible';
    else
        status = 'gate fail';
    end
    curve = struct('t',h.physicalTime(trusted), ...
        'cl',h.canonicalCL(trusted),'cw',h.canonicalCOmega(trusted), ...
        'rhoX',h.physicalRhoXInf(trusted), ...
        'label',sprintf('%s [%s]', ...
        strrep(run.cOmegaGauge,'_','\_'),status), ...
        'eligible',isfield(run,'productionEligible') && ...
        run.productionEligible,'negativeControl', ...
        strcmp(run.role,'negative_control'));
    curveCount = curveCount+1;
    curves(curveCount) = curve;
end
curves = curves(1:curveCount);
if isempty(curves)
    error('ipm:GaugeLabPlotEmpty','No successful tournament result exists.');
end
eligibleCurves = [curves.eligible];
if any(eligibleCurves)
    commonPhysicalEnd = min(arrayfun( ...
        @(curve)curve.t(end),curves(eligibleCurves)));
else
    commonPhysicalEnd = min(arrayfun(@(curve)curve.t(end),curves));
end

figureHandle = figure('Color','white','Position',[100,100,1200,820]);
try
    layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    colors = lines(numel(curves));
    labels = cell(numel(curves),1);
    for j = 1:numel(curves)
        curve = curves(j);
        keep = curve.t <= commonPhysicalEnd+ ...
            64*eps(max(1,abs(commonPhysicalEnd)));
        if curve.negativeControl
            lineStyle = '-.';
        elseif curve.eligible
            lineStyle = '-';
        else
            lineStyle = '--';
        end
        labels{j} = curve.label;

        nexttile(layout,1);
        plot(curve.t(keep),curve.cl(keep),'LineWidth',1.35, ...
            'LineStyle',lineStyle,'Color',colors(j,:)); hold on;
        nexttile(layout,2);
        plot(curve.t(keep),curve.cw(keep),'LineWidth',1.35, ...
            'LineStyle',lineStyle,'Color',colors(j,:)); hold on;
        nexttile(layout,3);
        plot(curve.t(keep),curve.cl(keep)-curve.cw(keep), ...
            'LineWidth',1.35,'LineStyle',lineStyle, ...
            'Color',colors(j,:)); hold on;
        nexttile(layout,4);
        semilogy(curve.t(keep),curve.rhoX(keep),'LineWidth',1.35, ...
            'LineStyle',lineStyle,'Color',colors(j,:)); hold on;
    end

    titles = {'canonical c_l','canonical c_\omega', ...
        '\kappa=c_l-c_\omega','physical max |\rho_{x_1}|'};
    ylabels = {'c_l','c_\omega','\kappa','maximum'};
    for tile = 1:4
        axisHandle = nexttile(layout,tile);
        grid(axisHandle,'on');
        xlabel(axisHandle,'physical time');
        ylabel(axisHandle,ylabels{tile});
        title(axisHandle,titles{tile});
    end
    legend(nexttile(layout,1),labels,'Location','best','Interpreter','tex');
    title(layout,sprintf( ...
        ['q256 exact-no-feedback amplitude-gauge tournament: ' ...
        'common trusted t <= %.6g'], ...
        commonPhysicalEnd));
    drawnow;
    exportgraphics(figureHandle,outputFile,'Resolution',240, ...
        'BackgroundColor','white');
catch exception
    if isgraphics(figureHandle)
        close(figureHandle);
    end
    rethrow(exception);
end
end
