function [figureHandle,analysis] = ipm_plot_lattice_spacing_balance(outputFile)
%IPM_PLOT_LATTICE_SPACING_BALANCE Compare q256 peak-lattice transitions.
%   The retained A=28 Lorentzian LAT is compared with the rounded-log F80
%   R16--R28 candidates and the denser F90/F92 R20 candidates.  Diagnostics
%   use the positive half-axis;
%   density jumps are |Delta log(dx)| and therefore do not depend on whether
%   they are described as cell-width or cell-density ratios.

analysisDirectory = fileparts(mfilename('fullpath'));
researchDirectory = fileparts(analysisDirectory);
projectDirectory = fileparts(researchDirectory);
experimentDirectory = fullfile(researchDirectory,'experiments');
addpath(projectDirectory,experimentDirectory,analysisDirectory);
if nargin < 1 || isempty(outputFile)
    outputFile = ipm.output.uniquePath(fullfile(projectDirectory, ...
        'result','verification', ...
        'fourth_order_q256_coordinated_lattice_v5_F80_F90_F92.png'), ...
        'spacing_balance');
end
validateattributes(outputFile,{'char','string'},{'scalartext'}, ...
    mfilename,'outputFile');
outputFile = char(outputFile);

caseNames = { ...
    'analytic_lorentzian_A28_w020_c114_quadrant_256', ...
    'coordinated_log_F80_R16_c114_quadrant_256', ...
    'coordinated_log_F80_R20_c114_quadrant_256', ...
    'coordinated_log_F80_R24_c114_quadrant_256', ...
    'coordinated_log_F80_R28_c114_quadrant_256', ...
    'coordinated_log_F90_R20_c114_quadrant_256', ...
    'coordinated_log_F92_R20_c114_quadrant_256'};
labels = {'A28 Lorentzian','coordinated R16','coordinated R20', ...
    'coordinated R24','coordinated R28','coordinated F90/R20', ...
    'coordinated F92/R20'};
colors = lines(numel(caseNames));
lineStyles = [{'--'},repmat({'-'},1,numel(caseNames)-1)];

template = struct('caseName','','label','','cellCenters',[], ...
    'gridInfo',struct(), ...
    'cellWidths',[],'interfaceX',[],'densityJumps',[], ...
    'activeDensityJumps',[],'activeInterfaceCount',NaN, ...
    'platformZeroInterfaceCount',NaN, ...
    'maximumAdjacentRatio',NaN,'maximumLogSpacingCurvature',NaN, ...
    'globalEqualizationEfficiency',NaN, ...
    'globalHalfVariationInterfaceFraction',NaN, ...
    'transitionEqualizationEfficiency',NaN, ...
    'transitionHalfVariationInterfaceFraction',NaN, ...
    'asymptoticSlopeBalance',NaN,'peakPathEffectiveCells',NaN, ...
    'peakPathMinimumSpacing',NaN,'peakPathMaximumSpacing',NaN, ...
    'fineCoverage',nan(1,2),'originSpacing',NaN, ...
    'remoteSpacing',NaN,'outerPeakToFarJumpRatio',NaN);
datasets = repmat(template,numel(caseNames),1);

for index = 1:numel(caseNames)
    opts = fourth_order_blowup_case(caseNames{index});
    gridInfo = opts.caseMetadata.smoothGrid;
    positive = opts.customX((opts.nx+1)/2:end);
    widths = diff(positive);
    centers = positive(1:end-1)+0.5*widths;
    interfaceX = positive(2:end-1);
    jumps = abs(diff(log(widths)));
    balance = jump_statistics(jumps);
    jumpTolerance = 1024*eps(max(1,max(abs(log(widths)))));
    activeJumps = jumps(jumps > jumpTolerance);
    path = linspace(0.97,1.28,1001);
    pathWidths = interp1(centers,widths,path,'pchip');
    fine = find(widths <= 1.5*min(widths));
    far = centers(2:end) >= 20;
    farJump = median(jumps(far));
    outer = interfaceX >= 1.30;

    data = template;
    data.caseName = caseNames{index};
    data.label = labels{index};
    data.gridInfo = gridInfo;
    data.cellCenters = centers;
    data.cellWidths = widths;
    data.interfaceX = interfaceX;
    data.densityJumps = jumps;
    data.activeDensityJumps = activeJumps;
    data.activeInterfaceCount = numel(activeJumps);
    data.platformZeroInterfaceCount = numel(jumps)-numel(activeJumps);
    data.maximumAdjacentRatio = gridInfo.maximumAdjacentCellRatioX;
    data.maximumLogSpacingCurvature = ...
        gridInfo.maximumLogSpacingCurvatureX;
    data.globalEqualizationEfficiency = balance.efficiency;
    data.globalHalfVariationInterfaceFraction = balance.halfFraction;
    if ~isempty(fieldnames(gridInfo.coordinatedDesign))
        transition = gridInfo.coordinatedDesign.transitionLogSpacing;
        data.transitionEqualizationEfficiency = ...
            transition.equalizationEfficiency;
        data.transitionHalfVariationInterfaceFraction = ...
            transition.halfVariationInterfaceFraction;
        data.asymptoticSlopeBalance = ...
            gridInfo.coordinatedDesign.asymptoticSlopeBalance;
    end
    data.peakPathEffectiveCells = trapz(path,1./pathWidths);
    data.peakPathMinimumSpacing = min(pathWidths);
    data.peakPathMaximumSpacing = max(pathWidths);
    data.fineCoverage = [positive(fine(1)),positive(fine(end)+1)];
    data.originSpacing = widths(1);
    data.remoteSpacing = widths(end);
    data.outerPeakToFarJumpRatio = max(jumps(outer))/farJump;
    datasets(index) = data;
end

figureHandle = figure('Color','w','Visible','off', ...
    'Position',[60,60,1680,1080],'Name','q256 lattice spacing balance');
layout = tiledlayout(figureHandle,2,2,'TileSpacing','compact', ...
    'Padding','compact');

axisHandle = nexttile(layout,1);
hold(axisHandle,'on');
for index = 1:numel(datasets)
    semilogx(axisHandle,datasets(index).cellCenters+1e-3, ...
        datasets(index).cellWidths,'Color',colors(index,:), ...
        'LineStyle',lineStyles{index},'LineWidth',1.7);
end
set(axisHandle,'XScale','log','YScale','log');
grid(axisHandle,'on');
xlabel(axisHandle,'positive-half cell center x (+10^{-3})');
ylabel(axisHandle,'cell width \Deltax');
title(axisHandle,'Positive-half spacing distribution');
legend(axisHandle,labels,'Location','northwest','Interpreter','none');

axisHandle = nexttile(layout,2);
hold(axisHandle,'on');
for index = 1:numel(datasets)
    semilogx(axisHandle,datasets(index).interfaceX+1e-3, ...
        exp(datasets(index).densityJumps),'Color',colors(index,:), ...
        'LineStyle',lineStyles{index},'LineWidth',1.7);
end
grid(axisHandle,'on');
xlabel(axisHandle,'cell interface x (+10^{-3})');
ylabel(axisHandle,'adjacent density ratio exp(|\Delta log \Deltax|)');
title(axisHandle,'Local density-ratio load');
set(axisHandle,'XScale','log');

axisHandle = nexttile(layout,3);
hold(axisHandle,'on');
plot(axisHandle,[0,1],[0,1],':','Color',[0.35,0.35,0.35], ...
    'LineWidth',1.2);
for index = 1:numel(datasets)
    sorted = sort(datasets(index).activeDensityJumps,'descend');
    fraction = (1:numel(sorted))/numel(sorted);
    plot(axisHandle,fraction,cumsum(sorted)/sum(sorted), ...
        'Color',colors(index,:),'LineStyle',lineStyles{index}, ...
        'LineWidth',1.7);
end
grid(axisHandle,'on');
xlabel(axisHandle,'fraction of active interfaces (largest jumps first)');
ylabel(axisHandle,'fraction of total |\Delta log \Deltax|');
title(axisHandle,{ ...
    'Active-transition spacing variation', ...
    'ratio=1 platform interfaces omitted'});
selected = datasets(end);
text(axisHandle,0.47,0.17,sprintf([ ...
    'F92/R20: active E=%.3f, global E=%.3f\n' ...
    '%d active + %d platform-zero interfaces'], ...
    selected.transitionEqualizationEfficiency, ...
    selected.globalEqualizationEfficiency, ...
    selected.activeInterfaceCount,selected.platformZeroInterfaceCount), ...
    'FontSize',9,'Color',[0.1,0.1,0.1], ...
    'BackgroundColor','w','Margin',4, ...
    'EdgeColor',[0.75,0.75,0.75]);

axisHandle = nexttile(layout,4);
hold(axisHandle,'on');
for index = 1:numel(datasets)
    plot(axisHandle,datasets(index).cellCenters, ...
        datasets(index).cellWidths,'Color',colors(index,:), ...
        'LineStyle',lineStyles{index},'LineWidth',1.7);
end
xline(axisHandle,0.95,':','Color',[0.25,0.25,0.25]);
xline(axisHandle,1.30,':','Color',[0.25,0.25,0.25]);
xline(axisHandle,0.97,'--','Color',[0.45,0.45,0.45]);
xline(axisHandle,1.28,'--','Color',[0.45,0.45,0.45]);
xlim(axisHandle,[0.82,1.43]);
grid(axisHandle,'on');
xlabel(axisHandle,'x');
ylabel(axisHandle,'cell width \Deltax');
title(axisHandle,'Solver/grid-x peak path and protected fine platform');

sgtitle(layout,{ ...
    'q256 coordinated LAT: adjacent density-ratio equalization', ...
    'fine platform x=0.95--1.30; measured peak path x=0.97--1.28'}, ...
    'FontWeight','bold');

outputDirectory = fileparts(outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle,outputFile,'Resolution',300);
[outputDirectory,outputStem] = fileparts(outputFile);
analysisFile = fullfile(outputDirectory,[outputStem,'_analysis.mat']);
analysis = struct('definition',struct( ...
    'densityJump','abs(diff(log(cellWidth)))', ...
    'activeInterface', ...
        'densityJump > 1024*eps(max(1,max(abs(log(cellWidth)))))', ...
    'equalizationEfficiency','sum(a)^2/(numel(a)*sum(a.^2))', ...
    'peakPath','[0.97,1.28]'), ...
    'datasets',datasets,'outputFile',outputFile, ...
    'analysisFile',analysisFile);
save(analysisFile,'analysis');
end

function statistics = jump_statistics(jumps)
jumps = jumps(:);
total = sum(jumps);
statistics = struct('efficiency',total^2/(numel(jumps)*sum(jumps.^2)), ...
    'halfFraction',find(cumsum(sort(jumps,'descend')) >= 0.5*total,1)/ ...
    numel(jumps));
end
