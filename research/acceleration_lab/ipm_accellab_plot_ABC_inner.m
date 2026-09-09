function ipm_accellab_plot_ABC_inner(out)
%IPM_ACCELLAB_PLOT_ABC_INNER Plot saved profiles without field recomputation.
d=load(fullfile(out,'profiles.mat'),'report','curves');r=d.report;c=d.curves;
f=figure('Visible','off','Color','w','Position',[0,0,1200,800]);cleanup=onCleanup(@()close(f));
layout=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
colors=[.13,.4,.67;.7,.094,.17;.137,.545,.27];
for direction=1:2
    ax=nexttile(layout,direction);hold(ax,'on');handles=gobjects(1,3);labels=cell(1,3);
    for k=1:3
        z=c{k,direction};handles(k)=plot(ax,z.coordinate,z.hermite,'Color',colors(k,:),'LineWidth',1.8);
        plot(ax,z.coordinate,z.linear,'Color',colors(k,:),'LineStyle','--','LineWidth',.6);
        labels{k}=sprintf('%s, tau=%.4f',r.cases(k).label,r.cases(k).canonicalTime);
    end
    ylabel(ax,'rho_x / continuous peak');legend(ax,handles,labels,'Interpreter','none','Location','best','FontSize',9);
    if direction==1,title(ax,'Wall: x = a + w_x xi, y = 0');xlabel(ax,'xi');
    else,title(ax,'Vertical: x = a, y = w_y eta');xlabel(ax,'eta');end
    grid(ax,'on');box(ax,'off');
    ax=nexttile(layout,direction+2);hold(ax,'on');pairs=[1,2;1,3;2,3];
    for k=1:3
        a=c{pairs(k,1),direction};b=c{pairs(k,2),direction};assert(isequal(a.coordinate,b.coordinate));
        handles(k)=plot(ax,a.coordinate,b.hermite-a.hermite,'Color',colors(k,:),'LineWidth',1.5);
        plot(ax,a.coordinate,b.linear-a.linear,'Color',colors(k,:),'LineStyle','--','LineWidth',.6);
    end
    yline(ax,0,'Color',[.5,.5,.5],'LineWidth',.5);legend(ax,handles,{'B - A','C - A','C - B'},'Location','best');
    ylabel(ax,'Aligned shape difference');grid(ax,'on');box(ax,'off');
    if direction==1,xlabel(ax,'xi');else,xlabel(ax,'eta');end
end
title(layout,sprintf(['Actual common physical time t = %.12f\n', ...
    'Own C1 peak + connected 90%% widths; solid Hermite / dashed linear'],r.samePhysicalEndpoint),'FontSize',14);
exportgraphics(f,fullfile(out,'ABC_aligned_inner_shapes.png'),'Resolution',180);
exportgraphics(f,fullfile(out,'ABC_aligned_inner_shapes.pdf'),'ContentType','vector');
end
