function plot_mesh_density_distribution(distributionFile,xSweepFile,ySweepFile,outputFile)
%PLOT_MESH_DENSITY_DISTRIBUTION Scientific fixed-node interface-ratio figure.
assert(isfile(distributionFile)&&isfile(xSweepFile)&&isfile(ySweepFile)&&~isfile(outputFile));
d=jsondecode(fileread(distributionFile));
xs=jsondecode(fileread(xSweepFile));ys=jsondecode(fileread(ySweepFile));
oldX=d.x.nodes(:);newX=xs.best.axis(:);
oldY=d.y.nodes(:);newY=ys.best.axis(:);
fig=figure('Visible','off','Color','w','Position',[100,100,1120,650]);
cleaner=onCleanup(@()close(fig));
subplot(2,1,1);
plot_ratio(oldX,newX,'X');
subplot(2,1,2);
plot_ratio(oldY,newY,'Y');
sgtitle(sprintf('Fixed-node mesh-density comparison at canonical \\tau = %.4f',d.tau));
exportgraphics(fig,outputFile,'Resolution',180);
end
function plot_ratio(oldAxis,newAxis,label)
a=ratio(oldAxis);b=ratio(newAxis);
plot(1:numel(a),a,'Color',[.75,.22,.17],'LineWidth',1.35);hold on;
plot(1:numel(b),b,'Color',[.12,.37,.70],'LineWidth',1.35);
grid on;xlim([1,numel(a)]);
ylabel(sprintf('%s adjacent cell ratio',label));
xlabel('interface index on positive axis');
legend(sprintf('existing: max %.5f',max(a)), ...
    sprintf('balanced density: max %.5f',max(b)),'Location','best');
end
function r=ratio(axis)
h=diff(axis);r=exp(abs(diff(log(h))));
end
