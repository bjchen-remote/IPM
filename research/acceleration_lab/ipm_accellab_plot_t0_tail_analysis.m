function ipm_accellab_plot_t0_tail_analysis(directory)
d=load(fullfile(directory,'run/analysis.mat'),'report','fieldData');r=d.report;H=r.registration.H;
cl=zeros(2,2);error=zeros(2,3);for k=1:2
 a=r.cases{k};cl(k,:)=[a.baselineRates(1),a.withTailRates(1)];
 error(k,:)=[abs(a.nativeBaselineMinusExactInsideCL),abs(a.discreteHarmonicAnchorResponseMinusExactTailCL),abs(a.correctedNativeMinusExactInfiniteCL)];
end
fig=figure('Visible','off','Color','w','Position',[60,60,1400,1000]);c=onCleanup(@()close(fig));
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;bar(H,cl);hold on;yline(r.cases{1}.infiniteExactInitialCL,'k--','LineWidth',1.3);
xlabel('H (each case retains its own initial axes)');ylabel('c_l at physical t = 0');title('Only the exterior Dirichlet trace changes');legend('Native interior source','With exact initial exterior tail','Exact initial full-space Green integral','Location','southwest');grid on;
nexttile;hold on;colors=lines(2);
for k=1:2
 f=d.fieldData{k};plot(f.wallX,f.wallU1Baseline,'--','Color',colors(k,:),'LineWidth',1.4);
 plot(f.wallX,f.wallU1WithTail,'-','Color',colors(k,:),'LineWidth',1.6);
end
xlabel('physical x on wall');ylabel('u_1');title('Actual fixed-operator wall velocities');legend('H8 native','H8 with tail','H32 native','H32 with tail','Location','best');grid on;
nexttile;hold on;
for branch={'baseline','withTail'}
 name=branch{1};a=d.fieldData{1}.(name);b=d.fieldData{2}.(name);v=hypot(b.u1-a.u1,b.u2-a.u2);v=reshape(v,11,21);
 plot(linspace(0,4,21),max(v,[],1),'-o','LineWidth',1.5);
end
set(gca,'YScale','log');xlabel('x (maximum over fixed y samples in [0,2])');ylabel('|u_{H32} - u_{H8}|');title('Cross-case difference remains after tail addition');legend('Native','With exact initial tail','Location','best');grid on;
nexttile;semilogy(H,error(:,1),'-o','LineWidth',2);hold on;
semilogy(H,error(:,2),'-s','LineWidth',1.5);semilogy(H,error(:,3),'--x','LineWidth',1.3);
xlabel('H');ylabel('absolute c_l discrepancy');title('The remaining error is the original interior-method gap');
legend('Native minus exact inside integral','Discrete tail response minus exact tail','With tail minus exact full integral','Location','best');grid on;
sgtitle('Original t0 k8 datum | one fixed LU per case, two Poisson solves, zero PDE steps | no evolved-tail closure');
ax=findall(fig,'Type','axes');set(ax,'Color','w','XColor',[.15,.15,.15],'YColor',[.15,.15,.15]);
set(findall(fig,'Type','text'),'Color',[.1,.1,.1]);set(findall(fig,'Type','legend'),'Color','w','TextColor',[.1,.1,.1]);
exportgraphics(fig,fullfile(directory,'t0_tail_operator_response.png'),'Resolution',150);
exportgraphics(fig,fullfile(directory,'t0_tail_operator_response.pdf'),'ContentType','vector');
end
