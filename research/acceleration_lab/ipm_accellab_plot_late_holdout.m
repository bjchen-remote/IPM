function audit=ipm_accellab_plot_late_holdout(rootDirectory)
% Saved-array postprocessing only, no native CP read or refit.
d=load(fullfile(rootDirectory,'run/evaluate/observations.mat'),'report','curveData');r=d.report;
z=load(fullfile(rootDirectory,'run/train/frozen_model.mat'),'frozen');f=z.frozen;reg=r.registration;
methods={'hermite','linear'};maxDifference=0;
for observer=1:2
 for direction=1:2
  q=linspace(reg.intervals(direction,1),reg.intervals(direction,2),reg.pointCounts(direction));a=zeros(8,numel(q));
  for k=1:8,a(k,:)=sample(f.trainingFrames{k}.traces{observer,direction},q,methods{observer});end
  w=ones(size(q));w([1,end])=.5;w=w/sum(w);W=a.*sqrt(w);delta=diff(W)';
  X=delta(:,1:end-1);Y=delta(:,2:end);[U,S,V]=svd(X,'econ');B=Y*V(:,1)/S(1,1);c=U(:,1)'*delta(:,end);A=U(:,1)'*B;
  v=W(end,:)';m=f.shapeModels{observer,direction};
  for k=1:reg.futureCount
   v=v+B*c;c=A*c;direct=v'./sqrt(w);comb=m.coefficients{3}(k,:)*a;
   maxDifference=max(maxDifference,max(abs(direct-comb)));
  end
 end
end
assert(maxDifference<=1e-11);
audit=struct('kind','frozen_forecast_linear_combination_vs_direct_DMD_audit','maximumAbsoluteDifference',maxDifference, ...
 'allFourShapeModelsAllHorizonsPassed',true,'modelsModified',false,'futureRefit',false,'noLU',true,'noPDE',true);
write_json(fullfile(rootDirectory,'direct_model_algebra_audit.json'),audit);
n=numel(r.futureCases);tau=zeros(n,1);L=zeros(n,2,3);P=zeros(n,3);X=P;Y=P;position=P;obs=zeros(n,2);
for k=1:n
 a=r.futureCases{k};tau(k)=a.record.equivalentTau;
 for j=1:2
  for m=1:3,L(k,j,m)=100*a.shape{1,j}.modelErrors{m}.relativeL2;end
  obs(k,j)=100*a.observerDifferences{j}.relativeL2;
 end
 for m=1:3
  P(k,m)=100*a.geometry{1}.models{m}.relativeError;position(k,m)=a.geometry{2}.models{m}.errorOverActualWallWidth;
  X(k,m)=100*a.geometry{3}.models{m}.relativeError;Y(k,m)=100*a.geometry{4}.models{m}.relativeError;
 end
end
fig=figure('Visible','off','Color','w','Position',[60,60,1600,1150]);clean=onCleanup(@()close(fig));
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');names={'Persistence','Linear trend','Fixed rank 1'};
nexttile;semilogy(tau,squeeze(L(:,1,:)),'-o','LineWidth',1.4);ylabel('relative L2 error (%)');title('Wall inner shape');grid on;legend(names,'Location','northwest');
nexttile;semilogy(tau,squeeze(L(:,2,:)),'-o','LineWidth',1.4);ylabel('relative L2 error (%)');title('Vertical inner shape');grid on;
nexttile;semilogy(tau,obs,'-o','LineWidth',1.4);ylabel('relative L2 difference (%)');title('Hermite vs bilinear observation');grid on;legend('wall','vertical','Location','best');
nexttile;semilogy(tau,P,'-o','LineWidth',1.4);ylabel('relative error (%)');title('Physical peak amplitude');grid on;
nexttile;semilogy(tau,X,'-o','LineWidth',1.4);ylabel('relative error (%)');title('Physical wall width');grid on;
nexttile;semilogy(tau,Y,'-o','LineWidth',1.4);ylabel('relative error (%)');title('Physical vertical width');grid on;
nexttile;semilogy(tau,position,'-o','LineWidth',1.4);ylabel('|position error| / actual wall width');title('Phase error at the inner scale');grid on;
for j=1:2
 nexttile;a=d.curveData{end,1,j};hold on
 plot(a.coordinate,a.actual,'k','LineWidth',2);
 for m=1:3,plot(a.coordinate,a.predictions{m},'LineWidth',1.2);end
 ylabel('normalized R_X');if j==1,xlabel('xi');title('Final wall shape (actual geometry)');else,xlabel('eta');title('Final vertical shape (actual geometry)');end
 grid on;if j==2,legend([{'Actual'},names],'Location','southwest');end
end
axesList=findall(fig,'Type','axes');set(axesList,'Color','w','XColor',[.15,.15,.15],'YColor',[.15,.15,.15]);
set(findall(fig,'Type','text'),'Color',[.1,.1,.1]);set(findall(fig,'Type','legend'),'Color','w','TextColor',[.1,.1,.1]);
sgtitle('Frozen late training tau 9.243-10.643; real future holdout 10.843-12.243 | rank-one screen rejected');
exportgraphics(fig,fullfile(rootDirectory,'late_inner_holdout.png'),'Resolution',150);
exportgraphics(fig,fullfile(rootDirectory,'late_inner_holdout.pdf'),'ContentType','vector');
end
function val=sample(t,z,method)
if strcmp(method,'linear'),val=interp1(t.x,t.v,z,'linear');return;end
j=discretize(z,t.x);j(z==t.x(end))=numel(t.x)-1;h=t.x(j+1)-t.x(j);q=(z-t.x(j))./h;
val=(2*q.^3-3*q.^2+1).*t.v(j)+(q.^3-2*q.^2+q).*h.*t.d(j)+(-2*q.^3+3*q.^2).*t.v(j+1)+(q.^3-q.^2).*h.*t.d(j+1);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
