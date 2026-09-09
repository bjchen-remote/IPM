function report=ipm_accellab_compare_smooth_cache_fields(smoothReportFile,outputRoot)
%IPM_ACCELLAB_COMPARE_SMOOTH_CACHE_FIELDS Fixed smooth-observer field errors.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
s=load(smoothReportFile,'report');source=s.report;assert(numel(source.cases)==3);
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['smooth_cache_field_comparison_',token]);mkdir(destination);
data=cell(1,3);old=cell(1,3);
for j=1:3
    q=load(source.cases{j}.sourceReport,'report');old{j}=q.report;
    data{j}=load(old{j}.registration.cacheFile,'Omega','FX','x','y','Dx','Dy','scale');
    assert(abs(old{j}.absolutePhysicalTime-old{1}.absolutePhysicalTime)<1e-12);
end
pairs=[1,2;2,3;1,3];comparisons=cell(3,3);
for k=1:3
    for row=1:3
        indices=pairs(row,:);axesX=[];axesY=[];
        for j=indices
            r=source.cases{j}.observers{k}.coordinatesUsedForResidual;
            axesX=[axesX;(data{j}.x(:)-r.peak.x)/r.wallCoreWidth]; %#ok<AGROW>
            axesY=[axesY;data{j}.y(:)/r.verticalCoreWidth]; %#ok<AGROW>
        end
        [xi,wx]=axis_gauss(axesX,[-2,-1,0,1,2]);[eta,wy]=axis_gauss(axesY,[0,1.5,3]);
        [XI,ETA]=meshgrid(xi,eta);core=abs(XI)<1 & ETA<1.5;weights=wy(:)*wx(:)';
        G=cell(1,2);U=cell(1,2);physicalWidths=zeros(2,2);betas=zeros(2,2);
        for jj=1:2
            j=indices(jj);d=data{j};r=source.cases{j}.observers{k};
            c=ipm_accellab_pullback_hermite(d.Omega,d.FX,d.x,d.y,d.Dx,d.Dy,r.coordinatesUsedForResidual,xi,eta);
            G{jj}=c.G*exp(d.scale.logC_l-d.scale.logC_omega);U{jj}=c.U;
            physicalWidths(jj,:)=[r.rawPhysicalWidthX,r.rawPhysicalWidthY];betas(jj,:)=[r.parentBeta,r.parentGamma];
        end
        regions=struct();
        for name={'core','holdout','fullLocal'}
            region=name{1};mask=true(size(core));if strcmp(region,'core'),mask=core;elseif strcmp(region,'holdout'),mask=~core;end
            w=weights(mask);g1=G{1}(mask);g2=G{2}(mask);u1=U{1}(mask);u2=U{2}(mask);
            regions.(region)=struct('firstResidualL2',normw(g1,w),'secondResidualL2',normw(g2,w), ...
                'secondOverFirstResidualNorm',normw(g2,w)/normw(g1,w), ...
                'residualDifferenceL2',normw(g2-g1,w),'residualDifferenceOverFirst',normw(g2-g1,w)/normw(g1,w), ...
                'profileDifferenceOverFirst',normw(u2-u1,w)/normw(u1,w), ...
                'residualDifferenceSampledInfinity',max(abs(g2-g1)), ...
                'residualCorrelation',sum(w.*g1.*g2)/sqrt(sum(w.*g1.^2)*sum(w.*g2.^2)));
        end
        comparisons{row,k}=struct('caseIndices',indices,'theta',source.cases{1}.observers{k}.theta, ...
            'isPrimary',k==1,'regions',regions,'rawPhysicalWidths',physicalWidths, ...
            'rawPhysicalWidthRelativeDifference',physicalWidths(2,:)./physicalWidths(1,:)-1, ...
            'parentBetaGamma',betas,'parentBetaGammaDifference',betas(2,:)-betas(1,:), ...
            'quadratureNodeCount',[numel(xi),numel(eta)]);
        save(fullfile(destination,sprintf('pair%d_theta%d.mat',row,round(source.cases{1}.observers{k}.theta*10))), ...
            'G','U','xi','eta','weights','core','-v7.3');
    end
end
report=struct('status','completed_fixed_smooth_observer_three_grid_fields','sourceSmoothReport',smoothReportFile, ...
    'outputDirectory',destination,'comparisons',{comparisons},'unitCalibration',source.unitCalibration, ...
    'primaryTheta',.5,'holdoutTheta',[.3,.7],'originalDecisionsChanged',false,'newAcceptanceGate',false, ...
    'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0, ...
    'interpretation','All fixed primary/holdout charts are compared at common physical time using union-native-cell Gauss quadrature. Raw uncalibrated physical widths and rate differences are retained. A more stable observer diagnoses sensitivity, not improved PDE accuracy or permission to accept old failed candidates. Nonuniform local refinements and shared transferred initial data do not establish a global convergence order.');
save(fullfile(destination,'report.mat'),'report');
fid=fopen(fullfile(destination,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('SMOOTH_FIELDS %s\n',destination);
for k=1:3
    for row=1:3
        c=comparisons{row,k};fprintf('SMOOTH_FIELD_PAIR %d-%d theta=%.1f coreDiff=%.8g holdDiff=%.8g betaDiff=%.3e\n', ...
            c.caseIndices,c.theta,c.regions.core.residualDifferenceOverFirst,c.regions.holdout.residualDifferenceOverFirst,c.parentBetaGammaDifference(1));
    end
end
end
function [q,w]=axis_gauss(native,breaks)
edges=unique([breaks(:);native(native>min(breaks)&native<max(breaks))]);
z=[-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7),sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
v=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
qq=(edges(1:end-1)+edges(2:end))/2+diff(edges)*z/2;ww=diff(edges)*v/2;
[q,~,map]=unique(qq(:));w=accumarray(map,ww(:));q=q(:)';w=w(:)';
end
function v=normw(f,w)
v=sqrt(sum(w.*f.^2)/sum(w));
end
