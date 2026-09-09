function report=ipm_accellab_compare_c1_cache_reports(reportFiles,outputRoot)
%IPM_ACCELLAB_COMPARE_C1_CACHE_REPORTS Paired physical-time C1 residual fields.
% Both charts use their own complete geometry; the common inner quadrature
% splits both native partitions, so no observer node density is selected.
assert(iscell(reportFiles) && numel(reportFiles)==2);
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['continuous_cache_pair_',token]);mkdir(destination);
registration=struct('reportFiles',{reportFiles},'outputDirectory',destination, ...
    'order',{{'baseline','finer_spatial_branch'}},'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0, ...
    'comparison','Actual U_t in each instantaneous normalized C1 chart, at the same physical time. No fitted alignment or rate adjustment.', ...
    'rectangles',struct('fullLocal',[-2,2,0,3],'core',[-1,1,0,1.5]));
save(fullfile(destination,'registration.mat'),'registration');
try
    sources=cell(1,2);data=cell(1,2);ax=[];ay=[];
    for k=1:2
        q=load(reportFiles{k},'report');sources{k}=q.report;
        assert(strcmp(sources{k}.status,'completed_audited_cache_no_LU'));
        data{k}=load(sources{k}.registration.cacheFile,'Omega','FX','x','y','Dx','Dy','scale','cacheAudit');
        assert(isequaln(data{k}.cacheAudit,sources{k}.cacheAudit) && isequaln(data{k}.scale,sources{k}.freshScale));
        r=sources{k}.coordinates;
        ax=[ax;(data{k}.x(:)-r.peak.x)/r.wallCoreWidth]; %#ok<AGROW>
        ay=[ay;data{k}.y(:)/r.verticalCoreWidth]; %#ok<AGROW>
    end
    timeError=abs(sources{1}.absolutePhysicalTime-sources{2}.absolutePhysicalTime);
    assert(timeError<1e-12,'ipm:C1CachePairTime','Different physical endpoints cannot be compared as a spatial pair.');
    assert(isequaln(data{1}.cacheAudit.lineage,data{2}.cacheAudit.lineage),'ipm:C1CachePairLineage','Fresh lineage differs.');
    [xi,wx]=axis_gauss(ax,[-2,-1,0,1,2]);[eta,wy]=axis_gauss(ay,[0,1.5,3]);
    w=wy(:)*wx(:)';[XI,ETA]=meshgrid(xi,eta);core=abs(XI)<1 & ETA<1.5;
    chart=cell(1,2);summary=cell(1,2);
    for k=1:2
        d=data{k};s=sources{k};r=s.coordinates;
        chart{k}=ipm_accellab_pullback_hermite(d.Omega,d.FX,d.x,d.y,d.Dx,d.Dy,r,xi,eta);
        chart{k}.GPhysical=chart{k}.G*s.physicalUnits.canonicalTimeRate;
        chart{k}.GParent=chart{k}.G/d.cacheAudit.lineage.canonicalCovarianceFactor;
        summary{k}=struct('sourceCheckpoint',s.cacheAudit.sourceCheckpoint,'gridSize',s.gridSize, ...
            'absolutePhysicalTime',s.absolutePhysicalTime,'parentEquivalentTau',s.parentEquivalentTau, ...
            'physicalUnits',s.physicalUnits,'parentUnits',s.parentUnits,'coordinates',r, ...
            'commonGaussPhysicalG',metric(chart{k}.GPhysical,core,w),'commonGaussU',metric(chart{k}.U,core,w));
        error=abs(summary{k}.commonGaussPhysicalG.coreL2-s.physicalUnits.UTimeDerivative.coreL2);
        summary{k}.independentPartitionCoreNormAbsoluteDifference=error;
        assert(error<1e-10*max(1,s.physicalUnits.UTimeDerivative.coreL2));
    end
    differenceG=chart{2}.GPhysical-chart{1}.GPhysical;differenceU=chart{2}.U-chart{1}.U;
    dG=metric(differenceG,core,w);dU=metric(differenceU,core,w);ratio=struct();correlation=struct();
    for key={'core','holdout','fullLocal'}
        name=key{1};metricName=[name,'L2'];
        ratio.(name)=struct('fineOverBaselineResidualNorm',summary{2}.commonGaussPhysicalG.(metricName)/summary{1}.commonGaussPhysicalG.(metricName), ...
            'residualFieldDifferenceOverBaseline',dG.(metricName)/summary{1}.commonGaussPhysicalG.(metricName), ...
            'profileFieldDifferenceOverBaseline',dU.(metricName)/summary{1}.commonGaussU.(metricName));
        mask=true(size(core));if strcmp(name,'core'),mask=core;elseif strcmp(name,'holdout'),mask=~core;end
        v1=chart{1}.GPhysical(mask);v2=chart{2}.GPhysical(mask);ww=w(mask);
        correlation.(name)=sum(ww.*v1.*v2)/sqrt(sum(ww.*v1.^2)*sum(ww.*v2.^2));
    end
    report=struct('status','completed_same_physical_time_spatial_residual_pair','registration',registration, ...
        'absolutePhysicalTimeDifference',timeError,'lineageExactlyPaired',true,'sources',{summary}, ...
        'residualDifference',dG,'profileDifference',dU,'ratios',ratio,'residualCorrelation',correlation, ...
        'commonGaussNodeCount',[numel(xi),numel(eta)],'coreArea',sum(w(core)),'holdoutArea',sum(w(~core)), ...
        'interpretation','Spatial plus transfer/evolution sensitivity between two registered branches. Norm and residual-field differences are both reported. This does not identify interpolation as the cause, establish continuum convergence, or accept an accelerator. Infinity norms are Gauss sampled maxima. Original mesh decisions remain unchanged.');
    save(fullfile(destination,'comparison_fields.mat'),'chart','xi','eta','w','core','differenceG','differenceU','-v7.3');
    save(fullfile(destination,'report.mat'),'report');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',registration,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('C1_CACHE_PAIR %s coreNormRatio=%.9g coreFieldDifference=%.9g holdNormRatio=%.9g holdFieldDifference=%.9g\n', ...
    destination,ratio.core.fineOverBaselineResidualNorm,ratio.core.residualFieldDifferenceOverBaseline, ...
    ratio.holdout.fineOverBaselineResidualNorm,ratio.holdout.residualFieldDifferenceOverBaseline);
end
function [q,w]=axis_gauss(native,breaks)
edges=unique([breaks(:);native(native>min(breaks)&native<max(breaks))]);
z=[-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7),sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
v=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
qq=(edges(1:end-1)+edges(2:end))/2+diff(edges)*z/2;ww=diff(edges)*v/2;
[q,~,map]=unique(qq(:));w=accumarray(map,ww(:));q=q(:)';w=w(:)';
end
function m=metric(v,core,w)
m=struct();
for names={'core','holdout','fullLocal'}
    name=names{1};mask=true(size(core));if strcmp(name,'core'),mask=core;elseif strcmp(name,'holdout'),mask=~core;end
    m.([name,'L2'])=sqrt(sum(w(mask).*v(mask).^2)/sum(w(mask)));
    m.([name,'SampledInfinity'])=max(abs(v(mask)));
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
