function report=ipm_accellab_smooth_cache_reports(reportFiles,outputRoot,calibrationFile)
%IPM_ACCELLAB_SMOOTH_CACHE_REPORTS New observer; original decisions retained.
% All theta values are fixed before data evaluation. Calibration is a single
% baseline unit convention, never recomputed for a finer case or later time.
if nargin<3,calibrationFile='';end
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['smooth_cache_observer_',token]);mkdir(destination);
p=struct('reportFiles',{reportFiles},'outputDirectory',destination,'primaryTheta',.5,'holdoutTheta',[.3,.7], ...
    'power',4,'calibrationFile',calibrationFile,'directionalSteps',[1e-5,5e-6,2.5e-6], ...
    'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0,'originalDecisionsChanged',false);
save(fullfile(destination,'registration.mat'),'p');write_json(fullfile(destination,'registration.json'),p);
thetas=[.5,.3,.7];cases=cell(1,numel(reportFiles));
unitCalibration=repmat(struct('theta',[],'horizontalFactor',[],'verticalFactor',[], ...
    'sourceCheckpoint','','sourceAbsolutePhysicalTime',[],'sourceOldWallWidth',[], ...
    'sourceOldVerticalWidth',[],'rule',''),1,3);
if ~isempty(calibrationFile)
    c=load(calibrationFile,'report');unitCalibration=c.report.unitCalibration;
    assert(isequal([unitCalibration.theta],thetas));
end
try
    for j=1:numel(reportFiles)
        old=load(reportFiles{j},'report');old=old.report;
        assert(strcmp(old.status,'completed_audited_cache_no_LU'));
        d=load(old.registration.cacheFile,'Omega','FX','x','y','Dx','Dy','scale','cacheAudit');
        assert(isequaln(d.cacheAudit,old.cacheAudit));rows=cell(1,3);
        for k=1:3
            theta=thetas(k);r=ipm_accellab_smooth_inner_rates(d.Omega,d.FX,d.x,d.y,d.Dx,d.Dy,old.localWindow,theta);
            assert(r.valid,'ipm:SmoothCachedGeometry','A registered smooth-width observer is invalid.');
            assert(r.peak.x==old.coordinates.peak.x && r.peak.value==old.coordinates.peak.value && ...
                r.peakPrime==old.coordinates.peakPrime && r.translationRate==old.coordinates.translationRate);
            if j==1 && isempty(calibrationFile)
                unitCalibration(k).theta=theta;
                unitCalibration(k).horizontalFactor=old.coordinates.wallCoreWidth/r.rawWallWidth;
                unitCalibration(k).verticalFactor=old.coordinates.verticalCoreWidth/r.rawVerticalWidth;
                unitCalibration(k).sourceCheckpoint=old.cacheAudit.sourceCheckpoint;
                unitCalibration(k).sourceAbsolutePhysicalTime=old.absolutePhysicalTime;
                unitCalibration(k).sourceOldWallWidth=old.coordinates.wallCoreWidth;
                unitCalibration(k).sourceOldVerticalWidth=old.coordinates.verticalCoreWidth;
                unitCalibration(k).rule='One fixed baseline unit factor matches the prior 90%-width rectangle; no recalibration by case, time, or error.';
            end
            c=unitCalibration(k);
            coordinates=struct('kind','new_smooth_width_with_fixed_baseline_units','valid',true,'peak',r.peak, ...
                'peakPrime',r.peakPrime,'translationRate',r.translationRate, ...
                'wallCoreWidth',c.horizontalFactor*r.rawWallWidth,'verticalCoreWidth',c.verticalFactor*r.rawVerticalWidth, ...
                'wallCoreWidthPrime',c.horizontalFactor*r.rawWallWidthPrime, ...
                'verticalCoreWidthPrime',c.verticalFactor*r.rawVerticalWidthPrime, ...
                'logScaleXRate',r.logScaleXRate,'logScaleYRate',r.logScaleYRate,'signature',r.signature);
            [residual,details]=ipm_accellab_continuous_residual(d.Omega,d.FX,d.x,d.y,d.Dx,d.Dy,coordinates);
            lambda=d.cacheAudit.lineage.canonicalCovarianceFactor;clock=exp(d.scale.logC_l-d.scale.logC_omega);
            directionErrors=zeros(size(p.directionalSteps));peakCellsStable=true(size(p.directionalSteps));
            for hIndex=1:numel(p.directionalSteps)
                h=p.directionalSteps(hIndex);zero=zeros(size(d.Omega));
                plus=ipm_accellab_smooth_inner_rates(d.Omega+h*d.FX,zero,d.x,d.y,d.Dx,d.Dy,old.localWindow,theta);
                minus=ipm_accellab_smooth_inner_rates(d.Omega-h*d.FX,zero,d.x,d.y,d.Dx,d.Dy,old.localWindow,theta);
                assert(plus.valid && minus.valid);
                derivative=([plus.rawWallWidth,plus.rawVerticalWidth]-[minus.rawWallWidth,minus.rawVerticalWidth])/(2*h);
                actual=[r.rawWallWidthPrime,r.rawVerticalWidthPrime];
                directionErrors(hIndex)=max(abs(derivative-actual)./max(abs(actual),1e-14));
                peakCellsStable(hIndex)=plus.peak.selectedCell==r.peak.selectedCell && minus.peak.selectedCell==r.peak.selectedCell;
            end
            rows{k}=struct('theta',theta,'isPrimary',theta==.5,'rawRates',r,'fixedUnitCalibration',c, ...
                'coordinatesUsedForResidual',coordinates,'parentBeta',r.logScaleXRate/lambda,'parentGamma',r.logScaleYRate/lambda, ...
                'rawPhysicalWidthX',r.rawWallWidth/exp(d.scale.logC_l),'rawPhysicalWidthY',r.rawVerticalWidth/exp(d.scale.logC_l), ...
                'residual',residual,'parentCoreResidualL2',residual.gauss.G.coreL2/lambda, ...
                'parentHoldoutResidualL2',residual.gauss.G.holdoutL2/lambda, ...
                'physicalCoreResidualL2',residual.gauss.G.coreL2*clock,'physicalHoldoutResidualL2',residual.gauss.G.holdoutL2*clock, ...
                'trueNativeDirectionalRelativeErrors',directionErrors,'directionalPeakCellsStable',peakCellsStable, ...
                'directionalMeaning','Central functional differences along Omega +/- h*trueCachedFX on the fixed native grid; not PDE time integration.', ...
                'directionalCheckPassed',all(peakCellsStable) && directionErrors(end)<1e-5);
            save(fullfile(destination,sprintf('case%d_theta%d_fields.mat',j,round(theta*10))),'details','-v7.3');
        end
        cases{j}=struct('sourceReport',reportFiles{j},'sourceCheckpoint',old.cacheAudit.sourceCheckpoint, ...
            'absolutePhysicalTime',old.absolutePhysicalTime,'parentEquivalentTau',old.parentEquivalentTau, ...
            'originalC1Residual',old.residual,'originalParentRates',old.parentUnits,'observers',{rows});
    end
    report=struct('status','completed_new_smooth_width_auxiliary_observation','registration',p,'cases',{cases}, ...
        'unitCalibration',unitCalibration,'allDirectionalChecksPassed',all(cellfun(@(c)all(cellfun(@(r)r.directionalCheckPassed,c.observers)),cases)), ...
        'originalDecisionsChanged',false,'interpretation','All primary/holdout thresholds are retained. Raw widths, rates and direction errors are reported independently of fixed baseline display units. This auxiliary observer cannot accept old candidates, replace the original 90%-width evidence, or establish continuum convergence.');
    save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
catch exception
    failure=struct('registration',p,'identifier',exception.identifier,'message',exception.message,'stack',exception.stack, ...
        'casesCompleted',{cases},'unitCalibration',unitCalibration);
    save(fullfile(destination,'failure.mat'),'failure');write_json(fullfile(destination,'failure.json'),failure);rethrow(exception);
end
fprintf('SMOOTH_CACHE_REPORT %s directions=%d\n',destination,report.allDirectionalChecksPassed);
for j=1:numel(cases)
    for k=1:3
        rr=cases{j}.observers{k};fprintf('SMOOTH_CASE %d theta=%.1f betaParent=%.10g gammaParent=%.10g coreParent=%.10g holdParent=%.10g fdError=%.3e\n', ...
            j,rr.theta,rr.parentBeta,rr.parentGamma,rr.parentCoreResidualL2,rr.parentHoldoutResidualL2,rr.trueNativeDirectionalRelativeErrors(end));
    end
end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
