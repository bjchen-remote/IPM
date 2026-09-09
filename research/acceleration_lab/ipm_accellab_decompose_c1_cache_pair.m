function report=ipm_accellab_decompose_c1_cache_pair(pairReportFile,outputRoot)
%IPM_ACCELLAB_DECOMPOSE_C1_CACHE_PAIR Auxiliary attribution, never a new gate.
% Compare each complete native inner derivative at common physical queries.
directory=fileparts(mfilename('fullpath'));addpath(directory,fileparts(fileparts(directory)));
s=load(pairReportFile,'report');pair=s.report;
assert(strcmp(pair.status,'completed_same_physical_time_spatial_residual_pair'));
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['continuous_pair_terms_',token]);mkdir(destination);
sources=cell(1,2);data=cell(1,2);
for k=1:2
    q=load(pair.registration.reportFiles{k},'report');sources{k}=q.report;
    data{k}=load(sources{k}.registration.cacheFile,'Omega','FX','x','y','Dx','Dy','scale','cacheAudit');
end
r0=sources{1}.coordinates;c0=exp(data{1}.scale.logC_l);
physicalA=(r0.peak.x-data{1}.scale.X_shift)/c0;physicalWX=r0.wallCoreWidth/c0;physicalWY=r0.verticalCoreWidth/c0;
ax=[];ay=[];
for k=1:2
    d=data{k};cx=exp(d.scale.logC_l);
    ax=[ax;((d.x(:)-d.scale.X_shift)/cx-physicalA)/physicalWX]; %#ok<AGROW>
    ay=[ay;d.y(:)/cx/physicalWY]; %#ok<AGROW>
end
[xi,wx]=axis_gauss(ax,[-2,-1,0,1,2]);[eta,wy]=axis_gauss(ay,[0,1.5,3]);
[XI,ETA]=meshgrid(xi,eta);weights=wy(:)*wx(:)';core=abs(XI)<1 & ETA<1.5;
names={'rawForcing','amplitude','translation','horizontal','vertical'};
terms=cell(1,2);kernels=cell(1,2);coefficients=zeros(2,5);totals=cell(1,2);
for k=1:2
    d=data{k};s=sources{k};r=s.coordinates;cx=exp(d.scale.logC_l);clock=s.physicalUnits.canonicalTimeRate;
    X=cx*(physicalA+physicalWX*XI)+d.scale.X_shift;Y=cx*physicalWY*ETA;
    o=ipm_accellab_tensor_hermite(d.Omega,d.x,d.y,d.Dx,d.Dy,X,Y);
    f=ipm_accellab_tensor_hermite(d.FX,d.x,d.y,d.Dx,d.Dy,X,Y);
    P=r.peak.value;
    kernels{k}=cat(3,f.value,o.value/P,o.derivativeX/P,(X-r.peak.x).*o.derivativeX/P,Y.*o.derivativeY/P);
    coefficients(k,:)=[clock/P,-clock*r.peakPrime/P,clock*r.translationRate,clock*r.logScaleXRate,clock*r.logScaleYRate];
    terms{k}=kernels{k}.*reshape(coefficients(k,:),1,1,5);totals{k}=sum(terms{k},3);
end
deltaTotal=totals{2}-totals{1};deltaTerms=terms{2}-terms{1};
ratePart=kernels{1}.*reshape(coefficients(2,:)-coefficients(1,:),1,1,5);
fieldPart=(kernels{2}-kernels{1}).*reshape(coefficients(2,:),1,1,5);
identityError=max(abs(sum(deltaTerms,3)-deltaTotal),[],'all');
splitError=max(abs(deltaTerms-ratePart-fieldPart),[],'all');
assert(identityError<1e-10 && splitError<1e-10);
regions=struct();
for label={'core','holdout','fullLocal'}
    name=label{1};mask=true(size(core));if strcmp(name,'core'),mask=core;elseif strcmp(name,'holdout'),mask=~core;end
    w=weights(mask);total=deltaTotal(mask);den=sum(w.*total.^2);part=struct();
    for k=1:5
        v=deltaTerms(:,:,k);vr=ratePart(:,:,k);vf=fieldPart(:,:,k);v0=terms{1}(:,:,k);v1=terms{2}(:,:,k);
        part.(names{k})=struct('baselineL2',l2(v0(mask),w),'fineL2',l2(v1(mask),w),'differenceL2',l2(v(mask),w), ...
            'signedFractionOfTotalDifferenceEnergy',sum(w.*v(mask).*total)/den, ...
            'coefficientDifferenceL2',l2(vr(mask),w),'kernelDifferenceL2',l2(vf(mask),w), ...
            'coefficientSignedFraction',sum(w.*vr(mask).*total)/den,'kernelSignedFraction',sum(w.*vf(mask).*total)/den);
    end
    regions.(name)=struct('baselineTotalL2',l2(totals{1}(mask),w),'fineTotalL2',l2(totals{2}(mask),w), ...
        'totalDifferenceL2',l2(total,w),'differenceOverBaseline',l2(total,w)/l2(totals{1}(mask),w),'terms',part);
end
report=struct('status','auxiliary_common_physical_query_component_attribution','sourcePairReport',pairReportFile, ...
    'outputDirectory',destination,'sourceMeshDecisionsUnchanged',true,'newAcceptanceGate',false,'fitPerformed',false, ...
    'commonGeometry',struct('physicalPeakPosition',physicalA,'physicalWallWidth',physicalWX,'physicalVerticalWidth',physicalWY), ...
    'coefficientNames',{names},'coefficients',coefficients,'regions',regions,'sumIdentityMaximumError',identityError, ...
    'rateKernelSplitMaximumError',splitError,'quadratureNodeCount',[numel(xi),numel(eta)], ...
    'timeUnits','Both derivatives and all contributions are per actual physical time.', ...
    'queryMeaning','The baseline physical peak and widths define identical physical query points. Each native complete inner derivative retains its own true PPrime, aPrime, beta, gamma and peak; it is evaluated at those common points. This evaluation is not a claim that the fixed baseline chart itself has those velocities.', ...
    'splitMeaning','For each coefficient times spatial kernel, DeltaTerm=DeltaCoefficient*BaselineKernel+FineCoefficient*DeltaKernel. Signed energy fractions are inner products with the total difference, can be negative, and sum to one; they are ordered algebraic attribution, not independent causal experiments.', ...
    'newRhsEvaluations',0,'poissonOperatorBuilds',0,'pdeSteps',0, ...
    'interpretation','This only locates cancellation and spatial sensitivity in two fixed saved states. It cannot identify the continuum error, change the source decisions, or accept a new gauge/accelerator.');
save(fullfile(destination,'component_fields.mat'),'xi','eta','weights','core','terms','totals','deltaTerms','ratePart','fieldPart','-v7.3');
save(fullfile(destination,'report.mat'),'report');
fid=fopen(fullfile(destination,'report.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('C1_PAIR_TERMS %s coreDiff=%.8g betaFraction=%.8g rawFraction=%.8g\n', ...
    destination,regions.core.differenceOverBaseline,regions.core.terms.horizontal.signedFractionOfTotalDifferenceEnergy, ...
    regions.core.terms.rawForcing.signedFractionOfTotalDifferenceEnergy);
end
function [q,w]=axis_gauss(native,breaks)
edges=unique([breaks(:);native(native>min(breaks)&native<max(breaks))]);
z=[-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7),sqrt((3-2*sqrt(6/5))/7),sqrt((3+2*sqrt(6/5))/7)];
v=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
qq=(edges(1:end-1)+edges(2:end))/2+diff(edges)*z/2;ww=diff(edges)*v/2;
[q,~,map]=unique(qq(:));w=accumarray(map,ww(:));q=q(:)';w=w(:)';
end
function v=l2(f,w)
v=sqrt(sum(w.*f.^2)/sum(w));
end
