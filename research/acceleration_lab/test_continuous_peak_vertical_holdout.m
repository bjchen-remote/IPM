function report=test_continuous_peak_vertical_holdout(checkpointFiles,outputDirectory)
%TEST_CONTINUOUS_PEAK_VERTICAL_HOLDOUT Compare vertical observation choices.
% Four fully trusted same-axis native frames; neither field nor PDE changes.
assert(iscell(checkpointFiles)&&numel(checkpointFiles)==4&& ...
    ~isfolder(outputDirectory));
eta=linspace(0,2,401);
nearest=zeros(4,numel(eta));continuous=nearest;
times=zeros(1,4);offsets=zeros(1,4);caseIds=cell(1,4);
firstX=[];firstY=[];
for k=1:4
    cp=ipm.output.readCheckpoint(checkpointFiles{k});
    s=cp.payload.state;h=cp.payload.log.history;
    assert(all(ipm.output.trustedMask(h,s.config)));
    x=s.x(:)';y=s.y(:);
    if k==1,firstX=x;firstY=y;
    else,assert(isequal(x,firstX)&&isequal(y,firstY), ...
        'ipm:ContinuousVerticalAxesChanged');end
    times(k)=s.scale.canonicalTime;
    caseIds{k}=char(s.runMetadata.caseId);
    peakX=h.gauge.omegaGaugeQuadraticPeakX(end);
    peak=h.gauge.omegaGaugeQuadraticPeakValue(end);
    wx=h.mesh.trackedWallCoreWidth(end);
    wy=h.mesh.trackedVerticalCoreWidth(end);
    assert(all(isfinite([peakX,peak,wx,wy]))&&peak>0&&wx>0&&wy>0);
    Dx=ipm.mesh.fdMatrix(x,1,7);
    source=s.rho*Dx';
    [~,ix]=min(abs(x-peakX));
    offsets(k)=(x(ix)-peakX)/wx;
    nodeTrace=source(:,ix);
    continuousTrace=interp1(x,source',peakX,'pchip')';
    assert(numel(continuousTrace)==numel(y)&& ...
        all(isfinite(continuousTrace)));
    nearest(k,:)=interp1(y/wy,nodeTrace/peak,eta,'pchip');
    continuous(k,:)=interp1(y/wy,continuousTrace/peak,eta,'pchip');
end
assert(all(diff(times)>0)&&isscalar(unique(caseIds))&& ...
    all(isfinite(nearest),'all')&&all(isfinite(continuous),'all'));
pointDifference=sqrt(sum((nearest-continuous).^2,2)./ ...
    sum(continuous.^2,2));
report=struct('kind','continuous_peak_vertical_observer_holdout_v1', ...
    'checkpoints',{checkpointFiles},'times',times, ...
    'sameAxes',true,'fullHistoriesTrusted',true, ...
    'nearestColumnOffsetsInWallCoreWidths',offsets, ...
    'nearestVsContinuousRelativeL2',pointDifference', ...
    'nearest',one_shape(nearest,times), ...
    'continuous',one_shape(continuous,times), ...
    'physicalFieldPredicted',false,'originalTrajectoryModified',false);
mkdir(outputDirectory);
fid=fopen(fullfile(outputDirectory,'report.json'),'w');assert(fid>=0);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('CONTINUOUS_VERTICAL_HOLDOUT tau=%.6f nearest=%.6g/%.6g/%.6g continuous=%.6g/%.6g/%.6g\n', ...
    times(end),report.nearest.carryForwardRelativeL2, ...
    report.nearest.linearRelativeL2,report.nearest.rankOneRelativeL2, ...
    report.continuous.carryForwardRelativeL2, ...
    report.continuous.linearRelativeL2,report.continuous.rankOneRelativeL2);
end

function result=one_shape(values,times)
first=values(1,:);second=values(2,:);
source=values(3,:);heldout=values(4,:);
ratio=(times(4)-times(3))/(times(3)-times(2));
d1=second-first;d2=source-second;
lambda=dot(d1,d2)/dot(d1,d1);
cosine=dot(d1,d2)/(norm(d1)*norm(d2));
rankValid=isfinite(lambda)&&lambda>0;
rankError=NaN;coefficient=NaN;
if rankValid
    if abs(lambda-1)<1e-8,coefficient=ratio;
    else,coefficient=lambda*(lambda^ratio-1)/(lambda-1);end
    rankError=relative_l2(source+coefficient*d2,heldout);
end
result=struct('carryForwardRelativeL2',relative_l2(source,heldout), ...
    'linearRelativeL2',relative_l2(source+ratio*d2,heldout), ...
    'rankOneRelativeL2',rankError,'rankOneValid',rankValid, ...
    'rankOneLambda',lambda,'rankOneForecastCoefficient',coefficient, ...
    'trainingIncrementCosine',cosine);
end

function value=relative_l2(candidate,reference)
value=sqrt(sum((candidate-reference).^2)/sum(reference.^2));
end
