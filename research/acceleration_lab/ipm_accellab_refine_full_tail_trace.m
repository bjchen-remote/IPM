function report=ipm_accellab_refine_full_tail_trace(registrationFile)
%IPM_ACCELLAB_REFINE_FULL_TAIL_TRACE One preregistered 512-point refinement.
% Keeps the failed 32/64/128/256 report immutable and uses its exact gates.
reg=jsondecode(fileread(registrationFile));out=reg.outputDirectory;
assert(maxNumCompThreads==10&&~isfolder(out));mkdir(out);
loaded=load(reg.previousReportMat,'report');old=loaded.report;
assert(~old.allPureChecksPassed&&isequal(old.orders,[32,64,128,256])&&reg.newOrder==512);
timer=tic;
[target,targetInfo]=ipm_accellab_green_full_tail_trace(old.query,old.H,method='target_polar',order=reg.newOrder);
[analytic,analyticInfo]=ipm_accellab_green_full_tail_trace(old.query,old.H,method='radial_analytic',order=reg.newOrder);
pert=ipm_accellab_green_full_tail_trace(old.query(1:old.boundaryQueryCount,:),old.H, ...
    method='target_polar',order=reg.newOrder,sourceOnly='perturbation');
scale=max(old.H,abs(analytic));nonzero=old.query(:,1)~=0&old.query(:,2)~=0;
qa=old.H*analytic(nonzero)./(old.query(nonzero,1).*old.query(nonzero,2));
qn=old.H*target(nonzero)./(old.query(nonzero,1).*old.query(nonzero,2));
report=reg;report.kind='complete_Green_exterior_trace_fixed_512_refinement_v1';
report.previousReportStillFailed=true;report.oldMaximumTargetLastChange=old.lastTargetPolarScaledChange;
report.target512=target;report.analytic512=analytic;report.perturbation512=pert;
report.targetInfo=targetInfo;report.analyticInfo=analyticInfo;
report.methodScaledError=max(abs(target-analytic)./scale);
report.target256To512ScaledChange=max(abs(target-old.targetCenteredPolar(:,end))./scale);
report.analytic256To512ScaledChange=max(abs(analytic-old.analyticRadial(:,end))./scale);
report.nearAxisQuotientScaledError=max(abs(qa-qn)./max(1,abs(qa)));
report.perturbation256To512ScaledChange=max(abs(pert-old.manufacturedPerturbation(:,end)))/max(1,max(abs(pert)));
report.originalTemplateVersus512ReferenceScaledError=max(abs(target(1:old.boundaryQueryCount)-old.analyticRadial(1:old.boundaryQueryCount,end)))/old.H;
report.sameQuadratureThresholds=struct('trace',old.maximumTraceScaledError, ...
    'nearAxisQuotient',old.maximumNearAxisQuotientScaledError);
report.refinedQuadraturePassed=max([report.methodScaledError,report.target256To512ScaledChange, ...
    report.analytic256To512ScaledChange,report.perturbation256To512ScaledChange])<=old.maximumTraceScaledError && ...
    report.nearAxisQuotientScaledError<=old.maximumNearAxisQuotientScaledError;
report.previouslyMeasuredOtherGatesPassed=old.allAxisTracesExactlyZero&&old.cornerDuplicatesExact&& ...
    old.leftRightOddSymmetryScaledError<=old.maximumSymmetryScaledError&& ...
    old.manufacturedSourceLinearityScaledError<=old.maximumTraceScaledError&& ...
    old.maximumCovarianceScaledError<=old.maximumCovarianceScaledError&& ...
    old.perturbationBoundaryHZeroScalingError<=old.maximumTraceScaledError&& ...
    old.maximumScalarReuseScaledError<=old.maximumTraceScaledError;
report.combinedPureEvidencePassed=report.refinedQuadraturePassed&&report.previouslyMeasuredOtherGatesPassed;
report.LUCount=0;report.PDECount=0;report.injectedIntoSolver=false;report.realEvolvedClosureValidated=false;
report.wallSeconds=toc(timer);
save(fullfile(out,'report.mat'),'report');write_json(fullfile(out,'report.json'),report);
disp(jsonencode(rmfield(report,{'target512','analytic512','perturbation512'})));
end
function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
