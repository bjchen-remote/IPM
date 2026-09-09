function report=ipm_accellab_recover_four_gauge(sourceDirectory,outputRoot)
%IPM_ACCELLAB_RECOVER_FOUR_GAUGE Read numerics saved before JSON export failed.
[~,token]=fileparts(tempname);destination=fullfile(outputRoot,['four_gauge_recovered_',token]);mkdir(destination);
registration=load(fullfile(sourceDirectory,'registration.mat'),'p');cases=cell(1,2);
for k=1:2
    grid=registration.p.grids(k,1);file=fullfile(sourceDirectory,sprintf('n%d',grid),'report.mat');
    data=load(file,'value');value=data.value;
    assert(strcmp(value.status,'completed') && value.pdeSteps==0);
    cases{k}=value;copyfile(file,fullfile(destination,sprintf('n%d_report_original.mat',grid)));
    write_json(fullfile(destination,sprintf('n%d_report.json',grid)),value);
    t=value.threeConstraintIsotropic;
    fprintf('RECOVER_GAUGES n%d fourResidual=%.3e fourRcond=%.8g seedResidual=%.3e threeResidual=%.3e gamma=%.8g covariance4=%.3e covariance3=%.3e\n', ...
        grid,norm(value.actualConstraintResidual,inf),value.finalJacobianRcond,norm(value.seedActualWenoResidual,inf), ...
        norm(t.actualConstraintResidual,inf),t.verticalLogWidthRateMeasured, ...
        value.physicalRhsCoordinateConsistency.wholeBoxRelativeL2,t.physicalRhsCoordinateConsistency.wholeBoxRelativeL2);
end
report=struct('status','recovered_completed_numerics_from_saved_mat','sourceDirectory',sourceDirectory, ...
    'sourceRegistration',registration.p,'cases',{cases},'fourConstraintPassed',all(cellfun(@(c)c.passed,cases)), ...
    'threeConstraintPassed',all(cellfun(@(c)c.threeConstraintIsotropic.passed,cases)), ...
    'sourceProcessExportFailed',true,'originalFailurePreserved',true,'newPoissonEvaluations',0,'newRhsEvaluations',0, ...
    'pdeSteps',0,'fullRhsFieldExportAvailable',false, ...
    'interpretation','Both instantaneous numerical studies completed and saved MAT reports before JSON serialization failed on an initial-condition function handle. This recovery only reads those exact reports; full RHS arrays were not exported before the original failure. No physical-time validation is claimed.');
save(fullfile(destination,'report.mat'),'report','-v7.3');write_json(fullfile(destination,'report.json'),report);
fprintf('FOUR_GAUGE_RECOVERED %s\n',destination);
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(ipm_accellab_json_view(value),'PrettyPrint',true));
end
