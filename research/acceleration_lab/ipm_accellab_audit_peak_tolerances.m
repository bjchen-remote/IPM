function report = ipm_accellab_audit_peak_tolerances(outputRoot,integrationDirectory)
%IPM_ACCELLAB_AUDIT_PEAK_TOLERANCES Preserve old false rejections explicitly.
directory = fileparts(mfilename('fullpath')); addpath(directory,fileparts(fileparts(directory)));
[~,token] = fileparts(tempname); destination = fullfile(outputRoot,['hermite_peak_tolerances_',token]);
mkdir(destination);
registration = struct('outputDirectory',destination,'centers',[.70123,.70710678,.734567,.713579], ...
    'widths',[.027,.0173,.043],'windowRight',[2,1e3,1e6], ...
    'coordinateUnitFactors',[1e-6,1,1e6],'amplitudeUnitFactors',[1e-12,1,1e12], ...
    'legacyHelper','ipm_accellab_hermite_peak_window_scaled_v0', ...
    'newPolicy','local_cell_polynomial_scale_and_local_coordinate_ulp_v1', ...
    'newPdeSteps',0,'poissonOperatorBuilds',0);
save(fullfile(destination,'registration.mat'),'registration');
x = [linspace(-1,2,601),logspace(log10(2.05),6,81)]; D = ipm.mesh.fdMatrix(x,1,7);
entries = {}; unitEntries = {}; amplitudeEntries = {}; legacyRejections = 0; maximumFieldDifference = 0;
try
    for a = registration.centers
        for w = registration.widths
            v = exp(-((x-a)/w).^2); forcing = .13*v+.002*(x-a)/w^2.*v;
            reference = ipm_accellab_hermite_peak(v,forcing,x,D,[0,2]);
            assert(reference.valid);
            for right = registration.windowRight
                old = ipm_accellab_hermite_peak_window_scaled_v0(v,forcing,x,D,[0,right]);
                current = ipm_accellab_hermite_peak(v,forcing,x,D,[0,right]);
                assert(current.valid && isequal(current.value,reference.value) && ...
                    isequal(current.PPrime,reference.PPrime) && isequal(current.x,reference.x), ...
                    'ipm:PeakToleranceFarEndpoint','A far search endpoint changed the local maximum or its validity.');
                legacyRejections = legacyRejections+double(~old.valid);
                entries{end+1} = struct('center',a,'width',w,'windowRight',right, ...
                    'legacy',compact(old),'current',compact(current), ...
                    'unchangedRawMaximumValue',isequal(old.value,current.value)); %#ok<AGROW>
            end
            for factor = registration.coordinateUnitFactors
                old = ipm_accellab_hermite_peak_window_scaled_v0(v,forcing,factor*x,D/factor,[0,factor*1e6]);
                current = ipm_accellab_hermite_peak(v,forcing,factor*x,D/factor,[0,factor*1e6]);
                errors = abs([current.value-reference.value,current.PPrime-reference.PPrime,current.x/factor-reference.x]);
                assert(current.valid && max(errors) < 1e-11, ...
                    'ipm:PeakToleranceCoordinateUnits','Changing coordinate units changed the resolved peak.');
                unitEntries{end+1} = struct('center',a,'width',w,'unitFactor',factor, ...
                    'valueForcingPositionErrors',errors,'legacy',compact(old),'current',compact(current)); %#ok<AGROW>
            end
            for factor = registration.amplitudeUnitFactors
                current = ipm_accellab_hermite_peak(factor*v,factor*forcing,x,D,[0,1e6]);
                errors = abs([current.value/factor-reference.value,current.PPrime/factor-reference.PPrime,current.x-reference.x]);
                assert(current.valid && max(errors) < 1e-11, ...
                    'ipm:PeakToleranceAmplitudeUnits','Changing amplitude units changed the resolved peak.');
                amplitudeEntries{end+1} = struct('center',a,'width',w,'unitFactor',factor, ...
                    'normalizedValueForcingPositionErrors',errors,'current',compact(current)); %#ok<AGROW>
            end
        end
    end
    assert(legacyRejections > 0,'ipm:PeakToleranceCounterexample','The registered old-policy false rejection was not reproduced.');
    regression = ipm_accellab_test_hermite_peak();
    auditedFrames = 0; bitwisePeakValues = true; bitwisePeakLocations = true; equalPeakCells = true;
    for nx = [49,65]
        for level = 1:3
            file = fullfile(integrationDirectory,sprintf('research_hermite_C1_n%d_level%d.mat',nx,level));
            data = load(file,'run'); r = data.run;
            for k = 1:numel(r.accepted)
                v = r.rhoSnapshots{k}(1,:)*r.Dx';
                p = ipm_accellab_hermite_peak(v,zeros(size(v)),r.x,r.Dx,[0,r.x(end)]);
                assert(p.valid);
                bitwisePeakValues = bitwisePeakValues && isequal(p.value,r.accepted{k}.H);
                bitwisePeakLocations = bitwisePeakLocations && isequal(p.x,r.accepted{k}.HLocation);
                equalPeakCells = equalPeakCells && p.selectedCell == r.accepted{k}.HCell;
                maximumFieldDifference = max(maximumFieldDifference,abs(p.value-r.accepted{k}.H));
                auditedFrames = auditedFrames+1;
            end
        end
    end
    assert(bitwisePeakValues && bitwisePeakLocations && equalPeakCells, ...
        'ipm:PeakToleranceArchivedValues','The new tolerance policy changed an already accepted tiny peak.');
    report = struct('status','passed_with_preserved_legacy_false_rejections','registration',registration, ...
        'farEndpointExamples',{entries},'coordinateUnitExamples',{unitEntries}, ...
        'amplitudeUnitExamples',{amplitudeEntries},'legacyFalseRejectionCount',legacyRejections, ...
        'mmsAndDegeneracyRegression',regression, ...
        'archivedTinyAudit',struct('sourceDirectory',integrationDirectory,'auditedFrames',auditedFrames, ...
        'bitwisePeakValues',bitwisePeakValues,'bitwisePeakLocations',bitwisePeakLocations, ...
        'equalPeakCells',equalPeakCells,'maximumPeakValueDifference',maximumFieldDifference), ...
        'interpretation','Only the validity metric uses local cell polynomial scales and local coordinate ULPs. The maximum functional and rate formulas are unchanged; weak/flat/tied/boundary guards remain explicit. No q512 or production rule is changed.');
    save(fullfile(destination,'report.mat'),'report'); write_json(fullfile(destination,'report.json'),report);
catch exception
    failure = struct('registration',registration,'farEndpointExamples',{entries},'coordinateUnitExamples',{unitEntries}, ...
        'amplitudeUnitExamples',{amplitudeEntries},'identifier',exception.identifier,'message',exception.message);
    save(fullfile(destination,'failure.mat'),'failure'); write_json(fullfile(destination,'failure.json'),failure); rethrow(exception);
end
fprintf('HERMITE_PEAK_TOLERANCE_AUDIT %s\n',destination);
end

function q = compact(p)
q = rmfield(p,{'candidatePositions','candidateValues','candidateCells'});
end

function write_json(file,value)
fid = fopen(file,'w'); assert(fid >= 0); cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
