function report = analyze_profile_series(checkpointFiles,outputDirectory)
%ANALYZE_PROFILE_SERIES Read-only profile/scale diagnostics on paired grids.
% Each native checkpoint is validated. No extrapolation is a time trajectory.
% Inner curves remove each sample's translation, width and amplitude drift;
% their agreement is a shape diagnostic, not convergence in the chosen gauge.
assert(iscell(checkpointFiles) && numel(checkpointFiles)>=2);
assert(~isfolder(outputDirectory),'Refusing to overwrite profile analysis.');
originalThreads = maxNumCompThreads(10);
threadCleanup = onCleanup(@() maxNumCompThreads(originalThreads));
count = numel(checkpointFiles);
profiles = cell(1,count); records = cell(1,count);
for k = 1:count
    cp = ipm.output.readCheckpoint(checkpointFiles{k});
    state = cp.payload.state; h = cp.payload.log.history;
    assert(strcmp(state.config.transport.spatialDiscretization,'high_order') && ...
        strcmp(state.config.scaling.dynamicScaleGeometry,'isotropic') && ...
        strcmp(state.config.scaling.cOmegaGauge,'wall_omega_quadratic_peak'), ...
        'This diagnostic requires the fourth-order isotropic quadratic-peak contract.');
    assert(all(ipm.output.trustedMask(h,state.config)));
    x = state.x(:)'; y = state.y(:);
    D = ipm.mesh.fdMatrix(x,1,7);
    source = state.rho*D';
    g = h.gauge; m = h.mesh; c = h.common;
    peakX = g.omegaGaugeQuadraticPeakX(end);
    peak = g.omegaGaugeQuadraticPeakValue(end);
    wx = m.trackedWallCoreWidth(end); wy = m.trackedVerticalCoreWidth(end);
    curvature = g.omegaGaugeQuadraticCurvature(end);
    assert(all(isfinite([peakX,peak,wx,wy,curvature])) && ...
        peak>0 && wx>0 && wy>0 && curvature<0, ...
        'Quadratic peak, negative curvature and both positive core widths are required.');
    [~,ix] = min(abs(x-peakX));
    profiles{k} = struct('x',x,'y',y,'wallRX',source(1,:), ...
        'wallRho',state.rho(1,:),'physicalX',(x-state.scale.X_shift)/exp(state.scale.logC_l), ...
        'physicalWallRho',state.rho(1,:)/exp(state.scale.logC_omega), ...
        'physicalWallRX',source(1,:)*exp(state.scale.logC_l-state.scale.logC_omega), ...
        'innerX',(x-peakX)/wx,'innerY',y/wy, ...
        'verticalRX',source(:,ix),'verticalSampleX',x(ix), ...
        'verticalSampleIndex',ix,'peak',peak);
    record = struct('checkpointFile',checkpointFiles{k}, ...
        'caseId',char(state.runMetadata.caseId),'nodeCount',[numel(x),numel(y)], ...
        'tau',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
        'physicalRhoXInf',c.physicalRhoXInf(end), ...
        'growthEFolds',log(c.physicalRhoXInf(end)/c.physicalRhoXInf(1)), ...
        'cL',c.c_l(end),'cOmega',c.c_omega(end), ...
        'quadraticPeak',peak,'peakX',peakX, ...
        'wallCurvature',curvature,'curvatureLength',sqrt(peak/abs(curvature)), ...
        'physicalPeakX',(peakX-state.scale.X_shift)/exp(state.scale.logC_l), ...
        'rescaledWallSecondDerivativeInf',max(abs(source(1,:)*D')), ...
        'strictWallWidth',wx,'strictVerticalWidth',wy, ...
        'physicalWallWidth',wx/exp(state.scale.logC_l), ...
        'physicalVerticalWidth',wy/exp(state.scale.logC_l), ...
        'verticalSampleX',x(ix),'verticalSampleOffset',x(ix)-peakX, ...
        'verticalSampleOffsetInCoreWidth',(x(ix)-peakX)/wx, ...
        'strictWallCells',m.coreGridPoints(end), ...
        'strictVerticalCells',m.verticalCoreGridPoints(end), ...
        'safety',m.safetyFactor(end), ...
        'fullHistoryTrusted',true);
    records{k} = record;
end
records = [records{:}];
tau = [records.tau];
assert(all(diff(tau)>0),'Choose increasing checkpoints, with no duplicate time.');
report = struct('kind','profile_and_inner_scale_diagnostics', ...
    'records',records,'interpretation', ...
    'Narrowing cores and growing curvature require space/time/box refinement; no singularity proof.', ...
    'checkpointValidationThreads',10, ...
    'innerNormalization', ...
    ['Each checkpoint uses its own quadratic peak position, full strict 0.9 widths ' ...
    'and peak amplitude. Agreement measures shape after removing drift and contraction, ' ...
    'not convergence of R in the original rescaled coordinates. The vertical trace ' ...
    'uses a nearest-node column whose position and offset are reported.'], ...
    'innerWallConsecutiveRelativeL2',NaN(1,count), ...
    'innerVerticalConsecutiveRelativeL2',NaN(1,count));
report.sameCaseId = isscalar(unique({records.caseId}));
xi = linspace(-1,1,401); eta = linspace(0,2,401);
for k = 2:count
    p = profiles{k-1}; q = profiles{k};
    a = interp1(p.innerX,p.wallRX/p.peak,xi,'pchip');
    b = interp1(q.innerX,q.wallRX/q.peak,xi,'pchip');
    report.innerWallConsecutiveRelativeL2(k) = shared_relative_l2(xi,a,b);
    a = interp1(p.innerY,p.verticalRX/p.peak,eta,'pchip');
    b = interp1(q.innerY,q.verticalRX/q.peak,eta,'pchip');
    report.innerVerticalConsecutiveRelativeL2(k) = shared_relative_l2(eta,a,b);
end
mkdir(outputDirectory);
f = figure('Visible','off','Color','w','Position',[100,100,1500,850]);
cleaner = onCleanup(@() close(f));
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
colors = parula(count); labels = compose('tau = %.3f',tau);
nexttile; hold on
for k=1:count
    p=profiles{k}; plot(p.x,p.wallRX/p.peak,'Color',colors(k,:),'LineWidth',1.4);
end
xlim([0.65,1.15]); xlabel('X'); ylabel('R_X(X,0) / P');
title('Wall profile in current rescaling'); grid on; legend(labels,'Location','best');
nexttile; hold on
for k=1:count
    p=profiles{k}; plot(p.innerX,p.wallRX/p.peak,'Color',colors(k,:),'LineWidth',1.4);
end
xlim([-1,1]); xlabel('(X-X_{peak}) / width_{0.9,x}'); ylabel('R_X / P');
title('Wall profile at its inner scale'); grid on
nexttile; hold on
for k=1:count
    p=profiles{k}; plot(p.innerY,p.verticalRX/p.peak,'Color',colors(k,:),'LineWidth',1.4);
end
xlim([0,3]); xlabel('Y / width_{0.9,y}'); ylabel('R_X(X_{node},Y) / P');
title('Vertical profile at nearest peak node'); grid on
nexttile; semilogy(tau,[records.strictWallWidth],'-o', ...
    tau,[records.strictVerticalWidth],'-s',tau,[records.curvatureLength],'-^');
xlabel('tau'); ylabel('Length in rescaled coordinates'); grid on
legend('wall 0.9 core','vertical 0.9 core','sqrt(P / |Q''''|)','Location','best');
title('Core widths continue to contract');
nexttile; plot(tau,[records.cL],'-o',tau,[records.cOmega],'-s');
xlabel('tau'); ylabel('Exact gauge rates'); grid on; legend('c_l','c_omega');
title('Rates retain their measured drift');
nexttile; semilogy([records.physicalTime],[records.physicalRhoXInf],'-o');
xlabel('Physical time'); ylabel('Physical max |rho_x|'); grid on
title(sprintf('Latest growth: %.3f e-folds',records(end).growthEFolds));
sgtitle('IPM: paired trusted profiles, finite-time evidence only');
exportgraphics(f,fullfile(outputDirectory,'profile_progress.png'),'Resolution',160);
% Actual amplitudes are retained in this companion plot. Only the spatial
% origin and additive density level are aligned independently at each frame.
f2=figure('Visible','off','Color','w','Position',[100,100,1500,440]);
cleaner2=onCleanup(@() close(f2));
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;hold on
for k=1:count
    p=profiles{k}; center=records(k).physicalPeakX;
    centerValue=interp1(p.physicalX,p.physicalWallRho,center,'pchip');
    plot(p.physicalX-center,p.physicalWallRho-centerValue, ...
        'Color',colors(k,:),'LineWidth',1.4);
end
xlim([-.15,.15]);xlabel('Physical x - current peak x');
ylabel('rho(x,0) - rho(current peak,0)');grid on
title('Physical wall density: amplitude retained');legend(labels,'Location','best');
nexttile;hold on
for k=1:count
    p=profiles{k};plot(p.physicalX-records(k).physicalPeakX,p.physicalWallRX, ...
        'Color',colors(k,:),'LineWidth',1.4);
end
xlim([-.15,.15]);xlabel('Physical x - current peak x');ylabel('Physical rho_x');
grid on;title('Resolved gradient sharpening');
nexttile;semilogy(tau,abs([records.wallCurvature])/abs(records(1).wallCurvature),'-o', ...
    tau,[records.rescaledWallSecondDerivativeInf]/records(1).rescaledWallSecondDerivativeInf,'-s');
xlabel('tau');ylabel('Growth relative to first displayed frame');grid on
legend('|quadratic curvature of R_X|','max |R_{XX}(X,0)|','Location','best');
title('Derivatives in the original rescaling');
sgtitle('Finite-time sharpening; no fitted singular time or limit exponent');
exportgraphics(f2,fullfile(outputDirectory,'density_sharpening.png'),'Resolution',160);
fid = fopen(fullfile(outputDirectory,'profile_metrics.json'),'w');
assert(fid>=0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
save(fullfile(outputDirectory,'profile_data.mat'),'report','profiles','-v7.3');
fprintf('PROFILE_ANALYSIS %s\n',jsonencode(report));
end

function value = shared_relative_l2(axis,reference,candidate)
assert(all(isfinite(reference)) && all(isfinite(candidate)), ...
    'The requested inner comparison window must be covered by both paired grids.');
denominator = trapz(axis,candidate.^2);
assert(isfinite(denominator) && denominator>0, ...
    'The inner comparison norm must be finite and nonzero.');
value = sqrt(trapz(axis,(candidate-reference).^2)/denominator);
end
