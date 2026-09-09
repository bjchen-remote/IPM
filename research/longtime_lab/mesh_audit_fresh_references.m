function report=mesh_audit_fresh_references(parentFile,childFile,outputFile)
%MESH_AUDIT_FRESH_REFERENCES Compare reference choices at one child state.
% Uses the child's one restored LU; no PDE, no checkpoint mutation. This
% audit does not install parent references into the fresh branch history.
assert(~isfile(outputFile));
originalThreads=maxNumCompThreads(10);cleanup=onCleanup(@()maxNumCompThreads(originalThreads));
parent=ipm.output.readCheckpoint(parentFile);p=parent.payload.state;
child=ipm.output.readCheckpoint(childFile);[state,~,~]=ipm.output.restoreCheckpoint(child,struct());
assert(strcmp(state.config.scaling.cOmegaGauge,'wall_omega_quadratic_peak') && ...
    strcmp(state.config.scaling.lengthGauge,'transport_anchor'));
Cx=exp(p.scale.logC_l);lambda=exp(p.scale.logC_l-p.scale.logC_omega);
transformed=p.rescaling;
transformed.pinX=transformed.pinX/Cx;
transformed.peakTrackingHalfWidth=transformed.peakTrackingHalfWidth/Cx;
transformed.strainTarget=transformed.strainTarget*lambda;
% This selected gauge has no frozen amplitude functional. A future gauge
% needs a reviewed transformation for every additional dimensional reference.
names=fieldnames(transformed);
assert(~any(startsWith(names,'referenceAnchor')), ...
    'Only the quadratic/transport reference transform is implemented.');
fresh=state.ops.rescaling;
alternative=ipm.output.restoreRuntimeReferences(fresh,transformed,state.ops.x,true);
[baseRhs,baseFlow]=ipm.evolve.flow(state.rho,state.ops,state.scale);
ops=state.ops;ops.rescaling=alternative;
[otherRhs,otherFlow]=ipm.evolve.flow(state.rho,ops,state.scale);
referenceNames=fieldnames(fresh);
referenceNames=referenceNames(startsWith(referenceNames,'reference') | ...
    startsWith(referenceNames,'adaptiveTarget') | startsWith(referenceNames,'safety') | ...
    ismember(referenceNames,{'pinX','strainTarget','peakTrackingHalfWidth'}));
freshReferences=struct();transformedReferences=struct();
for k=1:numel(referenceNames)
    name=referenceNames{k};freshReferences.(name)=fresh.(name);transformedReferences.(name)=alternative.(name);
end
report=struct('parentFile',parentFile,'childFile',childFile, ...
    'freshReferences',freshReferences,'transformedParentReferences',transformedReferences, ...
    'rhsRelativeInf',max(abs(baseRhs-otherRhs),[],'all')/max(max(abs(baseRhs),[],'all'),realmin), ...
    'cLAbsoluteDifference',abs(baseFlow.c_l-otherFlow.c_l), ...
    'cOmegaAbsoluteDifference',abs(baseFlow.c_omega-otherFlow.c_omega), ...
    'peakAbsoluteDifference',abs(baseFlow.omegaGaugeQuadraticPeakValue-otherFlow.omegaGaugeQuadraticPeakValue), ...
    'peakXAbsoluteDifference',abs(baseFlow.omegaGaugeQuadraticPeakX-otherFlow.omegaGaugeQuadraticPeakX), ...
    'pdeAdvanced',false,'referencesInstalled',false, ...
    'interpretation','Instantaneous reference-sensitivity audit at the child endpoint, not an evolved old-reference branch.');
save(outputFile,'report','-v7.3');
fprintf('FRESH_REFERENCE_AUDIT rhs=%.6g cL=%.6g cOmega=%.6g peak=%.6g\n', ...
    report.rhsRelativeInf,report.cLAbsoluteDifference,report.cOmegaAbsoluteDifference,report.peakAbsoluteDifference);
end
