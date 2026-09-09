function evaluation = ipm_accellab_peak_rhs(rho,z,ops,mode,nativeEvaluation)
%IPM_ACCELLAB_PEAK_RHS Independent tiny C1 rule or unchanged native RHS.
% z = [logCl; logComega; physicalTime; Xshift; canonicalTime].
assert(ops.nx <= 65 && ops.ny <= 33 && strcmp(ops.dynamicScaleGeometry,'isotropic') && ...
    ops.rescaling.enabled && ...
    strcmp(ops.symmetryMode,'double_odd_omega') && strcmp(ops.transportBoundaryMode,'open') && ...
    strcmp(ops.rescaling.lengthGauge,'transport_anchor') && ops.rescaling.transportAnchorX == 1, ...
    'ipm:PeakResearchScope','This independent RHS is limited to the registered tiny symmetric anchor case.');
assert(any(strcmp(mode,{'native_quadratic','research_hermite_C1'})));
if strcmp(mode,'native_quadratic')
    if nargin < 5
        [rhs,flow] = ipm.evolve.flow(rho,ops);
    else
        assert(isequaln(nativeEvaluation.rho,rho) && isequaln(nativeEvaluation.z,z), ...
            'ipm:PeakNativeObserverCache','Native observer must use the exact terminal field and scale.');
        rhs = nativeEvaluation.rhs; flow = nativeEvaluation.flow;
    end
    omega = flow.source; cl = flow.c_l; cw = flow.c_omega; cr = flow.c_r;
    poissonResidual = flow.poissonResidual;
    transportU1 = flow.transportU1; transportU2 = flow.transportU2;
    ratePeak = struct();
else
    [u1,u2,~,omega,poissonResidual] = ipm.field.velocity(rho,ops);
    anchor = ops.rescaling.transportAnchorX;
    cl = -interp1(ops.x,u1(1,:),anchor,'linear')/anchor;
    cr = 0;
    [baseRhs,~,transportU1,transportU2] = ipm.evolve.assembleRhs(rho,u1,u2,cl,cl,0,cr,ops,ops.transportBoundaryMode);
    baseFX = baseRhs(1,:)*ops.Dx';
    ratePeak = ipm_accellab_hermite_peak(omega(1,:),baseFX,ops.x,ops.Dx,[0,ops.x(end)]);
    if ~ratePeak.valid
        error('ipm:PeakResearchInvalid','Invalid C1 rate peak: %s',strjoin(ratePeak.invalidReasons,', '));
    end
    cw = -ratePeak.PPrime/ratePeak.value;
    rhs = baseRhs+cw*rho;
end
wallFX = rhs(1,:)*ops.Dx';
peak = ipm_accellab_hermite_peak(omega(1,:),wallFX,ops.x,ops.Dx,[0,ops.x(end)]);
quadValue = NaN; quadPrime = NaN; quadCenter = NaN; quadValid = true; quadError = '';
try
    [~,j] = max(omega(1,:));
    q = ipm.evolve.quadraticPeakFunctional(omega(1,:),ops.x,j);
    quadValue = q.value; quadPrime = sum(q.weights.*wallFX(q.indices)); quadCenter = ops.x(j);
catch exception
    quadValid = false; quadError = exception.message;
end
diagnostic = struct('canonicalTime',z(5),'physicalTime',z(3), ...
    'logC_l',z(1),'logC_omega',z(2),'X_shift',z(4), ...
    'c_l',cl,'c_omega',cw,'c_r',cr,'H',peak.value,'HPrime',peak.PPrime, ...
    'HLocation',peak.x,'HCell',peak.selectedCell,'HFraction',peak.selectedFraction, ...
    'HSecondDerivative',peak.selectedSecondDerivative, ...
    'HOneSidedSecondDerivatives',peak.activeOneSidedSecondDerivatives, ...
    'HStationarityResidual',peak.selectedStationarityResidual, ...
    'HActiveCount',peak.isolatedPeakCount,'HValid',peak.valid,'HInvalidReasons',{peak.invalidReasons}, ...
    'HActiveValueSpread',peak.activeValueSpread,'HActiveValueTolerance',peak.activeValueTolerance, ...
    'quadraticP',quadValue,'quadraticPPrime',quadPrime,'quadraticCenter',quadCenter, ...
    'quadraticValid',quadValid,'quadraticError',quadError,'poissonResidual',poissonResidual, ...
    'rhoInfinity',max(abs(rho),[],'all'),'rhsInfinity',max(abs(rhs),[],'all'), ...
    'transportRate',max(abs(transportU1)./ops.hx+abs(transportU2)./ops.hy,[],'all'));
if strcmp(mode,'research_hermite_C1')
    assert(peak.valid && isequal(peak.value,ratePeak.value), ...
        'ipm:PeakResearchObserver','The true full-RHS observer changed the field-only maximum.');
end
scaleRate = [cl;cw;exp(z(2)-z(1));cl*z(4)+cr;1];
assert(all(isfinite(rhs),'all') && all(isfinite(scaleRate)), ...
    'ipm:PeakResearchNonfinite','Nonfinite research field or scale RHS.');
evaluation = struct('rho',rho,'z',z,'rhoRate',rhs,'scaleRate',scaleRate, ...
    'diagnostic',diagnostic,'mode',mode);
end
