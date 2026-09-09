function report = ipm_accellab_hermite_amplitude(values,baseForcing,x,Dx,window,user)
%IPM_ACCELLAB_HERMITE_AMPLITUDE Pure instantaneous research amplitude rate.
%   Inputs are wall Omega and wall (F0)_X, with the independently assembled
%   decomposition F=F0+cOmega*rho. This never advances rho or restores P.
if nargin < 6, user = struct(); end
validateattributes(values,{'numeric'},{'vector','real','finite'});
validateattributes(baseForcing,{'numeric'},{'vector','real','finite','numel',numel(values)});
values = values(:)'; baseForcing = baseForcing(:)';
peak = ipm_accellab_hermite_peak(values,baseForcing,x,Dx,window,user);
rate = NaN; residual = NaN; normalizedResidual = NaN; correctedPeak = struct();
if peak.valid
    rate = -peak.PPrime/peak.value;
    correctedPeak = ipm_accellab_hermite_peak(values,baseForcing+rate*values,x,Dx,window,user);
    assert(correctedPeak.valid,'ipm:HermiteAmplitudeValidity','Forcing changed field-only peak validity.');
    residual = correctedPeak.PPrime;
    normalizedResidual = abs(residual)/max([abs(peak.PPrime),abs(rate*peak.value),realmin]);
end
report = struct('kind','independent_instantaneous_C1_peak_amplitude_rate', ...
    'valid',peak.valid,'basePeak',peak,'cOmega',rate,'correctedPeak',correctedPeak, ...
    'instantaneousPeakDerivative',residual,'relativeCancellationResidual',normalizedResidual, ...
    'usesTargetPeak',false,'usesFeedback',false,'pdeSteps',0,'isProductionGauge',false, ...
    'interpretation','An algebraic envelope derivative cancellation for F0+cOmega*rho; it is not a finite-step conservation or time-integration result.');
end
