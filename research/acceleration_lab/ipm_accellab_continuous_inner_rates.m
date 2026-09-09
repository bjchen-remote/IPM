function report = ipm_accellab_continuous_inner_rates(omega,forcing,x,y,Dx,Dy,window,user)
%IPM_ACCELLAB_CONTINUOUS_INNER_RATES True derivatives of the same C1 geometry.
% FORCING is the genuine Omega_tau, e.g. maintained rhs*Dx'. PPrime is
% measured, never set to zero. All rates are read-only coordinate rates.
if nargin < 8, user = struct(); end
opts = struct('level',.9,'minimumNormalizedRootSlope',1e-10, ...
    'maximumRelativeRootResidual',1e-11,'peakOptions',struct());
names = fieldnames(user); assert(all(ismember(names,fieldnames(opts))));
for k = 1:numel(names), opts.(names{k}) = user.(names{k}); end
x = x(:)'; y = y(:);
assert(isequal(size(omega),[numel(y),numel(x)]) && isequal(size(forcing),size(omega)) && y(1) == 0);
assert(opts.level > 0 && opts.level < 1);
peak = ipm_accellab_hermite_peak(omega(1,:),forcing(1,:),x,Dx,window,opts.peakOptions);
report = struct('kind','true_PDE_derivatives_of_continuous_C1_inner_geometry', ...
    'valid',false,'geometryValid',false,'invalidReasons',{{}},'peak',peak,'options',opts, ...
    'peakPrime',peak.PPrime,'isModifiedPde',false,'pdeSteps',0);
if ~peak.valid || ~peak.smoothUniquePeak
    report.invalidReasons = {'invalid_or_nonunique_continuous_peak'}; return;
end
P = peak.value; a = peak.x; PPrime = peak.PPrime;
wall = omega(1,:); wallJets = wall*Dx'; wallForcing = forcing(1,:); forcingJets = wallForcing*Dx';
wallRoots = ipm_accellab_hermite_level_roots(wall,wallJets,x,opts.level*P,window);
left = find(wallRoots.locations < a,1,'last'); right = find(wallRoots.locations > a,1,'first');
report.wallLevelRoots = wallRoots;
if isempty(left) || isempty(right)
    report.invalidReasons = {'missing_connected_wall_crossings'}; return;
end
indices = [left,right]; crosses = wallRoots.locations(indices); slopes = wallRoots.slopes(indices);
flat = wallRoots.flatCells;
if any(x(flat) <= crosses(2) & x(flat+1) >= crosses(1))
    report.invalidReasons = {'threshold_flat_cell_intersects_connected_wall_component'}; return;
end
if slopes(1) <= 0 || slopes(2) >= 0 || ...
        any(abs(wallRoots.normalizedSlopes(indices)) <= opts.minimumNormalizedRootSlope) || ...
        any(wallRoots.relativeResiduals(indices) > opts.maximumRelativeRootResidual)
    report.invalidReasons = {'unresolved_or_nontransverse_connected_wall_crossings'}; return;
end
trace = ipm_accellab_tensor_hermite(omega,x,y,Dx,Dy,a+zeros(size(y)),y);
verticalRoots = ipm_accellab_hermite_level_roots(trace.value,trace.derivativeY,y,opts.level*P,[0,y(end)]);
report.verticalLevelRoots = verticalRoots;
j = find(verticalRoots.locations > 0,1,'first');
flat = verticalRoots.flatCells;
if ~isempty(j) && any(y(flat) <= verticalRoots.locations(j) & y(flat+1) >= 0)
    report.invalidReasons = {'threshold_flat_cell_intersects_connected_vertical_component'}; return;
end
if isempty(j) || verticalRoots.slopes(j) >= 0 || ...
        abs(verticalRoots.normalizedSlopes(j)) <= opts.minimumNormalizedRootSlope || ...
        verticalRoots.relativeResiduals(j) > opts.maximumRelativeRootResidual
    report.invalidReasons = {'missing_or_unresolved_first_vertical_downcrossing'}; return;
end
wx = crosses(2)-crosses(1); wy = verticalRoots.locations(j);
report.geometryValid = true; report.wallCoreWidth = wx; report.verticalCoreWidth = wy;
report.wallCrossings = crosses; report.verticalCrossing = wy;
report.traceWallMinusPeak = trace.value(1)-P;
atPeak = ipm_accellab_hermite_line(wall,wallJets,x,a,peak.selectedCell);
forceAtPeak = ipm_accellab_hermite_line(wallForcing,forcingJets,x,a,peak.selectedCell);
[nodeDistance,node] = min(abs(x-a));
interface = nodeDistance <= peak.activePositionTolerance;
phase = struct('atOrNumericallyNearNativeInterface',interface,'nearestNode',node, ...
    'distanceToNearestNode',nodeDistance,'interfaceTolerance',peak.activePositionTolerance, ...
    'selectedCurvature',atPeak.secondDerivative,'selectedForcingDerivative',forceAtPeak.derivative, ...
    'selectedPeakStationarityResidual',atPeak.derivative, ...
    'forcingValueMinusPeakPrime',forceAtPeak.value-PPrime, ...
    'oneSidedNodeCurvatures',[],'oneSidedNodePhaseRates',[]);
if interface
    if node > 1 && node < numel(x)
        ql = ipm_accellab_hermite_line(wall,wallJets,x,x(node),node-1);
        qr = ipm_accellab_hermite_line(wall,wallJets,x,x(node),node);
        qf = ipm_accellab_hermite_line(wallForcing,forcingJets,x,x(node),node);
        phase.oneSidedNodeCurvatures = [ql.secondDerivative,qr.secondDerivative];
        phase.oneSidedNodePhaseRates = -qf.derivative./phase.oneSidedNodeCurvatures;
    end
    report.phase = phase;
    report.invalidReasons = {'peak_phase_derivative_at_native_interface_requires_one_sided_event_analysis'};
    return;
end
aPrime = -forceAtPeak.derivative/atPeak.secondDerivative;
phase.translationRate = aPrime;
phase.stationarityDerivativeResidual = forceAtPeak.derivative+aPrime*atPeak.secondDerivative;
report.phase = phase;
wallValues = ipm_accellab_hermite_line(wallForcing,forcingJets,x,crosses,wallRoots.cells(indices));
wallPrimes = (opts.level*PPrime-wallValues.value)./slopes;
wxPrime = wallPrimes(2)-wallPrimes(1);
vertical = ipm_accellab_tensor_hermite(omega,x,y,Dx,Dy,a,wy);
verticalForcing = ipm_accellab_tensor_hermite(forcing,x,y,Dx,Dy,a,wy);
movingTraceForcing = verticalForcing.value+aPrime*vertical.derivativeX;
wyPrime = (opts.level*PPrime-movingTraceForcing)/vertical.derivativeY;
report.valid = all(isfinite([aPrime,wxPrime,wyPrime])) && wx > 0 && wy > 0;
report.translationRate = aPrime; report.wallCoreWidthPrime = wxPrime; report.verticalCoreWidthPrime = wyPrime;
report.logScaleXRate = wxPrime/wx; report.logScaleYRate = wyPrime/wy;
report.wallCrossingPrimes = wallPrimes; report.wallCrossingForcing = wallValues.value;
report.wallCrossingSlopes = slopes; report.verticalCrossingSlope = vertical.derivativeY;
report.verticalFixedXForcing = verticalForcing.value;
report.verticalTranslationForcing = aPrime*vertical.derivativeX;
report.verticalMovingTraceForcing = movingTraceForcing;
report.verticalTensorMinusLevel = vertical.value-opts.level*P;
report.verticalTensorMinusRootSlope = vertical.derivativeY-verticalRoots.slopes(j);
report.wallImplicitDerivativeResiduals = slopes.*wallPrimes+wallValues.value-opts.level*PPrime;
report.verticalImplicitDerivativeResidual = vertical.derivativeY*wyPrime+movingTraceForcing-opts.level*PPrime;
report.signature = struct('peakCell',peak.selectedCell,'wallRootCells',wallRoots.cells(indices), ...
    'verticalRootCell',verticalRoots.cells(j),'level',opts.level);
report.interpretation = 'Read-only derivatives on the resolved connected-root/peak-cell branch. PPrime is genuine; peak interface phase derivatives are explicitly withheld. This is not a changed PDE or a singularity proof.';
end
