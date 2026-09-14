function result=test_outer_c_continuous_vertical_core()
%TEST_OUTER_C_CONTINUOUS_VERTICAL_CORE New C must keep Y shadow observable.
x=linspace(-4,4,129);
y=linspace(0,2,65);
[X,Y]=meshgrid(x,y);
source=exp(-((X-1)/0.22).^2).*exp(-(Y/0.35).^2);
flow=struct('source',source, ...
    'omegaGaugeQuadraticPeakX',NaN,'trackedPeakX',1);
levels=[0.25,0.5,0.75,0.9];
observed=ipm.remesh.continuousVerticalCore(flow,x,y,levels);
assert(abs(observed.peakX-1)<1e-12 && ...
    isfinite(observed.coreCells) && observed.coreCells>1 && ...
    isfinite(observed.coreWidth) && observed.coreWidth>0);
flow.omegaGaugeQuadraticPeakX=1.0625;
observedOldGauge=ipm.remesh.continuousVerticalCore(flow,x,y,levels);
assert(abs(observedOldGauge.peakX-1.0625)<1e-12);
result=struct('passed',true, ...
    'outerCGaugeCoreCells',observed.coreCells, ...
    'oldGaugePeakX',observedOldGauge.peakX);
end
