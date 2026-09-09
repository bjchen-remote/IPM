function rho = bessel_k1_two_scale_initial_data(x,y,parameters)
%BESSEL_K1_TWO_SCALE_INITIAL_DATA Regularized homogeneous-wall Bessel mode.
%   Away from its smooth core this is proportional to
%
%     exp(mu*x/2) K_1(mu*r/2) sin(theta).
%
%   The main component vanishes exactly on y=0.  An optional, explicitly
%   reported signed Gaussian wall seed can be added solely to satisfy the current
%   solver's wall-rho_x tracking initialization.  Setting monitorAmplitude
%   to zero recovers the exact zero wall trace.

if nargin < 3 || isempty(parameters)
    parameters = struct();
end
parameters = initial_parameters(parameters);
if ~isequal(size(x),size(y))
    error('ipm:BesselK1InitialGrid', ...
        'The x and y coordinate arrays must have the same size.');
end
if abs(parameters.perturbationAmplitude) >= 1
    error('ipm:BesselK1PerturbationAmplitude', ...
        'Use |perturbationAmplitude| < 1 to keep the multiplier positive.');
end

centeredX = x-parameters.centerX;
regularizedRadius = sqrt(centeredX.^2+y.^2+parameters.coreRadius^2);
argument = 0.5*parameters.mu*regularizedRadius;
exponential = exp(0.5*parameters.mu*centeredX-argument);
mainComponent = exponential.*besselk(1,argument,1).* ...
    (y./regularizedRadius);
mainComponent = parameters.amplitude*mainComponent / ...
    max(mainComponent,[],'all');

horizontalBump = exp(-((centeredX- ...
    parameters.perturbationCenterX)/parameters.perturbationScaleX).^2);
scaledY = y/parameters.perturbationScaleY;
verticalBump = scaledY.*exp(0.5*(1-scaledY.^2));
mainComponent = mainComponent.*(1+parameters.perturbationAmplitude* ...
    horizontalBump.*verticalBump);

monitorComponent = -exp(-(centeredX/parameters.monitorScaleX).^2- ...
    (y/parameters.monitorScaleY).^4);
rho = mainComponent+parameters.monitorAmplitude*parameters.amplitude* ...
    monitorComponent;
end

function parameters = initial_parameters(userParameters)
parameters = struct('mu',2,'coreRadius',0.45,'centerX',0, ...
    'amplitude',1,'perturbationAmplitude',0, ...
    'perturbationCenterX',5,'perturbationScaleX',2, ...
    'perturbationScaleY',1.5,'monitorAmplitude',0, ...
    'monitorScaleX',0.8,'monitorScaleY',2);
names = fieldnames(userParameters);
for index = 1:numel(names)
    parameters.(names{index}) = userParameters.(names{index});
end
positiveNames = {'mu','coreRadius','amplitude','perturbationScaleX', ...
    'perturbationScaleY','monitorScaleX','monitorScaleY'};
for index = 1:numel(positiveNames)
    validateattributes(parameters.(positiveNames{index}),{'numeric'}, ...
        {'scalar','real','finite','positive'});
end
finiteNames = {'centerX','perturbationAmplitude', ...
    'perturbationCenterX','monitorAmplitude'};
for index = 1:numel(finiteNames)
    validateattributes(parameters.(finiteNames{index}),{'numeric'}, ...
        {'scalar','real','finite'});
end
if parameters.monitorAmplitude < 0
    error('ipm:BesselK1MonitorAmplitude', ...
        'monitorAmplitude must be nonnegative.');
end
end
