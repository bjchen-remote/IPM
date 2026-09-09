function rho = bessel_two_scale_initial_data(x,y,parameters)
%BESSEL_TWO_SCALE_INITIAL_DATA Smooth right-focused Bessel-like density.
%   RHO = BESSEL_TWO_SCALE_INITIAL_DATA(X,Y,P) evaluates the K_0 member
%
%     exp(mu*(x-x0)/2) K_0(mu*sqrt((x-x0)^2+(y+h)^2)/2)
%
%   with its pole at (x0,-h), below the computational half-plane.  The
%   exponentially scaled Bessel function is used so the evaluation remains
%   safe when mu*r is large.  A bounded multiplicative bump, which vanishes
%   on y=0, supplies an optional controlled perturbation.

if nargin < 3 || isempty(parameters)
    parameters = struct();
end
parameters = initial_parameters(parameters);

if ~isequal(size(x),size(y))
    error('ipm:BesselInitialGrid', ...
        'The x and y coordinate arrays must have the same size.');
end
if abs(parameters.perturbationAmplitude) >= 1
    error('ipm:BesselPerturbationAmplitude', ...
        'Use |perturbationAmplitude| < 1 to keep the multiplier positive.');
end

centeredX = x-parameters.centerX;
radius = hypot(centeredX,y+parameters.poleDepth);
argument = 0.5*parameters.mu*radius;

% MATLAB's third BESSELK input returns exp(argument)*K_0(argument).
% Since radius >= centeredX, the remaining exponential never overflows.
exponential = exp(0.5*parameters.mu*centeredX-argument);
base = exponential.*besselk(0,argument,1);
base = parameters.amplitude*base/max(base,[],'all');

horizontalBump = exp(-((centeredX- ...
    parameters.perturbationCenterX)/parameters.perturbationScaleX).^2);
scaledY = y/parameters.perturbationScaleY;
verticalBump = scaledY.*exp(0.5*(1-scaledY.^2));
perturbation = horizontalBump.*verticalBump;
rho = base.*(1+parameters.perturbationAmplitude*perturbation);
end

function parameters = initial_parameters(userParameters)
parameters = struct('mu',2,'poleDepth',0.6,'centerX',0, ...
    'amplitude',1,'perturbationAmplitude',0, ...
    'perturbationCenterX',5,'perturbationScaleX',2, ...
    'perturbationScaleY',1.5);
names = fieldnames(userParameters);
for index = 1:numel(names)
    parameters.(names{index}) = userParameters.(names{index});
end
validateattributes(parameters.mu,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(parameters.poleDepth,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(parameters.centerX,{'numeric'}, ...
    {'scalar','real','finite'});
validateattributes(parameters.amplitude,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(parameters.perturbationAmplitude,{'numeric'}, ...
    {'scalar','real','finite'});
validateattributes(parameters.perturbationCenterX,{'numeric'}, ...
    {'scalar','real','finite'});
validateattributes(parameters.perturbationScaleX,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(parameters.perturbationScaleY,{'numeric'}, ...
    {'scalar','real','finite','positive'});
end
