function report = verify_ipm_profile_fit()
%VERIFY_IPM_PROFILE_FIT Synthetic accuracy check for polar profile fitting.

% Arbitrary synthetic exponent: this is a fitter regression, not a profile claim.
alphaExact = 1/2;
fanScale = 22*pi/180;
coreScale = 0.05;
sx = linspace(-1,1,1001);
sy = linspace(0,1,401)';
x = 100*sinh(4*sx)/sinh(4);
y = 100*sinh(4*sy)/sinh(4);
[X,Y] = meshgrid(x,y);
radius = hypot(X,Y);
theta = atan2(Y,X);
angular = exp(-(theta/fanScale).^8);
omega = sign(X).*radius./(coreScale+radius).^(1+alphaExact).*angular;
omega(:,X(1,:) == 0) = 0;
source = struct('x',x,'y',y,'omega',omega);
options = struct('radialWindow',[8,50],'makePlot',false, ...
    'minimumFitPoints',24);
fit = ipm_fit_polar_profile(source,options);

threshold = fit.polar.fanThresholds;
fanExact = fanScale*(-log(threshold)).^(1/8)*180/pi;
alphaError = abs(fit.wall.alpha-alphaExact);
fanError = max(abs(fit.polar.medianFanAngleDegrees-fanExact));
assert(alphaError < 0.025, ...
    'Wall exponent recovery failed: error %.3g.',alphaError);
assert(fit.wall.rSquared > 0.999, ...
    'Synthetic wall power fit has insufficient R^2.');
assert(fanError < 1.5, ...
    'Synthetic fan-angle recovery failed: error %.3g degrees.',fanError);

report = struct('alphaExact',alphaExact, ...
    'alphaFit',fit.wall.alpha,'alphaError',alphaError, ...
    'fanExactDegrees',fanExact, ...
    'fanFitDegrees',fit.polar.medianFanAngleDegrees, ...
    'fanMaximumErrorDegrees',fanError, ...
    'separabilityError',fit.polar.separabilityError);
fprintf(['Polar profile fit: alpha %.6f (error %.3g), ', ...
    'fan error %.3g deg, collapse %.3g.\n'],fit.wall.alpha, ...
    alphaError,fanError,fit.polar.separabilityError);
end
