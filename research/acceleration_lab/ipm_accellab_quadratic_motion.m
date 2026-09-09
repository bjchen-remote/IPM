function motion = ipm_accellab_quadratic_motion(values,forcing,x,index)
%IPM_ACCELLAB_QUADRATIC_MOTION Exact motion of the active discrete peak.
%   aPrime = -Q[forcing]'(a)/Q[values]''(a), with the same three nodes.
validateattributes(forcing,{'numeric'},{'vector','real','finite','numel',numel(values)});
peak = ipm.evolve.quadraticPeakFunctional(values,x,index);
indices = peak.indices;
hLeft = x(index)-x(index-1);
hRight = x(index+1)-x(index);
local = forcing(indices);
left = local(1)-local(2);
right = local(3)-local(2);
curvature = (left/hLeft+right/hRight)/(hLeft+hRight);
linear = right/hRight-curvature*hRight;
forcingDerivative = 2*curvature*(peak.x-x(index))+linear;
denominator = 2*peak.curvature;
rate = -forcingDerivative/denominator;
motion = struct('translationRate',rate,'forcingQuadraticDerivative',forcingDerivative, ...
    'sourceQuadraticSecondDerivative',denominator, ...
    'phaseResidual',forcingDerivative+rate*denominator, ...
    'peakX',peak.x,'stencilIndices',indices, ...
    'activeNodeMargin',peak.activeNodeMargin, ...
    'vertexInteriorFraction',peak.vertexInteriorFraction, ...
    'contract','Exact on the active fixed three-node polynomial; stencil switches remain nonsmooth events.');
end
