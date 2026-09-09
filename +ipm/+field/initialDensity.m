function rho = initialDensity(ops,physics)
%IPM.FIELD.INITIALDENSITY Physical density data for IPM.

if isa(physics.initialCondition,'function_handle')
    rho = physics.initialCondition(ops.X,ops.Y);
    if ~isequal(size(rho), [ops.ny,ops.nx])
        error('ipm:InitialSize', 'Initial-condition function returned the wrong size.');
    end
    return;
end

switch lower(string(physics.initialCondition))
    case "unstable_front"
        envelope = exp(-(ops.X/physics.frontWidth).^4);
        front = physics.frontHeight + physics.frontAmplitude*envelope.* ...
            cos(pi*ops.X/physics.frontWidth);
        % Larger rho lies above smaller rho: the Rayleigh-Taylor unstable
        % orientation for gravity in the negative y direction.
        rho = 0.5*(1+tanh((ops.Y-front)/physics.frontThickness));
    case "smooth_blob"
        rho = exp(-((ops.X/1.1).^2+((ops.Y-1.4)/0.7).^2));
    case "degenerate"
        k = physics.degeneratePower;
        if strcmp(ops.symmetryMode,'double_odd_omega')
            % The upper trace is even in x_1. Its extension across x_2=0
            % is odd and has the explicitly accepted jump at that axis.
            absoluteX = abs(ops.X);
            rho = absoluteX.^k ./ ...
                (1+absoluteX.^k+ops.Y.^k);
        else
            rho = zeros(ops.ny,ops.nx);
            positiveX = ops.X > 0;
            rho(positiveX) = ops.X(positiveX).^k ./ ...
                (1+ops.X(positiveX).^k+ops.Y(positiveX).^k);
        end
    case "degenerate_primitive"
        k = physics.degeneratePower;
        % Analytic primitive of x^(k-1)/(1+x^k+y^k), anchored by
        % rho(0,y)=0.
        if strcmp(ops.symmetryMode,'double_odd_omega')
            absoluteX = abs(ops.X);
            rho = log1p(absoluteX.^k./(1+ops.Y.^k))/k;
        else
            rho = zeros(ops.ny,ops.nx);
            positiveX = ops.X > 0;
            rho(positiveX) = log1p(ops.X(positiveX).^k ./ ...
                (1+ops.Y(positiveX).^k))/k;
        end
    case "ccf_heavy_tail"
        % A resolved approximation to the CCF boundary profile
        %   d_x rho = 1_{x>1}/sqrt(x-1).
        % For z=(|x|-a)_+, the primitive below has
        %   d_|x| rho = A*z/(z^2+epsilon^2)^(3/4),
        % which is finite initially and converges pointwise to the target
        % tail away from its onset as epsilon tends to zero.  The algebraic
        % wall-normal envelope keeps the source integrable in y without
        % changing its wall trace.
        if strcmp(ops.symmetryMode,'double_odd_omega')
            tailCoordinate = abs(ops.X);
        else
            tailCoordinate = ops.X;
        end
        z = max(tailCoordinate-physics.ccfTailOnset,0);
        unboundedPrimitive = 2*physics.ccfTailAmplitude * ...
            ((z.^2+physics.ccfTailRegularization^2).^(1/4) - ...
            sqrt(physics.ccfTailRegularization));
        if isfinite(physics.ccfTailCutoffScale)
            saturation = 2*physics.ccfTailAmplitude* ...
                sqrt(physics.ccfTailCutoffScale);
            radialPrimitive = saturation*tanh(unboundedPrimitive/saturation);
        else
            radialPrimitive = unboundedPrimitive;
        end
        verticalEnvelope = 1 ./ ...
            (1+(ops.Y/physics.ccfTailVerticalScale).^ ...
            physics.ccfTailVerticalPower);
        rho = radialPrimitive.*verticalEnvelope;
    otherwise
        error('ipm:UnknownInitialCondition', ...
            'Unknown initialCondition "%s".',physics.initialCondition);
end
end
