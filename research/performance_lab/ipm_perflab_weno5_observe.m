function [fluxDerivative,observation] = ipm_perflab_weno5_observe( ...
    q,a,J,computationalSpacing,options)
%IPM.FIELD.WENO5FLUXDERIVATIVE Mapped nodal WENO-Z5 flux derivative.
%   DFDX = IPM.FIELD.WENO5FLUXDERIVATIVE(Q,A,J,DXI) approximates
%
%       d(A.*Q)/dx = (1/J) d(A.*Q)/dxi,   J = dx/dxi,
%
%   on a uniform computational grid. Q and A are nLines-by-nNodes arrays;
%   each row is differentiated independently. J may have the same size as
%   Q or be a scalar, a 1-by-nNodes row shared by all lines, or an
%   nLines-by-1 column constant along each line. A accepts the same
%   broadcasting forms. A row vector is the one-line case.
%
%   The numerical flux is the finite-difference WENO-Z5 reconstruction of
%
%       F+ = 0.5*(A.*Q + alpha*J.*Q),
%       F- = 0.5*(A.*Q - alpha*J.*Q),
%
%   where alpha=max(abs(A./J)) on each line. Thus the Lax-Friedrichs split
%   acts on the mapped conserved variable J.*Q, while the returned result
%   is the physical derivative. The interface flux difference is divided
%   by J at the output nodes.
%
%   DFDX = IPM.FIELD.WENO5FLUXDERIVATIVE(...,OPTIONS) accepts:
%
%     lowerBoundary       'extrapolate' (default) or 'reflect'
%     upperBoundary       'extrapolate' (default) or 'reflect'
%     extrapolationDegree 5 (default) or 4
%     epsilon             floor for the normalized WENO-Z regularizer
%                         (default 1e-12). The actual regularizer is
%                         max(epsilon,1/(nNodes-1)^2), which preserves the
%                         optimal weights at high-order critical points.
%     alpha               optional scalar or one value per line. It must
%                         bound abs(A./J) on its line.
%
%   Extrapolation constructs three ghost nodes from the nearest six
%   (degree 5) or five (degree 4) physical nodes. Reflection extends Q and
%   J evenly and A oddly about the selected endpoint. For a genuinely odd
%   reflected velocity, its value at that endpoint should be zero.
%
%   This is a spatial operator only. The complete high-order path uses
%   SSPRK(5,4); its required multi-profile sweep currently recommends public
%   solver CFL 0.5 after an empirical margin. The usable limit still depends
%   on the multidimensional operator, grid mapping, and boundary treatment.

if nargin < 5 || isempty(options)
    options = struct();
end

validateattributes(q,{'single','double'}, ...
    {'2d','real','finite','nonempty'},mfilename,'q',1);
if size(q,2) == 1
    error('ipm:WenoFdOrientation', ...
        ['Q must store grid nodes in columns. Transpose a column vector ' ...
        'before calling this function.']);
end
validateattributes(computationalSpacing,{'single','double'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'computationalSpacing',4);
if ~isstruct(options) || ~isscalar(options)
    error('ipm:WenoFdOptions', 'OPTIONS must be a scalar structure.');
end

[numberOfLines,numberOfNodes] = size(q);
a = expand_to_size(a,numberOfLines,numberOfNodes,'a',q);
J = expand_to_size(J,numberOfLines,numberOfNodes,'J',q);
if any(J(:) <= 0)
    error('ipm:WenoFdMetric', ...
        'J must be finite and strictly positive at every node.');
end

lowerBoundary = boundary_option(options,'lowerBoundary','extrapolate');
upperBoundary = boundary_option(options,'upperBoundary','extrapolate');
extrapolationDegree = scalar_option( ...
    options,'extrapolationDegree',5);
if ~ismember(extrapolationDegree,[4,5])
    error('ipm:WenoFdExtrapolationDegree', ...
        'OPTIONS.extrapolationDegree must be 4 or 5.');
end
epsilonFloor = scalar_option(options,'epsilon',1e-12);
if ~isfinite(epsilonFloor) || epsilonFloor <= 0
    error('ipm:WenoFdEpsilon', ...
        'OPTIONS.epsilon must be a finite positive scalar.');
end

needsExtrapolation = strcmp(lowerBoundary,'extrapolate') || ...
    strcmp(upperBoundary,'extrapolate');
minimumNodes = 5;
if needsExtrapolation
    minimumNodes = max(minimumNodes,extrapolationDegree+1);
end
if numberOfNodes < minimumNodes
    error('ipm:WenoFdGridSize', ...
        ['At least %d nodes are required for the selected WENO5 ' ...
        'boundary extension.'],minimumNodes);
end
normalizedEpsilon = max(epsilonFloor,1/(numberOfNodes-1)^2);

minimumAlpha = max(abs(a./J),[],2);
alpha = alpha_option(options,minimumAlpha,numberOfLines,q);

extendedQ = extend_lines(q,lowerBoundary,upperBoundary, ...
    extrapolationDegree,1);
extendedA = extend_lines(a,lowerBoundary,upperBoundary, ...
    extrapolationDegree,-1);
extendedJ = extend_lines(J,lowerBoundary,upperBoundary, ...
    extrapolationDegree,1);

mappedState = extendedJ.*extendedQ;
physicalFlux = extendedA.*extendedQ;
positiveFlux = 0.5*(physicalFlux+alpha.*mappedState);
negativeFlux = 0.5*(physicalFlux-alpha.*mappedState);
physicalNodes = 4:numberOfNodes+3;
positiveScale = max(abs(positiveFlux(:,physicalNodes)),[],2);
negativeScale = max(abs(negativeFlux(:,physicalNodes)),[],2);
commonScale = max(positiveScale,negativeScale);
splitScaleFloor = sqrt(eps(class(q)))*commonScale;
positiveScale = max(positiveScale,splitScaleFloor);
negativeScale = max(negativeScale,splitScaleFloor);

faceFlux = zeros(numberOfLines,numberOfNodes+1,'like',q);
for face = 0:numberOfNodes
    positiveStencil = positiveFlux(:,face+(1:5));
    negativeStencil = negativeFlux(:,face+(6:-1:2));
    faceFlux(:,face+1) = weno_z_left( ...
        positiveStencil,positiveScale,normalizedEpsilon) + ...
        weno_z_left(negativeStencil,negativeScale,normalizedEpsilon);
end

computationalDerivative = diff(faceFlux,1,2)/computationalSpacing;
fluxDerivative = computationalDerivative./J;
% Research observation only: the original derivative above is unchanged.
lines=options.observeLines; faces=options.observeFaces;
positiveWeights=zeros(numel(lines),numel(faces),3);
negativeWeights=zeros(size(positiveWeights));
for k=1:numel(faces)
    face=faces(k);
    [~,wp]=weno_z_left(positiveFlux(lines,face+(1:5)),positiveScale(lines),normalizedEpsilon);
    [~,wm]=weno_z_left(negativeFlux(lines,face+(6:-1:2)),negativeScale(lines),normalizedEpsilon);
    positiveWeights(:,k,:)=reshape(wp,numel(lines),1,3);
    negativeWeights(:,k,:)=reshape(wm,numel(lines),1,3);
end
[~,maximumSpeedNodes]=max(abs(a./J),[],2);
observation=struct('alpha',alpha,'maximumSpeedNodes',maximumSpeedNodes, ...
    'positiveScale',positiveScale,'negativeScale',negativeScale, ...
    'normalizedEpsilon',normalizedEpsilon,'lines',lines,'faces',faces, ...
    'positiveWeights',positiveWeights,'negativeWeights',negativeWeights);
end

function values = expand_to_size(values,numberOfLines,numberOfNodes, ...
    argumentName,prototype)
validateattributes(values,{'single','double'}, ...
    {'2d','real','finite','nonempty'},mfilename,argumentName);
if isscalar(values)
    values = repmat(values,numberOfLines,numberOfNodes);
elseif isequal(size(values),[numberOfLines,numberOfNodes])
    % Already in the required form.
elseif isequal(size(values),[numberOfLines,1])
    values = repmat(values,1,numberOfNodes);
elseif isvector(values) && numel(values) == numberOfNodes
    values = repmat(reshape(values,1,numberOfNodes),numberOfLines,1);
else
    error('ipm:WenoFdSize', ...
        ['%s must be scalar, the size of Q, a 1-by-nNodes row, or an ' ...
        'nLines-by-1 column.'],upper(argumentName));
end
values = cast(values,'like',prototype);
end

function value = scalar_option(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('ipm:WenoFdOptions', ...
        'OPTIONS.%s must be a real numeric scalar.',name);
end
value = double(value);
end

function value = boundary_option(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('ipm:WenoFdBoundary', ...
        'OPTIONS.%s must be ''extrapolate'' or ''reflect''.',name);
end
value = lower(char(value));
if ~ismember(value,{'extrapolate','reflect'})
    error('ipm:WenoFdBoundary', ...
        'OPTIONS.%s must be ''extrapolate'' or ''reflect''.',name);
end
end

function alpha = alpha_option(options,minimumAlpha,numberOfLines,prototype)
if ~isfield(options,'alpha') || isempty(options.alpha)
    alpha = minimumAlpha;
    return
end
alpha = options.alpha;
validateattributes(alpha,{'single','double'}, ...
    {'vector','real','finite','nonnegative'},mfilename,'OPTIONS.alpha');
if isscalar(alpha)
    alpha = repmat(alpha,numberOfLines,1);
elseif numel(alpha) == numberOfLines
    alpha = reshape(alpha,numberOfLines,1);
else
    error('ipm:WenoFdAlphaSize', ...
        'OPTIONS.alpha must be scalar or have one value per line.');
end
alpha = cast(alpha,'like',prototype);
tolerance = 100*eps(class(prototype)).*max(1,minimumAlpha);
if any(alpha < minimumAlpha-tolerance)
    error('ipm:WenoFdAlphaBound', ...
        'OPTIONS.alpha must bound max(abs(A./J)) on every line.');
end
end

function extended = extend_lines(values,lowerBoundary,upperBoundary, ...
    degree,reflectionParity)
[numberOfLines,numberOfNodes] = size(values);
extended = zeros(numberOfLines,numberOfNodes+6,'like',values);
extended(:,4:numberOfNodes+3) = values;

if strcmp(lowerBoundary,'reflect')
    extended(:,1:3) = reflectionParity*values(:,[4,3,2]);
else
    nodes = 0:degree;
    targets = -3:-1;
    weights = lagrange_weights(nodes,targets);
    extended(:,1:3) = values(:,1:degree+1)*weights.';
end

if strcmp(upperBoundary,'reflect')
    extended(:,numberOfNodes+(4:6)) = reflectionParity* ...
        values(:,numberOfNodes-[1,2,3]);
else
    nodes = -degree:0;
    targets = 1:3;
    weights = lagrange_weights(nodes,targets);
    extended(:,numberOfNodes+(4:6)) = ...
        values(:,numberOfNodes-degree:numberOfNodes)*weights.';
end
end

function weights = lagrange_weights(nodes,targets)
numberOfNodes = numel(nodes);
numberOfTargets = numel(targets);
weights = ones(numberOfTargets,numberOfNodes);
for node = 1:numberOfNodes
    others = [1:node-1,node+1:numberOfNodes];
    weights(:,node) = prod((targets(:)-nodes(others)) ./ ...
        (nodes(node)-nodes(others)),2);
end
end

function [reconstructed,weights] = weno_z_left(values,lineScale,epsilon)
safeScale = max(lineScale,sqrt(realmin(class(values))));
values = values./safeScale;

v1 = values(:,1);
v2 = values(:,2);
v3 = values(:,3);
v4 = values(:,4);
v5 = values(:,5);

candidate0 = (1/3)*v1-(7/6)*v2+(11/6)*v3;
candidate1 = -(1/6)*v2+(5/6)*v3+(1/3)*v4;
candidate2 = (1/3)*v3+(5/6)*v4-(1/6)*v5;

beta0 = (13/12)*(v1-2*v2+v3).^2 + ...
    (1/4)*(v1-4*v2+3*v3).^2;
beta1 = (13/12)*(v2-2*v3+v4).^2 + ...
    (1/4)*(v2-v4).^2;
beta2 = (13/12)*(v3-2*v4+v5).^2 + ...
    (1/4)*(3*v3-4*v4+v5).^2;
tau5 = abs(beta0-beta2);

alpha0 = 0.1*(1+(tau5./(beta0+epsilon)).^2);
alpha1 = 0.6*(1+(tau5./(beta1+epsilon)).^2);
alpha2 = 0.3*(1+(tau5./(beta2+epsilon)).^2);
alphaSum = alpha0+alpha1+alpha2;
reconstructed = safeScale.*( ...
    alpha0.*candidate0+alpha1.*candidate1+alpha2.*candidate2) ./ ...
    alphaSum;
if nargout>1,weights=[alpha0./alphaSum,alpha1./alphaSum,alpha2./alphaSum];end
end
