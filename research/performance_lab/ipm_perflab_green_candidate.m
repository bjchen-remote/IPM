function boundary = ipm_perflab_green_candidate(source,ops,kappa)
%IPM_PERFLAB_GREEN_CANDIDATE Exact within-call distance reuse; research only.
%   DOUBLE_ODD_OMEGA integrates only the first-quadrant source and includes
%   its odd images across both axes. HALF_PLANE retains the single y image.

if nargin < 3
    kappa = 1;
end
validateattributes(kappa,{'numeric'},{'scalar','real','finite','positive'});

boundary = struct('left',zeros(ops.ny,1), ...
    'right',zeros(ops.ny,1),'bottom',zeros(1,ops.nx), ...
    'top',zeros(1,ops.nx));
if strcmp(ops.farBoundaryMode,'dirichlet_zero')
    return;
end

interiorSource = source(2:end-1,2:end-1);
sourceX = ops.X(2:end-1,2:end-1);
sourceY = ops.Y(2:end-1,2:end-1);
if strcmp(ops.symmetryMode,'double_odd_omega')
    admissible = sourceX > 0;
else
    admissible = true(size(interiorSource));
end
maskedSource = interiorSource;
maskedSource(~admissible) = 0;
sourceMaximum = max(abs(maskedSource),[],'all');
if sourceMaximum == 0
    return;
end
weights = ops.integrationWeights;
interiorWeights = weights(2:end-1,2:end-1);
weightedSource = interiorSource.*interiorWeights;
candidate = find(admissible & abs(interiorSource) >= ...
    ops.greenSourceTolerance*sourceMaximum);
[sourceX,sourceY,sourceStrength] = compress_sources(sourceX,sourceY, ...
    weightedSource,candidate,ops.greenMaxSources);

if kappa == 1 && strcmp(ops.symmetryMode,'double_odd_omega') && ops.x(1) == -ops.x(end)
    [boundary.left,boundary.right] = symmetric_sides(ops.x(end),ops.y, ...
        sourceX,sourceY,sourceStrength);
else
    boundary.left = evaluate_boundary(ops.x(1)*ones(ops.ny,1), ...
        ops.y,sourceX,sourceY,sourceStrength,ops.symmetryMode,kappa);
    boundary.right = evaluate_boundary(ops.x(end)*ones(ops.ny,1), ...
        ops.y,sourceX,sourceY,sourceStrength,ops.symmetryMode,kappa);
end
boundary.top = evaluate_boundary(ops.x', ...
    ops.y(end)*ones(ops.nx,1),sourceX,sourceY,sourceStrength, ...
    ops.symmetryMode,kappa)';
boundary.bottom(:) = 0;
end

function [sourceX,sourceY,sourceStrength] = compress_sources( ...
    sourceXGrid,sourceYGrid,weightedSource,candidate,maxSources)
% Preserve signed strength and centroid in uniform-reference grid blocks.

sourceX = sourceXGrid(candidate);
sourceY = sourceYGrid(candidate);
sourceStrength = weightedSource(candidate);
if numel(candidate) <= maxSources
    sourceX = sourceX';
    sourceY = sourceY';
    return;
end

[ny,nx] = size(weightedSource);
targetBins = floor(maxSources/2);
binsX = max(1,floor(sqrt(targetBins*nx/ny)));
binsY = max(1,floor(targetBins/binsX));
[row,column] = ind2sub([ny,nx],candidate);
binX = min(floor((column-1)*binsX/nx)+1,binsX);
binY = min(floor((row-1)*binsY/ny)+1,binsY);
bin = binY+(binX-1)*binsY;
numberOfBins = binsX*binsY;

[positiveX,positiveY,positiveStrength] = aggregate_sign( ...
    bin,sourceX,sourceY,max(sourceStrength,0),numberOfBins,1);
[negativeX,negativeY,negativeStrength] = aggregate_sign( ...
    bin,sourceX,sourceY,max(-sourceStrength,0),numberOfBins,-1);
sourceX = [positiveX;negativeX]';
sourceY = [positiveY;negativeY]';
sourceStrength = [positiveStrength;negativeStrength];
end

function [x,y,strength] = aggregate_sign(bin,sourceX,sourceY,magnitude, ...
    numberOfBins,signValue)
total = accumarray(bin,magnitude,[numberOfBins,1],@sum,0);
active = total > 0;
xMoment = accumarray(bin,magnitude.*sourceX,[numberOfBins,1],@sum,0);
yMoment = accumarray(bin,magnitude.*sourceY,[numberOfBins,1],@sum,0);
x = xMoment(active)./total(active);
y = yMoment(active)./total(active);
strength = signValue*total(active);
end

function values = evaluate_boundary( ...
        targetX,targetY,sourceX,sourceY,strength,symmetryMode,kappa)
targetX = targetX(:);
targetY = targetY(:);
values = zeros(size(targetX));
chunkSize = 128;
for first = 1:chunkSize:numel(targetX)
    last = min(first+chunkSize-1,numel(targetX));
    targetXChunk = targetX(first:last);
    targetYChunk = targetY(first:last);
    if kappa == 1
        % Keep the maintained isotropic arithmetic exactly unchanged.
        minusX = (targetXChunk-sourceX).^2;
        minusY = (targetYChunk-sourceY).^2;
        plusY = (targetYChunk+sourceY).^2;
        directSquared = minusX+minusY;
        yImageSquared = minusX+plusY;
        if strcmp(symmetryMode,'double_odd_omega')
            plusX = (targetXChunk+sourceX).^2;
            xImageSquared = plusX+minusY;
            xyImageSquared = plusX+plusY;
        end
    else
        directSquared = kappa*(targetXChunk-sourceX).^2+ ...
            (targetYChunk-sourceY).^2/kappa;
        yImageSquared = kappa*(targetXChunk-sourceX).^2+ ...
            (targetYChunk+sourceY).^2/kappa;
        if strcmp(symmetryMode,'double_odd_omega')
            xImageSquared = kappa*(targetXChunk+sourceX).^2+ ...
                (targetYChunk-sourceY).^2/kappa;
            xyImageSquared = kappa*(targetXChunk+sourceX).^2+ ...
                (targetYChunk+sourceY).^2/kappa;
        end
    end
    if strcmp(symmetryMode,'double_odd_omega')
        green = log((xImageSquared.*yImageSquared) ./ ...
            (directSquared.*xyImageSquared))/(4*pi);
    else
        green = log(yImageSquared./directSquared)/(4*pi);
    end
    values(first:last) = green*strength;
end
end

function [left,right] = symmetric_sides(H,y,sourceX,sourceY,strength)
% For exact x(1)=-x(end), squared offsets exchange bitwise under sign.
% Each side still evaluates the original reciprocal ratio, log and GEMV;
% do not replace one side by the negative of the other.
y=y(:);left=zeros(size(y));right=zeros(size(y));
chunkSize=128;
for first=1:chunkSize:numel(y)
    last=min(first+chunkSize-1,numel(y));yc=y(first:last);
    targetX=H*ones(size(yc));
    minusX=(targetX-sourceX).^2;plusX=(targetX+sourceX).^2;
    minusY=(yc-sourceY).^2;plusY=(yc+sourceY).^2;
    direct=minusX+minusY;yImage=minusX+plusY;
    xImage=plusX+minusY;xyImage=plusX+plusY;
    numerator=xImage.*yImage;denominator=direct.*xyImage;
    green=log(denominator./numerator)/(4*pi);left(first:last)=green*strength;
    green=log(numerator./denominator)/(4*pi);right(first:last)=green*strength;
end
end
