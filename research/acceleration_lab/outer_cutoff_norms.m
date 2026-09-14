function report = outer_cutoff_norms(state,radii,powers)
%OUTER_CUTOFF_NORMS Native semidiscrete R_{X tau} norms away from (1,0).
%   The coordinates are the current rescaled X,Y. A smooth half-disk cutoff
%   removes distance <= delta from (1,0), transitioning to one by 2*delta.
%   A fixed outer window is one on 0<=X<=3, 0<=Y<=1 and tapers to zero at
%   X=4 or Y=2. All finite-p norms use positive trapezoid weights. This is
%   a read-only diagnostic; candidate c_omega values are not applied.

if nargin < 2 || isempty(radii), radii = [0.05,0.1,0.2,0.4]; end
if nargin < 3 || isempty(powers), powers = [1,2,4,Inf]; end
assert(isstruct(state) && isscalar(state) && ...
    all(isfield(state,{'rho','ops','flow','rhsCache','scale'})), ...
    'ipm:OuterCutoffState','A native accepted solver state is required.');
assert(isnumeric(radii) && isvector(radii) && ...
    all(isfinite(radii)) && all(radii>0) && all(radii<0.5) && ...
    all(diff(radii)>0), ...
    'ipm:OuterCutoffRadii','Use strictly increasing radii in (0,0.5).');
assert(isnumeric(powers) && isvector(powers) && ...
    all(ismember(powers,[1,2,4,Inf])) && numel(unique(powers))==numel(powers), ...
    'ipm:OuterCutoffPowers','Allowed powers are 1, 2, 4, and Inf.');

ops = state.ops;
x = ops.x(:)'; y = ops.y(:);
R = state.rho; Rt = state.rhsCache.rhoRate;
assert(isequal(size(R),[numel(y),numel(x)]) && ...
    isequal(size(Rt),size(R)) && ...
    isequaln(state.rhsCache.rho,R) && ...
    isequaln(state.rhsCache.x,ops.x) && ...
    isequaln(state.rhsCache.y,ops.y) && ...
    state.rhsCache.remeshCount==ops.remeshCount && ...
    all(isfinite(Rt),'all'), ...
    'ipm:OuterCutoffCache','RHS cache must match the native accepted state.');
assert(x(1)<0 && x(end)>4 && y(1)==0 && y(end)>2, ...
    'ipm:OuterCutoffDomain','The fixed outer window must lie inside the grid.');

Rx = R*ops.Dx';
Rxt = Rt*ops.Dx';
cCurrent = state.flow.c_omega;
B = Rt-cCurrent*R;
Bx = B*ops.Dx';
assert(all(isfinite([cCurrent;Rx(:);Rxt(:);Bx(:)])), ...
    'ipm:OuterCutoffFinite','The mixed derivative must be finite.');

positive = find(x>=0);
assert(x(positive(1))==0, ...
    'ipm:OuterCutoffOrigin','The grid must contain the fixed X=0 wall point.');
wx = zeros(size(x));wx(positive) = trapezoid_weights(x(positive));
wy = trapezoid_weights(y);
area = wy*wx;
outerX = upper_taper(x,3,4).*(x>=0);
outerY = upper_taper(y,1,2);
outerWindow = outerY*outerX;
distance = hypot(y,x-1);

rows = struct([]);
for radius = radii(:)'
    chi = outerWindow.*core_cutoff(distance,radius);
    for domain = {'bulk','wall'}
        if strcmp(domain{1},'bulk')
            weights = area;
            mask = chi;
            derivative = Rxt;
            gradient = Rx;
            forcing = Bx;
        else
            weights = wx;
            mask = chi(1,:);
            derivative = Rxt(1,:);
            gradient = Rx(1,:);
            forcing = Bx(1,:);
        end
        assert(sum(weights.*mask,'all')>0 && nnz(mask)>10, ...
            'ipm:OuterCutoffSupport','The cutoff contains too few grid nodes.');
        for p = powers(:)'
            maximumLocation = [NaN,NaN];
            if isinf(p)
                absoluteDerivative = abs(mask.*derivative);
                [derivativeNorm,index] = max(absoluteDerivative(:));
                [iy,ix] = ind2sub(size(absoluteDerivative),index);
                maximumLocation = [x(ix),y(iy)];
                gradientNorm = max(abs(mask.*gradient),[],'all');
                pName = 'inf';
            else
                derivativeNorm = sum(weights.*abs(mask.*derivative).^p,'all')^(1/p);
                gradientNorm = sum(weights.*abs(mask.*gradient).^p,'all')^(1/p);
                pName = sprintf('%d',p);
            end
            candidateC = NaN;
            currentLogGradientNormRate = NaN;
            candidateRateResidual = NaN;
            if p==2 || p==4
                momentWeights = weights.*mask.^p;
                denominator = sum(momentWeights.*abs(gradient).^p,'all');
                numerator = sum(momentWeights.* ...
                    abs(gradient).^(p-2).*gradient.*forcing,'all');
                if isfinite(denominator) && denominator>100*realmin
                    candidateC = -numerator/denominator;
                    currentLogGradientNormRate = ...
                        (numerator+cCurrent*denominator)/denominator;
                    candidateRateResidual = numerator+candidateC*denominator;
                end
            end
            row = struct('radius',radius,'domain',domain{1},'p',pName, ...
                'effectiveArea',sum(weights.*mask,'all'), ...
                'supportNodes',nnz(mask), ...
                'mixedDerivativeNorm',derivativeNorm, ...
                'maximumMixedDerivativeLocation',maximumLocation, ...
                'gradientNorm',gradientNorm, ...
                'relativeMixedDerivativeNorm', ...
                    derivativeNorm/max(gradientNorm,realmin), ...
                'candidateOuterGradientGaugeC',candidateC, ...
                'currentLogGradientNormRate',currentLogGradientNormRate, ...
                'candidateGaugeRateResidual',candidateRateResidual);
            if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
        end
    end
end

report = struct('kind','outer_cutoff_native_mixed_derivative_v1', ...
    'tau',state.scale.canonicalTime,'step',state.step, ...
    'nodeCount',[numel(x),numel(y)],'remeshCount',ops.remeshCount, ...
    'currentCOmega',cCurrent, ...
    'derivative','R_X_tau_at_fixed_rescaled_XY_on_the_current_grid', ...
    'cutoffCenter',[1,0],'radii',radii(:)', ...
    'outerWindowPlateau',[0,3,0,1], ...
    'outerWindowSupport',[0,4,0,2], ...
    'powers',{arrayfun(@power_name,powers(:)','UniformOutput',false)}, ...
    'rows',rows,'pdeAdvanced',false, ...
    'interpretation', ...
    ['Vanishing outer R_X_tau is necessary for stationary outer gradients. ' ...
     'A convergence claim additionally needs a stable spatial/box limit ' ...
     'and control of the time integral. Candidate C rates are frozen-state ' ...
     'comparators, not a new trajectory.']);
end

function name = power_name(power)
if isinf(power),name='inf';else,name=sprintf('%d',power);end
end

function weights = trapezoid_weights(nodes)
vertical = iscolumn(nodes);
nodes = nodes(:)';
spacing = diff(nodes);
assert(numel(nodes)>=2 && all(spacing>0));
weights = [spacing(1),spacing(1:end-1)+spacing(2:end),spacing(end)]/2;
if vertical, weights = weights(:); end
end

function values = upper_taper(nodes,plateau,endPoint)
values = ones(size(nodes));
transition = nodes>plateau & nodes<endPoint;
values(transition) = cos(pi*(nodes(transition)-plateau)/ ...
    (2*(endPoint-plateau))).^2;
values(nodes>=endPoint) = 0;
end

function values = core_cutoff(distance,radius)
values = ones(size(distance));
values(distance<=radius) = 0;
transition = distance>radius & distance<2*radius;
values(transition) = sin(pi*(distance(transition)-radius)/(2*radius)).^2;
end
