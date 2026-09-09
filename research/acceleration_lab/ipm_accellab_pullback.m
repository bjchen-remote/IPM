function chart = ipm_accellab_pullback(omega,forcing,x,y,coordinates,xi,eta)
%IPM_ACCELLAB_PULLBACK Exact derivative of a specified bilinear observer.
%   The physical PDE field is untouched. The normalized local observation is
%   U=I_h[omega](a+wx*xi,wy*eta)/P. Its time derivative uses derivatives of
%   that SAME bilinear interpolant, not interpolated FD derivatives.
validateattributes(xi,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(eta,{'numeric'},{'vector','real','finite','increasing','nonnegative'});
assert(coordinates.valid,'ipm:PullbackCoordinates','Valid exact coordinates are required.');
x = x(:)'; y = y(:);
assert(isequal(size(omega),[numel(y),numel(x)]) && isequal(size(forcing),size(omega)) && ...
    all(isfinite(omega(:))) && all(isfinite(forcing(:))), ...
    'ipm:PullbackField','Finite paired native fields are required.');
[XI,ETA] = meshgrid(xi(:)',eta(:));
queryX = coordinates.peak.x+coordinates.wallCoreWidth*XI;
queryY = coordinates.verticalCoreWidth*ETA;
assert(all(queryX(:) > x(1) & queryX(:) < x(end)) && ...
    all(queryY(:) >= y(1) & queryY(:) < y(end)), ...
    'ipm:PullbackDomain','The complete observer rectangle must lie inside the native domain.');
column = discretize(queryX,x);
row = discretize(queryY,y);
dx = x(column+1)-x(column);
dy = y(row+1)-y(row);
fx = (queryX-x(column))./dx;
fy = (queryY-y(row))./dy;
shape = size(omega);
i00 = sub2ind(shape,row,column);
i10 = sub2ind(shape,row,column+1);
i01 = sub2ind(shape,row+1,column);
i11 = sub2ind(shape,row+1,column+1);
interpolated = bilinear(omega);
interpolatedForcing = bilinear(forcing);
derivativeX = ((1-fy).*(omega(i10)-omega(i00))+fy.*(omega(i11)-omega(i01)))./dx;
derivativeY = ((1-fx).*(omega(i01)-omega(i00))+fx.*(omega(i11)-omega(i10)))./dy;
queryXPrime = coordinates.translationRate+coordinates.wallCoreWidthPrime*XI;
queryYPrime = coordinates.verticalCoreWidthPrime*ETA;
P = coordinates.peak.value;
U = interpolated/P;
G = (interpolatedForcing+queryXPrime.*derivativeX+queryYPrime.*derivativeY)/P- ...
    (coordinates.peakPrime/P)*U;
fractionMargins = min(cat(3,fx,1-fx,fy,1-fy),[],3);
% Wall queries stay at Y=0, so that fixed boundary has no moving-cell event.
nonWall = queryY > y(1);
chart = struct('kind','read_only_bilinear_inner_pullback','U',U,'G',G, ...
    'xi',xi(:)','eta',eta(:),'queryX',queryX,'queryY',queryY, ...
    'cellRows',row,'cellColumns',column, ...
    'minimumMovingCellMargin',min(fractionMargins(nonWall)), ...
    'xFractionMargin',min(fx,1-fx),'yFractionMargin',min(fy,1-fy), ...
    'amplitudeDerivativeRetained',coordinates.peakPrime, ...
    'contract','G differentiates the same bilinear observation on its current cells; no PDE or C is changed.');

    function value = bilinear(field)
        value = (1-fx).*(1-fy).*field(i00)+fx.*(1-fy).*field(i10)+ ...
            (1-fx).*fy.*field(i01)+fx.*fy.*field(i11);
    end
end
