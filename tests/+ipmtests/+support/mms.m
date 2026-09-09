function state = mms(x,y,t)
%IPMTESTS.SUPPORT.MMS Smooth forced IPM solution on [-1,1] x [0,1].
%   STATE = IPMTESTS.SUPPORT.MMS(X,Y,T) returns a manufactured solution of
%
%       rho_t + u dot grad(rho) = forcing,
%       -Delta psi = rho_x,       u = grad^perp psi.
%
%   X and Y may be matching arrays or coordinate vectors.  The stream
%   function vanishes on all four sides, so homogeneous Dirichlet data
%   isolate the differential and elliptic discretizations from a Green
%   boundary approximation.

validateattributes(t,{'numeric'},{'scalar','real','finite'},mfilename,'t');
if isvector(x) && isvector(y)
    [X,Y] = meshgrid(x(:)',y(:));
elseif isequal(size(x),size(y))
    X = x;
    Y = y;
else
    error('ipm:HighOrderMmsGrid', ...
        'X and Y must be coordinate vectors or matching arrays.');
end
validateattributes(X,{'numeric'},{'2d','real','finite'},mfilename,'x');
validateattributes(Y,{'numeric'},{'2d','real','finite'},mfilename,'y');

g = Y.*(1-Y);
H = pi^2*g+2;
S = sin(pi*(X+1));
C = cos(pi*(X+1));
A = 1+0.1*sin(2*pi*t);
Aprime = 0.2*pi*cos(2*pi*t);
B = 2;

psi = A*S.*g;
rho = B-(A/pi)*H.*C;
rhoT = -(Aprime/pi)*H.*C;
rhoX = A*H.*S;
rhoY = -A*pi*(1-2*Y).*C;
u1 = -A*S.*(1-2*Y);
u2 = A*pi*C.*g;
transport = -A^2*(1-2*Y).*(H.*S.^2+pi^2*g.*C.^2);
forcing = -(Aprime/pi)*H.*C+transport;

% The MMS is only quadratic in y, which the seven-point first-derivative
% stencil differentiates to roundoff.  This smooth non-polynomial companion
% field makes the standalone D_x/D_y convergence gate informative.
probePhase = 1.7*X+0.6*Y;
derivativeProbe = exp(0.4*Y).*sin(probePhase);
derivativeProbeX = 1.7*exp(0.4*Y).*cos(probePhase);
derivativeProbeY = exp(0.4*Y).* ...
    (0.4*sin(probePhase)+0.6*cos(probePhase));

state = struct('A',A,'Aprime',Aprime,'B',B,'g',g,'H',H, ...
    'S',S,'C',C,'psi',psi,'rho',rho,'rhoT',rhoT, ...
    'rhoX',rhoX,'rhoY',rhoY,'source',rhoX,'u1',u1,'u2',u2, ...
    'transport',transport,'forcing',forcing, ...
    'derivativeProbe',derivativeProbe, ...
    'derivativeProbeX',derivativeProbeX, ...
    'derivativeProbeY',derivativeProbeY);
end
