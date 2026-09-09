function value = ipm_accellab_tensor_hermite(field,x,y,Dx,Dy,queryX,queryY)
%IPM_ACCELLAB_TENSOR_HERMITE Linear C1 bicubic observer with supplied FD jets.
%   Values, first derivatives, and mixed jets are fixed linear functions of
%   FIELD. Spatial derivatives differentiate the same tensor polynomial.
x = x(:)'; y = y(:);
validateattributes(x,{'numeric'},{'finite','real','increasing'});
validateattributes(y,{'numeric'},{'finite','real','increasing'});
assert(isequal(size(field),[numel(y),numel(x)]) && ...
    isequal(size(Dx),[numel(x),numel(x)]) && isequal(size(Dy),[numel(y),numel(y)]) && ...
    isequal(size(queryX),size(queryY)) && all(isfinite(field(:))) && isreal(field), ...
    'ipm:HermiteField','Finite field, paired query arrays and maintained derivative operators are required.');
assert(all(isfinite(queryX(:))) && all(isfinite(queryY(:))) && ...
    all(queryX(:) >= x(1) & queryX(:) <= x(end)) && ...
    all(queryY(:) >= y(1) & queryY(:) <= y(end)), ...
    'ipm:HermiteDomain','Queries must lie within the native tensor grid.');
column = discretize(queryX,x); row = discretize(queryY,y);
% Vector-axis indexing preserves the axis orientation for vector queries;
% reshape explicitly so a row query cannot broadcast against a column axis.
x0 = reshape(x(column),size(queryX)); y0 = reshape(y(row),size(queryY));
hx = reshape(x(column+1)-x(column),size(queryX));
hy = reshape(y(row+1)-y(row),size(queryY));
tx = (queryX-x0)./hx; ty = (queryY-y0)./hy;
[bx,dbx] = basis(tx,hx); [by,dby] = basis(ty,hy);
fieldX = field*Dx'; fieldY = Dy*field;
fieldXY = Dy*fieldX;
v = zeros(size(queryX)); vx = v; vy = v; vxy = v;
for ix = 1:2
    for iy = 1:2
        index = sub2ind(size(field),row+iy-1,column+ix-1);
        jets = {field(index),fieldX(index),fieldY(index),fieldXY(index)};
        jx = [ix,ix+2,ix,ix+2]; jy = [iy,iy,iy+2,iy+2];
        for k = 1:4
            v = v+bx{jx(k)}.*by{jy(k)}.*jets{k};
            vx = vx+dbx{jx(k)}.*by{jy(k)}.*jets{k};
            vy = vy+bx{jx(k)}.*dby{jy(k)}.*jets{k};
            vxy = vxy+dbx{jx(k)}.*dby{jy(k)}.*jets{k};
        end
    end
end
assert(isequal(size(v),size(queryX)) && all(isfinite(v(:))) && ...
    all(isfinite(vx(:))) && all(isfinite(vy(:))), ...
    'ipm:HermiteOutput','Observer value and derivatives must retain query shape and remain finite.');
value = struct('value',v,'derivativeX',vx,'derivativeY',vy,'derivativeXY',vxy, ...
    'cellRows',row,'cellColumns',column,'fractionX',tx,'fractionY',ty, ...
    'kind','linear_tensor_cubic_hermite_with_fixed_fd_jets', ...
    'mixedJetDefinition','Dy*(field*Dx_transpose)', ...
    'spatialContinuity','C1 across native cell interfaces from shared nodal jets');
end

function [b,db] = basis(t,h)
b = {2*t.^3-3*t.^2+1,-2*t.^3+3*t.^2, ...
    h.*(t.^3-2*t.^2+t),h.*(t.^3-t.^2)};
db = {(6*t.^2-6*t)./h,(-6*t.^2+6*t)./h, ...
    3*t.^2-4*t+1,3*t.^2-2*t};
end
