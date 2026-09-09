function [cache,info] = ipm_perflab_row_factor(A,method)
%IPM_PERFLAB_ROW_FACTOR Experimental left row equilibration for sparse LU.
%   This changes only the representation of A*x=b to D*A*x=D*b. It does
%   not replace the original A used for physical residual diagnostics.
%   METHOD is none, row_max, or power2_row_max. No warning is disabled.

if nargin < 2
    method = 'power2_row_max';
end
assert(issparse(A) && isreal(A) && isa(A,'double') && ...
    size(A,1) == size(A,2) && all(isfinite(nonzeros(A))), ...
    'Pass a finite real square sparse double matrix.');
timer = tic;
rowMaximum = full(max(abs(A),[],2));
assert(all(rowMaximum > 0),'Row scaling cannot repair an empty row.');
switch method
    case 'none'
        rowScale = ones(size(rowMaximum));
        scaledA = A;
    case 'row_max'
        rowScale = 1./rowMaximum;
        scaledA = spdiags(rowScale,0,size(A,1),size(A,1))*A;
    case 'power2_row_max'
        [~,exponent] = log2(rowMaximum);
        rowScale = pow2(-exponent);
        scaledA = spdiags(rowScale,0,size(A,1),size(A,1))*A;
    otherwise
        error('ipm:RowEquilibrationMethod','Unknown row-equilibration method.');
end
assert(all(isfinite(rowScale)) && all(rowScale > 0) && ...
    all(isfinite(nonzeros(scaledA))) && nnz(scaledA) == nnz(A), ...
    'ipm:RowEquilibrationRange','Scaling overflowed or lost nonzero entries.');
preparationSeconds = toc(timer);
timer = tic;
factor = decomposition(scaledA,'lu');
factorSeconds = toc(timer);
scaledRowMaximum = full(max(abs(scaledA),[],2));
cache = struct('factor',factor,'rowScale',rowScale,'method',method, ...
    'matrixSize',size(A));
info = struct('method',method,'matrixSize',size(A),'nonzeros',nnz(A), ...
    'rowMaximumMinimum',min(rowMaximum),'rowMaximumMaximum',max(rowMaximum), ...
    'rowScaleOrders',log10(max(rowMaximum))-log10(min(rowMaximum)), ...
    'scaledRowMaximumMinimum',min(scaledRowMaximum), ...
    'scaledRowMaximumMaximum',max(scaledRowMaximum), ...
    'rowScaleMinimum',min(rowScale),'rowScaleMaximum',max(rowScale), ...
    'rcond',rcond(factor),'preparationSeconds',preparationSeconds, ...
    'factorSeconds',factorSeconds,'factorClass',class(factor), ...
    'originalMatrixPreserved',true,'scaledMatrixCached',false);
if strcmp(method,'power2_row_max')
    restoredA = spdiags(1./rowScale,0,size(A,1),size(A,1))*scaledA;
    info.matrixRoundtripExact = isequaln(A,restoredA);
    assert(info.matrixRoundtripExact, ...
        'Power-of-two scaling changed matrix coefficients beyond exponent shifts.');
else
    info.matrixRoundtripExact = NaN;
end
end
