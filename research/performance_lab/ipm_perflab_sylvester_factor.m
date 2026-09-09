function [cache,info] = ipm_perflab_sylvester_factor(Ty,Tx,method)
%IPM_PERFLAB_SYLVESTER_FACTOR Research factor for Ty*X+X*Tx.' = B.
%   Only one-dimensional dense matrices are used by the Schur/eigen paths.
%   No PDE discretization, boundary assembly, or acceptance gate changes.

arguments
    Ty double
    Tx double
    method char = 'real_schur_cached'
end
assert(isreal(Ty) && isreal(Tx) && size(Ty,1) == size(Ty,2) && ...
    size(Tx,1) == size(Tx,2),'Pass real square one-dimensional operators.');
m = size(Ty,1);
n = size(Tx,1);
cache = struct('method',method,'matrixSize',[m,n]);
info = struct('method',method,'matrixSize',[m,n], ...
    'factorSeconds',NaN,'storedNumericBytes',NaN,'decompositionBytesOpaque',false, ...
    'leftRelativeReconstruction',NaN,'rightRelativeReconstruction',NaN, ...
    'leftComponentwiseReconstruction',NaN,'rightComponentwiseReconstruction',NaN, ...
    'leftVectorsRcond',NaN,'rightVectorsRcond',NaN, ...
    'minimumAbsEigenvalueSum',NaN,'directRcond',NaN);
timer = tic;
switch method
    case 'direct_cached'
        A = kron(speye(n),Ty)+kron(Tx,speye(m));
        cache.factor = decomposition(A,'lu');
        info.factorSeconds = toc(timer);
        info.directRcond = rcond(cache.factor);
        info.decompositionBytesOpaque = true;
    case 'real_schur_cached'
        [cache.leftVectors,cache.leftSchur] = schur(full(Ty),'real');
        [cache.rightVectors,cache.rightSchur] = schur(full(Tx.'),'real');
        info.factorSeconds = toc(timer);
        [info.leftRelativeReconstruction,info.leftComponentwiseReconstruction] = ...
            reconstruction(Ty,cache.leftVectors*cache.leftSchur*cache.leftVectors.');
        [info.rightRelativeReconstruction,info.rightComponentwiseReconstruction] = ...
            reconstruction(Tx.',cache.rightVectors*cache.rightSchur*cache.rightVectors.');
        info.storedNumericBytes = 16*(m^2+n^2);
    case 'fast_diagonalization'
        [cache.leftVectors,cache.leftValues] = eig(full(Ty),'vector');
        [cache.rightVectors,cache.rightValues] = eig(full(Tx.'),'vector');
        cache.leftFactor = decomposition(cache.leftVectors,'lu');
        cache.rightTransposeFactor = decomposition(cache.rightVectors.','lu');
        cache.denominator = cache.leftValues+cache.rightValues.';
        info.factorSeconds = toc(timer);
        info.leftVectorsRcond = rcond(cache.leftVectors);
        info.rightVectorsRcond = rcond(cache.rightVectors);
        info.minimumAbsEigenvalueSum = min(abs(cache.denominator),[],'all');
        [info.leftRelativeReconstruction,info.leftComponentwiseReconstruction] = ...
            reconstruction(Ty,cache.leftVectors*diag(cache.leftValues)/cache.leftVectors);
        [info.rightRelativeReconstruction,info.rightComponentwiseReconstruction] = ...
            reconstruction(Tx.',cache.rightVectors*diag(cache.rightValues)/cache.rightVectors);
        info.decompositionBytesOpaque = true;
    case 'matlab_sylvester'
        cache.leftMatrix = full(Ty);
        cache.rightMatrix = full(Tx.');
        info.factorSeconds = toc(timer);
        info.storedNumericBytes = 8*(m^2+n^2);
    otherwise
        error('ipm:SylvesterMethod','Unknown research Sylvester method.');
end
end

function [relative,componentwise] = reconstruction(original,reconstructed)
residual = full(original)-reconstructed;
relative = norm(residual,inf)/max(norm(original,inf),realmin);
% Row-relative defects include fill generated at structural zeros, without
% claiming relative accuracy for exact-zero coefficients.
rowMagnitude = full(max(abs(original),[],2));
componentwise = max(abs(residual)./max(rowMagnitude,realmin),[],'all');
end
