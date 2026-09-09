function X = ipm_perflab_sylvester_solve(cache,B)
%IPM_PERFLAB_SYLVESTER_SOLVE Solve the unchanged complete interior equation.
assert(isequal(size(B),cache.matrixSize),'Complete RHS has wrong size.');
switch cache.method
    case 'direct_cached'
        X = reshape(cache.factor\B(:),cache.matrixSize);
    case 'real_schur_cached'
        F = cache.leftVectors.'*B*cache.rightVectors;
        % Public sylvester accepts real quasi-triangular inputs. In this
        % MATLAB release it recognizes their Schur form and skips reduction.
        Z = sylvester(cache.leftSchur,cache.rightSchur,F);
        X = cache.leftVectors*Z*cache.rightVectors.';
    case 'fast_diagonalization'
        F = (cache.leftFactor\B)*cache.rightVectors;
        Z = F./cache.denominator;
        transformed = cache.leftVectors*Z;
        X = (cache.rightTransposeFactor\transformed.').';
    case 'matlab_sylvester'
        X = sylvester(cache.leftMatrix,cache.rightMatrix,B);
    otherwise
        error('ipm:SylvesterMethod','Unknown research Sylvester method.');
end
end
