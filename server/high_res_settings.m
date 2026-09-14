function settings = high_res_settings(projectRoot,n,maximumTotalNodes)
%HIGH_RES_SETTINGS Frozen 4n-by-n initial geometry with automatic REMESH.
validateattributes(n,{'numeric'}, ...
    {'scalar','real','finite','integer','>=',160});
assert(mod(n,4)==0,'ipm:HighResolutionAnchor', ...
    'n must be divisible by 4 so the uniform H8 root contains X=1 exactly.');
validateattributes(maximumTotalNodes,{'numeric'}, ...
    {'scalar','real','finite','integer','positive'});
settings = profile_settings(projectRoot);
settings.initialNodeCount = [4*n+1,n+1];
assert(maximumTotalNodes>=prod(settings.initialNodeCount), ...
    'ipm:HighResolutionNodeBudget', ...
    'The maximum total nodes must cover the initial 4n-by-n grid.');
settings.maximumTotalNodes = maximumTotalNodes;
settings.autonomousMeshVersion = 5;
settings.meshDensityVersion = 2;
settings.amplitudeGauge = 'outer_wall_density_window_l2';
end
